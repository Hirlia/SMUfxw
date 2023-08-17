#include "L0filter.h"

L0filter::L0filter(RichParameterSet* _para)
{
    cout << "L0 constructed!!" << endl;
    para = _para;
}

L0filter::~L0filter(void)
{
    cout << "L0 destroy!! " << endl;
}

void L0filter::setInput(DataMgr* pData)
{

    sample_ = pData->getCurrentSamples();
    cout<<"setInput"<<endl;
    //errpData = pData;
}

void L0filter::initial()
{
    total_iteration = para->getInt("Total Iteration");

    //knn
    normKnn = para->getInt("normKNN");
    norm_sparsity = para->getDouble("norm sparsity");
    norm_iteration = para->getInt("norm iteration");

    //vertex
    vertexKnn = para->getInt("vertKNN");
    vert_sparsity = para->getDouble("vert sparsity");
    vertex_iteration = para->getInt("vert iteration");
    vertex_stepsize_ = para->getDouble("Edge step size");

    //edge
    edgeKnn = para->getInt("EdgeKNN");
    edge_iteration = para->getInt("Edge Iteration");
    feature_scale_ = para->getDouble("Edge feature scale");
    edge_stepsize_ = para->getDouble("Edge step size");

    ptnum_ = sample_->VertexNumber();

    pnts_points_.clear();
    pnts_points_.resize(ptnum_);
    pnts_normals_.clear();
    pnts_normals_.resize(ptnum_);
    //nInit.resize(ptnum_,3);
    CMesh::VertexIterator vi;
    for (vi = sample_->vert.begin(); vi != sample_->vert.end(); ++vi)
    {
        vi->normal = vcg::Point3d(vi->N()[0],vi->N()[1],vi->N()[2]);
        pnts_points_[vi->m_index] = Eigen::Vector3d(vi->P()[0],vi->P()[1],vi->P()[2]);
        pnts_normals_[vi->m_index] = Eigen::Vector3d(vi->N()[0],vi->N()[1],vi->N()[2]);
        pnts_normals_[vi->m_index].normalize();

    }
    last_result_points_ = pnts_points_;


}

void L0filter::run()
{
    cout<<"run"<<endl;

    initial();


    for (int i = 0; i < total_iteration; ++i)
    {
        norm_auxiliary = para->getDouble("norm auxiliary");
        norm_increment = para->getDouble("norm increment");
        vert_auxiliary = para->getDouble("vert auxiliary");
        vert_increment = para->getDouble("vert increment");


        if (para->getBool("norm L0"))
        {
            //norm_auxiliary *=pow(1.2,(double)i);

            NormalFiltering(pnts_normals_);
            //getnewneighbor();
            //NormalFilteringplus();
            //NormalFiltering_plus(i);

            cout<<i<<endl;
        }
        else
        {
            cout<<"Do not Filter Normals"<<endl;
        }
        if (para->getBool("vert L0" ))
        {
            //vert_auxiliary *=pow(1.2,(double)i);
            //VertexUpdating(vertex_stepsize_,pnts_normals_,pnts_points_);
            vertexupdating_L0(pnts_normals_,pnts_points_);
        }
        else
        {
            cout<<"Do not Update Vertices"<<endl;
        }
        if (para->getBool("edge recovery" ))
        {
            recoveryEdge(feature_scale_,pnts_normals_,pnts_points_);
            //vertexupdating_L0(pnts_normals_,pnts_points_);
            //NormalFilteringplus(i);
        }
        else
        {
            cout<<"Do not Recovery Edges"<<endl;
        }
    }
}

void L0filter::NormalFiltering(std::vector<Eigen::Vector3d>& pnts_normals)
{
    int itr = 0;
    std::vector<Eigen::Vector3d> thetas(ptnum_*normKnn, Eigen::Vector3d::Zero());//n*k维
    //Eigen::MatrixXd thetas(ptnum_,normKnn);

    CMesh::VertexIterator vi;

    //input_pnts_normals_
    input_pnts_normals_.clear();
    for (int i=0; i < ptnum_; ++i)
    {
        Eigen::Vector3d n(pnts_normals[i][0], pnts_normals[i][1],
                pnts_normals[i][2]);
        input_pnts_normals_.push_back(n);
    }

    GlobalFun::computeAnnNeigbhors(sample_->vert, sample_->vert, normKnn, false, "L0 Algorithm");
    vertex_neighbors_.clear();
    vertex_neighbors_.resize(ptnum_);
    for (int i=0; i < ptnum_; ++i)
    {
            vertex_neighbors_[i] = sample_->vert[i].neighbors;
    }

    find_include_i_knn_idx(vertex_neighbors_,include_i_knn_idx,normKnn);

    do
    {
        //1.
        solve_theta_problem(pnts_normals, thetas);

        //2.
        solve_N_problem(thetas, pnts_normals);

        //3.
        norm_auxiliary *= norm_increment;

    }
    while (++itr < norm_iteration && norm_auxiliary < 5000);
    std::cout << "##  Finished " << itr <<" times Filtering Normals"<< "\n";


    for(vi = sample_->vert.begin(); vi != sample_->vert.end(); ++vi)
    {
        vi->N()[0] = pnts_normals[vi->m_index][0];
        vi->N()[1] = pnts_normals[vi->m_index][1];
        vi->N()[2] = pnts_normals[vi->m_index][2];

    }
}

void L0filter::find_include_i_knn_idx(const std::vector<std::vector<int>>& vertex_neighbors,											std::vector<std::vector<int>>& include_i_knn_idx,int Knn)
{
    include_i_knn_idx.clear();
    include_i_knn_idx.resize(ptnum_);

    int indx;
    for (size_t k = 0; k < ptnum_; ++k)
    {
        for (int j = 0; j < Knn; ++j)
        {
            indx = vertex_neighbors[k][j];
            include_i_knn_idx[indx].push_back(k);
        }
    }
}

void L0filter::solve_theta_problem(const vector<Eigen::Vector3d>& pnts_normals,std::vector<Eigen::Vector3d>& thetas)
{
    double theshold = norm_sparsity / norm_auxiliary;
    Eigen::Vector3d DN;
    for (size_t i = 0; i < ptnum_; ++i)
    {
        for(size_t j = 0; j < normKnn; ++j)
        {
            DN = pnts_normals[i] - pnts_normals[vertex_neighbors_[i][j]];

            if (DN.squaredNorm() >= theshold )
            {
                thetas[i*normKnn+j] = DN;
            }
            else
            {
                thetas[i*normKnn+j].setZero();
            }
        }
    }
}

void L0filter::solve_N_problem(const std::vector<Eigen::Vector3d>& thetas, vector<Eigen::Vector3d>& pnts_normals)
{
    std::vector<Eigen::Triplet<double>> coeff_triple; coeff_triple.clear();
    std::vector<Eigen::Triplet<double>> coeff_triple_ii; coeff_triple_ii.clear();
    std::vector<Eigen::Triplet<double>> coeff_triple_ij; coeff_triple_ij.clear();
    std::vector<Eigen::Triplet<double>> coeff_triple_ik; coeff_triple_ik.clear();
    Eigen::SparseMatrix<double> weight_matrix(ptnum_, ptnum_);
    Eigen::SparseMatrix<double> weight_matrix_ii(ptnum_, ptnum_);
    Eigen::SparseMatrix<double> weight_matrix_ij(ptnum_, ptnum_);
    Eigen::SparseMatrix<double> weight_matrix_ik(ptnum_, ptnum_);

    double coi;
    std::vector<int> vk_neighbors;
    std::vector<int>::iterator it;
    std::vector<int> vi_neighbors;
    std::vector<int> include_i_neighbors;
    int in_i_nei_size;

    for (size_t i = 0; i < ptnum_; ++i)
    {
        vi_neighbors = vertex_neighbors_[i];
        include_i_neighbors = include_i_knn_idx[i];
        in_i_nei_size = include_i_neighbors.size();
        int normKnn_size = normKnn;

        for (size_t j = 0; j < normKnn; ++j)
                coeff_triple_ij.push_back( Eigen::Triplet<double>(i, vi_neighbors[j], -norm_auxiliary) );
        for (size_t j = 0; j < include_i_neighbors.size(); ++j)
                coeff_triple_ik.push_back( Eigen::Triplet<double>(i, include_i_neighbors[j], -norm_auxiliary));

        coi = 1 + norm_auxiliary*normKnn_size + norm_auxiliary*in_i_nei_size;
        coeff_triple_ii.push_back( Eigen::Triplet<double>(i, i, coi));
    }

    weight_matrix_ii.setFromTriplets(coeff_triple_ii.begin(), coeff_triple_ii.end());
    weight_matrix_ij.setFromTriplets(coeff_triple_ij.begin(), coeff_triple_ij.end());
    weight_matrix_ik.setFromTriplets(coeff_triple_ik.begin(), coeff_triple_ik.end());
    weight_matrix = weight_matrix_ii + weight_matrix_ij + weight_matrix_ik;
    weight_matrix.makeCompressed();

    Eigen::Vector3d sum;
    Eigen::Vector3d theta_ij;
    Eigen::Vector3d theta_ki;
    Eigen::MatrixXd right_term(ptnum_, 3);

    for (size_t i = 0; i < ptnum_; ++i)
    {
        include_i_neighbors = include_i_knn_idx[i];
        in_i_nei_size = include_i_neighbors.size();
        sum = input_pnts_normals_[i];
        theta_ij = Eigen::Vector3d(0,0,0);
        theta_ki = Eigen::Vector3d(0,0,0);

        for (size_t j = 0; j < normKnn; ++j)
        {
            theta_ij += norm_auxiliary * thetas[i*normKnn+j];
        }
        for (size_t j = 0; j < in_i_nei_size; ++j)
        {
            vk_neighbors = vertex_neighbors_[include_i_neighbors[j]];
            it = std::find(vk_neighbors.begin(),vk_neighbors.end(),i);
            int k = std::distance(vk_neighbors.begin(),it);
            theta_ki -= norm_auxiliary * thetas[include_i_neighbors[j]*normKnn+k];
        }
        right_term.row(i) = sum + theta_ij + theta_ki;
    }

    //solve
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>   cgSolver_;
    cgSolver_.compute(weight_matrix);
    Eigen::MatrixX3d filtered_normals_matrix = cgSolver_.solve(right_term);
    filtered_normals_matrix.rowwise().normalize();

    for (size_t i = 0; i < ptnum_; ++i)
    {
        pnts_normals[i] = filtered_normals_matrix.row(i);
    }
}



void L0filter::VertexUpdating(const double stepsize,const vector<Eigen::Vector3d>& result_normals,vector<Eigen::Vector3d>& result_points)
{

    //存初始顶点
    vector<Eigen::Vector3d> oripoints;
    oripoints = result_points;

    //knn
    //std::vector<std::vector<int>> vertex_neighbors_point;
    GlobalFun::computeAnnNeigbhors(sample_->vert, sample_->vert, vertexKnn, false, "L0 verticesUpdate");
    neighbors_v_.clear();
    neighbors_v_.resize(ptnum_);
    for (int i=0; i < ptnum_; ++i)
    {
        neighbors_v_[i] = sample_->vert[i].neighbors;
    }

    vector<int> vertex_neighbor;
    Eigen::Vector3d fnormal,pij;
    //cout<<"***********Point Filtering***********"<<endl;
    for (int n = 0; n < vertex_iteration; ++n)
    {
        for (int i = 0; i < ptnum_; ++i)
        {
            Eigen::Vector3d tc = Eigen::Vector3d(0,0,0);
            vertex_neighbor.clear();
            vertex_neighbor = neighbors_v_[i];
            size_t num = vertex_neighbor.size();
            fnormal = result_normals[i];
            for (size_t j = 0; j < num; ++j)
            {
                tc += fnormal*(fnormal.dot(result_points[vertex_neighbor[j]] - result_points[i]));
            }
            pij = result_points[i] + stepsize*tc;
            result_points[i] = pij;
        }

    }
    cout<<"##  Finish "<<vertex_iteration<<"times Point filtering"<<endl;

    CMesh::VertexIterator vi;
    for(vi = sample_->vert.begin(); vi != sample_->vert.end(); ++vi)
    {
        vi->P()[0] = result_points[vi->m_index][0];
        vi->P()[1] = result_points[vi->m_index][1];
        vi->P()[2] = result_points[vi->m_index][2];
    }

}

//-------------------------------------------------------------------------------------------------
void L0filter::vertexupdating_L0(const vector<Eigen::Vector3d>& result_normals,vector<Eigen::Vector3d>& result_points)
{
    input_pnts_points_.clear();
    for (int i=0; i < ptnum_; ++i)
    {
        Eigen::Vector3d p(result_points[i][0], result_points[i][1],
               result_points[i][2]);
        input_pnts_points_.push_back(p);
    }

    GlobalFun::computeAnnNeigbhors(sample_->vert, sample_->vert, vertexKnn, false, "L0 verticesUpdate");
    neighbors_v_.clear();
    neighbors_v_.resize(ptnum_);
    for (int i=0; i < ptnum_; ++i)
    {
        neighbors_v_[i] = sample_->vert[i].neighbors;
    }
    find_include_i_knn_idx(neighbors_v_,neighbors_includeI_,vertexKnn);

    vector<Eigen::Vector3d> lastpoints = result_points; //初始化点为噪声点云
    vector<double> theta(ptnum_*vertexKnn,0);
    vector<double> alpha(ptnum_,0);

    int iter = 0;
    do
    {
        L0_vert_theta(result_normals,lastpoints,alpha,theta);

        L0_vert_alpha(result_normals,theta,lastpoints,alpha);

        //increase
        vert_auxiliary *= vert_increment;

    } while (++iter < vertex_iteration);

    CMesh::VertexIterator vi;
    for(vi = sample_->vert.begin(); vi != sample_->vert.end(); ++vi)
    {
        vi->P()[0] = lastpoints[vi->m_index][0];
        vi->P()[1] = lastpoints[vi->m_index][1];
        vi->P()[2] = lastpoints[vi->m_index][2];
    }
    result_points = lastpoints;
}

void L0filter::L0_vert_theta(const vector<Eigen::Vector3d>& result_normals,const vector<Eigen::Vector3d>& lastpoints,const vector<double>& alpha,vector<double>& theta)
{
    double theshold = vert_sparsity / vert_auxiliary;
    double temp, squtemp;
    for (size_t i = 0; i < ptnum_; ++i)
    {
        for(size_t j = 0; j < vertexKnn; ++j)
        {
            temp = (lastpoints[i] - lastpoints[neighbors_v_[i][j]]).dot(result_normals[i]);

            squtemp = pow(temp,2);

            if (squtemp >= theshold)
            {
                theta[i*vertexKnn+j] = temp;
            }
            else
            {
                theta[i*vertexKnn+j] = 0;
            }
        }
    }
}

void L0filter::L0_vert_alpha(const vector<Eigen::Vector3d>& result_normals,	const vector<double>& theta,vector<Eigen::Vector3d>& lastpoints,vector<double>& alpha)
{
    std::vector<Eigen::Triplet<double>> coeff_triple; coeff_triple.clear();
    std::vector<Eigen::Triplet<double>> coeff_triple_ii; coeff_triple_ii.clear();
    std::vector<Eigen::Triplet<double>> coeff_triple_ij; coeff_triple_ij.clear();
    std::vector<Eigen::Triplet<double>> coeff_triple_ik1; coeff_triple_ik1.clear();
    std::vector<Eigen::Triplet<double>> coeff_triple_ik2; coeff_triple_ik2.clear();
    Eigen::SparseMatrix<double> weight_matrix(ptnum_, ptnum_);
    Eigen::SparseMatrix<double> weight_matrix_ii(ptnum_, ptnum_);
    Eigen::SparseMatrix<double> weight_matrix_ij(ptnum_, ptnum_);
    Eigen::SparseMatrix<double> weight_matrix_ik1(ptnum_, ptnum_);
    Eigen::SparseMatrix<double> weight_matrix_ik2(ptnum_, ptnum_);
    size_t i,j,k;
    double coe_i, coe_ij, coe_ik1,coe_ik2;
    int include_i_size;

    for (i = 0; i < ptnum_; ++i)
    {
        include_i_size = neighbors_includeI_[i].size();
        coe_i = 1.0 + vertexKnn * vert_auxiliary;
        coeff_triple_ii.push_back( Eigen::Triplet<double>(i, i, coe_i));

        coe_ik1 = 0;
        coe_ik2 = 0;
        for (k = 0; k < include_i_size; ++k)
        {
            coe_ik1 += vert_auxiliary*pow(result_normals[neighbors_includeI_[i][k]].dot(result_normals[i]),2);

            coe_ik2 += vert_auxiliary*result_normals[neighbors_includeI_[i][k]].dot(result_normals[i]);
        }
        coeff_triple_ik1.push_back( Eigen::Triplet<double>(i, i, coe_ik1));
        coeff_triple_ik2.push_back( Eigen::Triplet<double>(i, i, -coe_ik2));
        for (j = 0; j < neighbors_v_[i].size(); ++j)
        {
            coe_ij = vert_auxiliary*result_normals[neighbors_v_[i][j]].dot(result_normals[i]);
            coeff_triple_ij.push_back( Eigen::Triplet<double>(i, neighbors_v_[i][j], -coe_ij));
        }
    }

    weight_matrix_ii.setFromTriplets(coeff_triple_ii.begin(),coeff_triple_ii.end());
    weight_matrix_ij.setFromTriplets(coeff_triple_ij.begin(),coeff_triple_ij.end());
    weight_matrix_ik1.setFromTriplets(coeff_triple_ik1.begin(),coeff_triple_ik1.end());
    weight_matrix_ik2.setFromTriplets(coeff_triple_ik2.begin(),coeff_triple_ik2.end());
    weight_matrix = weight_matrix_ii + weight_matrix_ij + weight_matrix_ik1 + weight_matrix_ik2;
    weight_matrix.makeCompressed();

    //right term
    Eigen::MatrixXd right_term(ptnum_, 1);
    std::vector<int> vk_neighbors;
    std::vector<int>::iterator it;

    double sumF, sumJ, sumK;
    Eigen::Vector3d	Dpik;
    for (i = 0; i < ptnum_; ++i)
    {
        sumF = (input_pnts_points_[i] - lastpoints[i]).dot(result_normals[i]);

        sumJ = 0;
        for (j = 0; j < neighbors_v_[i].size(); ++j)
        {
            sumJ += vert_auxiliary*((lastpoints[i] - lastpoints[neighbors_v_[i][j]]).dot(result_normals[i])-theta[i*vertexKnn+j]);
        }

        sumK = 0;
        for (k = 0; k < neighbors_includeI_[i].size(); ++k)
        {
            vk_neighbors = neighbors_v_[neighbors_includeI_[i][k]];
            it = std::find(vk_neighbors.begin(),vk_neighbors.end(),i);
            int kk = std::distance(vk_neighbors.begin(),it);

            Dpik = lastpoints[neighbors_includeI_[i][k]] - lastpoints[i];

            sumK += vert_auxiliary*(result_normals[neighbors_includeI_[i][k]].dot(result_normals[i]))*(-Dpik.dot(result_normals[neighbors_includeI_[i][k]])+theta[neighbors_includeI_[i][k]*vertexKnn+kk]);
        }
        right_term(i,0) = - sumF - sumJ - sumK;
    }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>   cgSolver_;
    cgSolver_.compute(weight_matrix);
    Eigen::MatrixXd result_alpha_matrix = cgSolver_.solve(right_term);
    for (i = 0; i < ptnum_; ++i)
    {
        alpha[i] = result_alpha_matrix(i);
        lastpoints[i] = lastpoints[i] + alpha[i]*result_normals[i];
    }
}

void L0filter::recoveryEdge(double feature_scale, vector<Eigen::Vector3d>& result_normals, vector<Eigen::Vector3d>& result_points)
{
    CMesh::VertexIterator vi;

    if (result_points.empty())
    {
        result_points.resize(ptnum_);
        for(vi = sample_->vert.begin(); vi != sample_->vert.end(); ++vi)
        {
            result_points[vi->m_index][0] = vi->P()[0];
            result_points[vi->m_index][1] = vi->P()[1];
            result_points[vi->m_index][2] = vi->P()[2];
        }
    }

    if (result_normals.empty())
    {
        result_normals.resize(ptnum_);
        for(vi = sample_->vert.begin(); vi != sample_->vert.end(); ++vi)
        {
            result_normals[vi->m_index][0] = vi->N()[0];
            result_normals[vi->m_index][1] = vi->N()[1];
            result_normals[vi->m_index][2] = vi->N()[2];
        }
    }

    vector<Eigen::Vector3d> putspoints;
    putspoints = result_points;

    vector<vector<int> > i_neighbor;i_neighbor.resize(ptnum_);
    vector<bool> feature_point;feature_point.resize(ptnum_);

    GlobalFun::computeAnnNeigbhors(sample_->vert, sample_->vert, edgeKnn, false, "Edge Recovery");
    for(size_t i = 0; i < ptnum_; ++i)
    {
        i_neighbor[i] = sample_->vert[i].neighbors;
        feature_point[i] = sample_->vert[i].feature;
    }

    Eigen::Vector3d sum;

    Eigen::MatrixXd right_term(ptnum_, 3);
    Eigen::Matrix3d A;
    for(int i=0;i<ptnum_;i++)
    {
        right_term.row(i) = putspoints[i] ;
        if(feature_point[i])
        {
            sum = Eigen::Vector3d(0,0,0);
            A.setZero();
            for(int j=0;j<edgeKnn;j++)
            {
                int index_j = i_neighbor[i][j];
                sum +=(putspoints[index_j].dot(result_normals[index_j]))*result_normals[index_j];
                for(int k=0;k<3;k++)
                    for(int t=0;t<3;t++)
                        A(k,t)+=result_normals[index_j][k]*result_normals[index_j][t];
            }
            A(0,0) += 1.0;
            A(1,1) += 1.0;
            A(2,2) += 1.0;
            right_term.row(i) += sum;
            right_term.row(i) *= A.inverse();
        }
        result_points[i] = right_term.row(i);
    }

    for(vi = sample_->vert.begin(); vi != sample_->vert.end(); ++vi)
    {
        vi->P()[0] = result_points[vi->m_index][0];
        vi->P()[1] = result_points[vi->m_index][1];
        vi->P()[2] = result_points[vi->m_index][2];
    }

    //cout<<"## finished"<<itre<<"times Edge iterations"<<endl;
}

//---------------------------------------------------------------------------------------------
void L0filter::vertexupdating_L0_plus(const vector<Eigen::Vector3d> &result_normals, vector<Eigen::Vector3d> &result_points)
{
    input_pnts_points_.clear();
    for (int i=0; i < ptnum_; ++i)
    {
        Eigen::Vector3d p(result_points[i][0], result_points[i][1],
               result_points[i][2]);
        input_pnts_points_.push_back(p);
    }

    GlobalFun::computeAnnNeigbhors(sample_->vert, sample_->vert, vertexKnn, false, "L0 verticesUpdate");
    neighbors_v_.clear();
    neighbors_v_.resize(ptnum_);
    for (int i=0; i < ptnum_; ++i)
    {
        neighbors_v_[i] = sample_->vert[i].neighbors;
    }
    find_include_i_knn_idx(neighbors_v_,neighbors_includeI_,vertexKnn);
    //find_include_i_knn_idx(vertex_neighbors_,include_i_knn_idx,normKnn);


    vector<Eigen::Vector3d> lastpoints = result_points; //初始化点为噪声点云

    std::vector<Eigen::Triplet<double>> coeff_triple; coeff_triple.clear();
    std::vector<Eigen::Triplet<double>> coeff_triple_ii; coeff_triple_ii.clear();
    std::vector<Eigen::Triplet<double>> coeff_triple_ij; coeff_triple_ij.clear();
    std::vector<Eigen::Triplet<double>> coeff_triple_ik; coeff_triple_ik.clear();
    Eigen::SparseMatrix<double> weight_matrix(ptnum_, ptnum_);
    Eigen::SparseMatrix<double> weight_matrix_ii(ptnum_, ptnum_);
    Eigen::SparseMatrix<double> weight_matrix_ij(ptnum_, ptnum_);
    Eigen::SparseMatrix<double> weight_matrix_ik(ptnum_, ptnum_);

    double coi;
    std::vector<int> vk_neighbors;
    std::vector<int>::iterator it;
    std::vector<int> vi_neighbors;
    std::vector<int> include_i_neighbors;
    int in_i_nei_size;

    for (size_t i = 0; i < ptnum_; ++i)
    {
        vi_neighbors = neighbors_v_[i];
        include_i_neighbors = neighbors_includeI_[i];
        in_i_nei_size = include_i_neighbors.size();
        int vert_size = vertexKnn;
        if(!sample_->vert[i].feature) vert_size = 0;
        else
        {
            for (size_t j = 0; j < vertexKnn; ++j)
            {
                    coeff_triple_ij.push_back( Eigen::Triplet<double>(i, vi_neighbors[j], -0.1) );
            }
        }
        for (size_t j = 0; j < include_i_neighbors.size(); ++j)
        {
            if(!sample_->vert[include_i_neighbors[j]].feature)
                in_i_nei_size--;
            else{
                    coeff_triple_ik.push_back( Eigen::Triplet<double>(i, include_i_neighbors[j], -0.1));
            }
        }
        coi = 1 + 0.1*vert_size + 0.1*in_i_nei_size;
        coeff_triple_ii.push_back( Eigen::Triplet<double>(i, i, coi));
    }

    weight_matrix_ii.setFromTriplets(coeff_triple_ii.begin(), coeff_triple_ii.end());
    weight_matrix_ij.setFromTriplets(coeff_triple_ij.begin(), coeff_triple_ij.end());
    weight_matrix_ik.setFromTriplets(coeff_triple_ik.begin(), coeff_triple_ik.end());
    weight_matrix = weight_matrix_ii + weight_matrix_ij + weight_matrix_ik;
    weight_matrix.makeCompressed();

    Eigen::Vector3d sum;
    Eigen::MatrixXd right_term(ptnum_, 3);

    for (size_t i = 0; i < ptnum_; ++i)
    {
        sum = input_pnts_points_[i];
        right_term.row(i) = sum ;
    }

    //solve
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>   cgSolver_;
    cgSolver_.compute(weight_matrix);
    Eigen::MatrixX3d filtered_normals_matrix = cgSolver_.solve(right_term);

    for (size_t i = 0; i < ptnum_; ++i)
    {
        lastpoints[i] = filtered_normals_matrix.row(i);
    }

    CMesh::VertexIterator vi;
    for(vi = sample_->vert.begin(); vi != sample_->vert.end(); ++vi)
    {
        vi->P()[0] = lastpoints[vi->m_index][0];
        vi->P()[1] = lastpoints[vi->m_index][1];
        vi->P()[2] = lastpoints[vi->m_index][2];
    }
    result_points = lastpoints;

}

void L0filter::getnewneighbor()
{
    for(int i=0;i<sample_->vn;++i){
        if(sample_->vert[i].feature){
            int k = sample_->vert[i].neighbors.size();
            vcg::Point3f ni = sample_->vert[i].N();
            vector<int> newneighbors;
            for(int j=0;j<k;++j){
                int index = sample_->vert[i].neighbors[j];
                vcg::Point3f nj = sample_->vert[index].N();
                if(sample_->vert[index].feature && ni * nj > 0.95){
                    newneighbors.push_back(index);
                }
               // sample_->vert[i].neighbors = newneighbors;
            }
            sample_->vert[i].neighbors = newneighbors;
        }
    }
}

void L0filter::NormalFilteringplus()
{
    double Rho = 0.95;
    double tao = 0.3;
    int iteration = 20;
    int vn = sample_->vn;
    vector<Point3f> n(vn);
    for(int i = 0; i < vn; ++i){
        n[i] = sample_->vert[i].N();
    }
    while(iteration--){
        for (int i = 0; i < vn; i++)
        {
            vcg::Matrix33<float> covariance_matrix;
            covariance_matrix.SetZero();
            double weightsum = 0.0;

            int k = sample_->vert[i].neighbors.size();
            for(int j = 0; j < k; ++j)
            {
                int nlnd = sample_->vert[i].neighbors[j];

                vcg::Point3f Nj = n[nlnd];

                double angle = n[i] * Nj;
                vcg::Matrix33<float> tmp;
                tmp.ExternalProduct(Nj , Nj);//Nj * Nj^T

                if(angle > Rho)weightsum += 1.0;
                else{
                    tmp *= 0.00001;
                    weightsum += 0.00001;
                }
                covariance_matrix += tmp;
            }

            if(weightsum != 0.0)
                covariance_matrix /= weightsum;
            vcg::Point3f eigenvalues;
            vcg::Matrix33<float> eigenvectors;
            int required_rotations;
            vcg::Jacobi<vcg::Matrix33<float>, vcg::Point3f>(covariance_matrix, eigenvalues, eigenvectors, required_rotations);
            vcg::SortEigenvaluesAndEigenvectors<vcg::Matrix33<float>, vcg::Point3f>(eigenvalues, eigenvectors);

            if(eigenvalues[2] < tao)
                eigenvalues[2] = 0.0;
            else eigenvalues[2] = 1.0;
            if(eigenvalues[1] < tao)
                eigenvalues[1] = 0.0;
            else eigenvalues[1] = 1.0;

            eigenvalues[0] = 1.0;

            vcg::Matrix33<float> Ti,tmp_;
            Ti.SetZero();
            vcg::Point3f x;//tmp的特征向量
            for(int j = 0; j < 3; j++)
            {
                x = eigenvectors.GetColumn(j);
                tmp_.ExternalProduct(x,x);
                Ti += tmp_ * eigenvalues[j];
            }

            n[i] = n[i] * 3 +  Ti * n[i];
            n[i].Normalize();
        }
    }
    for(int i = 0; i < vn; ++i){
        sample_->vert[i].N() = n[i];
    }
}

double L0filter::averageKNNDistance()
{
    double dist = 0;
    int k = 6;
    for(int i=0;i<ptnum_;i++){
        for(int j=0;j<k;j++){
            dist += (sample_->vert[i].P()-sample_->vert[sample_->vert[i].neighbors[j]].P()).Norm();
        }
    }
    dist /= k*ptnum_;
    return dist;
}

void L0filter::NormalFiltering_plus(int iteration)
{
    for(int i=0;i<ptnum_;i++)
    {
        vcg::Point3f ni = sample_->vert[i].N();
        double Rho;
        double tao = sample_->vert[i].sigma;
        //std::cout<<tao<<std::endl;
//        if(!sample_->vert[i].feature)
        Rho = 0.95 /*+ 0.5 * exp(-5 * tao)*/;
        //else
        //Rho = 0.55 + 0.4 * exp(-80.88 * tao * tao);
        //std::cout<<tao<<" "<<Rho<<std::endl;
        vcg::Point3f normal;
        vector<vcg::Point3f> locpoints;
        int k = sample_->vert[i].neighbors.size();
        for(int j=0;j<k;j++)
        {
            int indexj = sample_->vert[i].neighbors[j];
            vcg::Point3f nj = sample_->vert[indexj].N();
            if(ni*nj>Rho)
            {
                vcg::Point3f pj = sample_->vert[indexj].P();
                locpoints.push_back(pj);
            }
        }

        if(locpoints.size()!=0){
            vcg::Point3f center;
            vcg::Matrix33<float> covariance_matrix;
            covariance_matrix.SetZero();
            covariance_matrix.Covariance(locpoints,center);//计算协方差

            vcg::Point3f eigenvalues;
            vcg::Matrix33<float> eigenvectors;
            int required_rotations;
            vcg::Jacobi<vcg::Matrix33<float>, vcg::Point3f>(covariance_matrix, eigenvalues, eigenvectors, required_rotations);
            vcg::SortEigenvaluesAndEigenvectors<vcg::Matrix33<float>, vcg::Point3f>(eigenvalues, eigenvectors);

            normal = eigenvectors.GetColumn(2);

            normal.Normalize();
            if(normal*ni<0)
                normal *= -1;
        }
        normal += ni ;
        sample_->vert[i].N() = normal.Normalize();

    }
}
