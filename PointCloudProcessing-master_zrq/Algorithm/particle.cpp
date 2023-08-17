#include "particle.h"

#define LAMDA 20
#define myrandom(a) (2*a*rand()/double(RAND_MAX)-a)//生成-a~b的浮点随机数
typedef typename vcg::Matrix33<float> MatrixType;

Particle::Particle(RichParameterSet *_para)
{
    cout << "ParticleFilter constructed!!" << endl;
    para = _para;
}

Particle::~Particle(void)
{
    cout << "L0 destroy!! " << endl;
}

void Particle::setInput(DataMgr* pData)
{

    sample_ = pData->getCurrentSamples();
    origin_ = pData->getCurrentOriginal();
    cout<<"setInput"<<endl;
}

void Particle::initial()
{
    total_iteration = para->getInt("total iteration"); //总迭代次数
    ptnum_ = sample_->VertexNumber();//点云数目
    //计算球邻域（不包括自己）
    radius = global_paraMgr.norSmooth.getDouble("CGrid Radius");
    particlenum = para->getInt("particle number");
    double radius2 = radius * radius;
    for(CMesh::VertexIterator dest = sample_->vert.begin(); dest != sample_->vert.end(); dest++)
    {
        dest->neighbors.clear();
    }
    for(CMesh::VertexIterator dest = sample_->vert.begin(); dest != sample_->vert.end(); dest++)
    {
        Point3f p = dest->P();
        for(CMesh::VertexIterator origin = dest + 1; origin != sample_->vert.end(); origin++)
        {
            Point3f q = origin->P();
            Point3f diff = p - q;
            double dist2 = diff.SquaredNorm();
            if (dist2 < radius2)
            {
                dest->neighbors.push_back(origin->m_index);
                origin->neighbors.push_back(dest->m_index);
            }
        }
    }
//    for(int i=0;i<ptnum_;i++){
//        for(int j=0;j<sample_->vert[i].neighbors.size();j++){
//            sample_->vert[i].weight.push_back(1.0);
//        }
//    }
    //计算带宽
    bandwidth = calculate_bandwidth();
    bandwidth_.resize(ptnum_);
    for(int i=0;i<ptnum_;i++)
        bandwidth_[i] = bandwidth;
}

void Particle::run()
{
    cout<<"run"<<endl;

    initial();

    for(int i=0; i<total_iteration; i++)
    {

        if(para->getBool("edge"))
        {
            runparticleedge();
        }
        else if(para->getBool("normal"))
        {
            //runRANSAC();
            run_normal();
            //runnormalparticle();
        }
        else
        {
            beta = pow(1.2,i);
            runparticle();
        }

    }
}

bool cmp(pair<vcg::Point3f,double> a,pair<vcg::Point3f,double> b)
{
    return a.second < b.second;
}

void Particle::runparticle()
{
    srand((int)time(0));

    vector<vcg::Point3f> points;
    CMesh::VertexIterator iv;
    for(iv=sample_->vert.begin();iv!=sample_->vert.end();iv++)
    {
        //vcg::Point3f orip = origin_->vert[iv->m_index].P();
        vcg::Point3f pi = iv->P();
        vcg::Point3f ni = iv->N();
        int index = iv->m_index;

        vector<pair<vcg::Point3f,double> > particles(particlenum);
        for(int i=0;i<particlenum;i++)//随机生成particlenum个粒子
        {
            vcg::Point3f pm = pi + ni * myrandom(bandwidth_[index]);
            particles[i].first = pm;
            double cov = 0,weight;
            int k =iv->neighbors.size();
            for(int j=0;j<k;j++)//计算与邻域点的方差
            {
                vcg::Point3f pj = sample_->vert[iv->neighbors[j]].P();
                vcg::Point3f nj = sample_->vert[iv->neighbors[j]].N();

                weight = iv->weight[j]<0.3? 0.0:iv->weight[j] ;
                weight = iv->weight[j]>0.7? 1.0:iv->weight[j] ;
                cov += (pow((pm - pj).Normalize()*ni,2) + pow((pj - pm).Normalize()*nj,2))  * weight;
            }
            particles[i].second = cov/*beta * cov  + (orip - pm).Norm()*/;
        }

        sort(particles.begin(),particles.end(),cmp);//选择方差最小的粒子点作为先验
        points.push_back(particles[0].first);
        double band1 = (pi - particles[0].first).Norm();
        double band2 = bandwidth_[index]-(pi - particles[0].first).Norm();
        bandwidth_[index] = band1>band2?band2:band1;

    }

    int ii=0;
    for(iv=sample_->vert.begin();iv!=sample_->vert.end();iv++,ii++)
    {
        iv->P() = points[ii];
    }
}

void Particle::runnormalparticle()
{
    vcg::Point3f ni,pi;
    double omiga_f,omiga_d,omiga_n;
    vector<vcg::Point3f> nk;
    nk.resize(sample_->vn);
    for(int i=0;i<sample_->vn;i++)
        nk[i] = sample_->vert[i].N();//记录法向信息

    vector<vcg::Point3f> initpoint;
    initpoint.resize(sample_->vn);
    for(int i=0;i<sample_->vn;i++)
        initpoint[i] = sample_->vert[i].P();//记录点信息

    vector<vector<int> > neighbors;
    neighbors.resize(sample_->vn);
    for(int i=0;i<sample_->vn;i++)
        neighbors[i] = sample_->vert[i].neighbors;

    vector<double> feat;
    feat.resize(sample_->vn);
    for(int i=0;i<sample_->vn;i++)
        feat[i]=sample_->vert[i].sigma_plus;

    double feat_threshold = global_paraMgr.particle.getDouble("feat");
    for(int i=0;i<sample_->vn;i++)
    {
        ni = nk[i];
        pi = initpoint[i];
        int num = neighbors[i].size();
        vector<double> weight;

        double mean_n = 0;
        for(int j=0; j<num; j++)
        {
            vcg::Point3f nj = nk[neighbors[i][j]];
            mean_n += (ni-nj).Norm();
        }
        double sigma_n = mean_n / num;

        //cout<<sigma_n<<" "<<feat<<endl;
        double feature_scale_ = global_paraMgr.l0filter.getDouble("Edge feature scale");
        for(int j=0; j<num; j++)
        {
            vcg::Point3f nj = nk[neighbors[i][j]];
            vcg::Point3f pj = initpoint[neighbors[i][j]];
            //omiga_d = exp(-pow((pj-pi).Norm()/(radius),2));
            omiga_n = (nj-ni).Norm()>feature_scale_?1e-4:1;
            omiga_f = exp(-feat[neighbors[i][j]]/feat_threshold);
            //cout<<omiga_d<<endl;
            weight.push_back(omiga_f  /**omiga_d*/ * omiga_n);
        }

        double max_weight = *max_element(weight.begin(),weight.end());
        for(int j=0; j<num; j++)
            weight[j] /= max_weight;

        sample_->vert[i].weight = weight;
        for(int j=0;j<num;j++){
            vcg::Point3f nj = nk[neighbors[i][j]];
            ni += nj*weight[j];
        }
        ni.Normalize();
        if(nk[i]*ni<0)
            ni *= -1;
        sample_->vert[i].N() = ni;
    }

}

void Particle::runparticleedge()
{
    srand((int)time(0));

    int num; //取四分之邻域点的法向可以借用
    vector<vcg::Point3f> points;
    CMesh::VertexIterator iv;
    for(iv=sample_->vert.begin();iv!=sample_->vert.end();iv++)
    {
        vcg::Point3f orip = origin_->vert[iv->m_index].P();
        vcg::Point3f pi = iv->P();
        vcg::Point3f ni = iv->N();
        if(iv->feature)//对特征点进行粒子滤波
        {
            int k = iv->neighbors.size();
            num = k / 4;
            vector<pair<vcg::Point3f,double> > particles(particlenum*num*2);
            for(int t=0;t<num;t++)
            {
                vcg::Point3f nit;
                if(t==0)
                    nit = ni;
                else
                    nit = sample_->vert[iv->neighbors[t-1]].N();
                vcg::Point3f compose_nit = (ni + nit).Normalize();

                for(int i=0;i<particlenum;i++)//随机生成particlenum个粒子
                {
                    vcg::Point3f pm1 = pi + nit * myrandom(radius/6.0);
                    vcg::Point3f pm2 = pi + compose_nit * myrandom(radius/6.0);
                    particles[particlenum*t*2+i*2].first = pm1;
                    particles[particlenum*t*2+i*2+1].first = pm2;
                    double cov1 = 0,cov2 = 0,weight;

                    for(int j=0;j<k;j++)//计算与邻域点的方差
                    {
                        vcg::Point3f pj = sample_->vert[iv->neighbors[j]].P();
                        vcg::Point3f nj = sample_->vert[iv->neighbors[j]].N();
                        //weight = exp(-pow((pj-pi).Norm()/radius,2));
                        cov1 += (pm1 - pj).Norm() * fabs((pm1 - pj)*nj)/* * weight*/;
                        cov2 += (pm2 - pj).Norm() * fabs((pm2 - pj)*nj)/* * weight*/;
                    }
                    particles[particlenum*t*2+i*2].second = beta * cov1 + (orip - pm1).Norm();
                    particles[particlenum*t*2+i*2+1].second = beta * cov2 + (orip - pm2).Norm();
                }
            }
            sort(particles.begin(),particles.end(),cmp);//选择方差最小的粒子点作为先验
            points.push_back(particles[0].first);
        }
        else
            points.push_back(pi);
    }
    int ii=0;
    for(iv=sample_->vert.begin();iv!=sample_->vert.end();iv++,ii++)
    {
        iv->P() = points[ii];
    }
}

void Particle::runRANSAC()
{
    //计算球领域
//    GlobalFun::computeAnnNeigbhors(sample_->vert, sample_->vert, 6, false, "normal plus");
//    double radius = averageKNNDistance();radius *=2.5;
//    double radius2 = radius * radius;
//    for(CMesh::VertexIterator dest = sample_->vert.begin(); dest != sample_->vert.end(); dest++)
//        dest->neighbors.clear();
//    for(CMesh::VertexIterator dest = sample_->vert.begin(); dest != sample_->vert.end(); dest++)
//    {
//        Point3f p = dest->P();
//        for(CMesh::VertexIterator origin = dest + 1; origin != sample_->vert.end(); origin++)
//        {
//            Point3f q = origin->P();
//            Point3f diff = p - q;
//            double dist2 = diff.SquaredNorm();
//            if (dist2 < radius2)
//            {
//                dest->neighbors.push_back(origin->m_index);
//                origin->neighbors.push_back(dest->m_index);
//            }
//        }
//    }

    CMesh::VertexIterator iv;
    //计算带宽
    double beta = 1.5;
    double sum=0;int count=0;
    for(iv=sample_->vert.begin();iv!=sample_->vert.end();iv++){
        if(!iv->feature){
            sum += iv->errs;
            count++;
        }
    }
    double inner_threshold=beta*sum/count;

    int k = para->getInt("RANSAC iteration");

    for(iv=sample_->vert.begin();iv!=sample_->vert.end();iv++)
    {
        //if(!iv->feature) continue;
        int iterations = 0;
        vcg::Point3f ni = iv->N();
        vcg::Point3f pi = iv->P();
        vcg::Point3f this_normal;
        vcg::Point3f diff;
        double rj;
        double rou,omiga_d,omiga_n;
        vcg::Point3f best_normal = ni;
        vector<int> best_consensus_set;
        int this_num;
        int best_num = 0;
        double this_error;
        double best_error = -1;

        int num = iv->neighbors.size();

        //随机抽样
        while ( iterations < k )
        {
            iterations++;
            vector<int> randperm(num, 0);
            for (int j = 0; j < num; j++)
                randperm[j] = j;
            random_shuffle(randperm.begin(), randperm.begin() + num); //洗牌器,将领域点的序列顺序打乱

            vector<vcg::Point3f> points;//随机选三个点形成平面
            vcg::Point3f center;
            int index_;
            for(int i=0;i<3;i++)
            {
                index_ = iv->neighbors[randperm[i]];
                points.push_back(sample_->vert[index_].P());
            }

            double judge = ((points[0]-points[1]).Normalize() - (points[0]-points[2]).Normalize()).Norm();
            if(judge == 0.0) continue;//判断三点是否共线

            //计算这三个点所形成的平面
            vcg::Point3f normal;
            MatrixType covariance_matrix;
            vcg::Point3f eigenvalues;
            MatrixType eigenvectors;
            int required_rotations;
            covariance_matrix.Covariance(points,center);//计算协方差
            vcg::Jacobi<MatrixType, vcg::Point3f>(covariance_matrix, eigenvalues, eigenvectors, required_rotations);
            vcg::SortEigenvaluesAndEigenvectors<MatrixType, vcg::Point3f>(eigenvalues, eigenvectors);

            normal = eigenvectors.GetColumn(2);
            normal.Normalize();
            if(normal*ni<0)
                normal *= -1;
            //if(normal*ni<0.707) continue;//如果与原先法向夹角大于45度,则放弃该平面;

            //点自身要满足在平面带宽内
            if((pi-center).Norm() * fabs((pi-center)*normal)>inner_threshold) continue;


            double error;
            vcg::Point3f pj,nj;
            for(int i=3;i<num;i++) //对每个邻域集中不属于随机三点的邻域点,进行带宽判断
            {
                index_ = iv->neighbors[randperm[i]];
                pj = sample_->vert[index_].P();
                error = (pj-center).Norm() * fabs((pj-center)*normal);
                if(error>inner_threshold) continue;
                best_consensus_set.push_back(index_);
            }

            this_num = best_consensus_set.size();
            int num_threshold = num / 3;
            if(this_num < num_threshold) continue;//满足不了数量要求,继续

            //if(this_num < best_num) continue;//选数量最多中误差最小的

            this_normal = normal;
            this_error = 0;
            for(int i=3;i<num;i++)
            {
                index_ = iv->neighbors[randperm[i]];
                pj = sample_->vert[index_].P();
                nj = sample_->vert[index_].N();
                diff = pj-center;
                rj = diff.Norm() * fabs(diff*this_normal);//邻域点到平面的距离
                rou = exp(-pow(rj,2)/pow(inner_threshold,2));
                omiga_d = exp(-(pi-pj).Norm()/pow(radius,2));
                omiga_n = exp(-(nj-ni).Norm()/pow(0.707,2));
                this_error += rj * rou * omiga_d * omiga_n;
            }
            //this_error /= best_consensus_set.size();
            if ( this_error > best_error )
            {
                best_error = this_error;
                best_normal = this_normal;
            }
        }//end while
        iv->N() = best_normal;
    }//end for

}

void Particle::run_normal()
{
    double inner_threshold;
    double omiga_r,omiga_d,omiga_n;

    vcg::Point3f ni , pi;

    vector<vcg::Point3f> nk;
    nk.resize(sample_->vn);

    vector<vcg::Point3f> initpoint;
    initpoint.resize(sample_->vn);
    for(int i=0;i<sample_->vn;i++)
        initpoint[i] = sample_->vert[i].P();//记录点信息

    vector<vector<int> > neighbors;
    neighbors.resize(sample_->vn);
    for(int i=0;i<sample_->vn;i++)
        neighbors[i] = sample_->vert[i].neighbors;

    int k = para->getInt("particle number");
    int iterations = 0;
    double feat_threshold = global_paraMgr.particle.getDouble("feat");
    while(iterations < k)
    {
        iterations++;

        for(int i=0;i<sample_->vn;i++)
            nk[i] = sample_->vert[i].N();//记录法向信息

        for(int i=0;i<sample_->vn;i++)
        {
            //if(!sample_->vert[i].feature) continue;
            ni = nk[i];
            pi = initpoint[i];
            int num = neighbors[i].size();
            vcg::Point3f pm;
            pm.SetZero();
            vector<double> weight;
            double sum_weight = 0;

            double feat = sample_->vert[i].sigma_plus>feat_threshold? 1: sample_->vert[i].sigma_plus/feat_threshold;


            double max_n = 0;
            double max_r = 0;
            vcg::Point3f center = sample_->vert[i].center;
            for(int j=0; j<num; j++)
            {
                vcg::Point3f nj = nk[neighbors[i][j]];
                vcg::Point3f pj = initpoint[neighbors[i][j]];
                max_n = max_n > (ni-nj).Norm()? max_n: (ni-nj).Norm();
                max_r = max_r > fabs((pj-center)*ni)? max_r:  fabs((pj-center)*ni);
            }

            double sigma_n = max_n / (3*feat);
            double sigma_r = max_r / (3*feat);
            //cout<<sigma_n<<" "<<sigma_r<<" "<<feat<<endl;
            for(int j=0; j<num; j++)
            {
                vcg::Point3f nj = nk[neighbors[i][j]];
                vcg::Point3f pj = initpoint[neighbors[i][j]];
                omiga_d = exp(-pow((pj-pi).Norm()/(0.5*radius),2));
                omiga_n = exp(-pow((nj-ni).Norm()/sigma_n,2));
                omiga_r = exp(-pow((pj-center)*ni/sigma_r,2));
                //cout<<omiga_d<<endl;
                weight.push_back(omiga_r * omiga_d * omiga_n);
                sum_weight += weight[j];
                pm += pj * weight[j];
            }

            double max_weight = *max_element(weight.begin(),weight.end());
            for(int j=0; j<num; j++)
                weight[j] /= max_weight;

            sample_->vert[i].weight = weight;
            pm /= sum_weight;
            sample_->vert[i].center = pm;

            vcg::Point3f diff;
            MatrixType covariance_matrix;
            MatrixType cov;
            covariance_matrix.SetZero();
            for(int j=0;j<num;j++)
            {
                diff = initpoint[neighbors[i][j]] - pm;
                for (int m = 0; m < 3; m++)
                    for (int n = 0; n < 3; n++){
                        cov[m][n] = diff[m] * diff[n];
                    }
                covariance_matrix += cov * weight[j];
            }
            vcg::Point3f eigenvalues;
            MatrixType eigenvectors;
            int required_rotations;

            vcg::Jacobi<MatrixType, vcg::Point3f>(covariance_matrix, eigenvalues, eigenvectors, required_rotations);
            vcg::SortEigenvaluesAndEigenvectors<MatrixType, vcg::Point3f>(eigenvalues, eigenvectors);

            ni = eigenvectors.GetColumn(2);
            ni.Normalize();
            if(nk[i]*ni<0)
                ni *= -1;
            sample_->vert[i].N() = ni;
        }
    }

//        CMesh::VertexIterator iv;


//        vector<vcg::Point3f> nk;
//        nk.resize(sample_->vn);

//        vector<vcg::Point3f> initpoint;
//        initpoint.resize(sample_->vn);
//        for(int i=0;i<sample_->vn;i++)
//        {
//            nk[i] = sample_->vert[i].N();//记录法向
//            initpoint[i] = sample_->vert[i].P();//记录点信息
//        }

//        vector<vector<int> > neighbors;
//        neighbors.resize(sample_->vn);
//        for(int i=0;i<sample_->vn;i++)
//            neighbors[i] = sample_->vert[i].neighbors;

//        int k = para->getInt("particle number");
//        vcg::Point3f center;

//        for(int i=0;i<sample_->vn;i++)
//        {
//            vcg::Point3f pi = initpoint[i];
//            vcg::Point3f ni = nk[i];
//            center = sample_->vert[i].center;
//            int num = neighbors[i].size(); //邻域数量
//            vector<double> rj;
//            for(int j=0; j<num; j++)
//            {
//                vcg::Point3f nj = nk[neighbors[i][j]];
//                vcg::Point3f pj = initpoint[neighbors[i][j]];
//                double rj_dis = (pi - pj) * ni ;//每个邻域点到待计算点所在点平面的距离
//                rj.push_back(pow(rj_dis,2));
//            }
//            double mu = *max_element(rj.begin(),rj.end());
//            vector<double> weight;
//            int iterations = 0;
//            while(iterations<k)
//            {
//                iterations++;
//                rj.clear();
//                vcg::Point3f pm; //新平面的中点
//                double sum_weight = 0;
//                weight.clear();
//                for(int j=0; j<num; j++)
//                {
//                    vcg::Point3f nj = nk[neighbors[i][j]];
//                    vcg::Point3f pj = initpoint[neighbors[i][j]];
//                    double rj_dis = (pi - pj) * ni ;//每个邻域点到待计算点所在点平面的距离
//                    rj.push_back(pow(rj_dis,2));
//                    double weight_j = mu/(mu+rj[j]);
//                    weight.push_back(pow(weight_j,2));
//                    sum_weight += weight[j];
//                    pm += pj * weight[j];
//                }
//                pm /= sum_weight;
//                center = pm;
//                double max_weight = *max_element(weight.begin(),weight.end());
//                double min_weight = *min_element(weight.begin(),weight.end());
//                for(int j=0; j<num; j++)
//                    weight[j] = (weight[j]-min_weight)/(max_weight - min_weight);

//                vcg::Point3f diff;
//                MatrixType covariance_matrix;
//                MatrixType cov;
//                covariance_matrix.SetZero();
//                for(int j=0;j<num;j++)
//                {
//                    diff = initpoint[neighbors[i][j]] - pm;
//                    for (int m = 0; m < 3; m++)
//                        for (int n = 0; n < 3; n++){
//                            cov[m][n] = diff[m] * diff[n];
//                        }
//                    covariance_matrix += cov * weight[j];
//                }
//                vcg::Point3f eigenvalues;
//                MatrixType eigenvectors;
//                int required_rotations;

//                vcg::Jacobi<MatrixType, vcg::Point3f>(covariance_matrix, eigenvalues, eigenvectors, required_rotations);
//                vcg::SortEigenvaluesAndEigenvectors<MatrixType, vcg::Point3f>(eigenvalues, eigenvectors);

//                ni = eigenvectors.GetColumn(2);
//                ni.Normalize();
//                if(nk[i]*ni<0)
//                    ni *= -1;
//                mu /= 1.01;
//            }
//            sample_->vert[i].N() = ni;
//            sample_->vert[i].weight = weight;
//        }
}

double Particle::averageKNNDistance()
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

double Particle::calculate_bandwidth()
{
    CMesh::VertexIterator iv;
    double sum = 0;
    for(iv=sample_->vert.begin();iv!=sample_->vert.end();iv++)
    {
        if(!iv->feature){
            sum = sum > iv->errs? sum : iv->errs;//errs为pcvmn到所在平面的距离
        }
    }

    return sum;
}

double Particle::getPointtoFace(CMesh *sample,CMesh *original)
{
    double Projection_err = 0.0;
    vcg::tri::UpdateTopology<CMesh>::VertexFace(*original);
    for(CMesh::VertexIterator iter = original->vert.begin(),iter1 = sample->vert.begin(); iter != original->vert.end();iter++,iter1++)
    {
        CMesh::VertexPointer vp = &(*iter);
        CVertex *v = vp;
        Point3f P = iter1->P();
        vcg::face::VFIterator<CFace>vfi(v);
        double err = 0;
        for(;!vfi.End();++vfi)
        {
            Point3f quality(0,0,0);
            CFace* f = vfi.F();
            CVertex  *p0 = f->V(0);
            CVertex  *p1 = f->V(1);
            CVertex  *p2 = f->V(2);
            Point3f P0 = p0->P(); Point3f P1 = p1->P(); Point3f P2 = p2->P();
            quality = (P0 + P1 +P2)/3;
            Point3f n1 = (P0 - P1) ^ (P0 - P2).Normalize();
            Point3f edge1 = P1 - P0;
            Point3f edge2 = P1 - P2;
            double S = 0.5 * (edge1 ^ edge2).Norm();
            err =abs((P - quality).dot(n1)) * S;
            Projection_err += err;
        }
    }
     double Ak = 0.0;
    for (CMesh::FaceIterator f_it1 = original->face.begin();f_it1 != original->face.end();
         f_it1++)
    {
        CVertex  *p0 = f_it1->V(0);
        CVertex  *p1 = f_it1->V(1);
        CVertex  *p2 = f_it1->V(2);
        Point3f P0 = p0->P(); Point3f P1 = p1->P(); Point3f P2 = p2->P();
        Point3f edge1 = P1 - P0;
        Point3f edge2 = P1 - P2;
        Ak += (0.5 * (edge1 ^ edge2).Norm());
    }
    return Projection_err / (Ak*3);

}
