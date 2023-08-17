#include "pcvmn.h"
using namespace std;
using namespace Eigen;
using namespace vcg;

void PCVmn::compute_normal_NWR_EACH5_ExraFea(MatrixXd curPoints,MatrixXd curlocal_W,int ran_num,double inner_threshold,Point3f &normal_one,int &noncompute,MatrixXd &dis_vect,vector<Point3f>& plane_center)
{
    int num=curPoints.rows();
    double max_sort=0;
    MatrixXd x,y;
    x.setOnes(1,4);
    y.setOnes(num,1);
    dis_vect.setZero(num,1);
    MatrixXd poins_3(4,3);
    for(int i=0;i<ran_num;i++){
        /*STEP1:求得随机3个领域点及该点在内所构成的平面法向*/
        vector<int> randperm(num, 0);
        //srand(time(NULL));
        for (int j = 0; j < num; j++)
            randperm[j] = j;
        random_shuffle(randperm.begin()+1, randperm.begin() + num);
        for(int j=0;j<4;j++){
            for(int k=0;k<3;k++)
                poins_3(j,k)=curPoints(randperm[j],k);
        }
        MatrixXd mp=x*poins_3/4;
        MatrixXd points_center=curPoints-(y*mp);
        MatrixXd tmp(4,3);
        for(int j=0;j<4;j++)
            for(int k=0;k<3;k++)
                tmp(j,k)=points_center(randperm[j],k);
        MatrixXd C=tmp.transpose()*tmp/poins_3.rows();
        SelfAdjointEigenSolver<MatrixXd> eig(C);
        MatrixXd fix_normal = eig.eigenvectors().col(0);
        /*STEP2:计算得分*/
        MatrixXd temp_dis=points_center*fix_normal;
        temp_dis=temp_dis.array().pow(2);
        if(temp_dis(0,0)>pow(inner_threshold,2))continue;
        MatrixXd dis=-temp_dis/pow(inner_threshold,2);
        dis=dis.array().exp();
        MatrixXd dis_t=dis.transpose();
        MatrixXd cur_sort=dis_t*curlocal_W*dis;
        if(cur_sort(0,0)>max_sort){
            max_sort=cur_sort(0,0);
            normal_one=Point3f{fix_normal(0,0),fix_normal(1,0),fix_normal(2,0)};
            dis_vect=temp_dis;
            plane_center.push_back(Point3f{mp(0,0),mp(0,1),mp(0,2)});
        }
    }
    if(normal_one==vcg::Point3f(0,0,0)){
        normal_one=vcg::Point3f(1,0,0);
        noncompute=1;
    }
}

Point3f PCVmn::compute_normal_NWR_EACH4_ExraFea(MatrixXd points, MatrixXd local_W, int ran_num, double inner_threshold, double compart, MatrixXd local_density, int T,double sum_des,vector<Point3f>& plane_center)
{
    int num=points.rows();
    bool NonInner[num]={0};
    int computedNum=0;
    int noncompute=0;
    vcg::Point3f normal_one;normal_one.SetZero();
    bool flag=false;
    int num_curPoints=count(NonInner,NonInner+num,0),count(0);
    MatrixXd curPoints(num_curPoints,3);
    for(int i=0;i<num;i++){
        if(!NonInner[i]){
            for(int j=0;j<3;j++)
                curPoints(count,j)=points(i,j);
            count++;
        }
    }
    int rows(0),cols(0);
    MatrixXd curLocal_W(num_curPoints,num_curPoints);
    for(int i=0;i<num;i++){
        if(!NonInner[i]){
            for(int j=0;j<num;j++){
                if(!NonInner[j])curLocal_W(rows,cols++)=local_W(i,j);
            }
            rows++;
            cols=0;
        }
    }
    MatrixXd dis_vect;
    //单法向计算函数
    //feature_normal.clear();
    compute_normal_NWR_EACH5_ExraFea(curPoints,curLocal_W,ran_num,inner_threshold,normal_one,noncompute,dis_vect,plane_center);
    return normal_one;
}


void PCVmn::compute_normals(const VertexIterator &begin, const VertexIterator &end,int knn)
{
    int ran_num_min=150,ran_num_max=250;
    double beta=3.0;
    int alpha =6,T=2;
    double sum=0;int count=0;
    for(VertexIterator iter = begin; iter != end; iter++){
        if(!iter->feature){
            sum+=iter->errs;
            count++;
        }
    }
    double inner_threshold=beta*sum/count;
    for(VertexIterator iter = begin; iter != end; iter++){
        /*STEP1:判断是否为特征点*/
        if(!iter->feature){
            continue;
        }
        /*STEP2:计算带宽（由平滑区的方差决定）*/
        /*STEP3:计算法向权重（法向之间夹角越小，权重越大）*/
        vector<CVertex> points;
        MatrixXd p_points(knn,3);
        for(int i=0;i<knn;i++){
            points.push_back(*(begin+iter->neighbors[i]));
            for(int j=0;j<3;j++)
                p_points(i,j)=points[i].P()[j];
        }
        MatrixXd temp_nor(knn,3);
        vector<double> local_density_vector;
        MatrixXd local_density(knn,1);
        for(int i=0;i<knn;i++){
            local_density_vector.push_back(points[i].sigma_c);
            local_density(i,0)=local_density_vector[i];
            for(int j=0;j<3;j++)
                temp_nor(i,j)=points[i].normal_pcv[j];//knn×3
        }
        double sum_des=local_density.sum();
        MatrixXd temp_t=temp_nor.transpose();//3×knn
        MatrixXd temp_W=temp_nor*temp_t;//knn×knn
        temp_W=temp_W.array().abs();//knn×knn
        MatrixXd local_W=(temp_W.array().pow(alpha))/pow(0.707,alpha);
        local_W=local_W.array().exp();//knn×knn
        /*STEP4:计算密度权重（平均间隙越大权重越大）*/
        sort(local_density_vector.begin(),local_density_vector.end());
        double sum_maxd=0,sum_mind=0;
        for(int i=0;i<15;i++)
            sum_mind+=local_density_vector[i];
        for(int i=knn-15;i<knn;i++){
            sum_maxd+=local_density_vector[i];
        }
        sum_mind/=15;
        sum_maxd/=15;
        double density_r=sum_maxd/sum_mind;
        local_density=local_density.array().pow(2);//knn×1
        MatrixXd local_density_t=local_density.transpose();//1×knn
        MatrixXd density_w=local_density_t*local_density;//1×1
        /*STEP5:更新权重*/
        local_W= density_w(0,0)*local_W;//knn×knn
        int ran_num=density_r>2?ran_num_max:ran_num_min;
        double compart=0.03;
        vector<Point3f> plane_center;       //点的平面中心点；
        Point3f normal_pcvmn;
        normal_pcvmn=compute_normal_NWR_EACH4_ExraFea(p_points,local_W,ran_num,inner_threshold,compart,local_density,T,sum_des,plane_center);
        iter->N() = iter->N()*normal_pcvmn>0? normal_pcvmn:normal_pcvmn*(-1);

    }
}
