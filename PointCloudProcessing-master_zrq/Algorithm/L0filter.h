#pragma once
#include <iostream>
#include "GlobalFunction.h"
#include "normal_extrapolation.h"
#include "DataMgr.h"
#include "ParameterMgr.h"
#include "PointCloudAlgorithm.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <numeric>

using namespace std;
using namespace vcg;

class L0filter : public PointCloudAlgorithm
{
  public:
    L0filter(RichParameterSet* _para);
    ~L0filter(void);

    RichParameterSet* getParameterSet()
    {
        return para;
    }
    void setParameterSet(RichParameterSet* _para)
    {
        para = _para;
    }

    void initial();
    void run();
    void setInput(DataMgr* pData);
    void clear(){cout<<"clear"<<endl;}

    int total_iteration;
    int normKnn, norm_iteration;
    int vertexKnn, vertex_iteration;
    int edgeKnn, edge_iteration;
    double feature_scale_, edge_stepsize_;
    double norm_sparsity, norm_auxiliary, norm_increment ;
    double vert_sparsity, vert_auxiliary, vert_increment ;
    double vertex_stepsize_;
    int ptnum_;
    vector<Eigen::Vector3d> pnts_points_, last_result_points_, pnts_normals_,input_pnts_normals_,input_pnts_points_;
    vector<vector<int> > vertex_neighbors_, include_i_knn_idx, neighbors_v_, neighbors_includeI_;

    void NormalFiltering(std::vector<Eigen::Vector3d>& pnts_normals);
    void NormalFilteringplus();
    void getnewneighbor();
    void NormalFiltering_plus(int iteration);
    double averageKNNDistance();
    void find_include_i_knn_idx(const vector<vector<int>>& vertex_neighbors, vector<vector<int>>& include_i_knn_idx, int Knn);
    void solve_theta_problem(const vector<Eigen::Vector3d>& pnts_normals, std::vector<Eigen::Vector3d> &thetas);
    void solve_N_problem(const vector<Eigen::Vector3d>& thetas, vector<Eigen::Vector3d>& pnts_normals);

    void VertexUpdating(const double stepsize,const vector<Eigen::Vector3d>& result_normals,vector<Eigen::Vector3d>& result_points);
    void vertexupdating_L0(const vector<Eigen::Vector3d>& result_normals,vector<Eigen::Vector3d>& result_points);
    void L0_vert_theta(const vector<Eigen::Vector3d>& result_normals,const vector<Eigen::Vector3d>& lastpoints,const vector<double>& alpha,vector<double>& theta);
    void L0_vert_alpha(const vector<Eigen::Vector3d>& result_normals,	const vector<double>& theta,vector<Eigen::Vector3d>& lastpoints,vector<double>& alpha);

    void vertexupdating_L0_plus(const vector<Eigen::Vector3d>& result_normals,vector<Eigen::Vector3d>& result_points);
    void L0_vert_theta_plus(const vector<Eigen::Vector3d>& result_normals,const vector<Eigen::Vector3d>& lastpoints,const vector<double>& alpha,vector<double>& theta);
    void L0_vert_alpha_plus(const vector<Eigen::Vector3d>& result_normals,	const vector<double>& theta,vector<Eigen::Vector3d>& lastpoints,vector<double>& alpha);

    void recoveryEdge(double feature_scale, vector<Eigen::Vector3d>& result_normals, vector<Eigen::Vector3d>& result_points);

private:
    RichParameterSet* para;

public:
    CMesh* sample_;
    //DataMgr* errpData;
};
