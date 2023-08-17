#ifndef PCVMN_H
#define PCVMN_H
#include "ParameterMgr.h"
#include <iostream>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Core>

#include "Algorithm/normal_extrapolation.h"

using namespace std;
using namespace Eigen;
class PCVmn
{
    public:
    PCVmn(){}
    virtual ~PCVmn(){}
    typedef typename vector<CVertex>::iterator VertexIterator;

    static void compute_normals(const VertexIterator &begin, const VertexIterator &end, int knn);
    static Point3f compute_normal_NWR_EACH4_ExraFea(MatrixXd points, MatrixXd local_W, int ran_num, double inner_threshold, double compart, MatrixXd local_density, int T, double sum_des, vector<Point3f> &plane_center);
    static void compute_normal_NWR_EACH5_ExraFea(MatrixXd curPoints, MatrixXd curlocal_W, int ran_num, double inner_threshold, Point3f &normal_one, int &noncompute, MatrixXd &dis_vect, vector<Point3f> &plane_center);
    //static void compute_normal_NWR_EACH6_ExraFea(MatrixXd curPoints, MatrixXd curlocal_W, int ran_num, double inner_threshold, Point3f &normal_one, int &noncompute,Point3f &center,MatrixXd &dis_vect);
    private:

};

#endif // PCVMN_H
