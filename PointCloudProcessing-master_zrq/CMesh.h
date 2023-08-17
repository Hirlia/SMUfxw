
#pragma once
#include <vcg/simplex/edge/base.h>
#include <vcg/simplex/face/base.h>
#include <vcg/simplex/face/component_ocf.h>
#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/vertex/component_ocf.h>

#include <vcg/complex/trimesh/allocate.h>
#include <vcg/complex/trimesh/base.h>
#include <vcg/complex/used_types.h>

#include <vcg/simplex/face/topology.h>

#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/update/color.h>
#include <vcg/complex/trimesh/update/flag.h>
#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/update/position.h>
#include <vcg/complex/trimesh/update/quality.h>
#include <vcg/complex/trimesh/update/selection.h>
#include <vcg/complex/trimesh/update/topology.h>

#include <vcg/complex/trimesh/create/platonic.h>


#include <vcg/simplex/face/pos.h>

#include <vcg/space/point3.h>

#include <cstdlib>  //for rand()
#include <ctime>    //for time()

#include <vector>
using std::vector;
using namespace vcg;
//��vcg�ⶨ����ά�����ṹ

class CVertex;
class CFace;
class CEdge;

class CUsedTypes : public vcg::UsedTypes<vcg::Use<CVertex>::AsVertexType,           vcg::Use<CFace>::AsFaceType,vcg::Use<CEdge>::AsEdgeType>
{
};

class CVertex : public vcg::Vertex<CUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::Color4b,vcg::vertex::VFAdj>
{
  public:
    vector<int> neighbors;
    vector<int> big_neighbors;
    vector<int> original_neighbors;


    bool is_neighbor_sample;  //确定是不是所选的显示点

    bool bIsOriginal;
    int m_index;
    bool feature ;
    double sigma;
    double sigma_plus;
    double sigma_c;     //点密度

    Point3f center;
    Point3d normal;
    double errs;        //pcvmn到所在平面的距离
    vector<double> weight;  //邻域点的法向权重
    double weight_value;    //当该点为所选点的邻域点时,它的权重值
    Point3f normal_pcv;

    bool is_fixed_sample;  // feature points (blue color)
    bool is_skel_ignore;

    /* for skeletonization */
    double eigen_confidence;
    Point3f eigen_vector0;  // Associate with the biggest eigen value
    Point3f eigen_vector1;  // Also use for remember last better virtual point
    Point3f eigen_vector2; //The smallest eigen value : should be PCA normal N()
    float eigen_value0;
    float eigen_value1;
    float eigen_value2;

    Point3f norm_eigen_vector0;
    Point3f norm_eigen_vector1;
    Point3f norm_eigen_vector2;
    float norm_eigen_value0;
    float norm_eigen_value1;
    float norm_eigen_value2;


    bool is_skel_virtual;  // in our papaer, we said bridge point instead of virtual point
    bool is_skel_branch;
    bool is_fixed_original;

    double skel_radius;  // remember radius for branches

  public:
    operator Point3f &()
    {
        return P();
    }

    operator const Point3f &() const
    {
        return cP();
    }

    float &operator[](unsigned int i)
    {
        return P()[i];
    }

    CVertex()
        : m_index(0)
        , bIsOriginal(false)
        , is_fixed_sample(false)
        , eigen_confidence(0)
        , is_skel_branch(false)
        , is_skel_ignore(false)
        , is_skel_virtual(false)
        , is_fixed_original(false)
        , is_neighbor_sample(false)
        , eigen_vector0(Point3f(1, 0, 0))
        , eigen_vector1(Point3f(0, 1, 0))
        , skel_radius(-1.0)
    {
        N() = Point3f(0, 0, 0);
    }

    /* for skeletonization */
    void remove()  // important, some time we don't want to earse points, just remove them
    {
        neighbors.clear();
        original_neighbors.clear();
        is_skel_ignore = true;
        P() = Point3f(88888888888.8, 88888888888.8, 88888888888.8);
    }

    bool isSample_Moving()
    {
        return (!is_skel_ignore && !is_fixed_sample && !is_skel_branch);
    }

    bool isSample_JustMoving()
    {
        return (!is_skel_ignore && !is_fixed_sample && !is_skel_virtual && !is_skel_branch);
    }

    bool isSample_MovingAndVirtual()
    {
        return (!is_skel_ignore && !is_fixed_sample && is_skel_virtual && !is_skel_branch);
    }

    bool isSample_JustFixed()
    {
        return (!is_skel_ignore && is_fixed_sample && !is_skel_virtual && !is_skel_branch);
    }

    bool isSample_FixedAndBranched()
    {
        return (!is_skel_ignore && is_fixed_sample && !is_skel_virtual && is_skel_branch);
    }

    void setSample_JustMoving()
    {
        is_fixed_sample = false;
        is_skel_virtual = false;
        is_skel_branch = false;
    }

    void setSample_MovingAndVirtual()
    {
        is_fixed_sample = false;
        is_skel_virtual = true;
        is_skel_branch = false;
    }

    void setSample_JustFixed()
    {
        is_fixed_sample = true;
        is_skel_virtual = false;
        is_skel_branch = false;
    }

    void setSample_FixedAndBranched()
    {
        is_fixed_sample = true;
        is_skel_virtual = false;
        is_skel_branch = true;
    }

    void recompute_m_render()
    {
        srand(time(NULL));
        int x = rand() % 1000;
        int y = rand() % 1000;
        int z = rand() % 1000;

        Point3f normal = N();
        normal.Normalize();

        Point3f helper(x / 1000.0, y / 1000.0, z / 1000.0);
        Point3f new_m3_to_m5 = normal ^ helper;
        new_m3_to_m5.Normalize();
        Point3f new_m6_to_m8 = normal ^ new_m3_to_m5;

        eigen_vector0 = new_m3_to_m5;
        eigen_vector1 = new_m6_to_m8;

        // eigen_vector0 = N() ^ Point3f(0.12345, 0.12346, 0.12347);
        // eigen_vector0.Normalize();
        // eigen_vector1 = eigen_vector0 ^ N();
        // eigen_vector1.Normalize();
    }
};

class CFace : public vcg::Face<CUsedTypes, vcg::face::VertexRef, vcg::face::VFAdj>
{
};
class CEdge    : public Edge<CUsedTypes,edge::VertexRef, edge::EVAdj, edge::EEAdj, edge::BitFlags >{};
class CMesh : public vcg::tri::TriMesh<std::vector<CVertex>, std::vector<CFace>, std::vector<CEdge> >
{
public:
    CMesh &operator=(const CMesh &a)
    {
        if(this->vn == a.vn){
            for(int i = 0; i < a.vn; i++){
                this->vert[i].P() = a.vert[i].P();
                this->vert[i].N() = a.vert[i].cN();
            }
        }
        return *this;
    }
};
