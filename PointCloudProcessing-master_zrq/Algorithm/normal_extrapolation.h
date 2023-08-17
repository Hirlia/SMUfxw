/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef VCG_SPACE_NORMAL_EXTRAPOLATION_H
#define VCG_SPACE_NORMAL_EXTRAPOLATION_H

#include <vcg/math/disjoint_set.h>
#include <vcg/math/lin_algebra.h>
#include <vcg/math/linear.h>
#include <vcg/math/matrix33.h>
#include <vcg/space/box3.h>
#include <vcg/space/index/octree.h>
#include <vcg/space/point3.h>
#include <wrap/callback.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <CMesh.h>
#include "GlobalFunction.h"

#include <stdlib.h>
#include <algorithm>
#include <limits>
#include <queue>
#include <vector>

namespace vcg
{
/*!
        */
template <class VERTEX_CONTAINER>//template<class T>  和  template<typename T> 都可以用来定义函数模板和类模板//
class NormalExtrapolation
{
  public:
    typedef typename VERTEX_CONTAINER::value_type VertexType;//typedef:是用于定义类型用的,typename和class的意义完全一样。//
    typedef VertexType *VertexPointer;
    typedef typename VERTEX_CONTAINER::iterator VertexIterator;
    typedef typename VertexType::CoordType CoordType;
    typedef typename VertexType::NormalType NormalType;
    typedef typename VertexType::ScalarType ScalarType;
    typedef typename vcg::Box3<ScalarType> BoundingBoxType;
    typedef typename vcg::Matrix33<ScalarType> MatrixType;
    typedef Eigen::MatrixXd Matrix;
    typedef Eigen::SparseMatrix<double> SpMat;
    typedef Eigen::Triplet<double> Triplet;

    enum NormalOrientation//enum工具不仅能够创建符号常量，还能定义新的数据类型//
    {
        IsCorrect = 0,
        MustBeFlipped = 1
    };

  private:
    /*************************************************
                        *		Inner class definitions
                        **************************************************/
    // Dummy class: no object marker is needed
    class DummyObjectMarker
    {
    };

    // Object functor: compute the distance between a vertex and a point
    struct VertPointDistanceFunctor
    {//在c/c++中，为了解决一些频繁调用的小函数大量消耗栈空间（栈内存）的问题，特别的引入了inline修饰符，表示为内联函数，栈空间就是指放置程序的局部数据（也就是函数内数据）的内存空间。//
        inline bool operator()(const VertexType &v, const CoordType &p, ScalarType &d, CoordType &q) const
        {
            ScalarType distance = vcg::Distance(p, v.P());
            if (distance > d)
                return false;

            d = distance;
            q = v.P();
            return true;
        }
    };
    // Plane structure: identify a plain as a <center, normal> pair
    struct Plane
    {
        Plane()
        {
            center.SetZero();
            center_pcv.SetZero();
            normal.SetZero();
            normal_pcv.SetZero();
        };

        // Object functor: return the bounding-box enclosing a given plane
        inline void GetBBox(BoundingBoxType &bb)
        {
            bb.Set(center);
        };

        CoordType center;
        CoordType center_pcv;
        NormalType normal;
        NormalType normal_pcv;
        int index;
    };

    // Object functor: compute the distance between a point and the plane
    struct PlanePointDistanceFunctor
    {
        inline bool operator()(const Plane &plane, const CoordType &p, ScalarType &d, CoordType &q) const
        {
            ScalarType distance = vcg::Distance(p, plane.center);
            if (distance > d)
                return false;

            d = distance;
            q = plane.center;
            return true;
        }
    };

    // Represent an edge in the Riemannian graph
    struct RiemannianEdge
    {
        RiemannianEdge(Plane *p = NULL, ScalarType w = 1000000000. /* = std::numeric_limits<ScalarType>::max()*/)
        {
            plane = p;
            weight = w;
        }

        Plane *plane;
        ScalarType weight;
    };
    // Represent an edge in the MST tree
    struct MSTEdge
    {
        MSTEdge(Plane *p0 = NULL, Plane *p1 = NULL,
                ScalarType w = 1000000000. /*=std::numeric_limits<ScalarType>::max()*/)
        {
            u = p0;
            v = p1;
            weight = w;
        };
        inline bool operator<(const MSTEdge &e) const
        {
            return weight < e.weight;
        }

        Plane *u;
        Plane *v;
        ScalarType weight;
    };
    // Represent a node in the MST tree
    struct MSTNode
    {
        MSTNode(MSTNode *p = NULL)
        {
            parent = p;
        }

        MSTNode *parent;
        VertexPointer vertex;
        std::vector<MSTNode *> sons;
    };

    typedef std::vector<Plane> PlaneContainer;
    typedef typename PlaneContainer::iterator PlaneIterator;

  public:
    /*!
                */

    static void calculate_pointsigma(CMesh* sample_,const VertexIterator &begin, const VertexIterator &end,double radius,double feat)
    {
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
        //前面进行排序//
        CoordType center;
        int ii=0;
        for (VertexIterator iter = begin; iter != end; iter++,ii++)
        {
            std::vector<CoordType> nearest_points;
            center.SetZero();
            int k = iter->neighbors.size();
            for(int i=0;i<k;i++)
            {
                nearest_points.push_back(sample_->vert[iter->neighbors[i]].P());//获得k个最靠近的点//
            }
            center /= ScalarType(k);
            MatrixType covariance_matrix;
            covariance_matrix.Covariance(nearest_points,center);//计算协方差矩阵//

            CoordType eigenvalues;
            MatrixType eigenvectors;
            int required_rotations;
            vcg::Jacobi<MatrixType, CoordType>(covariance_matrix, eigenvalues, eigenvectors, required_rotations);//??//
            vcg::SortEigenvaluesAndEigenvectors<MatrixType, CoordType>(eigenvalues, eigenvectors);
            iter->sigma_plus = eigenvalues[2] / (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);//确认是否为特征点//
            iter->feature = false;
            iter->is_fixed_sample = false;
            if(iter->sigma_plus>feat){
                iter->is_fixed_sample = true;
                iter->feature = true;//更新特征属性
            }
        }
    }

    static void withoutshow_feature(const VertexIterator &begin, const VertexIterator &end)
    {
        for (VertexIterator iter = begin; iter != end; iter++)
            iter->is_fixed_sample = false;
    }


    static void calculate_normsigma(CMesh* sample_,const VertexIterator &begin, const VertexIterator &end,int k)
    {
        GlobalFun::computeAnnNeigbhors(sample_->vert, sample_->vert, k, false, "norm sigma");

        CoordType center;
        int ii=0;
        for (VertexIterator iter = begin; iter != end; iter++,ii++)
        {
            std::vector<CoordType> nearest_normals;
            std::vector<CoordType> nearest_points;
            center.SetZero();

            for(int i=0;i<k;i++)
            {
                nearest_normals.push_back(sample_->vert[iter->neighbors[i]].cN());
                nearest_points.push_back(sample_->vert[iter->neighbors[i]].P());
                center +=nearest_points[i];
            }
            center /= ScalarType(k);//中心点//
            MatrixType covariance_matrix1,covariance_matrix2;
            CoordType diff1,diff2;
            covariance_matrix1.SetZero();
            covariance_matrix2.SetZero();
            for (unsigned int n = 0; n < k; n++)
            {
                diff1 = nearest_normals[n] ;//法向差//
                diff2 = (nearest_points[n] - center);//点间坐标差//
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++){
                        covariance_matrix1[i][j] += diff1[i] * diff1[j];
                        covariance_matrix2[i][j] += diff2[i] * diff2[j];
                    }
            }

            CoordType eigenvalues;
            MatrixType eigenvectors;
            int required_rotations;
            vcg::Jacobi<MatrixType, CoordType>(covariance_matrix1, eigenvalues, eigenvectors, required_rotations);
            vcg::SortEigenvaluesAndEigenvectors<MatrixType, CoordType>(eigenvalues, eigenvectors);
            iter->norm_eigen_vector0=eigenvectors.GetColumn(2);
            iter->norm_eigen_vector1=eigenvectors.GetColumn(1);
            iter->norm_eigen_vector2=eigenvectors.GetColumn(0);
            iter->norm_eigen_value0=eigenvalues[2]/eigenvalues[0];
            iter->norm_eigen_value1=eigenvalues[1]/eigenvalues[0];
            iter->norm_eigen_value2=1.0;

            vcg::Jacobi<MatrixType, CoordType>(covariance_matrix2, eigenvalues, eigenvectors, required_rotations);
            vcg::SortEigenvaluesAndEigenvectors<MatrixType, CoordType>(eigenvalues, eigenvectors);
            iter->eigen_vector0=eigenvectors.GetColumn(2);
            iter->eigen_vector1=eigenvectors.GetColumn(1);
            iter->eigen_vector2=eigenvectors.GetColumn(0);
            iter->eigen_value0=eigenvalues[2]/eigenvalues[0];
            iter->eigen_value1=eigenvalues[1]/eigenvalues[0];
            iter->eigen_value2=1.0;

        }
    }

    static void calculate_feature(CMesh* sample_,const VertexIterator &begin, const VertexIterator &end,double sigma_n,double radius)
    {
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

        int ii=0;
        for (VertexIterator iter = begin; iter != end; iter++,ii++)
        {
            CoordType ni = iter->N();
            iter->feature = false;
            iter->is_fixed_sample = false;
            int k = iter->neighbors.size();
            for(int j=0;j<k;j++)
            {
                CoordType nj = sample_->vert[iter->neighbors[j]].N();
                if(fabs(ni*nj) <= sigma_n)
                {
                    iter->feature = true;
                    iter->is_fixed_sample = true;
                    break;
                }
            }
        }
    }

    static void feature_threshod_selection(const VertexIterator &begin, const VertexIterator &end)//阈值选择//
    {
        double feat;
        double sigma_threshold = 0.05;
        int length_taproot = 50;
        double lamda = 0.5;
        int taproot_hash[50] = {0};
        for(VertexIterator iter = begin; iter != end; iter++){
            taproot_hash[(int)floor(iter->sigma*length_taproot)]++;
        }
        int m_value=*std::max_element(taproot_hash,taproot_hash+50);                 //区域最大个数
        int m_id=std::max_element(taproot_hash,taproot_hash+50)-taproot_hash;        //第几个区域
        double taproot_standard_hash[50];                                       //标准化哈希表
        Matrix y(length_taproot,1),x(length_taproot,1);
        for(int i=0;i<length_taproot;i++){
            taproot_standard_hash[i]=(double)taproot_hash[i]/m_value;
            y(i,0) = taproot_standard_hash[i];
        }

        std::vector<Triplet> tripD;

        for(int i=0;i<length_taproot-1;i++)
        {
            tripD.push_back(Triplet(i, i, 1.0));
            tripD.push_back(Triplet(i, i+1, -1.0));
        }
        SpMat D(length_taproot-1, length_taproot);
        D.setFromTriplets(tripD.begin(), tripD.end());

        std::vector<Triplet> tripI(length_taproot);
        for(int i=0; i<length_taproot; i++) {
            tripI[i] = Triplet(i, i, 1.0);
        }
        SpMat I(length_taproot, length_taproot);
        I.setFromTriplets(tripI.begin(), tripI.end());

        SpMat  A = 0.5 * I +  lamda * D.transpose() * D;
        Matrix b = 0.5 * y;

        Eigen::ConjugateGradient<SpMat> solver;
        solver.compute(A);
        x = solver.solve(b);

        for(int i=m_id+1;i<length_taproot-1;i++){
            double diff_taproot=x(i-1,0)-x(i+1,0);
            if(diff_taproot<sigma_threshold&&x(i,0)<(x(m_id,0)-sigma_threshold)){
                feat=(double)i/length_taproot;
                break;
            }
        }
        for(VertexIterator iter = begin; iter != end; iter++){
            iter->feature = false;
            iter->is_fixed_sample = false;
            if(iter->sigma_plus>feat){
                iter->is_fixed_sample = true;
                iter->feature = true;//更新特征属性
            }
            //else
                //iter->N() = iter->eigen_vector0.Normalize();
        }
        std::cout<<"feature calculation finished"<<std::endl;
    }

    static void feature_threshod_selection_plus(const VertexIterator &begin, const VertexIterator &end,double feat)
    {
        for(VertexIterator iter = begin; iter != end; iter++)
        {
            iter->feature = false;
            iter->is_fixed_sample = false;
            if(iter->sigma>feat){
                iter->is_fixed_sample = true;
                iter->feature = true;//更新特征属性
            }
        }
        std::cout<<"feature calculation plus finished"<<std::endl;
    }

    static void ExtrapolateNormals(const VertexIterator &begin, const VertexIterator &end, const unsigned int k,
                                   const int root_index = -1, NormalOrientation orientation = IsCorrect,
                                   CallBackPos *callback = NULL)
    {
        BoundingBoxType dataset_bb;
        for (VertexIterator iter = begin; iter != end; iter++)
            dataset_bb.Add(iter->P());
        ScalarType max_distance = dataset_bb.Diag();

        // Step 1: identify the tangent planes used to locally approximate the surface
        int vertex_count = int(std::distance(begin, end));
        int step = int(vertex_count / 100) - 1;
        int progress = 0;
        int percentage;
        char message[128];
        sprintf(message, "Locating tangent planes...");
        std::vector<Plane> tangent_planes(vertex_count);
        vcg::Octree<VertexType, ScalarType> octree_for_planes;
        octree_for_planes.Set(begin, end);

        std::vector<VertexPointer> nearest_vertices;
        std::vector<CoordType> nearest_points;
        std::vector<ScalarType> distances;
        for (VertexIterator iter = begin; iter != end; iter++)
        {
            if (callback != NULL && (++progress % step) == 0 &&
                (percentage = int((progress * 100) / vertex_count)) < 100)
                (callback)(percentage, message);
            VertPointDistanceFunctor vpdf;
            DummyObjectMarker dom;
            octree_for_planes.GetKClosest(vpdf, dom, k, iter->P(), max_distance, nearest_vertices, distances,
                                          nearest_points);

            for(int j=0;j<k;j++)
                iter->neighbors.push_back(nearest_vertices[j]->m_index);

            double sum=0;
            for(int j=1;j<7;j++){  //此处 number=10
                sum+=distances[j];
            }

            iter->sigma_c=sum/6;
            // for each vertex *iter, compute the centroid as avarege of the k-nearest vertices of *iter
            Plane *plane = &tangent_planes[std::distance(begin, iter)];
            // then, identity the normal associated to the centroid
            MatrixType covariance_matrix,covariance_matrix_pcv;
            covariance_matrix.Covariance(nearest_points,plane->center);//计算协方差

            std::vector<CoordType> nearest_points_pcv(nearest_points.begin(),nearest_points.begin()+k/2);
            covariance_matrix_pcv.Covariance(nearest_points_pcv,plane->center_pcv);

            CoordType eigenvalues;
            MatrixType eigenvectors;
            int required_rotations;
            vcg::Jacobi<MatrixType, CoordType>(covariance_matrix, eigenvalues, eigenvectors, required_rotations);
            vcg::SortEigenvaluesAndEigenvectors<MatrixType, CoordType>(eigenvalues, eigenvectors);

            plane->normal = eigenvectors.GetColumn(2);
            plane->normal.Normalize();
            iter->N() = plane->normal;
            iter->center = plane->center;
            plane->index = int(std::distance(begin, iter));

            vcg::Jacobi<MatrixType, CoordType>(covariance_matrix_pcv, eigenvalues, eigenvectors, required_rotations);
            vcg::SortEigenvaluesAndEigenvectors<MatrixType, CoordType>(eigenvalues, eigenvectors);

            plane->normal_pcv = eigenvectors.GetColumn(2);
            plane->normal_pcv.Normalize();
            iter->normal_pcv = plane->normal_pcv;


            CoordType diff;
            diff = nearest_points[0] - plane->center;
            iter->errs = fabs(diff * plane->normal);
        }

        // Step 2: build the Riemannian graph, i.e. the graph where each point is connected to the k-nearest neigbours.
        dataset_bb.SetNull();
        PlaneIterator ePlane = tangent_planes.end();
        for (PlaneIterator iPlane = tangent_planes.begin(); iPlane != ePlane; iPlane++)
            dataset_bb.Add(iPlane->center);
        max_distance = dataset_bb.Diag();

        vcg::Octree<Plane, ScalarType> octree_for_plane;
        octree_for_plane.Set(tangent_planes.begin(), tangent_planes.end());
        std::vector<Plane *> nearest_planes(distances.size());
        std::vector<std::vector<RiemannianEdge> > riemannian_graph(
          vertex_count);  // it's probably that we are wasting the last position...
        progress = 0;
        sprintf(message, "Building Riemannian graph...");
        for (PlaneIterator iPlane = tangent_planes.begin(); iPlane != ePlane; iPlane++)
        {
            if (callback != NULL && (++progress % step) == 0 &&
                (percentage = int((progress * 100) / vertex_count)) < 100)
                (callback)(percentage, message);

            unsigned int kk = k;
            PlanePointDistanceFunctor ppdf;
            DummyObjectMarker dom;
            octree_for_plane.GetKClosest(ppdf, dom, kk, iPlane->center, max_distance, nearest_planes, distances,
                                         nearest_points, true, false);

            for (unsigned int n = 0; n < k; n++)
                if (iPlane->index < nearest_planes[n]->index)
                    riemannian_graph[iPlane->index].push_back(
                      RiemannianEdge(nearest_planes[n], 1.0f - fabs(iPlane->normal.dot(nearest_planes[n]->normal))));
        }

        // Step 3: compute the minimum spanning tree (MST) over the Riemannian graph (we use the Kruskal algorithm)
        std::vector<MSTEdge> E;
        typename std::vector<std::vector<RiemannianEdge> >::iterator iRiemannian = riemannian_graph.begin();
        typename std::vector<RiemannianEdge>::iterator iRiemannianEdge, eRiemannianEdge;
        for (int i = 0; i < vertex_count; i++, iRiemannian++)
            for (iRiemannianEdge = iRiemannian->begin(), eRiemannianEdge = iRiemannian->end();
                 iRiemannianEdge != eRiemannianEdge; iRiemannianEdge++)
                E.push_back(MSTEdge(&tangent_planes[i], iRiemannianEdge->plane, iRiemannianEdge->weight));

        std::sort(E.begin(), E.end());
        vcg::DisjointSet<Plane> planeset;

        for (typename std::vector<Plane>::iterator iPlane = tangent_planes.begin(); iPlane != ePlane; iPlane++)
            planeset.MakeSet(&*iPlane);

        typename std::vector<MSTEdge>::iterator iMSTEdge = E.begin();
        typename std::vector<MSTEdge>::iterator eMSTEdge = E.end();
        std::vector<MSTEdge> unoriented_tree;
        Plane *u, *v;
        for (; iMSTEdge != eMSTEdge; iMSTEdge++)
            if ((u = planeset.FindSet(iMSTEdge->u)) != (v = planeset.FindSet(iMSTEdge->v)))
                unoriented_tree.push_back(*iMSTEdge), planeset.Union(u, v);
        E.clear();

        // compute for each plane the list of sorting edges
        std::vector<std::vector<int> > incident_edges(vertex_count);
        iMSTEdge = unoriented_tree.begin();
        eMSTEdge = unoriented_tree.end();

        progress = 0;
        int mst_size = int(unoriented_tree.size());
        sprintf(message, "Building orieted graph...");
        for (; iMSTEdge != eMSTEdge; iMSTEdge++)
        {
            if (callback != NULL && (++progress % step) == 0 && (percentage = int((progress * 100) / mst_size)) < 100)
                (callback)(percentage, message);

            int u_index = int(iMSTEdge->u->index);
            int v_index = int(iMSTEdge->v->index);
            incident_edges[u_index].push_back(v_index), incident_edges[v_index].push_back(u_index);
        }

        // Traverse the incident_edges vector and build the MST
        VertexIterator iCurrentVertex, iSonVertex;
        std::vector<MSTNode> MST(vertex_count);

        typename std::vector<Plane>::iterator iFirstPlane = tangent_planes.begin();
        typename std::vector<Plane>::iterator iCurrentPlane, iSonPlane;

        MSTNode *mst_root;
        int r_index = (root_index != -1) ? root_index : rand() * vertex_count / RAND_MAX;
        mst_root = &MST[r_index];
        mst_root->parent = mst_root;  // the parent of the root is the root itself

        if (orientation == MustBeFlipped)
        {
            iCurrentVertex = begin;
            std::advance(iCurrentVertex, r_index);
            iCurrentVertex->N() = iCurrentVertex->N() * ScalarType(-1.0f);
        }

        {  // just to limit the scope of the variable border
            std::queue<int> border;
            border.push(r_index);
            int maxSize = 0;
            int queueSize = 0;
            progress = 0;
            sprintf(message, "Extracting the tree...");
            while ((queueSize = int(border.size())) > 0)
            {
                if (callback != NULL && ((++progress % step) == 0) &&
                    (percentage = int((maxSize - queueSize) * 100 / maxSize)) < 100)
                    (callback)(percentage, message);

                int current_node_index = border.front();
                border.pop();

                MSTNode *current_node = &MST[current_node_index];  // retrieve the pointer to the current MST node
                std::advance((iCurrentVertex = begin),
                             current_node_index);         // retrieve the pointer to the correspective vertex
                current_node->vertex = &*iCurrentVertex;  // and associate it to the MST node

                std::vector<int>::iterator iSon = incident_edges[current_node_index].begin();
                std::vector<int>::iterator eSon = incident_edges[current_node_index].end();
                for (; iSon != eSon; iSon++)
                {
                    MSTNode *son = &MST[*iSon];
                    if (son->parent == NULL)  // the node hasn't been visited
                    {
                        son->parent = current_node;  // Update the MST nodes
                        current_node->sons.push_back(son);
                        // std::advance((iSonVertex=begin), *iSon);//retrieve the pointer to the Vertex associated to
                        // son
                        border.push(*iSon);
                    }
                    maxSize = std::max<int>(maxSize, queueSize);
                }
            }
        }

        // and finally visit the MST tree in order to propagate the normals
        {
            std::queue<MSTNode *> border;
            border.push(mst_root);
            sprintf(message, "Orienting normals...");
            progress = 0;
            int maxSize = 0;
            int queueSize = 0;
            while ((queueSize = int(border.size())) > 0)
            {
                MSTNode *current_node = border.front();
                border.pop();
                // std::vector< MSTNode* >::iterator iMSTSon = current_node->sons.begin();
                // std::vector< MSTNode* >::iterator eMSTSon = current_node->sons.end();
                for (int s = 0; s < int(current_node->sons.size()); s++)
                {
                    if (callback != NULL && ((++progress % step) == 0) &&
                        (percentage = int((maxSize - queueSize) * 100 / maxSize)) < 100)
                        (callback)(percentage, message);

                    if (current_node->vertex->N().dot(current_node->sons[s]->vertex->N()) < ScalarType(0.0f))
                        current_node->sons[s]->vertex->N() *= ScalarType(-1.0f);
                    border.push(current_node->sons[s]);
                    maxSize = std::max<int>(maxSize, queueSize);
                }
            }
        }
        if (callback != NULL)
            (callback)(100, message);
    };
};

};  // end of namespace vcg

#endif  // end of VCG_SPACE_NORMAL_EXTRAPOLATION_H
