#ifndef NOISE_H
#define NOISE_H
#include "CMesh.h"
#include <vcg/math/disjoint_set.h>
#include <vcg/math/lin_algebra.h>
#include <vcg/math/linear.h>
#include <vcg/math/matrix33.h>
#include <vcg/space/box3.h>
#include <vcg/space/index/octree.h>
#include <vcg/space/point3.h>
#include <wrap/callback.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <stdlib.h>
#include <algorithm>
#include <limits>
#include <queue>
#include <vector>
#include <numeric>
#include <utility>
#include <QList>
#include "DataMgr.h"
#include "ParameterMgr.h"

using namespace Eigen;
using namespace std;
class Noise
{
public:
    Noise(CMesh &m);
public:
     void addNoise(CMesh *samples,const CMesh::VertexIterator &begin,const CMesh::VertexIterator  &end);
private:
     enum NoiseType{kGaussian, kImpulsive};
     enum NoiseDirection{kNormal, kRandom};

     double generateRandomGaussian(double mean, double StandardDerivation);
     vcg::Point3f generateRandomDirection();

     void randomGaussianNumbers(double mean, double StandardDerivation, int number, std::vector<double> &RandomNumbers);
     void randomImpulsiveNumbers(int min, int max, int number,
                                   double mean, double StandardDerivation,
                                   std::vector< std::pair<int, double> > &VertexListAndRandomNumbers);
     void randomDirections(int number, std::vector<vcg::Point3f> &RandomDirections);
     double averageKNNDistance(CMesh *sample_,double ptnum_);


};

#endif // NOISE_H
