#include "noise.h"

Noise::Noise(CMesh &m)
{

}
void Noise::addNoise(CMesh *samples,const CMesh::VertexIterator &begin,const CMesh::VertexIterator  &end)
{
    //设置参数
     double radius1 = global_paraMgr.norSmooth.getDouble("CGrid Radius");
     double noise_level = global_paraMgr.Noise.getDouble("Noise Level");
     double impulsive_level = global_paraMgr.Noise.getDouble("Impulsive Level");
     int noise_type_index,noise_direction_index;
     noise_type_index = global_paraMgr.Noise.getInt("Noise Type");
     noise_direction_index = global_paraMgr.Noise.getInt("Noise Direction");

     double radius2 = radius1 * radius1;
     for(CMesh::VertexIterator dest = samples->vert.begin(); dest != samples->vert.end(); dest++)
      dest->neighbors.clear();
     for(CMesh::VertexIterator dest = samples->vert.begin(); dest != samples->vert.end(); dest++)
     {
        Point3f p = dest->P();
         for(CMesh::VertexIterator origin = dest + 1; origin != samples->vert.end(); origin++)
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

     NoiseType noise_type = (noise_type_index == 0)? kGaussian : kImpulsive;
     NoiseDirection noise_direction = (noise_direction_index == 0)? kNormal : kRandom;
     // 计算平均边长

//     for(size_t i = 0; i < edge_number; ++i)
//     {
//        CVertex  *p0 = samples->edge[i].V(0);
//        CVertex  *p1 = samples->edge[i].V(1);
//        Point3f P0 = p0->P(); Point3f P1 = p1->P();
//        average_length += (P0 - P1).Norm();
//     }
//     average_length /= edge_number;

     //add noise

     std::vector<double> GaussianNumbers;
     std::vector<vcg::Point3f> RandomDirections;
     std::vector<std::pair<int, double> > VertexListAndGaussianNumbers;

     int ptnum_ = samples->vert.size();

     double average_length;
     average_length = averageKNNDistance(samples,ptnum_);

     double standard_derivation = average_length * noise_level;
     int impulsive_vertex_number = (ptnum_ * impulsive_level);


     if(noise_type == kGaussian){
         randomGaussianNumbers(0, standard_derivation, ptnum_, GaussianNumbers);
         if(noise_direction == kNormal){
             for(CMesh::VertexIterator v_it = begin; v_it != end; v_it++){

                 vcg::Point3f n = v_it->N();
                 vcg::Point3f p = v_it->P() + v_it->N() * GaussianNumbers[v_it->m_index];
                 v_it->P() = p;
             }
         }
         else if(noise_direction == kRandom){
             randomDirections(ptnum_, RandomDirections);
             for(CMesh::VertexIterator v_it = begin; v_it != end; v_it++){
                 int index = v_it->m_index;
                 vcg::Point3f p = v_it->P() + RandomDirections[index] * GaussianNumbers[index];
                 v_it->P() = p;
             }
         }
     }
     else if(noise_type == kImpulsive){
         randomImpulsiveNumbers(0,ptnum_ - 1, impulsive_vertex_number, 0, standard_derivation, VertexListAndGaussianNumbers);
         if(noise_direction == kNormal){
             for(int i = 0; i < (int)VertexListAndGaussianNumbers.size(); i++){
                 int index = VertexListAndGaussianNumbers[i].first;

                 vcg::Point3f p =  samples->vert[index].P() + samples->vert[index].N() * VertexListAndGaussianNumbers[i].second;
                 samples->vert[index].P() = p;
             }
         }
         else if(noise_direction == kRandom){
             randomDirections(impulsive_vertex_number, RandomDirections);
             for(int i = 0; i < (int)VertexListAndGaussianNumbers.size(); i++){
                 int index = VertexListAndGaussianNumbers[i].first;
                 vcg::Point3f p =  samples->vert[index].P() + RandomDirections[i] * VertexListAndGaussianNumbers[i].second;
                  samples->vert[index].P() = p;
             }
         }
     }


}
double Noise::generateRandomGaussian(double mean, double StandardDerivation)
{
    static double v1, v2, s;
    static int phase = 0;
    double x;

    if(phase == 0)
    {
        do
        {
            v1 = -1 + 2 * (double)rand() / (double) RAND_MAX;
            v2 = -1 + 2 * (double)rand() / (double) RAND_MAX;
            s = v1 * v1 + v2 * v2;
        }while (s >= 1 || s == 0);

        x = v1 * sqrt(-2 * log(s) / s);
    }
    else
        x = v2 * sqrt(-2 * log(s) / s);

    phase = 1 - phase;

    return x * StandardDerivation + mean;
}
vcg::Point3f Noise::generateRandomDirection()
{
    double x, y, z, length;
    do{
        x = -1 + 2 * (double)rand() / (double) RAND_MAX;
        y = -1 + 2 * (double)rand() / (double) RAND_MAX;
        length = x * x + y * y;
    }while(length>1);

    const double r = 2 * std::sqrt(1 - length);

    x *= r;
    y *= r;
    z = 1 - 2 * length;

    return Point3f(x, y, z);
}
void Noise::randomGaussianNumbers(double mean, double StandardDerivation, int number, std::vector<double> &RandomNumbers)
{
    RandomNumbers.resize(number, 0.0);

    srand((unsigned int)time(NULL));
    for(int i = 0; i < number; i++){
        RandomNumbers[i] = generateRandomGaussian(mean, StandardDerivation);
    }
}
void Noise::randomImpulsiveNumbers(int min, int max, int number, double mean, double StandardDerivation,
                                   std::vector<std::pair<int, double> > &VertexListAndRandomNumbers)
{
    int range = max - min + 1;
    if(number > range) return;

    VertexListAndRandomNumbers.resize(number, std::make_pair(0, 0.0));

    std::vector<double> randomNumbers;
    randomGaussianNumbers(mean, StandardDerivation, number, randomNumbers);

    srand((unsigned int)time(NULL));
    std::vector<int> rangeVector(range);
    for(int i = 0; i < range; i++)
        rangeVector[i] = min + i;

    srand((unsigned int)time(NULL));
    std::vector<int> vertexIndexList(number);
    for(int i = 0; i < number; i++){
        int pos = (int)((double)rand() / RAND_MAX * range);
        vertexIndexList[i] = rangeVector[pos];
        range--;
        std::swap(rangeVector[pos], rangeVector[range]);
    }

    for(int i = 0; i < number; i++)
        VertexListAndRandomNumbers[i] = std::make_pair(vertexIndexList[i], randomNumbers[i]);
}
void Noise::randomDirections(int number, std::vector<vcg::Point3f> &RandomDirections)
{
    RandomDirections.resize(number, Point3f(0.0, 0.0, 0.0));

    srand((unsigned int)time(NULL));
    for(int i = 0; i < number; i++){
        RandomDirections[i] = generateRandomDirection();
    }
}
double Noise::averageKNNDistance(CMesh *sample_,double ptnum_)
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
