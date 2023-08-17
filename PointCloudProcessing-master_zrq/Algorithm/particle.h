#ifndef PARTICLE_H
#define PARTICLE_H
#include "GlobalFunction.h"
#include "normal_extrapolation.h"
#include "DataMgr.h"
#include "ParameterMgr.h"
#include "PointCloudAlgorithm.h"
#include <cmath>

using namespace std;
using namespace vcg;

class Particle : public PointCloudAlgorithm
{
public:
    Particle(RichParameterSet* _para);
    ~Particle(void);

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
    int ptnum_;
    double radius;
    int particlenum;
    double beta;
    double bandwidth;
    vector<double> bandwidth_;

    void runparticle();
    void runnormalparticle();
    void runRANSAC();
    void run_normal();
    void runparticleedge();
    double averageKNNDistance();
    double calculate_bandwidth();
    double getPointtoFace(CMesh *sample,CMesh *original);

private:
    RichParameterSet* para;

public:
    CMesh* sample_;
    CMesh* origin_;
};

#endif // PARTICLE_H
