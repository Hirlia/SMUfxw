#pragma once
#include "CMesh.h"
#include "Parameter.h"

class ParameterMgr
{
  public:
    ParameterMgr(void);
    ~ParameterMgr(void);
    RichParameterSet* getDataParameterSet()
    {
        return &data;
    }
    RichParameterSet* getDrawerParameterSet()
    {
        return &drawer;
    }
    RichParameterSet* getGlareaParameterSet()
    {
        return &glarea;
    }
    RichParameterSet* getWLopParameterSet()
    {
        return &wLop;
    }
    RichParameterSet* getSkeletonParameterSet()
    {
        return &skeleton;
    }
    RichParameterSet* getNormalSmootherParameterSet()
    {
        return &norSmooth;
    }
    RichParameterSet* getUpsamplingParameterSet()
    {
        return &upsampling;
    }
    RichParameterSet* getL0filterParameterSet()
    {
        return &l0filter;
    }
    RichParameterSet* getBilateralParameterSet()
    {
        return &bilateral;
    }
    RichParameterSet* getParticleParameterSet()
    {
        return &particle;
    }
    RichParameterSet* getNoiseParameterSet()
    {
        return &Noise;
    }


    void setGlobalParameter(QString paraName, const Value& val);
    typedef enum { GLAREA, DATA, DRAWER, WLOP, NOR_SMOOTH, SKELETON, UPSAMPLING,L0FILTER,PARTICLE } ParaType;

  private:
    void initDataMgrParameter();
    void initDrawerParameter();
    void initGlareaParameter();
    void initWLopParameter();
    void initSkeletonParameter();
    void initNormalSmootherParameter();
    void initUpsamplingParameter();
    void initL0filterParameter();
    void initBilateralParameter();
    void initParticleParameter();
    void initNoise();

  public:
    RichParameterSet glarea;
    RichParameterSet data;
    RichParameterSet drawer;
    RichParameterSet wLop;
    RichParameterSet norSmooth;
    RichParameterSet skeleton;
    RichParameterSet upsampling;
    RichParameterSet l0filter;
    RichParameterSet bilateral;
    RichParameterSet particle;
    RichParameterSet Noise;

  private:
    static int init_time;
    double grid_r;
};

extern ParameterMgr global_paraMgr;
