#pragma once
#include <iostream>
#include "GlobalFunction.h"
#include "normal_extrapolation.h"
#include "DataMgr.h"
#include "ParameterMgr.h"
#include "PointCloudAlgorithm.h"

using namespace std;
using namespace vcg;

class Bilateral : public PointCloudAlgorithm
{
  public:
    Bilateral(RichParameterSet* para);
    ~Bilateral(void);

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

    int knn;
    int total;
    vector<CVertex> vert;

    void bilateralfilter(vector<CVertex> &vertices);
    void calculate_sigmac(vector<CVertex> &vertices);

private:
    RichParameterSet* para;

public:
    CMesh* sample_;
};
