#ifndef DLG_PARTICLE_PARA_H
#define DLG_PARTICLE_PARA_H

#include <QWidget>
#include "Algorithm/normal_extrapolation.h"
#include "GLArea.h"
#include "ParameterMgr.h"
#include "calculationthread.h"
#include "Algorithm/L0filter.h"

namespace Ui {
class dlg_particle_para;
}

class dlg_particle_para : public QWidget
{
    Q_OBJECT

public:
    explicit dlg_particle_para(QWidget *parent ,ParameterMgr *_paras, GLArea *_area);
    ~dlg_particle_para();
    void initConnects();

private slots:
    bool initWidgets();
    void gettotalnum(int _val);
    void getparticlenum(int _val);
    void getfeat(double _val);
    void runedge(bool _val);
    void runnormal(bool _val);

    void applyparticlefilter();
    void apply_Projection_err();


private:
    Ui::dlg_particle_para *ui;
    ParameterMgr *m_paras;
    GLArea *area;
    CalculationThread calculation_thread;
};

#endif // DLG_PARTICLE_PARA_H
