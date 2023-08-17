#pragma once

#include <QDialog>
#include <QtGui/QDockWidget>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QWidget>
#include <QtGui>
#include <iostream>

#include "ParameterMgr.h"
#include "UI/dlg_normal_para.h"
#include "UI/dlg_skeleton_para.h"
#include "UI/dlg_upsampling_para.h"
#include "UI/dlg_wlop_para.h"
#include "UI/dlg_l0filter_para.h"
#include "UI/dlg_particle_para.h"
#include "UI/dlg_noise_para.h"

#include "GLArea.h"

using namespace std;

class StdParaDlg : public QDockWidget
{
    Q_OBJECT
  public:
    StdParaDlg(ParameterMgr* _paras, GLArea* _area, QWidget* parent = 0);
    ~StdParaDlg();

    bool showWlopParaDialog();
    bool showNormalParaDlg();
    bool showSkeletonParaDlg();
    bool showUpsamplingParaDlg();
    bool showL0filterParaDlg();
    bool showParticleParaDlg();
    bool showNoiseParaDlg();

  private:
    void init();
    void createFrame();
    void loadWlopFrame();
    void loadNormalFrame();
    void loadSkeletonFrame();
    void loadUpsamplingFrame();
    void loadL0filterFrame();
    void loadParticleFrame();
    void loadNoiseFrame();

  private slots:
    void closeClick();

  private:
    WlopParaDlg* para_wlop;
    NormalParaDlg* para_normal;
    SkeletonParaDlg* para_skeleton;
    UpsamplingParaDlg* para_upsampling;
    dlg_l0filter_para* para_l0;
    dlg_particle_para* para_particle;
    NoiseParaDlg* para_noise;

    ParameterMgr* paras;

    QFrame* mainFrame;
    GLArea* gla;
};
