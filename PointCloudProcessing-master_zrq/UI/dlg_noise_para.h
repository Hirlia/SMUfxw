#pragma once

#include <QtGui/QFrame>
#include <QtGui/QWidget>
#include <QtGui>
#include <iostream>
#include "GLArea.h"
#include "ParameterMgr.h"
#include "Algorithm/PointCloudAlgorithm.h"
#include "Algorithm/noise.h"
#include <QWidget>
#include "ui_dlg_noise_para.h"
using namespace std;

class NoiseParaDlg : public QFrame
{
    Q_OBJECT

public:
    NoiseParaDlg(QWidget *p, ParameterMgr *_paras, GLArea *_area);
    ~NoiseParaDlg();
    void initConnects();


private slots:
  bool initWidgets();
  void getNoiseLevel(double _val);
  void getImpulsiveLevel(double _val);
  void getNoiseType(int _val);
  void getNoiseDirection(int _val);

  void apply_noise();
  void apply_reset();



private:
  Ui::dlg_noise_para *ui;
  ParameterMgr *m_paras;
  GLArea *area;

};

// DLG_NOISE_PARA_H
