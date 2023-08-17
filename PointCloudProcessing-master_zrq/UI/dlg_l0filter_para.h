#ifndef DLG_L0FILTER_PARA_H
#define DLG_L0FILTER_PARA_H

#include <QWidget>
#include "ParameterMgr.h"
#include "GLArea.h"
#include "Algorithm/L0filter.h"
#include "calculationthread.h"

namespace Ui {
class dlg_l0filter_para;
}

class dlg_l0filter_para : public QWidget
{
    Q_OBJECT

public:
    explicit dlg_l0filter_para(QWidget *parent , ParameterMgr *_paras, GLArea *_area);
    ~dlg_l0filter_para();
    void initConnects();
    //void initial();

private slots:
  bool initWidgets();
  void getnormKNN(int _val);
  void getvertKNN(int _val);
  void getedgeKNN(int _val);
  void getnormsparsity(double _val);
  void getvertsparsity(double _val);
  void getnormauxiliary(double _val);
  void getvertauxiliary(double _val);
  void getnormincrement(double _val);
  void getvertincrement(double _val);
  void getnormIterateNum(int _val);
  void getvertIterateNum(int _val);
  void gettotalIterateNum(int _val);
  void getedgeIterateNum(int _val);
  void getedgeStepsize(double _val);
  void isnormL0(bool _val);
  void isvertL0(bool _val);
  void isedgerecovery(bool _val);

  void applyL0filter();

private:
    Ui::dlg_l0filter_para *ui;
    ParameterMgr *m_paras;
    GLArea *area;
    CalculationThread calculation_thread;

};

#endif // DLG_L0FILTER_PARA_H
