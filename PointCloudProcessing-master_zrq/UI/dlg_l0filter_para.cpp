#include "dlg_l0filter_para.h"
#include "ui_dlg_l0filter_para.h"

dlg_l0filter_para::dlg_l0filter_para(QWidget *parent, ParameterMgr* _paras, GLArea* _area) :QWidget(parent)
{
    ui = new Ui::dlg_l0filter_para;
    ui->setupUi(this);
    area = _area;
    m_paras = _paras;

    if (!initWidgets())
    {
        cerr << " L0filterParaDlg::initWidgets failed." << endl;
        return;
    }
    initConnects();
}

dlg_l0filter_para::~dlg_l0filter_para()
{
    delete ui;
}

void dlg_l0filter_para::initConnects()
{
    connect(ui->normKnn, SIGNAL(valueChanged(int)), this, SLOT(getnormKNN(int)));
    connect(ui->vertKnn, SIGNAL(valueChanged(int)), this, SLOT(getvertKNN(int)));
    connect(ui->edgeKnn, SIGNAL(valueChanged(int)), this, SLOT(getedgeKNN(int)));
    connect(ui->norm_sparsity, SIGNAL(valueChanged(double)), this, SLOT(getnormsparsity(double)));
    connect(ui->vert_sparsity, SIGNAL(valueChanged(double)), this, SLOT(getvertsparsity(double)));
    connect(ui->norm_auxiliary, SIGNAL(valueChanged(double)), this, SLOT(getnormauxiliary(double)));
    connect(ui->vert_auxiliary, SIGNAL(valueChanged(double)), this, SLOT(getvertauxiliary(double)));
    connect(ui->norm_increment, SIGNAL(valueChanged(double)), this, SLOT(getnormincrement(double)));
    connect(ui->vert_increment, SIGNAL(valueChanged(double)), this, SLOT(getvertincrement(double)));
    connect(ui->norm_iteration, SIGNAL(valueChanged(int)), this, SLOT(getnormIterateNum(int)));
    connect(ui->vert_iteration, SIGNAL(valueChanged(int)), this, SLOT(getvertIterateNum(int)));
    connect(ui->total_iteration,SIGNAL(valueChanged(int)), this, SLOT(gettotalIterateNum(int)));
    connect(ui->edge_iteration,SIGNAL(valueChanged(int)), this, SLOT(getedgeIterateNum(int)));
    connect(ui->edge_stepsize, SIGNAL(valueChanged(double)), this, SLOT(getedgeStepsize(double)));

    connect(ui->norm_L0, SIGNAL(clicked(bool)), this, SLOT(isnormL0(bool)));
    connect(ui->vert_L0, SIGNAL(clicked(bool)), this, SLOT(isvertL0(bool)));
    connect(ui->edge_recovery, SIGNAL(clicked(bool)), this, SLOT(isedgerecovery(bool)));

    connect(ui->L0filter, SIGNAL(clicked()), this, SLOT(applyL0filter()));
}

bool dlg_l0filter_para::initWidgets()
{
    ui->normKnn->setValue(m_paras->l0filter.getInt("normKNN"));
    ui->vertKnn->setValue(m_paras->l0filter.getInt("vertKNN"));
    ui->norm_sparsity->setValue(m_paras->l0filter.getDouble("norm sparsity"));
    ui->vert_sparsity->setValue(m_paras->l0filter.getDouble("vert sparsity"));
    ui->norm_auxiliary->setValue(m_paras->l0filter.getDouble("norm auxiliary"));
    ui->vert_auxiliary->setValue(m_paras->l0filter.getDouble("vert auxiliary"));
    ui->norm_increment->setValue(m_paras->l0filter.getDouble("norm increment"));
    ui->vert_increment->setValue(m_paras->l0filter.getDouble("vert increment"));
    ui->norm_iteration->setValue(m_paras->l0filter.getInt("norm iteration"));
    ui->vert_iteration->setValue(m_paras->l0filter.getInt("vert iteration"));
    ui->total_iteration->setValue(m_paras->l0filter.getInt("Total Iteration"));
    ui->edgeKnn->setValue(m_paras->l0filter.getInt("EdgeKNN"));
    ui->edge_iteration->setValue(m_paras->l0filter.getInt("Edge Iteration"));
    ui->edge_stepsize->setValue(m_paras->l0filter.getDouble("Edge step size"));
    Qt::CheckState state1 =
    m_paras->l0filter.getBool("norm L0") ? (Qt::CheckState::Checked) : (Qt::CheckState::Unchecked);
    ui->norm_L0->setCheckState(state1);
    Qt::CheckState state2 =
    m_paras->l0filter.getBool("vert L0") ? (Qt::CheckState::Checked) : (Qt::CheckState::Unchecked);
    ui->vert_L0->setCheckState(state2);
    Qt::CheckState state3 =
    m_paras->l0filter.getBool("norm L0") ? (Qt::CheckState::Checked) : (Qt::CheckState::Unchecked);
    ui->edge_recovery->setCheckState(state3);

    return true;
}


void dlg_l0filter_para::getnormKNN(int _val)
{
    m_paras->l0filter.setValue("normKNN", IntValue(_val));
}

void dlg_l0filter_para::getvertKNN(int _val)
{
    m_paras->l0filter.setValue("vertKNN", IntValue(_val));
}

void dlg_l0filter_para::getedgeKNN(int _val)
{
    m_paras->l0filter.setValue("EdgeKNN", IntValue(_val));
}

void dlg_l0filter_para::getnormsparsity(double _val)
{
    m_paras->l0filter.setValue("norm sparsity", DoubleValue(_val));
}

void dlg_l0filter_para::getvertsparsity(double _val)
{
    m_paras->l0filter.setValue("vert sparsity", DoubleValue(_val));
}

void dlg_l0filter_para::getnormauxiliary(double _val)
{
    m_paras->l0filter.setValue("norm auxiliary", DoubleValue(_val));
}

void dlg_l0filter_para::getvertauxiliary(double _val)
{
    m_paras->l0filter.setValue("vert auxiliary", DoubleValue(_val));
}

void dlg_l0filter_para::getnormincrement(double _val)
{
    m_paras->l0filter.setValue("norm increment", DoubleValue(_val));
}

void dlg_l0filter_para::getvertincrement(double _val)
{
    m_paras->l0filter.setValue("vert increment", DoubleValue(_val));
}

void dlg_l0filter_para::getnormIterateNum(int _val)
{
    m_paras->l0filter.setValue("norm iteration", IntValue(_val));
}

void dlg_l0filter_para::getvertIterateNum(int _val)
{
    m_paras->l0filter.setValue("vert iteration", IntValue(_val));
}

void dlg_l0filter_para::gettotalIterateNum(int _val)
{
    m_paras->l0filter.setValue("Total Iteration", IntValue(_val));
}

void dlg_l0filter_para::getedgeIterateNum(int _val)
{
    m_paras->l0filter.setValue("Edge Iteration", IntValue(_val));
}

void dlg_l0filter_para::getedgeStepsize(double _val)
{
    m_paras->l0filter.setValue("Edge step size", DoubleValue(_val));
}

void dlg_l0filter_para::isnormL0(bool _val)
{
    m_paras->l0filter.setValue("norm L0", BoolValue(_val));
}

void dlg_l0filter_para::isvertL0(bool _val)
{
    m_paras->l0filter.setValue("vert L0", BoolValue(_val));
}

void dlg_l0filter_para::isedgerecovery(bool _val)
{
    m_paras->l0filter.setValue("edge recovery", BoolValue(_val));
}


void dlg_l0filter_para::applyL0filter()
{
    //m_paras->l0filter.setValue("Run Anisotropic LOP", BoolValue(true));
    global_paraMgr.glarea.setValue("Running Algorithm Name", StringValue("L0filter"));
    area->runL0filter();
    area->dataMgr.recomputeQuad();
    area->updateGL();
}
