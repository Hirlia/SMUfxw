#include "dlg_particle_para.h"
#include "ui_dlg_particle_para.h"

dlg_particle_para::dlg_particle_para(QWidget *parent, ParameterMgr *_paras, GLArea *_area) :
    QWidget(parent),
    ui(new Ui::dlg_particle_para)
{
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

dlg_particle_para::~dlg_particle_para()
{
    delete ui;
}

void dlg_particle_para::initConnects()
{
    connect(ui->totaliteration, SIGNAL(valueChanged(int)), this, SLOT(gettotalnum(int)));
    connect(ui->particlenum, SIGNAL(valueChanged(int)), this, SLOT(getparticlenum(int)));
    connect(ui->particle, SIGNAL(clicked()), this, SLOT(applyparticlefilter()));
    connect(ui->edge, SIGNAL(clicked(bool)), this, SLOT(runedge(bool)));
    connect(ui->normal, SIGNAL(clicked(bool)), this, SLOT(runnormal(bool)));
    connect(ui->feat, SIGNAL(valueChanged(double)), this, SLOT(getfeat(double)));
    connect(ui->pushButton_Projection_err,SIGNAL(clicked(bool)),this,SLOT(apply_Projection_err()));
}

bool dlg_particle_para::initWidgets()
{
    ui->totaliteration->setValue(m_paras->particle.getInt("total iteration"));
    ui->particlenum->setValue(m_paras->particle.getInt("particle number"));
    ui->feat->setValue(m_paras->particle.getDouble("feat"));
    Qt::CheckState state =
    m_paras->particle.getBool("edge") ? (Qt::CheckState::Checked) : (Qt::CheckState::Unchecked);
    ui->edge->setCheckState(state);
    Qt::CheckState state_normal =
    m_paras->particle.getBool("normal") ? (Qt::CheckState::Checked) : (Qt::CheckState::Unchecked);
    ui->normal->setCheckState(state_normal);
    return true;
}

void dlg_particle_para::gettotalnum(int _val)
{
    m_paras->particle.setValue("total iteration", IntValue(_val));
}

void dlg_particle_para::getparticlenum(int _val)
{
    m_paras->particle.setValue("particle number", IntValue(_val));
}

void dlg_particle_para::runedge(bool _val)
{
    m_paras->particle.setValue("edge", BoolValue(_val));
}

void dlg_particle_para::runnormal(bool _val)
{
    m_paras->particle.setValue("normal", BoolValue(_val));
}

void dlg_particle_para::getfeat(double _val)
{
    m_paras->particle.setValue("feat",DoubleValue(_val));
}

void dlg_particle_para::applyparticlefilter()
{
    global_paraMgr.glarea.setValue("Running Algorithm Name", StringValue("Particlefilter"));
    area->runParticle();
    area->dataMgr.recomputeQuad();
    area->updateGL();
}

void dlg_particle_para::apply_Projection_err()
{
    CMesh *samples = area->dataMgr.getCurrentSamples();
    CMesh *original = area->dataMgr.getCurrentOriginal();
    Particle *l0_minization1 ;
    double Projection_err = l0_minization1->getPointtoFace(samples,original);
    std::cout<< "Projection err = " << Projection_err<<std::endl;
}
