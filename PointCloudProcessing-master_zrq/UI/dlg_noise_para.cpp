#include "dlg_noise_para.h"

NoiseParaDlg::NoiseParaDlg(QWidget *p, ParameterMgr *_paras, GLArea *_area) : QFrame(p)
{
    ui = new Ui::dlg_noise_para;
    NoiseParaDlg::ui->setupUi(this);
    area = _area;
    m_paras = _paras;

    if (!initWidgets())
    {
        cerr << " dlg_lo_filter_para::initWidgets failed." << endl;
        return;
    }
    initConnects();

}
void NoiseParaDlg::initConnects()
{
    connect(ui->noise_level,SIGNAL(valueChanged(double)),this,SLOT(getNoiseLevel(double)));
    connect(ui->Impulsive_level,SIGNAL(valueChanged(double)),this,SLOT(getImpulsiveLevel(double)));
    connect(ui->Noise_type,SIGNAL(currentIndexChanged(int)),this,SLOT(getNoiseType(int)));
    connect(ui->Noise_direction,SIGNAL(currentIndexChanged(int)),this,SLOT(getNoiseDirection(int)));
    connect(ui->pushButton_noise,SIGNAL(clicked()),this,SLOT(apply_noise()));
    connect(ui->pushButton_reset,SIGNAL(clicked()),this,SLOT(apply_reset()));
}
bool NoiseParaDlg::initWidgets()
{
    ui->noise_level->setValue(m_paras->Noise.getDouble("Noise Level"));
    ui->Impulsive_level->setValue(m_paras->Noise.getDouble("Impulsive Level"));
    ui->Noise_type->setCurrentIndex(m_paras->Noise.getInt("Noise Type"));
    ui->Noise_direction->setCurrentIndex(m_paras->Noise.getInt("Noise Direction"));
    return true;
}



void NoiseParaDlg::getNoiseLevel(double _val)
{
     m_paras->setGlobalParameter("Noise Level", DoubleValue(_val));
}

void NoiseParaDlg::getImpulsiveLevel(double _val)
{
    m_paras->setGlobalParameter("Impulsive Level",DoubleValue(_val));
}

void NoiseParaDlg::getNoiseType(int _val)
{
    if(ui->Noise_type->currentIndex()==0)
       m_paras->setGlobalParameter("Gaussion",IntValue(_val));
    else if(ui->Noise_type->currentIndex()==1)
       m_paras->setGlobalParameter("Impulsive",IntValue(_val));
}

void NoiseParaDlg::getNoiseDirection(int _val)
{
    if(ui->Noise_direction->currentIndex()==0)
        m_paras->setGlobalParameter("Impulsive",IntValue(_val));
    else if(ui->Noise_direction->currentIndex()==1)
        m_paras->setGlobalParameter("Random",IntValue(_val));
}

void NoiseParaDlg::apply_noise()
{
    CMesh *samples = area->dataMgr.getCurrentNoise();
    CMesh *sample = area->dataMgr.getCurrentSamples();
    Noise *noise = new Noise(*samples);
    noise->addNoise(samples,samples->vert.begin(), samples->vert.end());
    *sample = *samples;
    area->initSetting();
    area->updateGL();
}
void NoiseParaDlg::apply_reset()
{
    CMesh *samples = area->dataMgr.getCurrentNoise();
    CMesh *sample = area->dataMgr.getCurrentSamples();
    *sample = *samples;
    area->initSetting();
    area->updateGL();
}

NoiseParaDlg::~NoiseParaDlg()
{
    delete ui;
}
