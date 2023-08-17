#include "ParameterMgr.h"
#include <iostream>

int ParameterMgr::init_time = 0;
ParameterMgr global_paraMgr;

ParameterMgr::ParameterMgr(void)
{
    init_time++;
    if (init_time > 1)
    {
        std::cout << "can not init ParameterMgr twice!" << endl;
        return;
    }

    grid_r = 0.2;

    initDataMgrParameter();
    initDrawerParameter();
    initGlareaParameter();
    initWLopParameter();
    initNormalSmootherParameter();
    initSkeletonParameter();
    initUpsamplingParameter();
    initL0filterParameter();
    initBilateralParameter();
    initParticleParameter();
    initNoise();
}

ParameterMgr::~ParameterMgr(void)
{
}

void ParameterMgr::setGlobalParameter(QString paraName, const Value& val)
{
    if (glarea.hasParameter(paraName))
        glarea.setValue(paraName, val);
    if (data.hasParameter(paraName))
        data.setValue(paraName, val);
    if (drawer.hasParameter(paraName))
        drawer.setValue(paraName, val);
    if (wLop.hasParameter(paraName))
        wLop.setValue(paraName, val);
    if (norSmooth.hasParameter(paraName))
        norSmooth.setValue(paraName, val);
    if (skeleton.hasParameter(paraName))
        skeleton.setValue(paraName, val);
    if (upsampling.hasParameter(paraName))
        upsampling.setValue(paraName, val);
    if (l0filter.hasParameter(paraName))
        l0filter.setValue(paraName, val);
    if (bilateral.hasParameter(paraName))
        bilateral.setValue(paraName,val);
    if (particle.hasParameter(paraName))
        particle.setValue(paraName,val);
    if (Noise.hasParameter(paraName))
        Noise.setValue(paraName, val);
}

void ParameterMgr::initDataMgrParameter()
{
    data.addParam(new RichDouble("Init Radius Para", 1.0));
    data.addParam(new RichDouble("Down Sample Num", 5000));
    data.addParam(new RichDouble("CGrid Radius", grid_r));
}

void ParameterMgr::initGlareaParameter()
{
    glarea.addParam(new RichString("Running Algorithm Name", ""));

    glarea.addParam(new RichBool("Light On or Off", false));

    glarea.addParam(new RichBool("Show Normal", false));
    glarea.addParam(new RichBool("Show Neighbor Normal", false));
    glarea.addParam(new RichBool("Show Face",true));
    glarea.addParam(new RichBool("Show EV",false));
    glarea.addParam(new RichBool("Show NEV",false));
    glarea.addParam(new RichBool("Show neighbors",false));

    glarea.addParam(new RichBool("Show Samples", true));
    glarea.addParam(new RichBool("Show Samples Quad", false));
    glarea.addParam(new RichBool("Show Samples Dot", true));
    glarea.addParam(new RichBool("Show Samples Circle", false));
    glarea.addParam(new RichBool("Show Samples Sphere", false));

    glarea.addParam(new RichBool("Show Original", true));
    glarea.addParam(new RichBool("Show Original Quad", false));
    glarea.addParam(new RichBool("Show Original Dot", true));
    glarea.addParam(new RichBool("Show Original Circle", false));
    glarea.addParam(new RichBool("Show Original Sphere", false));

    glarea.addParam(new RichBool("Show Noise", true));
    glarea.addParam(new RichBool("Show Noise Quad", false));
    glarea.addParam(new RichBool("Show Noise Dot", true));
    glarea.addParam(new RichBool("Show Noise Circle", false));
    glarea.addParam(new RichBool("Show Noise Sphere", false));

    glarea.addParam(new RichBool("Show Skeleton", true));

    glarea.addParam(new RichBool("Show Radius", true));
    glarea.addParam(new RichBool("Show All Radius", false));
    glarea.addParam(new RichBool("Show Radius Use Pick", true));
    glarea.addParam(new RichBool("Show Red Radius Line", true));
    glarea.addParam(new RichBool("Multiply Pick Point", true));

    glarea.addParam(new RichBool("GLarea Busying", false));
    glarea.addParam(new RichBool("Algorithom Stop", false));

    glarea.addParam(new RichPoint3f("Light Position", vcg::Point3f(-4.0, -4.0, -4.0)));
    glarea.addParam(new RichColor("Light Ambient Color", QColor(85, 85, 85)));
    glarea.addParam(new RichColor("Light Diffuse Color", QColor(218,165,32)));
    glarea.addParam(new RichColor("Light Specular Color", QColor(255, 255, 255)));

    // glarea.addParam(new RichPoint3f("Light Position", vcg::Point3f(4.0, 4.0, 4.0)));
    // glarea.addParam(new RichColor("Light Ambient Color", QColor(0.0, 0.0, 0.0)));
    // glarea.addParam(new RichColor("Light Diffuse Color", QColor(204, 204, 204)));
    // glarea.addParam(new RichColor("Light Specular Color", QColor(255, 255, 255)));

    glarea.addParam(new RichDouble("Snapshot Resolution", 2));
    glarea.addParam(new RichDouble("Snapshot Index", 1));
    glarea.addParam(new RichDouble("Radius Ball Transparency", 0.3));

    glarea.addParam(new RichBool("SnapShot Each Iteration", false));
    glarea.addParam(new RichBool("No Snap Radius", false));
}

void ParameterMgr::initDrawerParameter()
{
    drawer.addParam(new RichBool("Doing Pick", false));
    drawer.addParam(new RichBool("Need Cull Points", false));
    drawer.addParam(new RichBool("Use Pick Original", false));
    drawer.addParam(new RichBool("Use Pick Mode2", false));
    drawer.addParam(new RichBool("Skeleton Light", true));
    drawer.addParam(new RichBool("Show Individual Color", false));
    drawer.addParam(new RichBool("Use Color From Normal", false));
    drawer.addParam(new RichBool("Use Differ Branch Color", false));

    drawer.addParam(new RichDouble("Original Draw Width", 0.0015));
    drawer.addParam(new RichDouble("Sample Draw Width", 0.005));
    drawer.addParam(new RichDouble("Noise Draw Width", 0.005));
    drawer.addParam(new RichDouble("Sample Dot Size", 6));
    drawer.addParam(new RichDouble("Noise Dot Size", 6));
    drawer.addParam(new RichDouble("Original Dot Size", 1));
    drawer.addParam(new RichDouble("Normal Line Width", 6));
    drawer.addParam(new RichDouble("Normal Line Length", 1.5));

    drawer.addParam(new RichColor("Background Color", QColor(255, 255, 255)));
    drawer.addParam(new RichColor("Normal Line Color", QColor(0, 0, 255)));
    drawer.addParam(new RichColor("Sample Point Color", QColor(255, 0, 0)));
    drawer.addParam(new RichColor("Noise Point Color", QColor(48, 48, 48)));
    drawer.addParam(new RichColor("Original Point Color", QColor(48, 48, 48)));
    drawer.addParam(new RichColor("Feature Color", QColor(0, 0, 255)));
    drawer.addParam(new RichColor("Pick Point Color", QColor(128, 128, 0)));
    drawer.addParam(new RichColor("Pick Point DNN Color", QColor(0, 0, 155)));

    drawer.addParam(new RichColor("Skeleton Bone Color", QColor(200, 0, 0)));
    drawer.addParam(new RichColor("Skeleton Node Color", QColor(50, 250, 50)));
    drawer.addParam(new RichColor("Skeleton Branch Color", QColor(0, 0, 0)));
    drawer.addParam(new RichDouble("Skeleton Bone Width", 100));  // ./10000
    drawer.addParam(new RichDouble("Skeleton Node Size", 180));   // ./10000
    drawer.addParam(new RichDouble("Skeleton Branch Size", 30));  // abandoned
}

void ParameterMgr::initWLopParameter()
{
    wLop.addParam(new RichString("Algorithm Name", "WLOP"));
    wLop.addParam(new RichDouble("Num Of Iterate Time", 5));

    wLop.addParam(new RichDouble("CGrid Radius", grid_r));
    wLop.addParam(new RichDouble("H Gaussian Para", 1));
    wLop.addParam(new RichDouble("Repulsion Power", 1.0));
    wLop.addParam(new RichDouble("Average Power", 1.0));
    wLop.addParam(new RichBool("Need Compute Density", true));
    wLop.addParam(new RichBool("Need Compute PCA", false));
    wLop.addParam(new RichDouble("Repulsion Mu", 0.5));
    wLop.addParam(new RichDouble("Repulsion Mu2", 0.0));
    wLop.addParam(new RichBool("Run Anisotropic LOP", false));
    wLop.addParam(new RichDouble("Current Movement Error", 0.0));
}

void ParameterMgr::initSkeletonParameter()
{
    ///
    skeleton.addParam(new RichDouble("Repulsion Power", 1.0));
    skeleton.addParam(new RichDouble("Average Power", 2.0));

    ///
    skeleton.addParam(new RichDouble("Num Of Iterate Time", 1));
    skeleton.addParam(new RichString("Algorithm Name", "Skeletonization"));

    skeleton.addParam(new RichDouble("CGrid Radius", grid_r));
    skeleton.addParam(new RichDouble("H Gaussian Para", 4));
    skeleton.addParam(new RichBool("Need Compute Density", true));

    skeleton.addParam(new RichDouble("Current Movement Error", 0.0));
    skeleton.addParam(new RichBool("Run Auto Wlop One Step", false));
    skeleton.addParam(new RichBool("Run Auto Wlop One Stage", false));
    skeleton.addParam(new RichBool("The Skeletonlization Process Should Stop", false));

    skeleton.addParam(new RichBool("Step1 Detect Skeleton Feature", false));
    skeleton.addParam(new RichBool("Step2 Run Search New Branchs", false));
    skeleton.addParam(new RichBool("Step3 Clean And Update Radius", false));
    skeleton.addParam(new RichBool("Run Skeletonlization", false));

    // init
    skeleton.addParam(new RichDouble("Max Iterate Time", 55));
    skeleton.addParam(new RichDouble("Stop And Grow Error", 0.0005));
    skeleton.addParam(new RichDouble("Initial Radius", -1.));
    skeleton.addParam(new RichDouble("Radius Update Speed", 0.5));

    // step0
    skeleton.addParam(new RichDouble("Repulsion Mu", 0.35));
    skeleton.addParam(new RichDouble("Repulsion Mu2", 0.15));
    skeleton.addParam(new RichDouble("Follow Sample Radius", 0.33));
    skeleton.addParam(new RichDouble("Follow Sample Max Angle", 80));          // should add to UI
    skeleton.addParam(new RichDouble("Inactive And Keep Virtual Angle", 60));  // should add to UI
    skeleton.addParam(new RichDouble("Save Virtual Angle", 30));               // should add to UI

    skeleton.addParam(new RichDouble("Grow Accept Sigma", 0.8));  // should add to UI
    skeleton.addParam(new RichDouble("Bad Virtual Angle", 101));  // 2013-7-12

    // step1
    skeleton.addParam(new RichDouble("Combine Too Close Threshold", 0.01));
    skeleton.addParam(new RichDouble("Sigma KNN", 6));  // this one is hard to determine, should be small for narrow
                                                        // region, but will lead to unnecessary small branches
    skeleton.addParam(new RichDouble("Eigen Feature Identification Threshold", 0.901));

    // step2
    skeleton.addParam(new RichDouble("Branches Search Angle", 25));
    skeleton.addParam(new RichDouble("Virtual Head Accecpt Angle", 25));
    skeleton.addParam(new RichDouble("Snake Search Max Dist Blue", 0.4));
    skeleton.addParam(new RichDouble("Accept Branch Size", 6));  // important, and hard to determine
    skeleton.addParam(new RichDouble("Branch Search Max Dist Yellow", 0.1));

    skeleton.addParam(new RichDouble("Branches Merge Max Dist", 0.08));
    skeleton.addParam(new RichDouble("Branch Search KNN", 12));
    skeleton.addParam(new RichDouble("Combine Similar Angle", 140));
    skeleton.addParam(new RichDouble("Grow Search Radius", 0.15));
    skeleton.addParam(new RichDouble("Add Accept Branch Size", 1));

    // step3
    skeleton.addParam(new RichDouble("Clean Near Branches Dist", 0.05));
    skeleton.addParam(new RichDouble("Fix Original Weight", 0.91));
    skeleton.addParam(new RichDouble("Curve Segment Length", 0.051));
    skeleton.addParam(new RichInt("Fix Original Mode", 4));  // 1 for noisy , 4 for clean

    skeleton.addParam(new RichBool("Run ALL Segment", false));
    skeleton.addParam(new RichBool("Need Segment Right Away", true));
    skeleton.addParam(new RichDouble("Max Stop Radius", 1.99));

    // strategy...
    skeleton.addParam(new RichBool("Use Nearby Combine Strategy", true));
    skeleton.addParam(new RichBool("Use Go Through Strategy", false));
    skeleton.addParam(new RichBool("Use Aggresive Growth Strategy", false));
    skeleton.addParam(new RichBool("Use Clean Points When Following Strategy", true));
    skeleton.addParam(new RichBool("Use All Connect Strategy", true));
    skeleton.addParam(new RichBool("Use Plus Perpendicular Dist Strategy", false));
    skeleton.addParam(new RichBool("Use Kill Too Close Strategy", false));
    skeleton.addParam(new RichBool("Use Compute Eigen Ignore Branch Strategy", true));
    skeleton.addParam(new RichBool("Use Virtual Group Merge Strategy", false));
    skeleton.addParam(new RichBool("Use Final Merge Strategy", true));
    skeleton.addParam(new RichBool("Use Search New Twice Strategy", false));
    skeleton.addParam(new RichBool("Inactive Overlap Strategy", false));
    skeleton.addParam(new RichBool("Move Overlap Strategy", false));
    skeleton.addParam(new RichBool("Use Virtual Near Body Stop Strategy", true));
    skeleton.addParam(new RichBool("Need To Keep Big Bug", false));
    skeleton.addParam(new RichDouble("Change Strategy Radius", 0.45));
    skeleton.addParam(new RichBool("Need Recentering", true));
}

void ParameterMgr::initNormalSmootherParameter()
{
    norSmooth.addParam(new RichString("Algorithm Name", "NormalSmooth"));

    norSmooth.addParam(new RichInt("PCA KNN", 30));
    norSmooth.addParam(new RichDouble("CGrid Radius", grid_r));
    norSmooth.addParam(new RichDouble("Sharpe Feature Bandwidth Sigma", 30));
    norSmooth.addParam(new RichBool("Run Anistropic PCA", false));
    norSmooth.addParam(new RichBool("Run Init Samples Using Normal", false));

    norSmooth.addParam(new RichInt("Number Of Iterate", 1));
    norSmooth.addParam(new RichInt("Number of KNN", 400));

    norSmooth.addParam(new RichDouble("PCA Threshold", 0.8));
    norSmooth.addParam(new RichDouble("Feature threshold",0.03));
    norSmooth.addParam(new RichDouble("sigma n",0.3));

    norSmooth.addParam(new RichInt("Index",1));
}

void ParameterMgr::initUpsamplingParameter()
{
    upsampling.addParam(new RichString("Algorithm Name", "Upsampling"));

    upsampling.addParam(new RichDouble("CGrid Radius", 0.08));
    upsampling.addParam(new RichInt("Number of Add Point", 50000));
    upsampling.addParam(new RichDouble("Feature Sigma", 30));

    upsampling.addParam(new RichBool("Using Threshold Process", true));
    upsampling.addParam(new RichDouble("Dist Threshold", 0.02));
    upsampling.addParam(new RichDouble("Edge Parameter", 0.0));
    upsampling.addParam(new RichDouble("Z Parameter", 0.1));

    upsampling.addParam(new RichBool("Auto Recompute Radius For Dist", true));
    upsampling.addParam(new RichDouble("Min Dist Rate", 2.0));

    upsampling.addParam(new RichDouble("New Point Avg Sigma", 15));
    upsampling.addParam(new RichBool("Use Avg Normal Method", false));
    upsampling.addParam(new RichBool("Use Max Theta Psi Method", false));
    upsampling.addParam(new RichBool("Use Sigma Threshold Method", true));
    upsampling.addParam(new RichBool("Use No Psi Method", false));

    upsampling.addParam(new RichDouble("Upsample Radius", grid_r * 0.5));
    upsampling.addParam(new RichBool("Use Proj New Term", false));
    upsampling.addParam(new RichBool("Run Projection", false));
}

void ParameterMgr::initL0filterParameter()
{
    l0filter.addParam(new RichString("Algorithm Name", "L0filter"));
    l0filter.addParam(new RichInt("normKNN",15));
    l0filter.addParam(new RichDouble("norm sparsity",0.02));
    l0filter.addParam(new RichDouble("norm auxiliary",0.1));
    l0filter.addParam(new RichDouble("norm increment",2.0));
    l0filter.addParam(new RichInt("norm iteration",8));

    l0filter.addParam(new RichInt("vertKNN",10));
    l0filter.addParam(new RichDouble("vert sparsity",0.002));
    l0filter.addParam(new RichDouble("vert auxiliary",0.0005));
    l0filter.addParam(new RichDouble("vert increment",2.0));
    l0filter.addParam(new RichInt("vert iteration",6));


    l0filter.addParam(new RichInt("Total Iteration",6));
    l0filter.addParam(new RichBool("norm L0",true));
    l0filter.addParam(new RichBool("vert L0",true));
    l0filter.addParam(new RichBool("edge recovery",true));

    l0filter.addParam(new RichInt("EdgeKNN",10));
    l0filter.addParam(new RichInt("Edge Iteration",10));
    l0filter.addParam(new RichDouble("Edge feature scale",0.5));
    l0filter.addParam(new RichDouble("Edge step size",0.001));
}

void ParameterMgr::initBilateralParameter()
{
    bilateral.addParam(new RichString("Algorithm Name", "Bilateral"));
    bilateral.addParam(new RichInt("total iteration",3));
    bilateral.addParam(new RichInt("KNN",20));
}

void ParameterMgr::initParticleParameter()
{
    particle.addParam(new RichString("Algorithm Name", "Particle"));
    particle.addParam(new RichInt("total iteration",3));
    particle.addParam(new RichInt("particle number",10));
    particle.addParam(new RichInt("RANSAC iteration",150));
    particle.addParam(new RichInt("norm iteration",50));
    particle.addParam(new RichDouble("feat",0.05));
    //particle.addParam(new RichInt("miu",0.1));

    particle.addParam(new RichBool("edge",false));
    particle.addParam(new RichBool("normal",true));
}
void ParameterMgr::initNoise()
{
    Noise.addParam(new RichDouble("Noise Level",0.2));
    Noise.addParam(new RichDouble("Impulsive Level",0.2));
    Noise.addParam(new RichInt("Noise Type",0));
    Noise.addParam(new RichInt("Noise Direction",0));
}
