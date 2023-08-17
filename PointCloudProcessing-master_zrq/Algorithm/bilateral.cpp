#include "bilateral.h"

Bilateral::Bilateral(RichParameterSet* _para)
{
    cout << "L0 constructed!!" << endl;
    para = _para;
}

Bilateral::~Bilateral(void)
{
    cout << "L0 destroy!! " << endl;
}

void Bilateral::setInput(DataMgr* pData)
{

    sample_ = pData->getCurrentSamples();
    cout<<"setInput"<<endl;
}

void Bilateral::initial()
{
    total = para->getInt("total iteration");
    knn = para->getInt("KNN");
    vert = sample_->vert;

}
void Bilateral::run()
{
    cout<<"run"<<endl;
    initial();
    for(int i=0;i<total;i++){
        calculate_sigmac(vert);
        bilateralfilter(vert);
        vert = sample_->vert;
    }
}

void Bilateral::bilateralfilter(vector<CVertex> &vertices)
{
    double sigma_c, sigma_s;

    CMesh::VertexIterator vi;
    for(vi = sample_->vert.begin(); vi != sample_->vert.end(); ++vi)
    {

        // calculate sigma_c
        sigma_c = vi->sigma_c ;

        vcg::Point3f ni = vi->cN();
        vcg::Point3f pi = vi->P();

        double average_off_set = 0;
        std::vector<double> off_set_dis;
        off_set_dis.clear();
        for(int j = 0; j < knn; j++)
        {
            vcg::Point3f pj = vertices[vi->neighbors[j]].P();
            double t = (pj - pi) * ni;
            t = sqrt(t*t);
            average_off_set += t;
            off_set_dis.push_back(t);
        }
        //if(isnan(average_off_set)){std::cout<<i<<endl;break;}
        average_off_set = average_off_set / knn;
        double offset = 0;
        for(int j = 0; j < (int)off_set_dis.size(); j++)
            offset += (off_set_dis[j] - average_off_set) * (off_set_dis[j] - average_off_set);
        offset /= (double)off_set_dis.size();

        sigma_s = (sqrt(offset) < 1.0e-6) ? (sqrt(offset) + 1.0e-6) : sqrt(offset);//避免sigma_s等于0

        double sum = 0,normalizer = 0;
        for(int j = 0; j < knn; j++)
        {
            //vcg::Point3f nj = particles->vert[vert[i].neighbors[0][iv]].cN();
            vcg::Point3f pj = vertices[vi->neighbors[j]].P();
            double t = (pi - pj).Norm();
            double h = (pj - pi) * ni;
            double wc = std::exp(-0.25*t*t/(sigma_c *sigma_c));
            double ws = std::exp(-0.25*h*h/(sigma_s *sigma_s));
            sum += wc * ws * h;
            normalizer += wc * ws;
            //std::cout<<t<<" "<<h<<" "<<wc<<" "<<ws<<endl;
        }
        if(normalizer>1.0e-3)
        {
            vi->P() += ni * (sum / normalizer);
        }
    //std::cout<<vert[i].P()[0]<<" "<<vert[i].P()[1]<<" "<<vert[i].P()[2]<<endl;
    }
}

void Bilateral::calculate_sigmac(vector<CVertex> &vertices)
{
    GlobalFun::computeAnnNeigbhors(sample_->vert, sample_->vert, knn, false, "bilateral");
    CMesh::VertexIterator vi;
    for(vi = sample_->vert.begin(); vi != sample_->vert.end(); ++vi)
    {
        double sigmac = 0;
        for(int j = 0; j < knn; j++)
        {
            Point3f nearest_points = vertices[vi->neighbors[j]].P();
            sigmac += (vi->P()-nearest_points).Norm();
        }
        vi->sigma_c = (double)sigmac/(knn-1);
    }
}
