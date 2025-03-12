#include "PID_analyse.h"
#include "Global.h"
#include <iostream>
#include <TLinearFitter.h>
using namespace std;
PID_analyse::PID_analyse(vector<double>* _Hit_X,vector<double>* _Hit_Y,vector<double>* _Hit_Z, vector<double> *hit_energy,double a){
    // std::cout<<"Selection start"<<std::endl;
    Hit_X->assign(_Hit_X->begin(), _Hit_X->end());
    Hit_Y->assign(_Hit_Y->begin(), _Hit_Y->end());
    Hit_Z->assign(_Hit_Z->begin(), _Hit_Z->end());
    Hit_energy->assign(hit_energy->begin(), hit_energy->end());
    beam_energy = a;
    Init();
}
PID_analyse::~PID_analyse(){
    // cout<<"Selection done"<<endl;
}

void PID_analyse::Init()
{
    // cout << Hit_X->size() << endl;
    for(int i=0;i<Hit_X->size();i++){
        double e=Hit_energy->at(i);
        if(e<0.5*MIP_E){
            Hit_X->erase(Hit_X->begin()+i);
            Hit_Y->erase(Hit_Y->begin()+i);
            Hit_Z->erase(Hit_Z->begin()+i);
            Hit_energy->erase(Hit_energy->begin()+i);
            i--;
            continue;
        }
    }
    EnergyCenter();
    RMS_zeta();
    // MaxEnergy();
    E_Hit();
    FD();
}

void PID_analyse::EnergyCenter(){
    for (int i = 0; i < Hit_X->size(); i++)
    {
        int layer = Hit_Z->at(i) / 30;
        double e = Hit_energy->at(i);
        CenterX->at(layer) += Hit_X->at(i) * e;
        CenterY->at(layer) += Hit_Y->at(i) * e;
        hitno++;
        layerhitno[layer]++;
        lenergy[layer] += e;
    }
    for (int i = 0; i < 40; i++)
    {
        if (lenergy[i] == 0)
        {
            CenterX->at(i) = -500;
            CenterY->at(i) = -500;
            continue;
        }
        CenterX->at(i) /= lenergy[i];
        CenterY->at(i) /= lenergy[i];
        tenergy+=lenergy[i];
        if(lenergy[i]>0.5*MIP_E){
            ishit[i]=1;
            hitlayer++;
        }
    }
    shower_start=0;
    for(int i=0;i<40;i++)
    {
        // cout << i << "  " << ishit[i] << "  " << lenergy[i] << endl;
        if(layerhitno[i]>=4){
            shower_start=i;
            break;
        }
    }
    shower_end=0;
    for(int i=37;i>=0;i--){
        if(layerhitno[i]>=4){
            shower_end=i;
            break;
        }
    }
}
void PID_analyse::RMS_zeta(){
    for(int i=0;i<Hit_X->size();i++){
        int layer=Hit_Z->at(i)/30;
        double e=Hit_energy->at(i);
        rmsX[layer]+=(Hit_X->at(i)-CenterX->at(layer))*(Hit_X->at(i)-CenterX->at(layer))*e;
        rmsY[layer]+=(Hit_Y->at(i)-CenterY->at(layer))*(Hit_Y->at(i)-CenterY->at(layer))*e;
    }
    for(int i=0;i<40;i++){
        if(lenergy[i]==0){
            rmsX[i]=0;
            rmsY[i]=0;
            continue;
        }
        rmsX[i]/=lenergy[i];
        rmsY[i]/=lenergy[i];
        rms[i]=sqrt(rmsX[i]+rmsY[i]);
        t_rms+=rms[i];
    }
    // cout<<zeta<<"  ";
    zeta=lenergy[14]/tenergy*t_rms*t_rms*t_rms*t_rms;
    // cout<<lenergy[14]<<"  "<<lenergy[19]<<"  "<<tenergy<<"  "<<zeta<<endl;
}
void PID_analyse::MaxEnergy(){
    for(int i=0;i<Hit_X->size();i++){
        int layer=Hit_Z->at(i)/30;
        if(max_energy->at(layer) < Hit_energy->at(i)){
            max_energy->at(layer) = Hit_energy->at(i);
            max_x->at(layer) = Hit_X->at(i);
            max_y->at(layer) = Hit_Y->at(i);
        }
    }
}
void PID_analyse::E_Hit(){
    E_hit=0;
    for(int i=0;i<Hit_energy->size();i++){
        E_hit+=Hit_energy->at(i);
    }
    E_hit /= Hit_X->size();
}
void PID_analyse::FD(){
    // fitter->SetFormula("pol1");
    // fitter->StoreData(0);
    fd = 0;
    for (int i = 0; i < 9; i++)
    {
        alpha[i] = log(i+2);
        R[i] = log(R_alpha(i+2));
        fd+=R[i]/alpha[i];
    }
    alpha[9] = log(20);
    R[9] = log(R_alpha(20));
    fd+=R[9]/alpha[9];
    fd/=10;
    // fitter->AssignData(20, 1, alpha, R);
    // fitter->Eval();
    // fd = fitter->GetParameter(1);
}
double PID_analyse::R_alpha(int _alpha){
    double R_alpha=0;
    m_hit.clear();
    double x = 0, y = 0;
    if(_alpha%2==1){
        x = 20, y = 20;
    }
    for (int i = 0; i < Hit_X->size(); i++)
    {
        int a = floor((Hit_X->at(i) - x + 20 * (_alpha - 1)) / 40 / _alpha);
        int b = floor((Hit_Y->at(i) - x + 20 * (_alpha - 1)) / 40 / _alpha);
        m_hit[int(Hit_Z->at(i) / 30) * 10000 + a * 100 + b] = 1;
    }
    R_alpha=hitno/double(m_hit.size());
    return R_alpha;
}