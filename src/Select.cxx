#include "Select.h"
#include "Global.h"
#include <iostream>
using namespace std;
Select::Select(vector<double>* _Hit_X,vector<double>* _Hit_Y,vector<double>* _Hit_Z, vector<double> *hit_energy,double a){
    // std::cout<<"Selection start"<<std::endl;
    beam_energy = a;
    CellID = new vector<int>;
    Hit_X = new vector<double>;
    Hit_Y = new vector<double>;
    Hit_Z = new vector<double>;
    Hit_energy = new vector<double>;
    CenterX = new vector<double>(40, 0);
    CenterY = new vector<double>(40, 0);
    max_energy = new vector<double>(40, 0);
    max_x = new vector<double>(40, -500);
    max_y = new vector<double>(40, -500);

    Hit_X->assign(_Hit_X->begin(), _Hit_X->end());
    Hit_Y->assign(_Hit_Y->begin(), _Hit_Y->end());
    Hit_Z->assign(_Hit_Z->begin(), _Hit_Z->end());
    Hit_energy->assign(hit_energy->begin(), hit_energy->end());
    Init();
}
Select::~Select(){
    // cout<<"Selection done"<<endl;
}

void Select::Init()
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
            if(e<-99){
                big_hitno++;
            }
            continue;
        }
        int chip=0,channel=0;
        inverse(Hit_X->at(i),Hit_Y->at(i),chip,channel);
        int layer=Hit_Z->at(i)/30;
        int cellid = layer * 1e5 + chip * 1e4 + channel;
        CellID->push_back(cellid);
        true_hitno++;
        layerhitno[layer]++;
        lenergy[layer]+=e;
    }
    EnergyCenter();
    RMS_zeta();
    MaxEnergy();
    E_Hit();
    FD();
}

void Select::EnergyCenter(){
    vector<double> tmp_x(40,0);
    vector<double> tmp_y(40,0);
    for(int i=0;i<Hit_X->size();i++){
        int layer=Hit_Z->at(i)/30;
        double e = Hit_energy->at(i);
        tmp_x.at(layer)+=Hit_X->at(i)*e;
        tmp_y.at(layer)+=Hit_Y->at(i)*e;
    }
    for(int i=0;i<40;i++){
        if(lenergy[i]==0){
            tmp_x.at(i)=-500;
            tmp_y.at(i)=-500;
            continue;
        }
        tmp_x.at(i)/=lenergy[i];
        tmp_y.at(i)/=lenergy[i];
    }

    CenterX->resize(40,0);
    CenterY->resize(40,0);
    hitno = 0;
    for (int i = 0; i < 40;i++){
        layerhitno[i] = 0;
        lenergy[i] = 0;
    }
    for (int i = 0; i < Hit_X->size(); i++)
    {
        int layer = Hit_Z->at(i) / 30;
        double e = Hit_energy->at(i);
        double x = Hit_X->at(i) - tmp_x.at(layer);
        double y = Hit_Y->at(i) - tmp_y.at(layer);
        // if (x < 60 && x > -60 && y < 60 && y > -60)
        if (x < 100 && x > -100 && y < 100 && y > -100)
        {
            CenterX->at(layer) += Hit_X->at(i) * e;
            CenterY->at(layer) += Hit_Y->at(i) * e;
            int chip=0,channel=0;
            inverse(Hit_X->at(i),Hit_Y->at(i),chip,channel);
            int layer = Hit_Z->at(i) / 30;
            int cellid = layer * 1e5 + chip * 1e4 + channel;
            CellID->push_back(cellid);
            hitno++;
            layerhitno[layer]++;
            lenergy[layer] += e;
            // ishit[i] = 1;
        }
    }
    for (int i = 0; i < 40; i++)
    {
        tenergy += lenergy[i];
        if (lenergy[i] > 0)
        {
            ishit[i] = 1;
            hitlayer++;
            CenterX->at(i) /= lenergy[i];
            CenterY->at(i) /= lenergy[i];
        }
        else
        {
            CenterX->at(i) = -500;
            CenterY->at(i) = -500;
        }
    }
    for(int i=0;i<40;i++)
    {
        if(layerhitno[i]>=4){
            shower_start=i;
            break;
        }
    }
    for(int i=37;i>=0;i--){
        if(layerhitno[i]>=4){
            shower_end=i;
            break;
        }
    }
}
void Select::RMS_zeta(){
    for(int i=0;i<Hit_X->size();i++){
        int layer=Hit_Z->at(i)/30;
        double e=Hit_energy->at(i);
        rmsX[layer]+=(Hit_X->at(i)-CenterX->at(layer))*(Hit_X->at(i)-CenterX->at(layer))*e;
        rmsY[layer]+=(Hit_Y->at(i)-CenterY->at(layer))*(Hit_Y->at(i)-CenterY->at(layer))*e;
    }
    for(int i=0;i<40;i++){
        if(ishit[i]==0){
            rmsX[i]=0;
            rmsY[i]=0;
            continue;
            // break;
        }
        // cout<<i<<"  "<<ishit[i]<<"  "<<lenergy[i]<<endl;
        rmsX[i]/=lenergy[i];
        rmsY[i]/=lenergy[i];
        rms[i]=sqrt(rmsX[i]+rmsY[i]);
        t_rms+=rms[i];
    }
    // cout<<zeta<<"  ";
    zeta=lenergy[14]/tenergy*t_rms*t_rms*t_rms*t_rms;
    // cout<<lenergy[14]<<"  "<<lenergy[19]<<"  "<<tenergy<<"  "<<zeta<<endl;
}
void Select::MaxEnergy(){
    for(int i=0;i<Hit_X->size();i++){
        int layer=Hit_Z->at(i)/30;
        if(max_energy->at(layer) < Hit_energy->at(i)){
            max_energy->at(layer) = Hit_energy->at(i);
            max_x->at(layer) = Hit_X->at(i);
            max_y->at(layer) = Hit_Y->at(i);
        }
    }
    // cout<<maxenergyX[0]<<"  "<<maxenergyY[0]<<endl;
}

int Select::Result(double &total_energy,vector<double>* &layer_energy,vector<double>* &layer_hitno,int &_hitlayer,int &_hitno,int &_big_hitno)
{
    int flag_7 = 0;
    for (int i = 0; i < 20; i++)
    {
        if( ishit[i]&&(CenterX->at(i) > 120 || CenterX->at(i) < -120 || CenterY->at(i) > 120 || CenterY->at(i) < -120)){
            flag_7++;
        }
    }
    if(flag_7>2){
        return 7;
    }
    _hitlayer=hitlayer;
    _hitno=0;
    total_energy=0;
    _big_hitno = big_hitno;
    layer_energy->clear();
    layer_hitno->clear();
    for(int i=0;i<40;i++)
    {
        layer_energy->push_back(lenergy[i]);
        layer_hitno->push_back(layerhitno[i]);
    }
    for(int i=0;i<30;i++){
        _hitno += layer_hitno->at(i);
        total_energy += layer_energy->at(i);
    }
    if(total_energy<0.5*MIP_E){
        return 0;
    }
    if (hitno < hitnocut_down[beam_energy] || hitno > hitnocut_up[beam_energy])
    {
        return 1;
    }
    if (hitlayer < hitlayercut_down[beam_energy] ||hitlayer > hitlayercut_up[beam_energy])
    {
        return 2;
    }
    if(shower_start>=10){
		return 3;
    }
    if (shower_end>=25)
    {
        return 9;
    }
    if(Ismip()){
		return 4;
	}
	if(!Isstraight()){
		return 5;
	}
    if(layer_hitno->at(0)>1||layer_hitno->at(0)==0){
        return 8;
    }
    return 6;
}

bool Select::Isstraight(){
    double maxx=-360,minx=360,maxy=-360,miny=360;
    for(int i=0;i<30;i++){
        if(ishit[i]&&lenergy[i]/tenergy>0.05){
            // if(maxenergyX[i]>maxx)maxx=maxenergyX[i];
            // if(maxenergyX[i]<minx)minx=maxenergyX[i];
            // if(maxenergyY[i]>maxx)maxy=maxenergyY[i];
            // if(maxenergyY[i]<minx)miny=maxenergyY[i];
            if(CenterX->at(i)>maxx)maxx=CenterX->at(i);
            if(CenterX->at(i)<minx)minx=CenterX->at(i);
            if(CenterY->at(i)>maxx)maxy=CenterY->at(i);
            if(CenterY->at(i)<minx)miny=CenterY->at(i);
        }
    }
    // cout<<maxx<<"  "<<minx<<"  "<<maxy<<"  "<<miny<<endl;
    if((maxx-minx)<200&&(maxy-miny)<200){
        return true;
    }
    else{
        return false;
    }
}

bool Select::Ismip(){
    double lastenergy=0;
    for(int i=30;i<38;i++){
        lastenergy+=lenergy[i];
    }
    return (lastenergy/tenergy>0.05);
    // int vlayer[40]={0};
    // for(int i=0;i<Hit_Z->size();i++){
    //     int layer=Hit_Z->at(i)/30;
    //     if(vlayer[layer]>2){
    //         return true;
    //     }
    //     vlayer[layer]++;
    // }
    // return false;
}
bool Select::Ishadron(){
    for(int i=20;i<40;i++){
        if(layerhitno[i]>=3){
            return true;
        }
    }
    return false;
}
void Select::E_Hit(){
    E_hit=0;
    for(int i=0;i<Hit_energy->size();i++){
        E_hit+=Hit_energy->at(i);
    }
    E_hit /= Hit_X->size();
}
void Select::FD(){
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
double Select::R_alpha(int _alpha){
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
    R_alpha=true_hitno/double(m_hit.size());
    return R_alpha;
}
