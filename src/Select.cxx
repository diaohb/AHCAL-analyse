#include "Select.h"
// #include "Global.h"
#include <iostream>
using namespace std;
Select::Select(vector<double>* _Hit_X,vector<double>* _Hit_Y,vector<double>* _Hit_Z, vector<double> *hit_energy,double a){
    // std::cout<<"Selection start"<<std::endl;
    // CellID=cellID;
    // copy(_Hit_X->begin(),_Hit_X->end(),Hit_X->begin());
    // copy(_Hit_Y->begin(),_Hit_Y->end(),Hit_Y->begin());
    // copy(_Hit_Z->begin(),_Hit_Z->end(),Hit_Z->begin());
    // copy(hit_energy->begin(),hit_energy->end(),Hit_energy->begin());
    Hit_X=_Hit_X;
    Hit_Y=_Hit_Y;
    Hit_Z=_Hit_Z;
    Hit_energy=hit_energy;
    beam_energy=a;
    Init();
}
Select::~Select(){
    // cout<<"Selection done"<<endl;
}

void Select::Init(){
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
        int layer=Hit_Z->at(i)/30;
        hitno++;
        layerhitno[layer]++;
        lenergy[layer]+=e;
    }

    // for(int i=0;i<Hit_X->size();i++){
    //     int layer=Hit_Z->at(i)/30;
    //     double e=Hit_energy->at(i);
    //     lenergy[layer]+=e;
    //     if(e>0.5*MIP_E){
    //         hitno++;
    //         layerhitno[layer]++;
    //     }
    // }
    for(int i=0;i<40;i++){
        tenergy+=lenergy[i];
        if(lenergy[i]>0){//0.3*MIP_E){
            ishit[i]=1;
            hitlayer++;
        }
    }
    for(int i=0;i<40;i++){
        if(layerhitno[i]>4){
            shower_start==i;
            break;
        }
    }
    for(int i=39;i>=0;i--){
        if(layerhitno[i]>4){
            shower_end==i;
            break;
        }
    }
}

void Select::EnergyCenter(){
    for(int i=0;i<Hit_X->size();i++){
        int layer=Hit_Z->at(i)/30;
        double e=Hit_energy->at(i);
        CenterX[layer]=Hit_X->at(i)*e;
        CenterY[layer]=Hit_Y->at(i)*e;
    }
    for(int i=0;i<40;i++){
        if(ishit[i]==0){
            // noenergy=i;
            noenergylayer.push_back(i);
            CenterX[i]=0;
            CenterY[i]=0;
            continue;
            // break;
        }
        CenterX[i]/=lenergy[i];
        CenterY[i]/=lenergy[i];
        // cout<<i<<"  "<<lenergy[i]<<" "<<CenterX[i]<<" "<<CenterY[i]<<endl;
    }
}
void Select::RMS_zeta(){
    for(int i=0;i<Hit_X->size();i++){
        int layer=Hit_Z->at(i)/30;
        double e=Hit_energy->at(i);
        rmsX[layer]+=(Hit_X->at(i)-CenterX[layer])*(Hit_X->at(i)-CenterX[layer])*e;
        rmsY[layer]+=(Hit_Y->at(i)-CenterY[layer])*(Hit_Y->at(i)-CenterY[layer])*e;
    }
    for(int i=0;i<40;i++){
        if(ishit[i]==0){
            rmsX[i]=0;
            rmsY[i]=0;
            continue;
            // break;
        }
        rmsX[i]/=lenergy[i];
        rmsY[i]/=lenergy[i];
        rms[i]=sqrt(rmsX[i]+rmsY[i]);
        t_rms+=rms[i];
    }
    // cout<<zeta<<"  ";
    zeta=lenergy[14]/tenergy;//*t_rms*t_rms*t_rms*t_rms/8E6;
    // cout<<lenergy[14]<<"  "<<lenergy[19]<<"  "<<tenergy<<"  "<<zeta<<endl;
}
void Select::MaxEnergy(){
    for(int i=0;i<Hit_X->size();i++){
        int layer=Hit_Z->at(i)/30;
        if(maxenergy[layer]<Hit_energy->at(i)){
            maxenergy[layer]=Hit_energy->at(i);
            maxenergyX[layer]=Hit_X->at(i);
            maxenergyY[layer]=Hit_Y->at(i);
        }
    }
    // cout<<maxenergyX[0]<<"  "<<maxenergyY[0]<<endl;
}

int Select::Result(double &total_energy,vector<double>* &layer_energy,int &_hitlayer,int &_hitno){
    EnergyCenter();
    RMS_zeta();
    // MaxEnergy();
    for(int i=0;i<20;i++){
        if( ishit[i]&&(CenterX[i] > 280 || CenterX[i] < -280 || CenterY[i] > 280 || CenterY[i] < -280)){
            return 7;
        }
    }
    _hitlayer=hitlayer;
    _hitno=hitno;
    total_energy=0;
    layer_energy->clear();
    for(int i=0;i<40;i++){
        layer_energy->push_back(0);
    }
    for(int i=0;i<Hit_X->size();i++){
        int layer=Hit_Z->at(i)/30;
        // if(layer>=noenergy)continue;
        double x=Hit_X->at(i)-CenterX[layer];
        double y=Hit_Y->at(i)-CenterY[layer];
        if(x<100&&x>-100&&y<100&&y>-100){
            layer_energy->at(layer)+=Hit_energy->at(i);
        }
    }
    for(int i=0;i<40;i++){
        total_energy+=layer_energy->at(i);
    }
    if(tenergy<0.5*MIP_E){
        return 0;
    }
    if(hitno<fhitnocut->Eval(beam_energy)){//hitnocut[beam_energy]){
        return 1;
    }
    if(hitlayer<=10){//hitlayercut[beam_energy]){
        return 2;
    }
    if(shower_start>=5){
        return 3;
    }
    if(Ismip()){
        return 4;
    }
    if(!Isstraight()){
        return 5;
    }
    return 6;
}

bool Select::Isstraight(){
    double maxx=-360,minx=360,maxy=-360,miny=360;
    for(int i=0;i<40;i++){
        if(ishit[i]&&lenergy[i]/tenergy>0.05){
            // if(maxenergyX[i]>maxx)maxx=maxenergyX[i];
            // if(maxenergyX[i]<minx)minx=maxenergyX[i];
            // if(maxenergyY[i]>maxx)maxy=maxenergyY[i];
            // if(maxenergyY[i]<minx)miny=maxenergyY[i];
            if(CenterX[i]>maxx)maxx=CenterX[i];
            if(CenterX[i]<minx)minx=CenterX[i];
            if(CenterY[i]>maxx)maxy=CenterY[i];
            if(CenterY[i]<minx)miny=CenterY[i];
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
    for(int i=35;i<40;i++){
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

void Select::Layerd(){
    
}