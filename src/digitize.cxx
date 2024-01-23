#include "RawtoRoot.h"
#include "Global.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <string>
#include <TF1.h>
#include <TH1.h>
#include <TStyle.h>
#include "TMath.h"
#include "TRandom3.h"
#include <TSpectrum.h>
#include "TLegend.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGaxis.h"
#include </afs/ihep.ac.cn/users/s/shiyk/YukunToolBox/Root.h>
//#include "Event.h"
int layer,chip,channel;
using namespace std;
char char_tmp[200];
double ped_high[Layer_No][chip_No][channel_No];
double ped_low[Layer_No][chip_No][channel_No];
double rms_high[Layer_No][chip_No][channel_No];
double rms_low[Layer_No][chip_No][channel_No];
double gain_ratio[Layer_No][chip_No][channel_No];
double gain_plat[Layer_No][chip_No][channel_No];
double lowgain_plat[Layer_No][chip_No][channel_No];
double _MIP[Layer_No][chip_No][channel_No];
double SPE[Layer_No][chip_No][channel_No];
double chi2_ndf[Layer_No][chip_No][channel_No];
double NDF[Layer_No][chip_No][channel_No];
double Sigma[Layer_No][chip_No][channel_No];
double MIP_E=0.461;//MeV
int pixel=7284;
double pde=0.32;
Int_t main(int argc,char *argv[])
{
    double start = clock();
    raw2Root tw;
    tw.Digitize(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
    double end = clock();
    cout<<"end of RawToRoot : Time : "<<(end-start)/CLOCKS_PER_SEC<<endl;
    return 0;
}
int raw2Root::Digitize(string str_dat,string str_ped,string str_dac,string str_MIP,string str_SPE,string output_file){
    //string str_root=find_datname(str_in);
    //string str_out=outputDir+"/"+"cos_ana.root";
    string str_out=output_file;
    TFile *fin,*fout;
    TTree *tree_in,*tree_out;
    double hitE,hitE_layer[Layer_No],hitNo_layer[Layer_No];
    double Edep=0;
    double SwitchPoint=500;
    const double ref_ped_high=374;
    const double ref_rms_high=3.5;
    const double ref_ped_low=371;
    const double ref_rms_low=2.3;
    const double ref_MIP=344.3;
    const double ref_gain_ratio=26;
    const double ref_gain_plat=2927;
    const double ref_spe=27.7;
    const double ref_sigma=7;
    int Select_EventNo=0;
    int HitNo=0;
    int CellID=0;
    double MPV=0;
    float slope=0;
    float plat=0;
    float low_plat=0;
    double pedestal_high=0,sigma_high=0;
    double pedestal_low=0,sigma_low=0;
    double spe=0,sigma=0,chi2=0,ndf=0;
    double Layer_E[Layer_No]={0};
    TH1D *h_Edep;
    TH2D *h2_HitE_Layer;
    TH2D *h2_HitMap;
    sprintf(char_tmp,"Energy deposition");
    // h_Edep = new TH1D(char_tmp,char_tmp,1000,0,3);
    // sprintf(char_tmp,"Layer-HitE");
    // h2_HitE_Layer = new TH2D(char_tmp,char_tmp,100,0,100,40,0,40);
    sprintf(char_tmp,"Hit Map");
    h2_HitMap = new TH2D(char_tmp,char_tmp,18,-360,360,18,-360.,360.);
    //read ped
    fin = TFile::Open(str_ped.c_str(),"READ");
    if (!fin){
        cout<<"cant open "<<str_ped<<endl;
        return 0;
    }
    tree_in = (TTree*)fin->Get("pedestal");
    if(!tree_in){
        cout<<"cant get tree pedestal"<<endl;
        return 0;
    }
    tree_in->SetBranchAddress("cellid",&CellID);
    tree_in->SetBranchAddress("highgain_peak",&pedestal_high);
    tree_in->SetBranchAddress("lowgain_peak",&pedestal_low);
    tree_in->SetBranchAddress("highgain_rms",&sigma_high);
    tree_in->SetBranchAddress("lowgain_rms",&sigma_low);
    for (int i_layer = 0; i_layer < Layer_No; ++i_layer){
        for (int i_chip = 0; i_chip < chip_No; ++i_chip){
            for (int i_chan = 0; i_chan < channel_No; ++i_chan){
                ped_high[i_layer][i_chip][i_chan]=-1;
                ped_low[i_layer][i_chip][i_chan]=-1;
                rms_high[i_layer][i_chip][i_chan]=-1;
                rms_low[i_layer][i_chip][i_chan]=-1;
                _MIP[i_layer][i_chip][i_chan]=-1;
                gain_ratio[i_layer][i_chip][i_chan]=-1;
                SPE[i_layer][i_chip][i_chan]=-1;
                Sigma[i_layer][i_chip][i_chan]=-1;
            }
        }
    }
    for (int i = 0; i < tree_in->GetEntries(); ++i){
        tree_in->GetEntry(i);
        decode_cellid(CellID,layer,chip,channel);
        //cout<<layer<<" "<<chip<<" "<<channel<<" "<<pedestal_time<<" "<<pedestal_charge<<endl;
        ped_high[layer][chip][channel]=pedestal_high;
        rms_high[layer][chip][channel]=sigma_high;
        ped_low[layer][chip][channel]=pedestal_low;
        rms_low[layer][chip][channel]=sigma_low;
    }
    tree_in->Delete();
    fin->Close();
    //read dac
    fin = TFile::Open(str_dac.c_str(),"READ");
    if (!fin){
        cout<<"cant open "<<str_dac<<endl;
        return 0;
    }
    //tree_in = (TTree*)fin->Get("calib");
    tree_in = (TTree*)fin->Get("dac");
    if(!tree_in){
        cout<<"cant get tree dac"<<endl;
        return 0;
    }
    tree_in->SetBranchAddress("cellid",&CellID);
    tree_in->SetBranchAddress("slope",&slope);
    tree_in->SetBranchAddress("plat",&plat);
    tree_in->SetBranchAddress("lowgain_satu_point",&low_plat);
    for (int i = 0; i < tree_in->GetEntries(); ++i){
        tree_in->GetEntry(i);
        decode_cellid(CellID,layer,chip,channel);
        //cout<<layer<<" "<<chip<<" "<<channel<<" "<<slope<<" "<<endl;
        gain_ratio[layer][chip][channel]=slope;
        gain_plat[layer][chip][channel]=plat;
        lowgain_plat[layer][chip][channel]=low_plat;
        if(gain_ratio[layer][chip][channel]<10 ||gain_ratio[layer][chip][channel]>50)cout<<CellID<<" abnormal gain ratio "<<layer<<" "<<chip<<" "<<channel<<" "<<slope<<endl;
    }
    tree_in->Delete();
    fin->Close();
    //read MIP
    fin = TFile::Open(str_MIP.c_str(),"READ");
    if (!fin){
        cout<<"cant open "<<str_MIP<<endl;
        return 0;
    }
    tree_in = (TTree*)fin->Get("MIP_Calibration");
    if(!tree_in){
        cout<<"cant get MIP tree"<<endl;
        return 0;
    }
    tree_in->SetBranchAddress("CellID",&CellID);
    tree_in->SetBranchAddress("MPV",&MPV);
    //tree_in->SetBranchAddress("ID",&CellID);
    //tree_in->SetBranchAddress("mpv",&MPV);
    for (int i = 0; i < tree_in->GetEntries(); ++i){
        tree_in->GetEntry(i);
        decode_cellid(CellID,layer,chip,channel);
        //MIP[layer][chip][channel]=MPV - ped_time[layer][chip][channel];
        _MIP[layer][chip][channel]=MPV;
        //cout<<layer<<" "<<chip<<" "<<channel<<" "<<MPV<<endl;
        // if(_MIP[layer][chip][channel]<200)cout<<"abnormal MIP "<<layer<<" "<<chip<<" "<<channel<<" "<<MPV<<endl;
    }
    tree_in->Delete();
    fin->Close();
    //read SPE
    fin = TFile::Open(str_SPE.c_str(),"READ");
    if (!fin){
        cout<<"cant open "<<str_MIP<<endl;
        return 0;
    }
    tree_in = (TTree*)fin->Get("spe");
    if(!tree_in){
        cout<<"cant get sipm gain tree"<<endl;
        return 0;
    }
    tree_in->SetBranchAddress("cellid",&CellID);
    tree_in->SetBranchAddress("spe",&spe);
    tree_in->SetBranchAddress("sigma",&sigma);
    tree_in->SetBranchAddress("chi2",&chi2);
    tree_in->SetBranchAddress("ndf",&ndf);
    for (int i = 0; i < tree_in->GetEntries(); ++i){
        tree_in->GetEntry(i);
        decode_cellid(CellID,layer,chip,channel);
        SPE[layer][chip][channel]=spe;
        chi2_ndf[layer][chip][channel]=chi2/ndf;
        NDF[layer][chip][channel]=ndf;
        Sigma[layer][chip][channel]=sigma;
    }
    tree_in->Delete();
    fin->Close();
    // smooth parameters
    for (int i_layer = 0; i_layer < Layer_No; ++i_layer){
        for (int i_chip = 0; i_chip < chip_No; ++i_chip){
            double mean_gain=0;
            int n_good=0;
            for (int i_chan = 0; i_chan < channel_No; ++i_chan){
                if(ped_high[i_layer][i_chip][i_chan]<0) ped_high[i_layer][i_chip][i_chan]=ref_ped_high;
                if(rms_high[i_layer][i_chip][i_chan]<0) rms_high[i_layer][i_chip][i_chan]=ref_rms_high;
                if(ped_low[i_layer][i_chip][i_chan]<0) ped_low[i_layer][i_chip][i_chan]=ref_ped_low;
                if(rms_low[i_layer][i_chip][i_chan]<0) rms_low[i_layer][i_chip][i_chan]=ref_rms_high;
                if(_MIP[i_layer][i_chip][i_chan]<100) _MIP[i_layer][i_chip][i_chan]=ref_MIP;
                if(gain_ratio[i_layer][i_chip][i_chan]<0) gain_ratio[i_layer][i_chip][i_chan]=ref_gain_ratio;
                if(lowgain_plat[i_layer][i_chip][i_chan]<2000) lowgain_plat[i_layer][i_chip][i_chan]=2000;
                if(gain_plat[i_layer][i_chip][i_chan]<0) gain_plat[i_layer][i_chip][i_chan]=ref_gain_plat;
                if(Sigma[i_layer][i_chip][i_chan]<0) Sigma[i_layer][i_chip][i_chan]=ref_sigma;
                if(chi2_ndf[layer][chip][channel]<20||NDF[layer][chip][channel]<=3||SPE[layer][chip][channel]<15){
                    mean_gain+=SPE[layer][chip][channel];
                    n_good++;
                }
            }
            mean_gain/=n_good;
            for (int i_chan = 0; i_chan < channel_No; ++i_chan){
                if(chi2_ndf[layer][chip][channel]>=20||SPE[layer][chip][channel]<mean_gain-5||SPE[layer][chip][channel]>mean_gain+5||NDF[layer][chip][channel]<=3||SPE[layer][chip][channel]<15){
                    SPE[layer][chip][channel]=mean_gain;
                }
            }
        }
    } 
    //read dat
    fin = TFile::Open(str_dat.c_str(),"READ");
    if (!fin){
        cout<<"cant open "<<str_dat<<endl;
        return 0;
    }
    cout<<"Read TTree "<<endl;
    tree_in = (TTree*)fin->Get("EventTree");
    ReadMCTree(tree_in);
    cout<<"Read TTree Over"<<endl;
    fout = TFile::Open(str_out.c_str(),"RECREATE");
    if (!fout){
        cout<<"cant open "<<str_out<<endl;
        return 0;
    }
    cout<<"digitizing "<<str_dat<<endl;
    tree_out = new TTree("Raw_Hit","Hits afer digitization");
    SetDataTree(tree_out);
    for (int i = 0; i < tree_in->GetEntries(); ++i){
        BranchClear();
        if((i%1000)==0)cout<<i<<" out of "<<tree_in->GetEntries()<<endl;
        // if(i>39000)cout<<i<<endl;
        tree_in->GetEntry(i);
        for(int i_hit=0;i_hit<cellID->size();i_hit++){
            // decode_cellid(cellID->at(i_hit),layer,chip,channel);
            hitTag->push_back(1);
            int cid=cellID->at(i_hit);
            decode_cellid(cid,layer,chip,channel);
            // layer=Hit_Z->at(i_hit)/30;
            // inverse(Hit_X->at(i_hit),Hit_Y->at(i_hit),chip,channel);
            // cout<<cid<<" "<<layer<<" "<<chip<<" "<<channel<<endl;
            data_cellid->push_back(cid);
            double HG=0,LG=0;
            digi(Hit_E->at(i_hit),HG,LG,i);
            HG_Charge->push_back(HG);
            LG_Charge->push_back(LG);
            h2_HitMap->Fill(Hit_X->at(i_hit)/40.3,Hit_Y->at(i_hit)/40.3);
            // Hit_Time->push
        }
        tree_out->Fill();
    }
    fout->cd();
    tree_out->Write();
    h2_HitMap->Write();
    // h_Edep->Write();
    fout->Close();
    return 1;
}
int raw2Root::digi(double energy,double &HG,double &LG,int i){
    int seed = 0;//time(0)%4294967296;
    TRandom3 *ran=new TRandom3(seed);
    // nonuniform of Scintillators
    energy=ran->Gaus(energy,energy*0.048);
    // energy to photon
    int n_photon=energy/MIP_E * _MIP[layer][chip][channel]/SPE[layer][chip][channel]/pde;
    n_photon=ran->Poisson(n_photon);
    // photon to phoelectron
    int n_phoelectron=ran->Binomial(n_photon,pde);
    // SiPM model
    // int n_fired=0;
    // bool isa[pixel];
    // std::fill(isa,isa+pixel,0);
    // for(int n=0;n<n_phoelectron;n++){
    //     int npx=ran->Integer(pixel);
    //     if(isa[npx])continue;
    //     isa[npx]=1;
    //     n_fired++;
    // }
    double n_fired_mean=pixel * (1.0 - TMath::Exp(-n_phoelectron * 1.0 / pixel));
    double n_fired_sigma=32.62 * (1.0 - TMath::Exp(-n_phoelectron / 3847));
    int n_fired = round(ran->Gaus(n_fired_mean , n_fired_sigma));

    // spe error
    double adc=n_fired*SPE[layer][chip][channel];
    double adc_sigma=sqrt(n_fired)*3;//Sigma[layer][chip][channel];
    // TODO
    adc=ran->Gaus(adc,adc_sigma);
    // ADC to HG LG
    HG=adc+ran->Gaus(ped_high[layer][chip][channel],rms_high[layer][chip][channel]);
    LG=adc/gain_ratio[layer][chip][channel]+ran->Gaus(ped_low[layer][chip][channel],rms_low[layer][chip][channel]);
    if(HG>gain_plat[layer][chip][channel]){
        HG=gain_plat[layer][chip][channel];
    }
    // if(LG>lowgain_plat[layer][chip][channel]){
    //     LG=lowgain_plat[layer][chip][channel];
    // }
    if(LG>2000){
        LG=2000;
    }
    HG=round(HG);
    LG=round(LG);
    delete ran;
    // cout<<n_phoelectron<<"  "<<n_fired<<"  "<<HG<<endl;
    return 1;
}