#include "RawtoRoot.h"
#include "Global.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <algorithm>
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
#include "TRandom.h"
#include <TSpectrum.h>
#include "TLegend.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGaxis.h"
// #include </afs/ihep.ac.cn/users/s/shiyk/YukunToolBox/Root.h>
//#include "Event.h"
int layer,chip,channel;
using namespace std;
char char_tmp[200];
Int_t main(int argc,char *argv[])
{
    double start = clock();
    raw2Root tw;
    tw.neEnergyCalib(argv[1],argv[2],argv[3],argv[4]);
    double end = clock();
    cout<<"end of RawToRoot : Time : "<<(end-start)/CLOCKS_PER_SEC<<endl;
    return 0;
}
int raw2Root::neEnergyCalib(string str_dat,string str_MIP,string output_file,string mode){
    //string str_root=find_datname(str_in);
    //string str_out=outputDir+"/"+"cos_ana.root";
    string str_out=output_file;
    TFile *fin,*fout;
    TTree *tree_in,*tree_out;
    double ped_time[Layer_No][chip_No][channel_No];
    double ped_charge[Layer_No][chip_No][channel_No];
    double gain_ratio[Layer_No][chip_No][channel_No];
    double gain_plat[Layer_No][chip_No][channel_No];
    double MIP[Layer_No][chip_No][channel_No];
    double hitE,hitE_layer[Layer_No],hitNo_layer[Layer_No];
    int nphoton=0;
    int totaln=0;
    double Edep=0;
    double MIP_E=0.461;//MeV
    const double ref_ped_time=390;
    const double ref_ped_charge=384;
    const double ref_MIP=900;
    const double ref_gain_ratio=26;
    int Select_EventNo=0;
    int HitNo=0;
    int CellID=0;
    float slope=0;
    double MPV=0;
    float plat=0;
    double pedestal_time=0;
    double pedestal_charge=0;
    double Layer_E[Layer_No]={0};
    TH1D *h_Edep;
    TH2D *h2_HitE_Layer;
    TH2D *h2_HitMap;
    sprintf(char_tmp,"Energy deposition");
    h_Edep = new TH1D(char_tmp,char_tmp,1000,0,3);
    sprintf(char_tmp,"Layer-HitE");
    h2_HitE_Layer = new TH2D(char_tmp,char_tmp,100,0,100,40,0,40);
    sprintf(char_tmp,"Hit Map");
    h2_HitMap = new TH2D(char_tmp,char_tmp,18,-HBU_X*3/2.,HBU_X*3/2.,18,-HBU_Y/2.,HBU_Y/2.);
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
        MIP[layer][chip][channel]=MPV;
        //cout<<layer<<" "<<chip<<" "<<channel<<" "<<MPV<<endl;
        if(MIP[layer][chip][channel]<200)cout<<"abnormal MIP "<<layer<<" "<<chip<<" "<<channel<<" "<<MPV<<endl;
    }
    tree_in->Delete();
    fin->Close();
    //read dat
    fin = TFile::Open(str_dat.c_str(),"READ");
    if (!fin){
        cout<<"cant open "<<str_dat<<endl;
        return 0;
    }
    cout<<"Read TTree "<<endl;
    tree_in = (TTree*)fin->Get("EventTree");
    if(mode=="Digi"){
        ReadCalibTree(tree_in);
    }
    if(mode=="MC"){
        ReadMCTree(tree_in);
    }
    cout<<"Read TTree Over"<<endl;
    for (int i_layer = 0; i_layer < Layer_No; ++i_layer){
        for (int i_chip = 0; i_chip < chip_No; ++i_chip){
            for (int i_chan = 0; i_chan < channel_No; ++i_chan){
                if(MIP[i_layer][i_chip][i_chan]<100) MIP[i_layer][i_chip][i_chan]=ref_MIP - ref_ped_time;
                //cout<<i_layer<<" "<<i_chip<<" "<<i_chan<<endl;
                // cout<<MIP[i_layer][i_chip][i_chan]<<endl;
            }
        }
    }
    cout<<"Create outfile"<<endl;
    fout = TFile::Open(str_out.c_str(),"RECREATE");
    if (!fout){
        cout<<"cant open "<<str_out<<endl;
        return 0;
    }
    tree_out = new TTree("EventTree","Hits afer n calibration");
    int t_Run_No=0,t_Event_Time=0,t_Event_No=0,t_Detector_ID=0;
    double t_Digi_Energy=0;
    vector<int> *t_cellID=0;
    vector<int> *t_cherenkov=0;
    vector<double> *t_Hit_E=0;
    vector<double> *t_Hit_X=0;
    vector<double> *t_Hit_Y=0;
    vector<double> *t_Hit_Z=0;
    vector<double> *t_Hit_Time=0;

    // tree_out ->Branch("Run_Num",&t_Run_No);
    // tree_out ->Branch("Event_Time",&t_Event_Time);
    tree_out ->Branch("Event_Num",&t_Event_No);
    // tree_out ->Branch("DetectorID",&t_Detector_ID);
    tree_out ->Branch("CellID",&t_cellID);
    tree_out ->Branch("Digi_Energy_HCAL",&t_Digi_Energy);
    tree_out ->Branch("Digi_Hit_Energy",&t_Hit_E);
    tree_out ->Branch("Hit_X",&t_Hit_X);
    tree_out ->Branch("Hit_Y",&t_Hit_Y);
    tree_out ->Branch("Hit_Z",&t_Hit_Z);    
    tree_out ->Branch("Hit_Time",&t_Hit_Time);
    // tree_out ->Branch("Cherenkov",&t_cherenkov);
    
    cout<<"Create outfile done"<<endl;

    for (int i = 0; i < tree_in->GetEntries(); ++i){
        totaln=0;
        t_Event_No=0;t_Digi_Energy=0;
        t_cellID->clear();t_Hit_E->clear();t_Hit_X->clear();t_Hit_Y->clear();t_Hit_Z->clear();t_Hit_Time->clear();
        if((i%1000)==0)cout<<i<<" out of "<<tree_in->GetEntries()<<endl;
        tree_in->GetEntry(i);
        // t_Run_No=_Run_No;
        // t_Event_Time=_Event_Time;
        t_Event_No=_Event_No;
        // t_Detector_ID=_Detector_ID;
        t_cellID=cellID;
        t_Hit_X=Hit_X;
        t_Hit_Y=Hit_Y;
        t_Hit_Z=Hit_Z;
        t_Hit_Time=Hit_Time;
        // t_cherenkov=cherenkov;
        for(int i_hit=0;i_hit<cellID->size();i_hit++){
            // decode_cellid(cellID->at(i_hit),layer,chip,channel);
            layer=Hit_Z->at(i_hit)/30;
            inverse(Hit_X->at(i_hit),Hit_Y->at(i_hit),chip,channel);
            nphoton=Hit_E->at(i_hit)/MIP_E * MIP[layer][chip][channel]/30;
            cout<<Hit_E->at(i_hit)<<"  "<<layer<<" "<<chip<<" "<<channel<<"  "<<MIP[layer][chip][channel]<<"  "<<nphoton<<endl;
            t_Hit_E->push_back(nphoton);
            totaln+=nphoton;
        }
        t_Digi_Energy=totaln;
        tree_out->Fill();
    }
    fout->cd();
    tree_out->Write();
    fout->Close();
    return 1;
}
// void raw2Root::ReadMCTree(TTree *tree){
//     Reset();
//     // tree ->SetBranchAddress("Run_Num",&_Run_No);
//     // tree ->SetBranchAddress("Event_Time",&_Event_Time);
//     tree ->SetBranchAddress("EventNum",&_Event_No);
//     // tree ->SetBranchAddress("DetectorID",&_Detector_ID);
//     tree ->SetBranchAddress("CellID",&cellID);
//     tree ->SetBranchAddress("Hit_Energy",&Hit_E);
//     tree ->SetBranchAddress("Hit_X",&Hit_X);
//     tree ->SetBranchAddress("Hit_Y",&Hit_Y);
//     tree ->SetBranchAddress("Hit_Z",&Hit_Z);    
//     // tree ->SetBranchAddress("Digi_Energy_HCAL",&Digi_Energy);
//     tree ->SetBranchAddress("Hit_Time",&Hit_Time);
//     // tree ->SetBranchAddress("Cherenkov",&cherenkov);
// }
