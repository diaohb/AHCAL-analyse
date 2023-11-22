#include "HBase.h"
#include <TGraph.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include "Select.h"
using namespace std;
//analyse list digi/mc select gaus/cb
//default digi yes gaus
int main(int argc,char* argv[]){
    HBase *hbase=new HBase();
    hbase->ReadList(argv[1]);
    TString s=argv[1];
    string mode="digi";
    bool isselect=true;
    string func="gaus";
    for(int i=2;i<argc;i++){
        string s=argv[i];
        if(s=="mc"){mode="mc";continue;}
        if(s=="no"){isselect=false;continue;}
        if(s=="cb"){func="cb";}
    }
    cout<<"\e[44mmode: "<<mode<<" , selection: "<<isselect<<" , function: "<<func<<"\e[0m"<<endl;
    TFile *fout=TFile::Open(s+"_analyse.root","RECREATE");
    TTree *tzprofile=new TTree("tzprofile","tzprofile");
    TTree *tfitEn=new TTree("tfitEn","tfitEn");
    double a;
    double total_energy=0;
    int event_no=0;
    int hitlayer=0;
    int hitno=0;
    double zeta=0;
    double rms=0;
    double f_15=0; //15th, index=14
    double f_20=0;
    int select_flag_l=0;
    vector<double> *layer_energy=0;
    vector<double> *layer_energy_origin=0;
    vector<double> *layer_energy_threshold=0;
    vector<double> *layer_energy_hitno=0;
    vector<double> *layer_energy_hitlayer=0;
    vector<double> *layer_energy_showerstart=0;
    vector<double> *layer_energy_lastratio=0;
    vector<double> *layer_energy_straight=0;
    vector<int> *layer=0;
    
    tzprofile->Branch("beam_energy",&a);
    tzprofile->Branch("total_energy",&total_energy);
    tzprofile->Branch("layer",&layer);
    tzprofile->Branch("event_no",&event_no);
    tzprofile->Branch("hitlayer",&hitlayer);
    tzprofile->Branch("hitno",&hitno);
    tzprofile->Branch("zeta",&zeta);
    tzprofile->Branch("rms",&rms);
    tzprofile->Branch("f_15",&f_15);
    tzprofile->Branch("f_20",&f_20);
    tzprofile->Branch("layer_energy",&layer_energy);
    // tzprofile->Branch("layer_energy_origin",&layer_energy_origin);
    // tzprofile->Branch("layer_energy_threshold",&layer_energy_threshold);
    // tzprofile->Branch("layer_energy_hitno",&layer_energy_hitno);
    // tzprofile->Branch("layer_energy_hitlayer",&layer_energy_hitlayer);
    // tzprofile->Branch("layer_energy_showerstart",&layer_energy_showerstart);
    // tzprofile->Branch("layer_energy_lastratio",&layer_energy_lastratio);
    // tzprofile->Branch("layer_energy_straight",&layer_energy_straight);
    tzprofile->Branch("select_flag",&select_flag_l);

    double mean=0;
    double sigma=0;
    int select_flag=0;
    tfitEn->Branch("beam_energy",&a);
    tfitEn->Branch("mean",&mean);
    tfitEn->Branch("sigma",&sigma);
    tfitEn->Branch("select_flag",&select_flag);
    TF1* f;
    if(func=="gaus"){
        f=new TF1("f","gaus",0,10000);
    }
    else if(func=="cb"){
        f=new TF1("f","crystalball",0,10000);
    }
    vector<double> be;
    unordered_map<double ,TH1D*> h_origin;
    unordered_map<double ,TH1D*> h_threshold;
    unordered_map<double ,TH1D*> h_hitno;
    unordered_map<double ,TH1D*> h_hitlayer;
    unordered_map<double ,TH1D*> h_showerstart;
    unordered_map<double ,TH1D*> h_lastratio;
    unordered_map<double ,TH1D*> h_straight;
    for(int l=0;l<40;l++){
        layer->push_back(l);
    }
    for_each(hbase->list.begin(),hbase->list.end(),[&](string tmp){
        string beamenergy = tmp;
		// beamenergy = beamenergy.substr(0,beamenergy.find_last_of('/'));
		beamenergy = beamenergy.substr(0,beamenergy.find("eV")-1);
		beamenergy = beamenergy.substr(beamenergy.rfind('/')+1);
        beamenergy = beamenergy.substr(beamenergy.rfind('_')+1);
        a=stoi(beamenergy);
        if(a==500)a=0.5;
        vector<double>::iterator it=find(be.begin(),be.end(),a);
        if(it==be.end()){
            be.push_back(a);
            TString hname="energy_"+beamenergy+"GeV";
            if(a==0.5)hname="energy_500MeV";
            h_origin[a]=new TH1D(hname+"_origin",hname+"_origin",10000,0,10000);
            h_threshold[a]=new TH1D(hname+"_threshold",hname+"_threshold",10000,0,10000);
            h_hitno[a]=new TH1D(hname+"_hitno",hname+"_hitno",10000,0,10000);
            h_hitlayer[a]=new TH1D(hname+"_hitlayer",hname+"_hitlayer",10000,0,10000);
            h_showerstart[a]=new TH1D(hname+"_showerstart",hname+"_showerstart",10000,0,10000);
            h_lastratio[a]=new TH1D(hname+"_lastratio",hname+"_lastratio",10000,0,10000);
            h_straight[a]=new TH1D(hname+"_straight",hname+"_straight",10000,0,10000);
        }
		hbase->ReadTree(tmp,"EventTree",mode);
        int flag[8]={0};
        for(int ientry=0;ientry<hbase->tin->GetEntries();ientry++){
            hbase->tin->GetEntry(ientry);
            event_no=ientry;
            double totalenergy=0;
            // layer->clear();
            layer_energy->clear();
            // layer_energy_origin->clear();
            // layer_energy_threshold->clear();
            // layer_energy_hitno->clear();
            // layer_energy_hitlayer->clear();
            // layer_energy_showerstart->clear();
            // layer_energy_lastratio->clear();
            // layer_energy_straight->clear();
            // vector<double> vtmp={1};
            // vector<double>* layer_energy=&vtmp;
            if(isselect==false){
		        double l_e[40]={0};
                for(int icell=0;icell<hbase->_cellID->size();icell++){
                    int lay=hbase->_Hit_Z->at(icell)/30;
                    double e=hbase->_Digi_Hit_Energy->at(icell);
                    l_e[lay]+=e;
                    totalenergy+=e;
                    // std::cout<<l_e[lay]<<" "<<lay<<std::endl;
                }
                // cout<<totalenergy<<endl;
                for(int l=0;l<40;l++){
                    layer_energy_origin->push_back(l_e[l]);
                }
                h_origin[a]->Fill(totalenergy);
            }
            else{
                Select select(hbase->_Hit_X,hbase->_Hit_Y,hbase->_Hit_Z,hbase->_Digi_Hit_Energy,a);
                int result=select.Result(total_energy,layer_energy,hitlayer,hitno);
                zeta=select.GetZeta();
                rms=select.GetRMS();
                flag[result]++;
                select_flag_l=result;
                if(result==7){
                    continue;
                }
                f_15=layer_energy->at(14)/total_energy;
                f_20=layer_energy->at(19)/total_energy;
                h_origin[a]->Fill(total_energy);
                if(result>=1){
                    h_threshold[a]->Fill(total_energy);
                }
                if(result>=2){
                    h_hitno[a]->Fill(total_energy);
                }
                if(result>=3){
                    h_hitlayer[a]->Fill(total_energy);
                }
                if(result>=4){
                    h_showerstart[a]->Fill(total_energy);
                }
                if(result>=5){
                    h_lastratio[a]->Fill(total_energy);
                }
                if(result==6){
                    h_straight[a]->Fill(total_energy);
                }
            }
            tzprofile->Fill();         
        }
        cout<<"----------------------------------------------\e[45m"<<a<<"GeV\e[0m----------------------------------------------------"<<endl;
        if(isselect==true)cout<<"\e[44mcorner: "<<flag[7]<<"  < 0.5MIP: "<<flag[0]<<"  hitno: "<<flag[1]<<"  hitlayer: "<<flag[2]<<"  shower start: "<<flag[3]<<"  lastratio: "<<flag[4]<<"  not straight: "<<flag[5]<<"  ramain: "<<flag[6]<<"  total eneries: "<<hbase->tin->GetEntries()<<"\e[0m"<<endl;
        cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
        hbase->fin->Close();
    });
    auto fsave=[&](unordered_map<double,TH1D*> h){for(int i=0;i<be.size();i++){
        a=be.at(i);
        TH1D *htmp=h[a];
        mean=htmp->GetMean();
        sigma=htmp->GetRMS();
        double fs=mean-2*sigma;
        // if(fs<0)fs==5
        if(mean>4*fs)fs=mean/4;
        f->SetParameter("Mean",mean);
        f->SetParameter("Sigma",sigma);
        if(func=="cb"){
            f->SetParameter("Alpha",1);
            f->SetParameter("N",1);
            double y=htmp->GetBinContent(htmp->FindBin(mean));
            f->SetParLimits(0,y-100,y*5);       //constant
            f->SetParLimits(1,mean-sigma,mean+3*sigma);//mean
            f->SetParLimits(2,mean/120,sigma);   //sigma
            f->SetParLimits(3,-2,2);   //alpha
            f->SetParLimits(4,0.5,100);  //N
            fs=mean/2;
        }
        if(a<=2&&fs<0){fs=0;}
        htmp->Fit(f,"q","",fs,mean+2*sigma);
        mean=f->GetParameter("Mean");
        sigma=f->GetParameter("Sigma");
        fs=mean-2*sigma;
        if(mean>4*fs)fs=mean/4;
        if(func=="cb"){
            fs=mean/2;
        }
        if(a<=2&&fs<0){fs=0;}
        htmp->Fit(f,"q","",fs,mean+2*sigma);
        mean=f->GetParameter("Mean");
        sigma=f->GetParameter("Sigma");
        fs=mean-2*sigma;
        if(mean>4*fs)fs=mean/4;
        if(func=="cb"){
            fs=mean/2;
        }
        if(a<=2&&fs<0){fs=0;}
        htmp->Fit(f,"q","",fs,mean+2*sigma);
        mean=f->GetParameter("Mean");
        sigma=f->GetParameter("Sigma");
        // if(func=="cb"){
            // double n=f->GetParameter("N");
            // double alpha=f->GetParameter("Alpha");
            // double constant=f->GetParameter("Constant");
            // if(alpha<1.177){
            //     // double fwhm=sigma+1;
            // }
        // }
        tfitEn->Fill();
        htmp->Write();
    }};
    fout->mkdir("origin");
    fout->cd("origin");
    select_flag=0;
    fsave(h_origin);
    if(isselect==true){
        fout->mkdir("threshold");
        fout->cd("threshold");
        select_flag=1;
        fsave(h_threshold);
        fout->mkdir("hitno");
        fout->cd("hitno");
        select_flag=2;
        fsave(h_hitno);
        fout->mkdir("hitlayer");
        fout->cd("hitlayer");
        select_flag=3;
        fsave(h_hitlayer);
        fout->mkdir("showerstart");
        fout->cd("showerstart");
        select_flag=4;
        fsave(h_showerstart);
        fout->mkdir("lastratio");
        fout->cd("lastratio");
        select_flag=5;
        fsave(h_lastratio);
        fout->mkdir("straight");
        fout->cd("straight");
        select_flag=6;
        fsave(h_straight);
    }
    fout->cd();
    tzprofile->Write();
    tfitEn->Write();
    // fout->Write();
}