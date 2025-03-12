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
int main(int argc,char* argv[]){
    HBase *hbase=new HBase();
    hbase->ReadList(argv[1]);
    TString s=argv[1];
    string mode;
    if(argc==3){
        mode=argv[2];
    }
    else{
        mode="digi";
    }
    TFile *fout=TFile::Open(s+"_analyse.root","RECREATE");
    TTree *tzprofile=new TTree("tzprofile","tzprofile");
    TTree *tfitEn=new TTree("tfitEn","tfitEn");
    double a;
    int event_no=0;
    vector<double> *layer_energy=0;
    vector<int> *layer=0;
    tzprofile->Branch("beam_energy",&a);
    tzprofile->Branch("layer",&layer);
    tzprofile->Branch("layer_energy",&layer_energy);
    tzprofile->Branch("event_no",&event_no);

    double mean=0;
    double sigma=0;
    tfitEn->Branch("beam_energy",&a);
    tfitEn->Branch("mean",&mean);
    tfitEn->Branch("sigma",&sigma);
    TF1* f=new TF1("f","gaus",0,10000);
    vector<double> be;
    unordered_map<double ,TH1D*> h;
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
            h[a]=new TH1D(hname,hname,10000,0,10000);
        }
		hbase->ReadTree(tmp,"EventTree");
        int flag=0;
        for(int ientry=0;ientry<hbase->tin->GetEntries();ientry++){
            hbase->tin->GetEntry(ientry);
            layer->clear();
            layer_energy->clear();
            event_no=ientry;
            double totalenergy=0;

            Select select(hbase->_Hit_X,hbase->_Hit_Y,hbase->_Hit_Z,hbase->_Digi_Hit_Energy);
            if(select.Result(totalenergy,layer_energy)==0){
                flag++;
                continue;
            }

            
            for(int l=0;l<40;l++){
                layer->push_back(l);
            }
            tzprofile->Fill();
            h[a]->Fill(totalenergy);          
        }
        cout<<flag<<" / "<<hbase->tin->GetEntries()<<endl;
        hbase->fin->Close();
    });
    // std::cout<<"Read done"<<std::endl;
    // std::unordered_map<int,TH1D*>::iterator it;
    // for(it=h.begin();it!=h.end();it++){
    fout->cd();
    for(int i=0;i<be.size();i++){
        a=be.at(i);
        TH1D *htmp=h[a];
        // mean=a/100;
        // sigma=mean/10;
        mean=htmp->GetMean();
        // std::cout<<htmp->GetEntries()<<std::endl;
        sigma=htmp->GetRMS();
        double fs=mean-2*sigma;
        if(mean<4*sigma)fs=mean/2;
        htmp->Fit(f,"q","",fs,mean+2*sigma);
        mean=f->GetParameter(1);
        sigma=f->GetParameter(2);
        fs=mean-2*sigma;
        if(mean<4*sigma)fs=mean/2;
        htmp->Fit(f,"q","",fs,mean+2*sigma);
        mean=f->GetParameter(1);
        sigma=f->GetParameter(2);
        tfitEn->Fill();
        htmp->Write();
    }
    // for(auto i:h){
    //     i.second->Write();
    // }
    tzprofile->Write();
    tfitEn->Write();
    // fout->Write();
}