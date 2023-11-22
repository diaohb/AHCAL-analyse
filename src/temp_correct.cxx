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
using namespace std;
int main(int argc,char* argv[]){
    HBase *hbase=new HBase();
    hbase->ReadList(argv[1]);
    ifstream in;
    string ftxt=argv[2];
    in.open(ftxt,ios::in);
    double frac;
    double f[18][18];
    int i=0,j=0;
    while(!in.eof()){
        in>>frac;
        if(in.fail()){    
		    break;    
		}
        f[i][j]=frac;
        cout<<"("<<f[i][j]<<")";
        i=i+(j+1)/18;
        j=(j+1)%18;
        if(j==0)cout<<endl;
    }

    for_each(hbase->list.begin(),hbase->list.end(),[&](string tmp){
		hbase->ReadTree(tmp,"EventTree","digi");
        string s=tmp;
        s=s.substr(s.find_last_of('/')+1);
        TString fname="correct_"+TString(s);
        TFile *fout=TFile::Open(fname,"RECREATE");
        TTree *tout=new TTree("EventTree","EventTree");
        int   _Event_No=0;
	    double _Digi_Energy_ECAL=0;
	    double _Digi_Energy_HCAL=0;
	    vector< int > *_detectorID=0;
	    vector< int > *_cellID=0;
	    vector< double > *_Digi_Hit_Energy=0;
	    vector< double > *_Hit_Time=0;
	    vector< double > *_Hit_X=0;
	    vector< double > *_Hit_Y=0;
	    vector< double > *_Hit_Z=0;
        // tout->Branch("EventNum",&_Event_No);
	    // tout->Branch("DetectorID",&_detectorID);
	    // tout->Branch("CellID",&_cellID);
	    tout->Branch("Digi_Hit_Energy",&_Digi_Hit_Energy);
	    // tout->Branch("Digi_Energy_ECAL",&_Digi_Energy_ECAL);
	    tout->Branch("Digi_Energy_HCAL",&_Digi_Energy_HCAL);
	    // tout->Branch("Hit_Time",&_Hit_Time);
	    tout->Branch("Hit_X",&_Hit_X);
	    tout->Branch("Hit_Y",&_Hit_Y);
	    tout->Branch("Hit_Z",&_Hit_Z);
        for(int ientry=0;ientry<hbase->tin->GetEntries();ientry++){
            hbase->tin->GetEntry(ientry);
            _Event_No=0;
            _Digi_Energy_ECAL=0;
            // _detectorID->clear();
            // _cellID->clear();
            _Digi_Hit_Energy->clear();
            // _Hit_Time->clear();
            _Hit_X->clear();
            _Hit_Y->clear();
            _Hit_Z->clear();
            _Digi_Energy_HCAL=0;      
            for(int icell=0;icell<hbase->_Hit_X->size();icell++){
                // _cellID->push_back(hbase->_cellID->at(icell));
                // _Hit_Time->push_back(hbase->_Hit_Time->at(icell));
                double x=hbase->_Hit_X->at(icell);
                double y=hbase->_Hit_Y->at(icell);
                _Hit_X->push_back(x);
                _Hit_Y->push_back(y);
                _Hit_Z->push_back(hbase->_Hit_Z->at(icell));
                int j=(x+360)/40;
                int i=(360-y)/40;
                double e=hbase->_Digi_Hit_Energy->at(icell);
                // cout<<x<<" , "<<j<<"    "<<y<<" , "<<i<<endl;
                double ee=e*f[i][j];
                if(hbase->_Hit_Z->at(icell)<30&&x>-40&&x<0&&y>-40&&y>0)cout<<e<<"  "<<ee<<endl;
                _Digi_Hit_Energy->push_back(ee);
                _Digi_Energy_HCAL+=ee;
                // std::cout<<e<<" "<<_Digi_Energy_HCAL<<" "<<f[i][j]<<std::endl;
            }
            if(ientry%10000==0)std::cout<<ientry<<" / "<<hbase->tin->GetEntries()<<std::endl;
            tout->Fill();         
        }
        hbase->fin->Close();
        fout->cd();
        tout->Write();
    });
}