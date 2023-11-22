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
using namespace std;
int main(int argc,char *argv[]){
    double start = clock();
    raw2Root tw;
    tw.MIPlist(argv[1]);
    double end = clock();
    cout<<"end of mip : Time : "<<(end-start)/CLOCKS_PER_SEC<<endl;
    return 0;
}
int raw2Root::MIPlist(const string _list){
    ReadList(_list);
    TFile *fin,*fout;
    TTree *tin,*tout;
    fout=TFile::Open(TString(_list)+"_mip.root","recreate");
    tout=new TTree("mip","mip");
    double mip=0;
    int cellid=0;
    tout->Branch("mip",&mip);
    tout->Branch("cellid",&cellid);
    for_each(list.begin(),list.end(),[&](string tmp){
        cout<<"Reading: "<<tmp<<endl;
        fin=TFile::Open(TString(tmp),"read");
        tin=(TTree*)fin->Get("Raw_Hit");
        ReadTreeBranch(tin);
        for(int n=0;n<tin->GetEntries();n++){
            tin->GetEntry(n);
            if(MIP(cellID,HG_Charge,hitTag)){
                for(int i=0;i<cellID->size();i++){
                    if(hitTag->at(i)==1){
                        mip=HG_Charge->at(i);
                        cellid=cellID->at(i);
                        tout->Fill();
                    }
                }
            }
        }    
        fin->Close();
    });
    fout->cd();
    tout->Write();
    return 1;
}
int raw2Root::ReadList(const string _list){
    ifstream data(_list);
	while(!data.eof())
	{
			string temp;
			data>>temp;
			if(temp=="")continue;
			list.push_back(temp);
	}
    return 1;
}
int raw2Root::MIP(vector<int> *_cellid,vector<double> *_HG_Charge,vector<int> *_hitTag){
    int layerlen[40]={0};
    for(int i=0;i<_cellid->size();i++){
        if(_hitTag->at(i)==1){
            int layer=_cellid->at(i)/100000;
            layerlen[layer]++;
        }
    }
    if(layerlen[5]!=1)return 0;
    for(int i=0;i<40;i++){
        if(layerlen[i]>2){
            return 0;
        }
    }
    return 1;
}