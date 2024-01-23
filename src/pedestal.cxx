#include "RawtoRoot.h"
#include "Global.h"
#include "langaus.h"
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
#include <unordered_map>
using namespace std;
int main(int argc,char *argv[]){
    double start = clock();
    raw2Root tw;
    tw.PEDlist(argv[1]);
    double end = clock();
    cout<<"end of mip : Time : "<<(end-start)/CLOCKS_PER_SEC<<endl;
    return 0;
}
int raw2Root::PEDlist(const string _list){
    ReadList(_list);
    TFile *fin,*fout;
    TTree *tin,*tout;
    
    fout=TFile::Open(TString(_list)+"_pedestal.root","recreate");
    unordered_map<int, TH1D*> mped_high;
    unordered_map<int, TH1D*> mped_low;
    for(int layer=0;layer<40;layer++){
        TString slayer="layer"+TString(to_string(layer).c_str());
        for(int chip=0;chip<9;chip++){
            TString schip="chip"+TString(to_string(chip).c_str());
            TString name_chip ="MIP Spectrum "+slayer+" "+schip;
            for(int channel=0;channel<36;channel++){
                TString schannel="channel"+TString(to_string(channel).c_str());
                TString name ="highgain pedestal "+slayer+" "+schip+" "+schannel;
                int cellid=layer*1e5+chip*1e4+channel;
                mped_high[cellid]=new TH1D(name,name,500,0,1000);
                name ="lowgain pedestal "+slayer+" "+schip+" "+schannel;
                mped_low[cellid]=new TH1D(name,name,500,0,1000);
            }
        }
    }
    tout=new TTree("pedestal","pedestal");
    double highpeak=0,highrms=0,lowpeak=0,lowrms=0;
    int cellid=0;
    tout->Branch("cellid",&cellid);
    tout->Branch("highgain_peak",&highpeak);
    tout->Branch("highgain_rms",&highrms);
    tout->Branch("lowgain_peak",&lowpeak);
    tout->Branch("lowgain_rms",&lowrms);
    for_each(list.begin(),list.end(),[&](string tmp){
        cout<<"Reading: "<<tmp<<endl;
        fin=TFile::Open(TString(tmp),"read");
        tin=(TTree*)fin->Get("Raw_Hit");
        ReadTreeBranch(tin);
        for(int n=0;n<tin->GetEntries();n++){
            tin->GetEntry(n);
            for(int i=0;i<cellID->size();i++){
                mped_high[cellID->at(i)]->Fill(HG_Charge->at(i));
                mped_low[cellID->at(i)]->Fill(LG_Charge->at(i));
            }
        }    
        fin->Close();
    });
    fout->cd();
    TString dir="histogram";
    fout->mkdir(dir);
    fout->mkdir(dir+"/highgain");
    fout->mkdir(dir+"/lowgain");
    TF1 *ffit=new TF1("ffit","gaus",0,1000);
    for(int layer=0;layer<40;layer++){
        TString slayer="layer"+TString(to_string(layer).c_str());
        cout<<"fitting "<<slayer<<" ..."<<endl;
        fout->mkdir(dir+"/highgain/"+slayer);
        fout->mkdir(dir+"/lowgain/"+slayer);
        for(int chip=0;chip<9;chip++){
            TString schip="chip"+TString(to_string(chip).c_str());
            fout->mkdir(dir+"/highgain/"+slayer+"/"+schip);
            fout->mkdir(dir+"/lowgain/"+slayer+"/"+schip);
            for(int channel=0;channel<36;channel++){
                cellid=layer*1e5+chip*1e4+channel;
                fout->cd(dir+"/highgain/"+slayer+"/"+schip);
                double mean=mped_high[cellid]->GetMean();
                double sigma=mped_high[cellid]->GetRMS();
                mped_high[cellid]->Fit(ffit,"q","",mean-2*sigma,mean+2*sigma);
                highpeak=ffit->GetParameter(1);
                highrms=ffit->GetParameter(2);
                mped_high[cellid]->Write();

                fout->cd(dir+"/lowgain/"+slayer+"/"+schip);
                mean=mped_low[cellid]->GetMean();
                sigma=mped_low[cellid]->GetRMS();
                mped_low[cellid]->Fit(ffit,"q","",mean-2*sigma,mean+2*sigma);
                lowpeak=ffit->GetParameter(1);
                lowrms=ffit->GetParameter(2);
                mped_low[cellid]->Write();

                tout->Fill();
            }
        }
    }
    fout->cd();
    tout->Write();
    fout->Close();
    return 1;
}
// int raw2Root::ReadList(const string _list){
//     ifstream data(_list);
// 	while(!data.eof())
// 	{
// 			string temp;
// 			data>>temp;
// 			if(temp=="")continue;
//             if (data.fail())break;
// 			list.push_back(temp);
// 	}
//     return 1;
// }