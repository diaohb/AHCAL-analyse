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
#include <thread>
#include <mutex>
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
    unordered_map<int, TH1D*> mmip;
    for(int layer=0;layer<40;layer++){
        TString slayer="layer"+TString(to_string(layer).c_str());
        for(int chip=0;chip<9;chip++){
            TString schip="chip"+TString(to_string(chip).c_str());
            for(int channel=0;channel<36;channel++){
                TString schannel="channel"+TString(to_string(channel).c_str());
                TString name ="MIP Spectrum "+slayer+" "+schip+" "+schannel;
                int cellid=layer*1e5+chip*1e4+channel;
                mmip[cellid]=new TH1D(name,name,500,0,3);
            }
        }
    }
    tout=new TTree("mip","mip");
    TH2I *hitmap = new TH2I("hitmap", "hitmap", 18, -360, 360, 18, -360, 360);
    double MPV=0,width=0,gaus_sigma=0,max_x=0,FWHM=0;
    double _MPV[40][9][36],_width[40][9][36],_gaus_sigma[40][9][36],_max_x[40][9][36],_FWHM[40][9][36];
    std::fill(_MPV[0][0],_MPV[0][0]+12960,0);
    std::fill(_width[0][0],_width[0][0]+12960,0);
    std::fill(_max_x[0][0],_max_x[0][0]+12960,0);
    std::fill(_FWHM[0][0],_FWHM[0][0]+12960,0);
    std::fill(_gaus_sigma[0][0],_gaus_sigma[0][0]+12960,0);
    int cellid=0,entries=0;
    tout->Branch("MPV",&MPV);
    tout->Branch("width",&width);
    tout->Branch("gaus_sigma",&gaus_sigma);
    tout->Branch("cellid",&cellid);
    tout->Branch("entries",&entries);
    tout->Branch("max_x",&max_x);
    tout->Branch("FWHM",&FWHM);
    for_each(list.begin(),list.end(),[&](string tmp){
        cout<<"Reading: "<<tmp<<endl;
        fin=TFile::Open(TString(tmp),"read");
        tin=(TTree*)fin->Get("EventTree");
        ReadCalibTree(tin);
        for(int n=0;n<tin->GetEntries();n++){
            tin->GetEntry(n);
            // cout << MIP(cellID, Hit_E) << endl;
            if (MIP(cellID, Hit_E))
            {
                for(int i=0;i<cellID->size();i++){
                    int cid=cellID->at(i);
                    if(cid/100%100==0){
                        mmip[cid]->Fill(Hit_E->at(i));
                        double x=Pos_X_1(cid);
                        double y=Pos_Y_1(cid);
                        hitmap->Fill(x,y);
                    }
                }
            }
        }    
        fin->Close();
    });
    mutex g_mutex;
    auto ffit=[&](int layer){
        double fr[2];
        double sv[4], pllo[4], plhi[4], fps[4], fpe[4];
        double chisqr;
        int ndf;
        printf("fitting layer %d ...\n",layer);
        for(int chip=0;chip<9;chip++){
            for(int channel=0;channel<36;channel++){
                int cellid=layer*1e5+chip*1e4+channel;
                int entries=mmip[cellid]->GetEntries();
                if(entries<100){
                    // tout->Fill();
                    // mmip[cellid]->Write();
                    continue;
                }
                fr[0] = 0.25;fr[1] = 1;
                pllo[0] = 0;  pllo[1] = 0.4; pllo[2] = 0.5; pllo[3] = 0;
                plhi[0] = 0.5;plhi[1] = 0.5;plhi[2] = 5000;plhi[3] = 0.5;
                sv[0] = 0.02; sv[1] = 0.461; sv[2] = 100;   sv[3] = 0.01;
        // g_mutex.lock();
                langaufit(mmip[cellid], fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf);
        // g_mutex.unlock();
                double maxx=0,fwhm;
                langaupro(fps,maxx,fwhm);
                _max_x[layer][chip][channel]=maxx;
                _FWHM[layer][chip][channel]=fwhm;
                _MPV[layer][chip][channel]=fps[1];
                _width[layer][chip][channel]=fps[0];
                _gaus_sigma[layer][chip][channel]=fps[3];
            }
        }
        printf("fitted layer %d\n",layer);
    };
    for(int i = 0; i < 40; i++){
        ffit(i);
    }
        // thread th[40];
        // for(int i=0;i<40;i++){
        //     th[i]=thread(ffit,i);
        // }
        // for(int i=0;i<40;i++){
        //     th[i].join();
        // }
    fout->cd();
    hitmap->Write();
    TString dir="histogram";
    fout->mkdir(dir);
    for(int layer=0;layer<40;layer++){
        TString slayer="layer"+TString(to_string(layer).c_str());
        fout->mkdir(dir+"/"+slayer);
        for(int chip=0;chip<9;chip++){
            TString schip="chip"+TString(to_string(chip).c_str());
            fout->mkdir(dir+"/"+slayer+"/"+schip);
            fout->cd(dir+"/"+slayer+"/"+schip);
            for(int channel=0;channel<36;channel++){
                cellid=layer*1e5+chip*1e4+channel;
                entries=mmip[cellid]->GetEntries();
                if(entries<100){
                    tout->Fill();
                    mmip[cellid]->Write();
                    continue;
                }
                mmip[cellid]->Write();
                MPV=_MPV[layer][chip][channel];
                width=_width[layer][chip][channel];
                gaus_sigma=_gaus_sigma[layer][chip][channel];
                max_x=_max_x[layer][chip][channel];
                FWHM=_FWHM[layer][chip][channel];
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
int raw2Root::MIP(vector<int> *_cellid,vector<double> *Hit_E){
    int layerlen[40]={0};
    for(int i=0;i<_cellid->size();i++){
        int layer=_cellid->at(i)/100000;
        layerlen[layer]++;
    }
    for(int i=0;i<40;i++){
        if(layerlen[i]>2){
            return 0;
        }
    }
    return 1;
}