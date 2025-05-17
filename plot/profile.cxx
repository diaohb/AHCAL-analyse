#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TStyle.h"
#include <vector>
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include <iostream>
using namespace std;
int profile()
{
    gStyle->SetOptStat(0);
    // TString s[6] = {"/cefs/higgs/diaohb/SIM/v2code/select_run/data/PS/e-list_analyse.root", "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0116_25000126_1010_trigger0/list_PS/calib_analyse.root", "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0116_25000126_1010_trigger40/list_PS/calib_analyse.root", "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0116_25000126_1010_trigger60/list_PS/calib_analyse.root", "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0115_25000126_1010_trigger100/list_PS/calib_analyse.root", "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0116_25000126_1010_trigger40/list_nobirks_non/truth_analyse.root"};
    // TString title[6] = {"data", "MC digi0", "MC digi40", "MC digi60", "MC digi100", "MC truth40"};
    TString pdfname="layer_profile_hitno.pdf";
    TString s[3] = {"/cefs/higgs/diaohb/SIM/v2code/select_run/data/PS_newmip/e-list_analyse55.root", "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0119_25000126_trigger40/list_PS_newmip/calib_analyse55.root", "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0119_25000126_trigger40/list_PS/truth_analyse55.root"};
    int color[3] = {2, 1, 3};
    TString title[6] = {"data", "MC digi", "MC truth"};
    double energy[6] = {0.5, 1, 2, 3, 4, 5};
    double range[6] = {2,3,4,4,5,6};
    // double range[6] = {3, 4.5, 8, 12, 15, 22};
    int rms_limit[6] = {100, 200, 300, 400, 400, 500};
    TFile *fin;
    TTree* tin;
    TH2D *hlayer[6][6];
    for (int i = 0; i < 3;i++){
        for(int n=0;n<6;n++){
            hlayer[i][n]=new TH2D(title[i],title[i],40,0,40,8000,0,200);
        }
        fin=TFile::Open(s[i],"READ");
        tin=(TTree*)fin->Get("tzprofile");
        double beam_energy=0;
        vector<int> *layer=0;
        vector<double> *layer_energy=0;
        vector<double> *layer_hitno=0;
        int flag=0;
        double rms = 0;
        tin->SetBranchAddress("beam_energy", &beam_energy);
        tin->SetBranchAddress("layer",&layer);
        tin->SetBranchAddress("layer_energy",&layer_energy);
        tin->SetBranchAddress("layer_hitno",&layer_hitno);
        tin->SetBranchAddress("select_flag", &flag);
        tin->SetBranchAddress("rms", &rms);
        // cout<<"flag"<<endl;
        for (int nentries = 0; nentries < tin->GetEntries(); nentries++)
        {
            tin->GetEntry(nentries);
            for (int n = 0; n < 6; n++)
            {
                if (beam_energy == energy[n] && flag == 6 && rms < rms_limit[n])
                {
                    // cout<<energy[n]<<endl;
                    for(int i_cell=0;i_cell<layer->size();i_cell++){
                        hlayer[i][n]->Fill(layer->at(i_cell), layer_hitno->at(i_cell));
                    }
                    break;
                }
            }
        }
        fin->Close();
    }
    cout<<"fill done"<<endl;
    TCanvas *c=new TCanvas("c","c");
    TH2D *htitle;
    c->Print(pdfname+"[");
    for (int n = 0; n < 6 ; n++){
        c->Clear();
        TLegend *tl=new TLegend(0.5,0.6,0.9,0.9);
        TString s;
        s.Form("%.1lfGeV e-",energy[n]);
        htitle=new TH2D(s,s,40,0,40,8000,0,200);
        htitle->GetXaxis()->SetTitle("layer");
        htitle->GetYaxis()->SetTitle("hitno");
        // htitle->GetYaxis()->SetTitle("energy/MeV");
        htitle->GetYaxis()->SetRangeUser(0,range[n]);
        htitle->Draw();
        for (int i = 0; i < 3;i++){
            hlayer[i][n]->SetMarkerColor(color[i]);
            hlayer[i][n]->SetLineColor(color[i]);
            hlayer[i][n]->SetMarkerSize(0.7);
            hlayer[i][n]->SetMarkerStyle(20);
            hlayer[i][n]->SetTitle(title[i]);
            TProfile *p = hlayer[i][n]->ProfileX("p");
            p->DrawCopy("same");
            tl->AddEntry(hlayer[i][n],title[i],"pl");
        }
        tl->Draw();
        c->Print(pdfname);
    }
    c->Print(pdfname+"]");
    return 1;
}
