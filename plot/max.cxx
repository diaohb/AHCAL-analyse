#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include "TH1.h"
#include "TH2.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TF1.h"
#include "/cefs/higgs/diaohb/plotstyle/plotstyle.c"
using namespace std;
int max()
{
    SetStyle();
    SetPrelimStyle();
    int nbins[6] = {200, 200, 200, 200, 200, 200};
    double energy[9]={10,20,30,40,50,60,70,80,100};
    double ranges[6] = {100, 250, 400, 550, 700, 700};
    double rangee[9] = {40, 60, 80, 100, 130, 150,200,200,250};
    // double ranges[6] = {0, 40, 60, 80, 100, 100};
    // double rangee[6] = {100, 140, 160, 180, 200, 220};
    // double energy[6] = {0.5, 1, 2, 3, 4, 5};
    // double ranges[6] = {0, 0, 0, 0, 0, 0};
    // double rangee[6] = {30, 50, 80, 100, 150, 200};
    // double rangee[6] = {20, 30, 50, 50, 60, 80};
    TString filename[3] = {"/cefs/higgs/diaohb/SIM/v2code/select_run/data/low800/e-list_analyse.root",
                           "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0115_25000126_1010_trigger100/list_800/calib_analyse.root",
                           "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0115_25000126_1010_trigger100/list/truth_analyse.root"};
    TString filename_a[2]={"/cefs/higgs/diaohb/SIM/v2code/select_run/data/origin/e-list_analyse.root","/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0115_25000126_1010_trigger100/list/truth_analyse.root"};
    TString histoname[3] = {"data","MC digi", "MC truth"};
    int color[3] = {2,1,3};
    int n_total = 9;
    TH2D *h[3][n_total];
    for (int i = 0; i < 3; i++)
    {
        for (int n = 0; n < n_total; n++){
            TString sn = TString::Format("%d GeV  ", int(energy[n]));
            h[i][n] = new TH2D(sn, sn, 40, 0, 40, 300, 0, 300);
        }
    }
    TFile *fin,*_fin;
    TTree *tin,*_tin;
    TCanvas *c = new TCanvas();
    for (int i_file = 0; i_file < 3; i_file++){
        cout << histoname[i_file] << endl;
        fin = TFile::Open(filename[i_file], "READ");
        tin = (TTree *)fin->Get("tzprofile");
        double beam_energy = 0;
        int select_flag = 0;
        vector<int> *layer = 0;
        vector<double> *layer_hitno = 0;
        vector<double> *layer_energy = 0;
        vector<double> *layer_hit_x = 0;
        vector<double> *layer_hit_y = 0;
        vector<double> *layer_max_energy = 0;
        vector<double> *layer_max_x = 0;
        vector<double> *layer_max_y = 0;
        double total_energy = 0;
        int hitno = 0;
        tin->SetBranchAddress("beam_energy", &beam_energy);
        tin->SetBranchAddress("select_flag", &select_flag);
        tin->SetBranchAddress("layer", &layer);
        tin->SetBranchAddress("layer_energy", &layer_energy);
        tin->SetBranchAddress("layer_hitno", &layer_hitno);
        tin->SetBranchAddress("layer_hit_x", &layer_hit_x);
        tin->SetBranchAddress("layer_hit_y", &layer_hit_y);
        tin->SetBranchAddress("layer_max_energy", &layer_max_energy);
        tin->SetBranchAddress("layer_max_x", &layer_max_x);
        tin->SetBranchAddress("layer_max_y", &layer_max_y);
        tin->SetBranchAddress("total_energy", &total_energy);
        tin->SetBranchAddress("hitno", &hitno);
        double _beam_energy = 0;
            int _select_flag = 0;
            vector<int> *_layer = 0;
            vector<double> *_layer_hitno = 0;
            vector<double> *_layer_energy = 0;
            vector<double> *_layer_hit_x = 0;
            vector<double> *_layer_hit_y = 0;
            double _total_energy = 0;
            int _hitno = 0;
        if(i_file<2){
            _fin = TFile::Open(filename_a[i_file], "READ");
            _tin = (TTree *)_fin->Get("tzprofile");
            _tin->SetBranchAddress("beam_energy", &_beam_energy);
            _tin->SetBranchAddress("select_flag", &_select_flag);
            _tin->SetBranchAddress("layer", &_layer);
            _tin->SetBranchAddress("layer_energy", &_layer_energy);
            _tin->SetBranchAddress("layer_hitno", &_layer_hitno);
            _tin->SetBranchAddress("layer_hit_x", &_layer_hit_x);
            _tin->SetBranchAddress("layer_hit_y", &_layer_hit_y);
            _tin->SetBranchAddress("total_energy", &_total_energy);
            _tin->SetBranchAddress("hitno", &_hitno);
        }
        for (int n_entries = 0; n_entries < _tin->GetEntries(); n_entries++)
        {
             _tin->GetEntry(n_entries);
            for (int n_e = 0; n_e < n_total;n_e++){
                if (_beam_energy == energy[n_e] && _select_flag == 6)
                {
                    int a = 0;
                    for (int i = 0; i < _layer->size(); i++)
                    {
                        if (_layer->at(i) >= 4 && _layer->at(i) < 5 && _layer_hit_x->at(i) != -500 && _layer_hit_y->at(i) != -500 && (abs(_layer_hit_x->at(i) - 20) > 10 || abs(_layer_hit_y->at(i) - 20) > 10))
                        {
                            a = 1;
                            break;
                        }
                    }
                    if (a == 1)
                        break;
                    tin->GetEntry(n_entries);
                    if(layer_energy->size()==0)
                        break;
                    for (int i = 0; i < 40; i++)
                    {
                        // if(i==6)continue;
                        if(layer_max_energy->at(i)!=0)
                            h[i_file][n_e]->Fill(i,layer_max_energy->at(i));
                    }
                    break;
                }
            }
        }
    }
    cout << "fill done" << endl;
    gStyle->SetOptTitle(1);
    gStyle->SetFrameLineWidth(3);
    gStyle->SetPadTopMargin(0.15);
    for (int n_e = 0; n_e < n_total; n_e++)
    {
        cout << n_e << endl;
        TLegend *tl = new TLegend(0.62, 0.24, 0.8, 0.44);
        TPad *pad1 = new TPad("p1", "p1", 0, 0, 1, 1);
        pad1->Draw();
        pad1->cd();
        TH1D *hframe = new TH1D("hframe", "hframe", 40,0,40);
        hframe->SetTitle(TString::Format("%d GeV e-", int(energy[n_e])));
        NameAxis(hframe, "layer", "Max Energy");
        FormatData(hframe);
        hframe->GetYaxis()->SetRangeUser(40, 100);
        hframe->GetYaxis()->SetMaxDigits(3);
        hframe->DrawCopy();
        for (int i = 0; i < 3; i++)
        {
            h[i][n_e]->SetTitle(histoname[i]);
            FormatData(h[i][n_e], color[i], 20);
            TProfile *p = h[i][n_e]->ProfileX("p");
            if (histoname[i] == "data")
            {
                p->DrawCopy("same p");
                tl->AddEntry(h[i][n_e], "", "p");
            }
            else
            {
                p->DrawCopy("same p");
                tl->AddEntry(h[i][n_e], "", "p");
            }
        }
        tl->Draw();
        tl->SetFillColor(0);
        tl->SetTextSize(0.04);
        TString s_suffix = "";
        if (n_e == 0)
            s_suffix = "(";
        if (n_e == n_total - 1)
            s_suffix = ")";
        c->SaveAs("trigger100_SPS_max.pdf" + s_suffix);
    }
    return 1;
}
