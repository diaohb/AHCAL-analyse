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
#include "TStyle.h"
#include "TF1.h"
#include "/cefs/higgs/diaohb/plotstyle/plotstyle.c"
using namespace std;
int violin()
{
    SetStyle();
    SetPrelimStyle();
    // int nbins[9] = {200, 200, 200, 200, 200, 200, 200, 200, 200};
    // double energy[9]={10,20,30,40,50,60,70,80,100};
    // double ranges[9] = {100, 250, 100, 500, 200, 200,300,300,300};
    // double rangee[9] = {300, 550, 800, 1000, 1300, 1500,2000,2500,3000};
    // double ranges[6] = {0, 40, 60, 80, 100, 100};
    // double rangee[6] = {100, 140, 160, 180, 200, 220};
    int nbins[6] = {200, 200, 200, 200, 200, 200};
    // int nbins[6] = {20, 20, 20, 20, 20, 20};
    double energy[6] = {0.5, 1, 2, 3, 4, 5};
    double rangee[6] = {20, 30, 40, 50, 60, 70};
    // double rangee[6] = {20, 20, 20, 20, 20, 20};
    int rms_limit[6]={100,200,300,400,400,500};
    TString yaxis="energy";
    TString pdfname="layer_"+yaxis;
    TString filename[3] = {"/cefs/higgs/diaohb/SIM/v2code/select_run/data/PS_newmip/e-list_analyse.root",
                           "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0119_25000126_trigger40/list_PS_newmip/calib_analyse.root",
                           "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0119_25000126_trigger40/list_PS/truth_analyse.root"};
    TString histoname[3] = {"data","MC_digi", "MC_truth"};
    int n_total = 6;
    TH2D *hen[3][n_total];
    for (int i = 0; i < 3; i++)
    {
        for (int n = 0; n < n_total; n++)
        {
            TString sn = TString::Format("%d GeV  ", int(energy[n]))+histoname[i];
            hen[i][n] = new TH2D(sn , sn ,40,0,40, nbins[n], 0, rangee[n]);
        }
    }
    TFile *fin;
    TTree *tin;
    TCanvas *c = new TCanvas();
    for (int i = 0; i < 3; i++){
        cout << histoname[i] << endl;
        fin = TFile::Open(filename[i], "READ");
        tin = (TTree *)fin->Get("tzprofile");
        double beam_energy = 0;
        int select_flag = 0;
        vector<int> *layer = 0;
        vector<double> *layer_hitno = 0;
        vector<double> *layer_energy = 0;
        vector<double> *layer_hit_x = 0;
        vector<double> *layer_hit_y = 0;
        double total_energy = 0;
        double rms=0;
        int hitno = 0;
        tin->SetBranchAddress("beam_energy", &beam_energy);
        tin->SetBranchAddress("select_flag", &select_flag);
        tin->SetBranchAddress("layer", &layer);
        tin->SetBranchAddress("layer_energy", &layer_energy);
        tin->SetBranchAddress("layer_hitno", &layer_hitno);
        tin->SetBranchAddress("layer_hit_x", &layer_hit_x);
        tin->SetBranchAddress("layer_hit_y", &layer_hit_y);
        tin->SetBranchAddress("total_energy", &total_energy);
        tin->SetBranchAddress("hitno", &hitno);
        tin->SetBranchAddress("rms", &rms);
        for (int n_entries = 0; n_entries < tin->GetEntries(); n_entries++)
        {
            tin->GetEntry(n_entries);
            for (int n_e = 0; n_e < n_total;n_e++)
            {
                if (beam_energy == energy[n_e] && select_flag == 6 && rms<rms_limit[n_e])
                {
                    // int a = 0;
                    // for (int i = 0; i < layer->size(); i++)
                    // {
                    //     if (layer->at(i) >= 4 && layer->at(i) < 7 && layer_hit_x->at(i) != -500 && layer_hit_y->at(i) != -500 && (abs(layer_hit_x->at(i) - 20) > 10 || abs(layer_hit_y->at(i) - 20) > 10))
                    //     {
                    //         a = 1;
                    //         break;
                    //     }
                    // }
                    // if (a == 1)
                    //     continue;
                    // double tote = 0;
                    // for (int layer = 0; layer < 20; layer++)
                    // {
                    //     // if (layer == 6)
                    //     //     continue;
                    //     tote += layer_energy->at(layer);
                    // }
                    // if(tote>70)
                    //     break;
                    for (int layer = 0; layer < 40; layer++)
                    {
                        hen[i][n_e]->Fill(layer,layer_energy->at(layer));
                    }
                    break;
                }
            }
        }
    }
        gStyle->SetOptTitle(1);
        gStyle->SetFrameLineWidth(3);
        gStyle->SetPadTopMargin(0.15);
    for(int i=0;i<3;i++){   
        for (int n_e = 0; n_e < n_total; n_e++)
        {
            FormatData(hen[i][n_e], 1, 20);
            // hen[i][n_e] ->SetMaximum() ;
            c->SetLogz();
            NameAxis(hen[i][n_e],"Layer",yaxis);
            hen[i][n_e]->Draw("colz");
            TString s_suffix = "";
            if (n_e == 0)
                s_suffix = "(";
            if (n_e == n_total - 1)
                s_suffix = ")";
            c->SaveAs(pdfname+"_"+histoname[i]+".pdf" + s_suffix);
        }
    }
    return 1;
}
