#include "/home/diaohb/plotstyle/plotstyle.c"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include <iostream>
#include <vector>
using namespace std;
int violin() {
    SetStyle();
    SetPrelimStyle();
    // int nbins[9] = {200, 200, 200, 200, 200, 200, 200, 200, 200};
    // double energy[9]={10,20,30,40,50,60,70,80,100};
    // double ranges[9] = {100, 250, 100, 500, 200, 200,300,300,300};
    // double rangee[9] = {300, 550, 800, 1000, 1300, 1500,2000,2500,3000};
    // double ranges[6] = {0, 40, 60, 80, 100, 100};
    // double rangee[6] = {100, 140, 160, 180, 200, 220};
    // int nbins[6] = {200, 200, 200, 200, 200, 200};
    int nbins[6] = {20, 20, 20, 20, 20, 20};
    double energy[6] = {0.5, 1, 2, 3, 4, 5};
    // double rangee[6] = {20, 30, 40, 50, 60, 70};
    double rangee[6] = {20, 20, 20, 20, 20, 20};
    int rms_limit[6] = {100, 200, 300, 400, 400, 500};
    // int nbins[11] = {200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200};
    // double energy[11] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15};
    // double ranges[11] = {0, 0, 0, 00, 0, 0, 0, 0, 0, 0, 0};
    // double rangee[11] = {40, 80, 100, 100, 120, 120, 120, 120, 130, 150, 150};
    // int rms_limit[11] = {300, 600, 1000, 1200, 1500, 1700, 2000, 2200, 2500, 2700, 3000};
    // double limit_2345[11] = {7, 15, 20, 25, 35, 40, 45, 50, 60, 70, 90};
    TString yaxis = "hitno";
    TString pdfname = "figure/layer_" + yaxis;
    // TString filename[3] = {"/home/diaohb/CEPC/ihep/PS/analyse/intercept/pi-list_analyse.root",
    //                        "/home/diaohb/CEPC/AHCAL-simulation/run/0422_30080126_trigger40/list_PS_birks/calib_analyse.root",
    //                        "/home/diaohb/CEPC/AHCAL-simulation/run/0422_30080126_trigger40/list_PS_birks/truth_analyse.root"};
    TString filename[3] = {"/home/diaohb/CEPC/ihep/PS/analyse/intercept/e-list_analyse.root",
                           "/home/diaohb/CEPC/AHCAL-simulation/run/0405_30080126_trigger40/list_PS9/calib_analyse.root",
                           "/home/diaohb/CEPC/AHCAL-simulation/run/0405_30080126_trigger40/list_PS/truth_analyse.root"};
    TString histoname[3] = {"data", "MC_digi", "MC_truth"};
    const int n_total = 6;
    TH2D *hen[3][n_total];
    for (int i = 0; i < 3; i++) {
        for (int n = 0; n < n_total; n++) {
            TString sn = TString::Format("%d GeV  ", int(energy[n])) + histoname[i];
            hen[i][n] = new TH2D(sn, sn, 40, 0, 40, nbins[n], 0, rangee[n]);
        }
    }
    TFile *fin;
    TTree *tin;
    TCanvas *c = new TCanvas();
    for (int i = 0; i < 3; i++) {
        cout << histoname[i] << endl;
        fin = TFile::Open(filename[i], "READ");
        tin = (TTree *) fin->Get("tanalyse");
        double beam_energy = 0;
        int select_flag = 0;
        vector<int> *layer = 0;
        vector<double> *layer_hitno = 0;
        vector<double> *layer_energy = 0;
        vector<double> *layer_hit_x = 0;
        vector<double> *layer_hit_y = 0;
        double total_energy = 0;
        double rms = 0;
        int cherenkov = 0;
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
        tin->SetBranchAddress("cherenkov", &cherenkov);
        tin->SetBranchAddress("rms", &rms);
        for (int n_entries = 0; n_entries < tin->GetEntries(); n_entries++) {
            tin->GetEntry(n_entries);
            for (int n_e = 0; n_e < n_total; n_e++) {
                if (beam_energy == energy[n_e] && select_flag == 6 && rms < rms_limit[n_e]) {
                    // if (cherenkov == 1) continue;
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
                    // double energy2345 = 0;
                    // // int flag = 1;
                    // for (int layer = 2; layer <= 5; layer++) {
                    //     energy2345 += layer_energy->at(layer);
                    // }
                    // if (energy2345 > limit_2345[n_e]) continue;
                    for (int layer = 0; layer < 40; layer++) {
                        hen[i][n_e]->Fill(layer, layer_hitno->at(layer));
                    }
                    break;
                }
            }
        }
    }
    gStyle->SetOptTitle(1);
    gStyle->SetFrameLineWidth(3);
    gStyle->SetPadTopMargin(0.15);
    for (int i = 0; i < 3; i++) {
        for (int n_e = 0; n_e < n_total; n_e++) {
            FormatData(hen[i][n_e], 1, 20);
            // hen[i][n_e] ->SetMaximum() ;
            c->SetLogz();
            NameAxis(hen[i][n_e], "Layer", yaxis);
            hen[i][n_e]->Draw("colz");
            TString s_suffix = "";
            if (n_e == 0)
                s_suffix = "(";
            if (n_e == n_total - 1)
                s_suffix = ")";
            c->SaveAs(pdfname + "_" + histoname[i] + ".pdf" + s_suffix);
        }
    }
    return 1;
}
