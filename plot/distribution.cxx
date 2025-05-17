#include "/home/diaohb/plotstyle/plotstyle.c"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;
int distribution() {
    SetStyle();
    SetPrelimStyle();
    // int nbins[9] = {200, 200, 200, 200, 200, 200, 200, 200, 200};
    // double energy[9]={10,20,30,40,50,60,70,80,100};
    // double ranges[9] = {0, 0, 0, 0, 0, 0,0,0,0};
    // double rangee[9] = {300, 550, 800, 1000, 1300, 1500,2000,2500,3000};
    // double ranges[6] = {0, 40, 60, 80, 100, 100};
    // double rangee[6] = {100, 140, 160, 180, 200, 220};
    // int nbins[6] = {200, 200, 200, 200, 200, 200};
    // double energy[6] = {0.5, 1, 2, 3, 4, 5};
    // double ranges[6] = {0, 0, 10, 30, 40, 50};
    // double rangee[6] = {25, 40, 70, 80, 120, 130};
    // int rms_limit[6] = {100, 200, 300, 400, 400, 500};
    int nbins[11] = {200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200};
    double energy[11] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15};
    double ranges[11] = {0, 0, 0, 00, 0, 0, 0, 0, 0, 0, 0};
    double rangee[11] = {40, 80, 100, 150, 180, 200, 200, 250, 300, 450, 500};
    int rms_limit[11] = {300, 600, 1000, 1200, 1500, 1700, 2000, 2200, 2500, 2700, 3000};
    double limit_2345[11] = {7, 15, 20, 25, 35, 40, 45, 50, 60, 70, 90};
    TString pdf_name = "pi-figure/trigger40_PS_energy_3008_cherenkov2345.pdf";
    string txtname = "pi-txt/energy_PS_3008_cherenkov2345.txt";
    TString filename[7] = {"/home/diaohb/CEPC/ihep/PS/analyse/intercept/pi-list_analyse.root",
                           //    "/home/diaohb/CEPC/ihep/PS/analyse/split/list2_analyse.root",
                           "/home/diaohb/CEPC/AHCAL-simulation/run/0422_30080126_trigger40/list_PS_birks/calib_analyse.root",
                           "/home/diaohb/CEPC/AHCAL-simulation/run/0422_30080126_trigger40/list_PS_birks/truth_analyse.root",
                           "/home/diaohb/CEPC/AHCAL-simulation/run/0427_30080126_trigger40/list_PS/truth_analyse.root",
                           "/home/diaohb/CEPC/AHCAL-simulation/run/0428_30080126_trigger40/list_PS/truth_analyse.root",
                           "/home/diaohb/CEPC/AHCAL-simulation/run/0429_30080126_trigger40/list_PS/truth_analyse.root",
                           "/home/diaohb/CEPC/AHCAL-simulation/run/0422_30080126_trigger40/list_PS/truth_analyse.root"};
    // TString filename_a[3] = {"/cefs/higgs/diaohb/SIM/v2code/select_run/data/PS_newmip/e-list_analyse.root",
    //                         "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0312_30050126_trigger40/list_PS_lrms_spe_mpv/calib_analyse.root",
    //                         "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0312_30050126_trigger40/list_PS_lrms/truth_analyse.root"};
    TString histoname[7] = {"data", "MC digi", "MC truth", "MC truth w/ birks 0.07943", "MC truth w/ birks 0.04", "MC truth w/ birks 0.02", "MC truth w/o birks"};
    double parameters[12][2] = {{0.2, -0.3}, {0.35, -0.3}, {0.45, -0.25}, {0.5, -0.2}, {0.6, -0.2}, {0.6, -0.17}, {0.65, -0.15}, {0.8, -0.2}, {0.8, -0.1778}, {0.8, -0.14545}, {0.85, -0.14166}, {0.8, -0.1}};
    TF1 *fun = new TF1("fun", "[0]+[1]*x", 0, 10);
    int color[7] = {2, 1, 3, 6, 7, 8, 9};
    const int n_total = 11;
    TH1D *hen[7][n_total];
    for (int i = 0; i < 7; i++) {
        for (int n = 0; n < n_total; n++) {
            TString sn = TString::Format("%d GeV  ", int(energy[n]));
            hen[i][n] = new TH1D(sn + "hen", sn + "hen", nbins[n], ranges[n], rangee[n]);
        }
    }
    TFile *fin, *_fin;
    TTree *tin, *_tin;
    TCanvas *c = new TCanvas();
    for (int i = 0; i < 7; i++) {
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
        int hitno = 0;
        double fd = 0;
        double e_hit = 0;
        int cherenkov = 0;
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
        tin->SetBranchAddress("FD", &fd);
        tin->SetBranchAddress("E_Hit", &e_hit);
        tin->SetBranchAddress("cherenkov", &cherenkov);
        double _beam_energy = 0;
        int _select_flag = 0;
        vector<int> *_layer = 0;
        vector<double> *_layer_hitno = 0;
        vector<double> *_layer_energy = 0;
        vector<double> *_layer_hit_x = 0;
        vector<double> *_layer_hit_y = 0;
        double _total_energy = 0;
        int _hitno = 0;
        double _rms = 0;
        double _fd = 0;
        double _e_hit = 0;
        int _cherenkov = 0;
        _fin = TFile::Open(filename[i], "READ");
        _tin = (TTree *) _fin->Get("tanalyse");
        _tin->SetBranchAddress("beam_energy", &_beam_energy);
        _tin->SetBranchAddress("select_flag", &_select_flag);
        _tin->SetBranchAddress("layer", &_layer);
        _tin->SetBranchAddress("layer_energy", &_layer_energy);
        _tin->SetBranchAddress("layer_hitno", &_layer_hitno);
        _tin->SetBranchAddress("layer_hit_x", &_layer_hit_x);
        _tin->SetBranchAddress("layer_hit_y", &_layer_hit_y);
        _tin->SetBranchAddress("total_energy", &_total_energy);
        _tin->SetBranchAddress("hitno", &_hitno);
        _tin->SetBranchAddress("rms", &_rms);
        _tin->SetBranchAddress("FD", &_fd);
        _tin->SetBranchAddress("cherenkov", &_cherenkov);
        _tin->SetBranchAddress("E_Hit", &_e_hit);
        for (int n_entries = 0; n_entries < tin->GetEntries(); n_entries++) {
            _tin->GetEntry(n_entries);
            tin->GetEntry(n_entries);
            for (int n_e = 0; n_e < n_total; n_e++) {
                // if (n_e < 5) {
                fun->SetParameters(parameters[n_e + 1]);
                // }
                if (_beam_energy == energy[n_e] && _select_flag == 6 && _rms < rms_limit[n_e]) {
                    if (_cherenkov == 1) continue;
                    // if (fd > fun->Eval(e_hit)) break;
                    // int a = 0;
                    // for (int i = 0; i < _layer->size(); i++)
                    // {
                    //     if (_layer->at(i) >= 2 && _layer->at(i) < 3 && _layer_hit_x->at(i) != -500 && _layer_hit_y->at(i) != -500 && (abs(_layer_hit_x->at(i) + 20) > 15 || abs(_layer_hit_y->at(i) - 20) > 15))
                    //     {
                    //         a = 1;
                    //         break;
                    //     }
                    // }
                    // if (a == 1)
                    //     break;
                    if (layer_energy->size() == 0)
                        break;

                    double energy2345 = 0;
                    // int flag = 1;
                    for (int layer = 2; layer <= 5; layer++) {
                        energy2345 += layer_energy->at(layer);
                    }
                    if (energy2345 > limit_2345[n_e]) continue;
                    double tote = 0;
                    for (int layer = 0; layer < 38; layer++) {
                        // if (layer == 6)
                        //     continue;
                        tote += layer_energy->at(layer);
                    }
                    hen[i][n_e]->Fill(tote);
                    break;
                }
            }
        }
    }
    ofstream result(txtname, ios::out);
    gStyle->SetOptTitle(1);
    gStyle->SetFrameLineWidth(3);
    gStyle->SetPadTopMargin(0.15);
    for (int n_e = 0; n_e < n_total; n_e++) {
        int n_entries = hen[0][n_e]->GetEntries();
        cout << n_entries << endl;
        TF1 *fun = new TF1("fun", "gaus", ranges[n_e], rangee[n_e]);
        TLegend *tl = new TLegend(0.55, 0.64, 0.84, 0.84);
        TPad *pad1 = new TPad("p1", "p1", 0, 0, 1, 1);
        pad1->Draw();
        pad1->cd();
        TH1D *hframe = new TH1D("hframe", "hframe", nbins[n_e], ranges[n_e], rangee[n_e]);
        hframe->SetTitle(TString::Format("%d GeV pi-", int(energy[n_e])));
        NameAxis(hframe, "Energy [MeV]", "Event Number");
        FormatData(hframe);
        double ymax = hen[0][n_e]->GetBinContent(hen[0][n_e]->GetMaximumBin());
        hframe->GetYaxis()->SetRangeUser(0, 1.5 * ymax);
        hframe->GetYaxis()->SetMaxDigits(3);
        hframe->DrawCopy();
        for (int i = 0; i < 3; i++) {
            double mean = hen[i][n_e]->GetMean();
            double rms = hen[i][n_e]->GetRMS();
            // cout<<mean<<"  "<<rms<<endl;
            fun->SetRange(mean - 1.5 * rms, mean + 1.5 * rms);
            hen[i][n_e]->SetTitle(histoname[i]);
            FormatData(hen[i][n_e], color[i], 20);
            hen[i][n_e]->Scale(n_entries / double(hen[i][n_e]->GetEntries()));
            cout << histoname[i] + " entries: " << hen[i][n_e]->GetEntries() << endl;
            auto f = hen[i][n_e]->Fit(fun, "0rS", "");
            if (histoname[i] == "data") {
                hen[i][n_e]->Draw("sames hist p");
                tl->AddEntry(hen[i][n_e], "", "p");
            } else {
                hen[i][n_e]->Draw("sames hist");
                tl->AddEntry(hen[i][n_e], "", "l");
            }
            // fun->DrawCopy("same");
            mean = f->Parameter(1);
            rms = f->Parameter(2);
            result << i << "  " << int(energy[n_e]) << "  " << mean << "  " << rms << "\n";
        }
        tl->Draw();
        tl->SetFillColor(0);
        tl->SetTextSize(0.04);
        TString s_suffix = "";
        if (n_e == 0)
            s_suffix = "(";
        if (n_e == n_total - 1)
            s_suffix = ")";
        c->SaveAs(pdf_name + s_suffix);
    }
    result.close();
    return 1;
}
