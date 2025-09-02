#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"
#include "iostream"
using namespace std;
int select_flow() {
    const int number = 11;
    // double energy[6] = {0.5, 1, 2, 3, 4, 5};
    // int rms_limit[6] = {100, 200, 300, 400, 400, 500};
    // int hitno_limit[6] = {15, 25, 35, 45, 55, 60};
    double energy[11] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15};
    double rangee[11] = {40, 80, 100, 150, 180, 200, 200, 250, 300, 450, 500};
    // double energy[number] = {0.5, 1, 2, 3, 4, 5};
    // double rangee[number] = {40, 80, 100, 150, 180, 200};
    double limit_2345[11] = {7, 15, 20, 25, 35, 40, 45, 50, 60, 70, 90};

    // int rms_limit[number] = {300, 600, 1000, 1200, 1500, 1700, 2000, 2200, 2500, 2700, 3000};
    // int hitno_limit[number] = {40, 50, 60, 70, 80, 100, 110, 120, 140, 160, 200};
    // double range2345[number] = {30, 60, 80, 90, 110, 120, 130, 150, 170, 190, 220};
    // ROOT::RDataFrame df("tanalyse", "/home/diaohb/CEPC/ihep/PS/analyse/intercept/pi-list_analyse.root");
    ROOT::RDataFrame df("tanalyse", "/home/diaohb/CEPC/AHCAL-simulation/run/0422_30080126_trigger40/list_PS_birks/calib_analyse.root");
    // ROOT::RDataFrame df("tanalyse", "/home/diaohb/CEPC/AHCAL-simulation/run/0405_30080126_trigger40/list_PS9/calib_analyse.root");
    TString pdf_name = "pi-figure/select_flow_MC.pdf";
    TString h_name[11] = {"origin", "0.5 mip", "start", "end", "hitno", "hitlayer", "mip", "1st layer 1", "che", "2345"};
    TH1D *h1[number][11];
    for (int i = 0; i < number; i++) {
        for (int j = 0; j < 10; j++) {
            h1[i][j] = new TH1D("energy", TString::Format("%d GeV;energy [MeV];Events;", int(energy[i])), 200, 0, rangee[i]);
        }
        // cout << h6[i] << endl;
    }
    df.Foreach([&](double beam_energy, int select_flag, vector<double> layer_energy, int cherenkov) {
        double *pos = find(energy, energy + number, beam_energy);
        int i = pos - energy;
        if (i < number && layer_energy.size() != 0) {
            double total_energy = 0;
            for (int layer = 0; layer < 38; layer++) {
                total_energy += layer_energy.at(layer);
            }
            h1[i][0]->Fill(total_energy);
            if (select_flag != 0) {
                h1[i][1]->Fill(total_energy);
                if (select_flag != 3) {
                    h1[i][2]->Fill(total_energy);
                    if (select_flag != 9) {
                        h1[i][3]->Fill(total_energy);
                        if (select_flag != 1) {
                            h1[i][4]->Fill(total_energy);
                            if (select_flag != 2) {
                                h1[i][5]->Fill(total_energy);
                                if (select_flag != 4) {
                                    h1[i][6]->Fill(total_energy);
                                    if (select_flag != 8) {
                                        h1[i][7]->Fill(total_energy);
                                        if (cherenkov != 1) {
                                            h1[i][8]->Fill(total_energy);
                                            double energy2345 = 0;
                                            for (int j = 2; j <= 5; j++) {
                                                energy2345 += layer_energy.at(j);
                                            }
                                            if (energy2345 < limit_2345[i]) {
                                                h1[i][9]->Fill(total_energy);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    },
               {"beam_energy", "select_flag", "layer_energy", "cherenkov"});
    TCanvas *c = new TCanvas("c", "c", 2400, 2400);
    TF1 *fun = new TF1("fun", "gaus", 0, 1000);
    fun->SetNpx(1000);
    gStyle->SetOptStat(0);
    c->Print(pdf_name + "[");
    for (int i = 0; i < number; i++) {
        c->Clear();
        TLegend *tl = new TLegend(0.55, 0.64, 0.94, 0.94);
        h1[i][0]->GetYaxis()->SetRangeUser(1, 1e5);
        h1[i][0]->Draw("plc");
        double mean = h1[i][0]->GetMean();
        double rms = h1[i][0]->GetRMS();
        fun->SetRange(mean - 2. * rms, mean + 2. * rms);
        h1[i][0]->Fit(fun, "rq", "");
        mean = fun->GetParameter(1);
        rms = fun->GetParameter(2);
        fun->SetRange(mean - 1. * rms, mean + 1. * rms);
        h1[i][0]->Fit(fun, "rq", "");
        mean = fun->GetParameter(1);
        rms = fun->GetParameter(2);
        cout << 0 << "  " << h1[i][0]->GetEntries() << "  " << mean << "  " << rms << "  " << rms / mean << endl;
        for (int j = 1; j < 10; j++) {
            h1[i][j]->Draw("same plc");
            tl->AddEntry(h1[i][j], TString::Format("%dk  ", int(h1[i][j]->GetEntries() / 1000)) + h_name[j], "pl");
            c->SetLogy(1);

            // mean = h1[i][j]->GetMean();
            // rms = h1[i][j]->GetRMS();
            fun->SetRange(mean - 1. * rms, mean + 1. * rms);
            h1[i][j]->Fit(fun, "qr", "");
            mean = fun->GetParameter(1);
            rms = fun->GetParameter(2);
            cout << j << "  " << h1[i][j]->GetEntries() << "  " << mean << "  " << rms << "  " << rms / mean << endl;
        }
        cout << endl;
        tl->Draw();
        tl->SetFillColor(0);
        tl->SetTextSize(0.04);
        c->Print(pdf_name, TString::Format("Title:%d GeV", int(energy[i])));
    }
    c->Print(pdf_name + "]");
    return 1;
}