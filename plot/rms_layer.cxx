#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH2.h"
#include "iostream"
using namespace std;
int rms_layer() {
    const int number = 6;
    double energy[6] = {0.5, 1, 2, 3, 4, 5};
    int rms_limit[6] = {100, 200, 300, 400, 400, 500};
    // int rms_limit[number] = {300, 600, 1000, 1200, 1500, 1700, 2000, 2200, 2500, 2700, 3000};
    // double energy[number] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15};
    ROOT::RDataFrame df("tanalyse", "/home/diaohb/CEPC/AHCAL-simulation/run/0405_30080126_trigger40/list_PS9/calib_analyse.root");
    // ROOT::RDataFrame df("tanalyse", "/home/diaohb/CEPC/ihep/PS/analyse/intercept/pi-list_analyse.root");
    TString pdf_name = "figure/MC_RMS.pdf";
    TH2D *h[number][40];
    for (int i = 0; i < number; i++) {
        for (int j = 0; j < 40; j++) {
            h[i][j] = new TH2D("rms lf", TString::Format("%d GeV %d layer;RMS;lf", int(energy[i]), j), 200, 0, rms_limit[i] * 0.5, 200, 0, 1);
        }
    }
    df.Filter("select_flag==6").Foreach([&](double beam_energy, double total_energy, vector<double> &layer_energy, vector<double> &layer_rms) {
        double *pos = find(energy, energy + number, beam_energy);
        int i = pos - energy;
        for (int k = 0; k < 40; k++) {
            double f = layer_energy.at(k) / total_energy;
            h[i][k]->Fill(layer_rms.at(k), f);
        }
    },
                                        {"beam_energy", "total_energy", "layer_energy", "layer_rms"});
    TCanvas *c = new TCanvas("c", "c", 4096, 6144);
    c->Print(pdf_name + "[");
    for (int i = 0; i < number; i++) {
        c->Clear();
        c->Divide(5, 6);
        for (int j = 0; j < 30; j++) {
            c->cd(j + 1);
            h[i][j]->Draw("colz");
        }
        c->Print(pdf_name, TString::Format("Title:%d GeV", int(energy[i])));
    }
    c->Print(pdf_name + "]");
    return 1;
}