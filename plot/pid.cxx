#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TH2.h"
#include "iostream"
using namespace std;
int pid() {
    const int number = 11;
    // double energy[6] = {0.5, 1, 2, 3, 4, 5};
    // int rms_limit[6] = {100, 200, 300, 400, 400, 500};
    // int hitno_limit[6] = {15, 25, 35, 45, 55, 60};
    int rms_limit[number] = {300, 600, 1000, 1200, 1500, 1700, 2000, 2200, 2500, 2700, 3000};
    double energy[number] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15};
    int hitno_limit[number] = {40, 50, 60, 70, 80, 100, 110, 120, 140, 160, 200};
    double range2345[number] = {30, 60, 80, 90, 110, 120, 130, 150, 170, 190, 220};
    ROOT::RDataFrame df("tanalyse", "/home/diaohb/CEPC/ihep/PS/analyse/intercept/pi-list_analyse.root");
    // ROOT::RDataFrame df("tanalyse", "/home/diaohb/CEPC/AHCAL-simulation/run/0423_30080126_trigger40/list_PS/calib_analyse.root");
    // ROOT::RDataFrame df("tanalyse", "/home/diaohb/CEPC/AHCAL-simulation/run/0422_30080126_trigger40/list_PS_birks/calib_analyse.root");
    TString pdf_name = "pi-figure/PID_cherenkov.pdf";
    TH2D *h1[number];
    TH2D *h8[number];
    TH2D *h9[number];
    TH1D *h2[number];
    TH1D *h3[number];
    TH1D *h4[number];
    TH1D *h5[number];
    TH1D *h6[number];
    TH1D *h7[number];
    for (int i = 0; i < number; i++) {
        h1[i] = new TH2D("pid", TString::Format("%d GeV;E_Hit;FD", int(energy[i])), 200, 0, 10, 200, 0, 1);
        h8[i] = new TH2D("pid", TString::Format("%d GeV;E_Hit;layer2345", int(energy[i])), 200, 0, 10, 200, 0, range2345[i]);
        h9[i] = new TH2D("pid", TString::Format("%d GeV;hitno;layer2345", int(energy[i])), hitno_limit[i], 0, hitno_limit[i], 200, 0, range2345[i]);
        h2[i] = new TH1D("FD", TString::Format("%d GeV;FD;", int(energy[i])), 200, 0, 1);
        h3[i] = new TH1D("E_Hit", TString::Format("%d GeV;E_Hit;", int(energy[i])), 200, 0, 10);
        h4[i] = new TH1D("RMS", TString::Format("%d GeV;RMS;", int(energy[i])), 200, 0, rms_limit[i] * 1.2);
        h5[i] = new TH1D("hitno", TString::Format("%d GeV;hitno;", int(energy[i])), hitno_limit[i], 0, hitno_limit[i]);
        h6[i] = new TH1D("hitlayer", TString::Format("%d GeV;hitlayer;", int(energy[i])), 40, 0, 40);
        h7[i] = new TH1D("2345", TString::Format("%d GeV;layer 2345;", int(energy[i])), 200, 0, range2345[i]);
        // cout << h6[i] << endl;
    }
    df.Filter("select_flag==6&&cherenkov!=1").Foreach([&](double beam_energy, double FD, double E_Hit, double rms, int hitno, int hitlayer, vector<double> layer_energy, double total_energy) {
        double *pos = find(energy, energy + number, beam_energy);
        int i = pos - energy;
        if (i < number) {
            h1[i]->Fill(E_Hit, FD);
            h2[i]->Fill(FD);
            h3[i]->Fill(E_Hit);
            h4[i]->Fill(rms);
            h5[i]->Fill(hitno);
            h6[i]->Fill(hitlayer);
            double energy2345 = 0;
            for (int j = 2; j <= 5; j++) {
                energy2345 += layer_energy.at(j);
            }
            h7[i]->Fill(energy2345);
            h8[i]->Fill(E_Hit, energy2345);
            h9[i]->Fill(hitno, energy2345);
        }
    },
                                                      {"beam_energy", "FD", "E_Hit", "rms", "hitno", "hitlayer", "layer_energy", "total_energy"});
    TCanvas *c = new TCanvas("c", "c", 2400, 2400);
    c->Print(pdf_name + "[");
    for (int i = 0; i < number; i++) {
        c->Clear();
        c->Divide(3, 3);
        c->cd(1);
        h2[i]->Draw();
        c->cd(2);
        gPad->SetLogz(1);
        h1[i]->Draw("colz");
        c->cd(3);
        h5[i]->Draw();
        c->cd(4);
        gPad->SetLogy(1);
        h4[i]->Draw();
        c->cd(5);
        h3[i]->Draw();
        c->cd(6);
        h6[i]->Draw();
        c->cd(7);
        gPad->SetLogy(1);
        h7[i]->Draw();
        c->cd(8);
        gPad->SetLogz(1);
        h8[i]->Draw("colz");
        c->cd(9);
        gPad->SetLogz(1);
        h9[i]->Draw("colz");
        c->Print(pdf_name, TString::Format("Title:%d GeV", int(energy[i])));
    }
    c->Print(pdf_name + "]");
    return 1;
}