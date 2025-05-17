#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStyle.h"
#include </home/diaohb/plotstyle/plotstyle.c>
#include <fstream>
#include <iostream>
using namespace std;
int energy() {
    SetStyle();
    gStyle->SetOptTitle(1);
    // double energy[11]={0.5,1,2,3,4,5,10,20,30,40,50};
    // double energy_data[11]={7.135,15.68,34.43,52.55,70.72,88.68,190.3,370.6,538.5,690.5,827.6};
    // double sigma_data[11]={2.977,4.206,5.693,7.123,8.33,9.734,13.55,19.23,23.89,27.9,31.39};
    // double energy_MC[11]={8.74,18.22,37.32,56.33,75.3,94.21,186.9,367.8,542.4,710.4,872.7};
    // double sigma_MC[11]={2.553,3.816,5.55,6.862,8.008,8.945,12.92,17.82,22.72,26.58,32.37};
    TCanvas *c = new TCanvas("c", "c", 800, 1200);
    TGraph *gr1 = new TGraph();
    TGraph *gr2 = new TGraph();
    TGraph *gr3 = new TGraph();
    TGraph *grframe = new TGraph();

    ifstream in("pi-txt/energy_PS_3008_ehit.txt", ios::in);
    TString outfile = "pi-figure/energy_PS_3008_ehit.pdf";
    while (!in.eof()) {
        double energy = 0, mean = 0, sigma = 0;
        // string mode;
        int mode;
        in >> mode >> energy >> mean >> sigma;
        if (in.fail()) break;
        if (energy == 0) energy = 0.5;
        if (mode == 0) gr1->AddPoint(energy, mean);
        if (mode == 1) gr2->AddPoint(energy, mean);
        if (mode == 2) gr3->AddPoint(energy, mean);
    }
    // for(int i=0;i<11;i++){
    //     double tmp=sigma_data[i]/energy_data[i];
    //     cout<<tmp<<endl;
    //     gr1->AddPoint(energy[i],100*tmp);
    //     tmp=sigma_MC[i]/energy_MC[i];
    //     cout<<tmp<<endl;
    //     gr2->AddPoint(energy[i],100*tmp);
    // }
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1);
    pad1->SetBottomMargin(0.4);
    pad1->Draw();
    pad1->SetFillStyle(0);
    pad1->cd();

    grframe->AddPoint(0, 0);
    grframe->AddPoint(100, 1000);
    grframe->GetXaxis()->SetTitle("Beam Energy [GeV]");
    grframe->GetYaxis()->SetTitle("mean [MeV]");
    grframe->GetXaxis()->CenterTitle();
    grframe->GetYaxis()->CenterTitle();
    grframe->SetTitle("Energy Linearity");
    grframe->GetXaxis()->SetRangeUser(0, 15);
    grframe->GetYaxis()->SetRangeUser(0, 300);
    grframe->GetYaxis()->SetMaxDigits(3);
    grframe->GetXaxis()->SetLabelSize(0);
    grframe->GetXaxis()->SetTitleSize(0);
    grframe->DrawClone("ap");

    TF1 *fun1 = new TF1("fun1", "[0]+[1]*x", 0, 100);
    // TF1* fun1=new TF1("fun1","(([0])^2/x+[1]^2)^0.5",0,100);
    // TF1* fun12=new TF1("fun12","(([0])^2/x+[1]^2)^0.5",0,100);
    fun1->SetNpx(10000);
    // fun1->SetParameters(0.23,0.002);
    // fun12->SetParameters(0.25,0.01);
    gr1->Fit(fun1, "0", "", 2, 12);
    TString s1, s2, s3;
    s1.Form("data: %.2f x+%.2f", fun1->GetParameter(1), fun1->GetParameter(0));
    // s1.Form("data: #frac{%.3f}{#sqrt{E}}#oplus%.3f",fun1->GetParameter(0),fun1->GetParameter(1));
    FormatData(gr1, 2, 20);
    gr1->SetTitle(s1);
    // gr1->SetMarkerSize(0.8);
    gr1->DrawClone("p");
    fun1->DrawCopy("lsame");
    for (Int_t i = 0; i < gr1->GetN(); i++) {
        gr1->SetPointY(i, -100 + 100 * gr1->GetPointY(i) / fun1->Eval(gr1->GetPointX(i)));
    }

    gr2->Fit(fun1, "0", "", 2, 12);
    s2.Form("digi: %.2f x+%.2f", fun1->GetParameter(1), fun1->GetParameter(0));
    // s2.Form("MC: #frac{%.3f}{#sqrt{E}}#oplus%.3f",fun12->GetParameter(0),fun12->GetParameter(1));
    FormatData(gr2, 1, 21);
    gr2->SetTitle(s2);
    // gr2->SetMarkerSize(0.8);
    gr2->DrawClone("p");
    fun1->SetLineColor(1);
    fun1->DrawCopy("lsame");
    for (Int_t i = 0; i < gr2->GetN(); i++) {
        gr2->SetPointY(i, -100 + 100 * gr2->GetPointY(i) / fun1->Eval(gr2->GetPointX(i)));
    }

    gr3->Fit(fun1, "0", "", 2, 12);
    s3.Form("truth: %.2f x+%.2f", fun1->GetParameter(1), fun1->GetParameter(0));
    // s2.Form("MC: #frac{%.3f}{#sqrt{E}}#oplus%.3f",fun12->GetParameter(0),fun12->GetParameter(1));
    FormatData(gr3, 3, 21);
    gr3->SetTitle(s3);
    // gr2->SetMarkerSize(0.8);
    gr3->DrawClone("p");
    fun1->SetLineColor(3);
    fun1->DrawCopy("lsame");
    for (Int_t i = 0; i < gr3->GetN(); i++) {
        gr3->SetPointY(i, -100 + 100 * gr3->GetPointY(i) / fun1->Eval(gr3->GetPointX(i)));
    }

    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 1);
    pad2->SetTopMargin(0.6);
    pad2->Draw();
    pad2->SetFillStyle(0);
    pad2->cd();
    grframe->GetXaxis()->SetLabelSize(0.06);
    grframe->GetYaxis()->SetNdivisions(505);
    grframe->GetYaxis()->SetTitle("residual [%]");
    grframe->GetYaxis()->SetRangeUser(-10, 10);
    grframe->GetXaxis()->SetTitleSize(0.07);
    grframe->Draw("ap");
    gr1->Draw("p");
    gr2->Draw("p");
    gr3->Draw("p");

    TLegend *tl = new TLegend(0.21, 0.65, 0.5, 0.84);
    tl->AddEntry(gr1);
    tl->AddEntry(gr2);
    tl->AddEntry(gr3);
    tl->Draw();
    gStyle->SetLegendTextSize(0.04);
    c->SaveAs(outfile);

    return 1;
}