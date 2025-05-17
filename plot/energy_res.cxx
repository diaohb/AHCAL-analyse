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
int energy_res() {
    SetStyle();
    SetPrelimStyle();
    gStyle->SetLabelSize(0.04, "xyz");
    gStyle->SetTitleSize(0.04, "xyz");
    gStyle->SetFrameLineWidth(3);
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetTitleOffset(1.2, "xyz");
    gStyle->SetOptTitle(1);
    double fit_range[2] = {2, 12};
    TCanvas *c = new TCanvas("c", "c", 900, 800);
    TGraph *gr1 = new TGraph();
    TGraph *gr2 = new TGraph();
    TGraph *gr3 = new TGraph();
    TGraph *gr11 = new TGraph();
    TGraph *gr21 = new TGraph();
    TGraph *gr31 = new TGraph();
    TGraph *grframe = new TGraph();

    ifstream in("pi-txt/energy_PS_3008_cherenkov.txt", ios::in);
    TString outfile = "pi-figure/energy_res_PS_3008_cherenkov.pdf";
    while (!in.eof()) {
        double energy = 0, mean = 0, sigma = 0;
        // string mode;
        int mode;
        in >> mode >> energy >> mean >> sigma;
        if (in.fail()) break;
        if (energy == 0)
            energy = 0.5;
        if (mode == 0)
            gr1->AddPoint(energy, 100 * sigma / mean);
        if (mode == 1) gr2->AddPoint(energy, 100 * sigma / mean);
        if (mode == 2) gr3->AddPoint(energy, 100 * sigma / mean);
    }

    // TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1);
    // pad1->SetBottomMargin(0.4);
    // pad1->Draw();
    // pad1->SetFillStyle(0);
    // pad1->SetGridx();
    // pad1->SetGridy();
    // pad1->cd();
    grframe->AddPoint(0, 0);
    grframe->AddPoint(120, 5000);
    NameAxis(grframe, "Beam Energy [GeV]", "Resolution [%]");
    FormatData(grframe, 0, 21);
    grframe->SetMarkerSize(0);
    // grframe->GetYaxis()->SetTitle("#sigma(E)/E(%)");
    // grframe->SetTitle("Energy Linearity");
    // grframe->SetTitle("Energy Resolution");
    grframe->GetXaxis()->SetRangeUser(0, 20);
    // grframe->GetXaxis()->SetLabelSize(0);
    grframe->GetYaxis()->SetRangeUser(0.00, 100);
    grframe->GetYaxis()->SetMaxDigits(3);
    grframe->Draw("ap");

    TF1 *fun1 = new TF1("fun1", "(([0])^2/x+[1]^2)^0.5", 0, 100);
    // TF1 *fun12 = new TF1("fun12", "(([0])^2/x+[1]^2)^0.5", 0, 100);
    fun1->SetNpx(10000);
    fun1->SetParameters(0.23, 0.002);
    gr1->Fit(fun1, "0", "", fit_range[0], fit_range[1]);
    TString s1, s2, s3;
    s1.Form("data: #frac{%.2f%%}{#sqrt{E}} #oplus %.2f%%", abs(fun1->GetParameter(0)), abs(fun1->GetParameter(1)));
    // s1.Form("data");
    FormatData(gr1, 2, 20);
    gr1->SetTitle(s1);
    gr1->Draw("p");
    fun1->DrawCopy("lsame");

    gr2->Fit(fun1, "0", "", fit_range[0], fit_range[1]);
    s2.Form("digi: #frac{%.2f%%}{#sqrt{E}} #oplus %.2f%%", abs(fun1->GetParameter(0)), abs(fun1->GetParameter(1)));
    // s2.Form("digi");
    FormatData(gr2, 1, 21);
    gr2->SetTitle(s2);
    gr2->Draw("p");
    fun1->SetLineColor(1);
    fun1->DrawCopy("lsame");

    gr3->Fit(fun1, "0", "", fit_range[0], fit_range[1]);
    s3.Form("truth: #frac{%.2f%%}{#sqrt{E}} #oplus %.2f%%", abs(fun1->GetParameter(0)), abs(fun1->GetParameter(1)));
    // s3.Form("truth");
    FormatData(gr3, 3, 21);
    gr3->SetTitle(s3);
    gr3->Draw("p");
    fun1->SetLineColor(3);
    fun1->DrawCopy("lsame");

    TLegend *tl = new TLegend(0.4, 0.55, 0.84, 0.84);
    tl->AddEntry(gr1);
    tl->AddEntry(gr2);
    tl->AddEntry(gr3);
    tl->Draw();
    gStyle->SetLegendTextSize(0.04);
    // TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 1);
    // pad2->SetGridx();
    // pad2->SetGridy();
    // pad2->SetTopMargin(0.6);
    // pad2->Draw();
    // pad2->SetFillStyle(0);
    // pad2->cd();
    // for (int i = 0; i < gr1->GetN(); i++)
    // {
    //     gr11->AddPoint(gr1->GetPointX(i), -100 + 100 * gr1->GetPointY(i) / gr2->GetPointY(i));
    // }
    // TGraph *grframe1 = new TGraph();
    // grframe1->AddPoint(0, -30);
    // grframe1->AddPoint(100, 30);
    // grframe1->GetXaxis()->SetTitle("Beam Energy(GeV)");
    // grframe1->GetYaxis()->SetTitle("#frac{mean-fit}{mean} (%)");
    // // TString ss;
    // // ss.Form("#mu {#pi}");
    // // cout << ss << endl;
    // grframe1->GetXaxis()->CenterTitle();
    // grframe1->GetYaxis()->CenterTitle();
    // // grframe1->GetYaxis()->SetNdivisions(510);

    // grframe1->SetTitle("Energy Linearity");
    // grframe1->GetXaxis()->SetRangeUser(0, 120);
    // grframe1->GetYaxis()->SetRangeUser(-20, 20);
    // grframe1->Draw("ap");

    // gr11->SetMarkerStyle(20);
    // gr11->SetMarkerSize(0.8);
    // gr11->SetMarkerColor(2);
    // // gr11->SetLineColor(kOrange);
    // gr11->SetLineColor(2);
    // gr11->Draw("pl");

    c->SaveAs(outfile);
    return 1;
}
