#include "/home/diaohb/plotstyle/plotstyle.c"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <ctime>
#include <fstream>
#include <iostream>
#include <time.h>
using namespace std;
int temperature() {
    SetStyle();
    gStyle->SetCanvasPreferGL(true);// 启用OpenGL渲染
    ifstream in("temp.txt");
    string time, str1, str2;
    double temp = 0, humi = 0;
    tm tm;
    time_t tmp_t;
    TGraph *gr1 = new TGraph();
    TGraph *gr2 = new TGraph();
    while (1) {
        in >> str1 >> str2 >> temp >> humi;
        if (in.eof()) break;
        time = str1 + " " + str2;
        strptime(time.c_str(), "%Y-%m-%d %H:%M:%S", &tm);
        // cout << tm.tm_year << endl;
        tmp_t = mktime(&tm);
        // cout << tmp_t << "  " << ctime(&tmp_t) << endl;
        // cout << tmp_t << "  " << temp << "  " << humi << endl;
        gr1->AddPoint(tmp_t, temp);
        gr2->AddPoint(tmp_t, humi);
    }
    TCanvas *c = new TCanvas();
    // c->Divide(1, 2);
    TPad *pad1 = new TPad();
    pad1->SetFillStyle(0);
    pad1->SetFillColorAlpha(0, 0);
    pad1->Draw();
    pad1->SetGrid(1, 1);
    pad1->cd();
    FormatData(gr1);
    NameAxis(gr1, "time", "temperature #circ C");
    gr1->GetXaxis()->SetTimeDisplay(true);
    gr1->GetXaxis()->SetTimeFormat("%d-%H");
    gr1->Draw("ap");
    TPad *pad2 = new TPad();
    pad2->SetFillColorAlpha(0, 0);
    pad2->SetFillStyle(0);
    pad2->Draw();
    pad2->cd();
    FormatData(gr2, kRed, 20);
    NameAxis(gr2, "time", "humidity");
    gr2->GetXaxis()->SetTimeDisplay(true);
    gr2->GetXaxis()->SetTimeFormat("%d-%H");
    gr2->GetYaxis()->SetAxisColor(kRed);
    gr2->GetYaxis()->SetLabelColor(kRed);
    gr2->Draw("ap Y+ ");
    c->Print("temperatur.pdf");
    return 1;
}