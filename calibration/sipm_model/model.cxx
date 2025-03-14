#define info(X) std::cout << #X << "  " << (X) << std::endl;
#include <vector>
#include <unordered_map>
#include "TRandom3.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
using namespace std;
class Sim
{
public:
    Sim(int pixel, double p = 1, double x = 0.05); //{n_pixel=pixel;pde=p;xt=x;}
    ~Sim();
    void generate_pixel();
    void run(int photon, int entries);
    int ava(int npx, int npy);
    int n_pixel = 0;
    int n_photon = 15000;
    double pde = 1;
    double xt = 0;
    int n_fired = 0;
    int length = 0;
    ;
    int area = 0;
    vector<int> n_fired_series;
    TRandom *ran = new TRandom3(0);
    vector<vector<bool>> isp, isa;
    unordered_map<int, TH1I *> h_fired;
    TF1 *fcor = new TF1("fcor", "[0]*(1-e^(-x*[1]/[0]))", 0, 15000);
};
Sim::Sim(int pixel, double p = 1, double x = 0.05)
{
    n_pixel = pixel;
    pde = p;
    xt = x;
    fcor->SetParameters(n_pixel, pde);
    length = ceil(sqrt(n_pixel));
    area = length * length;
    vector<bool> tmp(length, 0);
    isp.resize(length, tmp);
    isa.resize(length, tmp);
}
Sim::~Sim()
{
}
void Sim::generate_pixel()
{
    double p;
    int n = n_pixel, m = area;
    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < length; j++)
        {
            p = gRandom->Uniform(0, 1);
            if (p < n / double(m))
            {
                isp[i][j] = 1;
                m--;
                n--;
            }
            else
            {
                isp[i][j] = 0;
                m--;
            }
        }
    }
}
void Sim::run(int photon, int entries)
{
    n_photon = photon;
    vector<bool> tmp(length, 0);
    for (int n = 0; n < n_photon; n++)
    {
        TString name = "incident_" + TString(std::to_string(n + 1).c_str()) + "_photons";
        double m = fcor->Eval(n + 1);
        int fs = m - 5 * sqrt(m), fe = ceil(m + 5 * sqrt(m));
        if (fs < 0)
            fs = 0;
        h_fired[n] = new TH1I(name, name, fe - fs, fs, fe);
    }
    for (int n = 0; n < entries; n++)
    {
        if ((n + 1) % 10000 == 0)
            printf("fill histogram %d\n", n + 1);
        std::fill(isa.begin(), isa.end(), tmp);
        n_fired = 0;
        for (int i = 0; i < n_photon; i++)
        {
            // n_fired_series.push_back(n_fired);
            if (pde < 1 && ran->Uniform(0, 1) > pde)
                continue;
            int npx, npy;
            do
            {
                npx = ran->Integer(86);
                npy = ran->Integer(86);
            } while (!isp[npx][npy]);
            ava(npx, npy);
            h_fired[i]->Fill(n_fired);
        }
    }
}
int Sim::ava(int npx, int npy)
{
    if (isa[npx][npy])
        return 0;
    isa[npx][npy] = 1;
    n_fired++;
    if (xt == 0)
        return 1;
    if (npy >= 1 && isp[npx][npy - 1] && ran->Uniform(0, 1) < xt)
    {
        ava(npx, npy - 1);
    }
    if (npx >= 1 && isp[npx - 1][npy] && ran->Uniform(0, 1) < xt)
    {
        ava(npx - 1, npy);
    }
    if (npy <= length - 2 && isp[npx][npy + 1] && ran->Uniform(0, 1) < xt)
    {
        ava(npx, npy + 1);
    }
    if (npx <= length - 2 && isp[npx + 1][npy] && ran->Uniform(0, 1) < xt)
    {
        ava(npx + 1, npy);
    }
    return 1;
}
int model()
{
    int n_pixel = 7284;
    double pde = 1;
    double xt = 0.00;
    int n_photons = 15000;
    int entries = 100000;
    Sim sim(n_pixel, pde, xt);
    sim.generate_pixel();
    sim.run(n_photons, entries);
    TFile *fout = TFile::Open("sipm_model_0.0xt.root", "RECREATE");
    TTree *tout = new TTree("sipm_fire", "sipm_fire");
    int n_photon = 0;
    double mean = 0, sigma = 0;
    tout->Branch("n_photon", &n_photon);
    tout->Branch("mean", &mean);
    tout->Branch("sigma", &sigma);
    fout->mkdir("hist");
    fout->cd("hist");
    unordered_map<int, TH1I *> h = sim.h_fired;
    TH1I *htmp;
    for (int i = 0; i < n_photons; i++)
    {
        htmp = (TH1I *)h[i]->Clone();
        htmp->Write();
        n_photon = i;
        mean = htmp->GetMean();
        sigma = htmp->GetRMS();
        tout->Fill();
    }
    fout->cd();
    tout->Write();
    fout->Close();
    return 1;
}
