#ifndef SELECT_HH
#define SELECT_HH
#include <TF1.h>
#include <TLinearFitter.h>
#include <unordered_map>
#include <vector>
using namespace std;
class Select {
public:
    Select();
    ~Select();

    void Delete();
    int Result(double &total_energy, vector<double> *&layer_energy, vector<double> *&layer_hitno, int &_hitlayer, int &_hitno, int &_big_hitno);
    double GetZeta() { return zeta; }
    double GetRMS() { return t_rms; }
    double GetShower_Start() { return shower_start; }
    double GetShower_End() { return shower_end; }
    double GetE_Hit() { return E_hit; }
    double GetFD() { return fd; }
    void GetR_alpha(double *_R) {
        for (int i = 0; i < 10; i++)
            _R[i] = R[i];
    }
    void GetHitEnergy(vector<int> *&cellid, vector<double> *&hitenergy) {
        cellid->assign(CellID.begin(), CellID.end());
        // cout << cellid->size() << endl;
        // cout << Hit_Energy->size()<<"  "<<hitenergy->capacity() << endl;
        hitenergy->assign(Hit_Energy.begin(), Hit_Energy.end());
    }
    void GetCenter(vector<double> *&layer_hit_X, vector<double> *&layer_hit_Y) {
        layer_hit_X->assign(CenterX.begin(), CenterX.end());
        layer_hit_Y->assign(CenterY.begin(), CenterY.end());
    }
    void GetMax(vector<double> *&layer_max_energy, vector<double> *&layer_max_x, vector<double> *&layer_max_y, double &_shower_max) {
        layer_max_energy->assign(max_energy.begin(), max_energy.end());
        layer_max_x->assign(max_x.begin(), max_x.end());
        layer_max_y->assign(max_y.begin(), max_y.end());
        _shower_max = shower_max;
    }
    void GetLayerRMS(vector<double> *&layer_rms) {
        layer_rms->assign(rms.begin(), rms.end());
    }
    void ResetEvent(vector<double> *_Hit_X, vector<double> *_Hit_Y, vector<double> *_Hit_Z, vector<double> *_Hit_Energy, double a, int _particleID = 1);

private:
    void Init();
    void EnergyCenter();
    void RMS_zeta();
    void MaxEnergy();
    bool Isstraight();
    bool Ismip();
    bool Ishadron();
    void E_Hit();
    void FD();
    double R_alpha(int alpha);

    const double MIP_E = 0.461;//MeV
    int particleID = 1;
    int range[2] = {100, 360};
    int shower_start_ref[2] = {10, 10};
    int shower_end_ref[2] = {25, 35};

    vector<int> CellID;
    vector<double> Hit_X;
    vector<double> Hit_Y;
    vector<double> Hit_Z;
    vector<double> Hit_Energy;
    vector<double> CenterX = vector<double>(40, 0);
    vector<double> CenterY = vector<double>(40, 0);
    vector<double> max_energy = vector<double>(40, 0);
    vector<double> max_x = vector<double>(40, 0);
    vector<double> max_y = vector<double>(40, 0);
    vector<double> rms = vector<double>(40, 0);
    double lenergy[40] = {0};
    double rmsX[40] = {0};
    double rmsY[40] = {0};
    double zeta = 0;
    double t_rms = 0;
    double E_hit = 0;
    double alpha[10] = {0};
    double R[10] = {0};
    double fd = 0;
    bool ishit[40] = {0};
    int layerhitno[40] = {0};
    int hitlayer = 0;
    int hitno = 0;
    int true_hitno = 0;
    int big_hitno = 0;
    int shower_start = 0;
    int shower_end = 0;
    double shower_max = 0;
    double tenergy = 0;
    double beam_energy = 0;
    // TLinearFitter *fitter = new TLinearFitter();
    unordered_map<int, bool> m_hit;
    TF1 *fhitnocut = new TF1("fhitnocut", "372.285*(1-exp(-0.0064584*x))", 0, 350);
    // TF1 *fhitlayercut=new TF1("fhitlayercut","19.5*(1-exp(-0.0576*x))",0,350);
    unordered_map<double, int> hitnocut_down[2] = {
            {{0.5, 4}, {1, 4}, {2, 10}, {3, 12}, {4, 15}, {5, 19}, {6, 22}, {7, 25}, {8, 27}, {10, 30}, {12, 36}, {15, 42}, {20, 60}, {30, 80}, {40, 100}, {50, 110}, {60, 120}, {70, 140}, {80, 150}, {100, 170}, {120, 200}, {150, 230}, {250, 320}},
            {{1, 3}, {2, 8}, {3, 15}, {4, 20}, {5, 28}, {6, 33}, {7, 38}, {8, 40}, {10, 40}, {12, 50}, {15, 60}}};
    unordered_map<double, int> hitnocut_up[2] = {
            {{0.5, 14}, {1, 20}, {2, 30}, {3, 40}, {4, 44}, {5, 50}, {6, 55}, {7, 58}, {8, 60}, {10, 70}, {12, 74}, {15, 82}, {20, 60}, {30, 80}, {40, 100}, {50, 110}, {60, 120}, {70, 140}, {80, 150}, {100, 170}, {120, 200}, {150, 230}, {250, 320}},
            {{1, 31}, {2, 37}, {3, 50}, {4, 60}, {5, 70}, {6, 85}, {7, 95}, {8, 105}, {10, 120}, {12, 140}, {15, 160}}};
    unordered_map<double, int> hitlayercut_down[2] = {
            {{0.5, 4}, {1, 4}, {2, 7}, {3, 8}, {4, 9}, {5, 10}, {6, 11}, {7, 11}, {8, 11}, {10, 12}, {12, 13}, {15, 13}, {20, 15}, {30, 15}, {40, 16}, {50, 17}, {60, 18}, {70, 18}, {80, 20}, {100, 20}, {120, 20}, {150, 20}, {250, 22}},
            {{1, 2}, {2, 4}, {3, 6}, {4, 8}, {5, 10}, {6, 11}, {7, 12}, {8, 13}, {10, 14}, {12, 15}, {15, 16}}};
    unordered_map<double, int> hitlayercut_up[2] = {
            {{0.5, 10}, {1, 13}, {2, 16}, {3, 18}, {4, 18}, {5, 20}, {6, 20}, {7, 21}, {8, 22}, {10, 22}, {12, 23}, {15, 24}, {20, 15}, {30, 15}, {40, 16}, {50, 17}, {60, 18}, {70, 18}, {80, 20}, {100, 20}, {120, 20}, {150, 20}, {250, 22}},
            {{1, 37}, {2, 37}, {3, 37}, {4, 37}, {5, 37}, {6, 37}, {7, 37}, {8, 37}, {10, 37}, {12, 37}, {15, 37}}};
};


#endif