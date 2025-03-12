#ifndef PID_ANALYSE_HH
#define PID_ANALYSE_HH
#include <vector>
#include <map>
#include <unordered_map>
#include <TLinearFitter.h>
#include <TF1.h>
using namespace std;
class PID_analyse{
    public:
        PID_analyse(vector<double>* _Hit_X,vector<double>* _Hit_Y,vector<double>* _Hit_Z, vector<double> *hit_energy,double a);
        ~PID_analyse();

        double GetZeta(){return zeta;}
        double GetRMS(){return t_rms;}
        double GetShower_Start(){return shower_start;}
        double GetE_Hit(){return E_hit;}
        double GetFD(){return fd;}
        void GetR_alpha(double* _R){for(int i=0;i<10;i++)_R[i]=R[i];}
        void GetCenter(vector<double>* &layer_hit_X,vector<double>* &layer_hit_Y){layer_hit_X=CenterX;layer_hit_Y=CenterY;}
        void GetMax(vector<double>* &layer_max_energy,vector<double>* &layer_max_x,vector<double>* &layer_max_y){layer_max_energy=max_energy;layer_max_x=max_x;layer_max_y=max_y;}

    private:
        void Init();
        void EnergyCenter();
        void RMS_zeta();
        void MaxEnergy();
        void E_Hit();
        void FD();
        double R_alpha(int alpha);

        const double MIP_E=0.461; //MeV
        vector<int> *CellID=new vector <int>;
        vector<double>* Hit_X=new vector<double>;
        vector<double> *Hit_Y = new vector<double>;
        vector<double> *Hit_Z = new vector<double>;
        vector<double> *Hit_energy = new vector<double>;
        vector<double> *CenterX = new vector<double>(40,0);
        vector<double> *CenterY = new vector<double>(40, 0);
        vector<double> *max_energy = new vector<double>(40, 0);
        vector<double> *max_x = new vector<double>(40, -500);
        vector<double> *max_y = new vector<double>(40, -500);
        double lenergy[40]={0};
        double rmsX[40]={0};
        double rmsY[40]={0};
        double rms[40]={0};
        double zeta=0;
        double t_rms=0;
        double E_hit=0;
        double alpha[10]={0};
        double R[10] = {0};
        bool ishit[40] = {0};
        int layerhitno[40]={0};
        int hitlayer=0;
        int hitno=0;
        int big_hitno=0;
        int shower_start=0;
        int shower_end=0;
        double tenergy=0;
        double beam_energy=0;
        double fd = 0;
        TLinearFitter *fitter = new TLinearFitter();
        unordered_map<int, bool> m_hit;
};


#endif