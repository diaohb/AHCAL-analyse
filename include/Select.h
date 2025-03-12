#ifndef SELECT_HH
#define SELECT_HH
#include <vector>
#include <map>
#include <unordered_map>
#include <TLinearFitter.h>
#include <TF1.h>
using namespace std;
class Select{
    public:
        Select(vector<double>* _Hit_X,vector<double>* _Hit_Y,vector<double>* _Hit_Z, vector<double> *hit_energy,double a);
        ~Select();

        int Result(double &total_energy,vector<double>* &layer_energy,vector<double>* &layer_hitno,int &_hitlayer,int &_hitno,int &_big_hitno);
        double GetZeta(){return zeta;}
        double GetRMS(){return t_rms;}
        double GetShower_Start() { return shower_start; }
        double GetShower_End() { return shower_end; }
        double GetE_Hit() { return E_hit; }
        double GetFD() { return fd; }
        void GetR_alpha(double *_R)
        {
            for (int i = 0; i < 10; i++)
                _R[i] = R[i];
        }
        void GetHitEnergy(vector<int>* &cellid,vector<double>* &hitenergy)
        {
            cellid = CellID;
            // cout << cellid->size() << endl;
            // cout << Hit_energy->size()<<"  "<<hitenergy->capacity() << endl;
            hitenergy = Hit_energy;
            // hitenergy->assign(Hit_energy->begin(), Hit_energy->end());
        }
        void GetCenter(vector<double>* &layer_hit_X,vector<double>* &layer_hit_Y){layer_hit_X=CenterX;layer_hit_Y=CenterY;}
        void GetMax(vector<double>* &layer_max_energy,vector<double>* &layer_max_x,vector<double>* &layer_max_y){layer_max_energy=max_energy;layer_max_x=max_x;layer_max_y=max_y;}

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

        const double MIP_E=0.461; //MeV

        vector<int> *CellID=nullptr;
        vector<double> *Hit_X=nullptr;
        vector<double> *Hit_Y=nullptr;
        vector<double> *Hit_Z=nullptr;
        vector<double> *Hit_energy=nullptr;
        vector<double> *CenterX=nullptr;
        vector<double> *CenterY=nullptr;
        vector<double> *max_energy=nullptr;
        vector<double> *max_x=nullptr;
        vector<double> *max_y=nullptr;
        double lenergy[40]={0};
        double rmsX[40]={0};
        double rmsY[40]={0};
        double rms[40]={0};
        double zeta=0;
        double t_rms=0;
        double E_hit=0;
        double alpha[10]={0};
        double R[10] = {0};
        double fd = 0;
        bool ishit[40] = {0};
        int layerhitno[40]={0};
        int hitlayer=0;
        int hitno=0;
        int true_hitno=0;
        int big_hitno=0;
        int shower_start=0;
        int shower_end=0;
        double tenergy=0;
        double beam_energy = 0;
        // TLinearFitter *fitter = new TLinearFitter();
        unordered_map<int, bool> m_hit;
        TF1 *fhitnocut=new TF1("fhitnocut","372.285*(1-exp(-0.0064584*x))",0,350);
        // TF1 *fhitlayercut=new TF1("fhitlayercut","19.5*(1-exp(-0.0576*x))",0,350);
        unordered_map<double,int> hitnocut_down={{0.5,4},{1,7},{2,10},{3,12},{4,15},{5,20},{10,20},{20,60},{30,80},{40,100},{50,110},{60,120},{70,140},{80,150},{100,170},{120,200},{150,230},{250,320}};
        unordered_map<double,int> hitnocut_up={{0.5,14},{1,22},{2,33},{3,40},{4,48},{5,50},{10,20},{20,60},{30,80},{40,100},{50,110},{60,120},{70,140},{80,150},{100,170},{120,200},{150,230},{250,320}};
        unordered_map<double,int> hitlayercut_down={{0.5,4},{1,6},{2,8},{3,9},{4,10},{5,11},{10,10},{20,15},{30,15},{40,16},{50,17},{60,18},{70,18},{80,20},{100,20},{120,20},{150,20},{250,22}};
        unordered_map<double,int> hitlayercut_up={{0.5,10},{1,13},{2,16},{3,18},{4,19},{5,20},{10,10},{20,15},{30,15},{40,16},{50,17},{60,18},{70,18},{80,20},{100,20},{120,20},{150,20},{250,22}};
        
};


#endif