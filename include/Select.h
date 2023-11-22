#ifndef SELECT_HH
#define SELECT_HH
#include <vector>
#include <map>
#include <unordered_map>
#include <TF1.h>
using namespace std;
class Select{
    public:
        Select(vector<double>* _Hit_X,vector<double>* _Hit_Y,vector<double>* _Hit_Z, vector<double> *hit_energy,double a);
        ~Select();

        int Result(double &total_energy,vector<double>* &layer_energy,int &_hitlayer,int &_hitno);
        double GetZeta(){return zeta;}
        double GetRMS(){return t_rms;}

    private:
        void Init();
        void Layerd();
        void EnergyCenter();
        void RMS_zeta();
        void MaxEnergy();
        bool Isstraight();
        bool Ismip();

        const double MIP_E=0.461; //MeV
        int noenergy=0;
        vector<int> *CellID;
        vector<double>* Hit_X;
        vector<double>* Hit_Y;
        vector<double>* Hit_Z;
        vector<double> *Hit_energy;
        vector<vector<int>> layerd_cellid=vector<vector<int>>(40);
        vector<vector<double>> layerd_hit_energy=vector<vector<double>>(40);
        double CenterX[40]={0};
        double CenterY[40]={0};
        int maxenergycellid[40];
        int maxchannel[40];
        int maxchip[40];
        double maxenergyX[40]={0};
        double maxenergyY[40]={0};
        double maxenergy[40]={0};
        double lenergy[40]={0};
        double rmsX[40]={0};
        double rmsY[40]={0};
        double rms[40]={0};
        double zeta=0;
        double t_rms=0;

        bool ishit[40]={0};
        int layerhitno[40]={0};
        int hitlayer=0;
        int hitno=0;
        int shower_start=0;
        int shower_end=39;
        int last_layer=39;
        vector<int> noenergylayer;
        vector<double> layerEnergy;
        double tenergy=0;
        double beam_energy=0;
        TF1 *fhitnocut=new TF1("fhitnocut","372.285*(1-exp(-0.0064584*x))",0,350);
        // TF1 *fhitlayercut=new TF1("fhitlayercut","19.5*(1-exp(-0.0576*x))",0,350);
        // unordered_map<double,int> hitnocut={{10,20},{20,60},{30,80},{40,100},{50,110},{60,120},{70,140},{80,150},{100,170},{120,200},{150,230},{250,320}};
        // unordered_map<double,int> hitlayercut={{10,10},{20,15},{30,15},{40,16},{50,17},{60,18},{70,18},{80,20},{100,20},{120,20},{150,20},{250,22}};
        
};


#endif