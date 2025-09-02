#include "Select.h"
#include "Global.h"
#include "TGraph.h"
#include <algorithm>
using namespace std;
Select::Select() {
    // std::cout<<"Selection start"<<std::endl;

    // Init();
}
Select::~Select() {

    // cout<<"Selection done"<<endl;
}
void Select::ResetEvent(vector<double> *_Hit_X, vector<double> *_Hit_Y, vector<double> *_Hit_Z, vector<double> *_Hit_Energy, double a, int _particleID) {
    Hit_X = *_Hit_X;
    Hit_Y = *_Hit_Y;
    Hit_Z = *_Hit_Z;
    Hit_Energy = *_Hit_Energy;
    beam_energy = a;
    particleID = _particleID;

    fill(CenterX.begin(), CenterX.end(), 0);
    fill(CenterY.begin(), CenterY.end(), 0);
    fill(max_energy.begin(), max_energy.end(), 0);
    fill(max_x.begin(), max_x.end(), 0);
    fill(max_y.begin(), max_y.end(), 0);
    fill(rms.begin(), rms.end(), 0);
    fill(lenergy, lenergy + 40, 0);
    fill(rmsX, rmsX + 40, 0);
    fill(rmsY, rmsY + 40, 0);
    fill(alpha, alpha + 10, 0);
    fill(R, R + 10, 0);
    fill(ishit, ishit + 40, 0);
    fill(layerhitno, layerhitno + 40, 0);
    // fill(iscluster, iscluster + 5 * 5 * 25, 0);
    is_hit.assign(40, vector<vector<int>>(9, vector<int>(36, 0)));
    is_visited.assign(40, vector<vector<bool>>(9, vector<bool>(36, false)));
    clusters.clear();
    cluster_size.clear();
    hitno = 0;
    hitlayer = 0;
    big_hitno = 0;
    true_hitno = 0;
    tenergy = 0;
    t_rms = 0;
    zeta = 0;
    fd = 0;
    E_hit = 0;
    shower_start = 0;
    shower_end = 0;
    shower_max = 0;
    track_chi2x = 0;
    track_chi2y = 0;
    cont_layers = 0;
    grx->Set(0);
    gry->Set(0);
    fitterx->ClearPoints();
    fittery->ClearPoints();
    Init();
}

void Select::Delete() {
}

void Select::Init() {
    // cout << Hit_X->size() << endl;
    for (int i = 0; i < Hit_X.size(); i++) {
        double e = Hit_Energy.at(i);
        if (e < 0.5 * MIP_E) {
            Hit_X.erase(Hit_X.begin() + i);
            Hit_Y.erase(Hit_Y.begin() + i);
            Hit_Z.erase(Hit_Z.begin() + i);
            Hit_Energy.erase(Hit_Energy.begin() + i);
            i--;
            if (e < -99) {
                big_hitno++;
            }
            continue;
        }
        // int chip = 0, channel = 0;
        // inverse(Hit_X.at(i), Hit_Y.at(i), chip, channel);
        int layer = Hit_Z.at(i) / 30;
        // int cellid = layer * 1e5 + chip * 1e4 + channel;
        // CellID.push_back(cellid);
        true_hitno++;
        layerhitno[layer]++;
        lenergy[layer] += e;
    }
    EnergyCenter();
    RMS_zeta();
    MaxEnergy();
    E_Hit();
    FD();
    FindTrack();
    // FindCluster();
    // EvalEnergy();
}

void Select::EnergyCenter() {
    vector<double> tmp_x(40, 0);
    vector<double> tmp_y(40, 0);
    for (int i = 0; i < Hit_X.size(); i++) {
        int layer = Hit_Z.at(i) / 30;
        double e = Hit_Energy.at(i);
        tmp_x.at(layer) += Hit_X.at(i) * e;
        tmp_y.at(layer) += Hit_Y.at(i) * e;
    }
    for (int i = 0; i < 40; i++) {
        if (lenergy[i] == 0) {
            tmp_x.at(i) = -500;
            tmp_y.at(i) = -500;
            continue;
        }
        tmp_x.at(i) /= lenergy[i];
        tmp_y.at(i) /= lenergy[i];
    }
    for (int k = 0; k < 1; k++) {
        hitno = 0;
        hitlayer = 0;
        fill(CenterX.begin(), CenterX.end(), 0);
        fill(CenterY.begin(), CenterY.end(), 0);
        fill(layerhitno, layerhitno + 40, 0);
        fill(lenergy, lenergy + 40, 0);
        fill(ishit, ishit + 40, 0);
        // double erasee = 0;
        for (int i = 0; i < Hit_X.size(); i++) {
            int layer = Hit_Z.at(i) / 30;
            int chip = 0, channel = 0;
            // int layer = Hit_Z.at(i) / 30;
            inverse(Hit_X.at(i), Hit_Y.at(i), chip, channel);
            int cellid = layer * 1e5 + chip * 1e4 + channel;
            double e = Hit_Energy.at(i);
            double x = Hit_X.at(i) - tmp_x.at(layer);
            double y = Hit_Y.at(i) - tmp_y.at(layer);
            // if (x < 60 && x > -60 && y < 60 && y > -60)
            int range_tmp = range[particleID - 1];
            if (x < range_tmp && x > -range_tmp && y < range_tmp && y > -range_tmp) {
                // cout << CenterX.at(layer) << endl;
                CenterX.at(layer) += Hit_X.at(i) * e;
                CenterY.at(layer) += Hit_Y.at(i) * e;
                CellID.push_back(cellid);
                hitno++;
                layerhitno[layer]++;
                lenergy[layer] += e;
            }
            is_hit[layer][chip][channel] = i + 1;
        }
        // cout << lenergy[0] << endl;
        for (int i = 0; i < 40; i++) {
            if (lenergy[i] > 0) {
                tenergy += lenergy[i];
                ishit[i] = 1;
                hitlayer++;
                CenterX.at(i) /= lenergy[i];
                CenterY.at(i) /= lenergy[i];
            } else {
                ishit[i] = 0;
                CenterX.at(i) = -500;
                CenterY.at(i) = -500;
            }
        }
        tmp_x = CenterX;
        tmp_y = CenterY;
    }
    for (int i = 0; i < 40; i++) {
        if (layerhitno[i] >= 4) {
            shower_start = i;
            break;
        }
    }
    for (int i = 37; i >= 0; i--) {
        if (layerhitno[i] >= 4) {
            shower_end = i;
            break;
        }
    }
    for (int i = 0; i < 40; i++) {
        if (layerhitno[i] == 0) {
            cont_layers = i;
            break;
        }
    }
}
void Select::RMS_zeta() {
    for (int i = 0; i < Hit_X.size(); i++) {
        int layer = Hit_Z.at(i) / 30;
        double e = Hit_Energy.at(i);
        rmsX[layer] +=
                (Hit_X.at(i) - CenterX.at(layer)) * (Hit_X.at(i) - CenterX.at(layer)) * e;
        rmsY[layer] +=
                (Hit_Y.at(i) - CenterY.at(layer)) * (Hit_Y.at(i) - CenterY.at(layer)) * e;
    }
    for (int i = 0; i < 40; i++) {
        if (ishit[i] == 0) {
            rmsX[i] = 0;
            rmsY[i] = 0;
            continue;
            // break;
        }
        // cout<<i<<"  "<<ishit[i]<<"  "<<lenergy[i]<<endl;
        rmsX[i] /= lenergy[i];
        rmsY[i] /= lenergy[i];
        rms[i] = sqrt(rmsX[i] + rmsY[i]);
        t_rms += rms[i];
    }
    // cout<<zeta<<"  ";
    zeta = lenergy[14] / tenergy * t_rms * t_rms * t_rms * t_rms;
    // cout<<lenergy[14]<<"  "<<lenergy[19]<<"  "<<tenergy<<"  "<<zeta<<endl;
}
void Select::MaxEnergy() {
    for (int i = 0; i < Hit_X.size(); i++) {
        int layer = Hit_Z.at(i) / 30;
        if (max_energy.at(layer) < Hit_Energy.at(i)) {
            max_energy.at(layer) = Hit_Energy.at(i);
            max_x.at(layer) = Hit_X.at(i);
            max_y.at(layer) = Hit_Y.at(i);
        }
    }
    shower_max = *max_element(max_energy.begin(), max_energy.end() - 2);
    // cout<<maxenergyX[0]<<"  "<<maxenergyY[0]<<endl;
}

int Select::Result(double &total_energy, vector<double> *&layer_energy,
                   vector<double> *&layer_hitno, int &_hitlayer, int &_hitno, int &_big_hitno) {
    int flag_7 = 0;
    for (int i = 0; i < 20; i++) {
        // cout << CenterX.at(i) << "  ";
        if (ishit[i] && (CenterX.at(i) > 180 || CenterX.at(i) < -180 || CenterY.at(i) > 180 || CenterY.at(i) < -180)) {
            flag_7++;
        }
    }
    // cout << flag_7 << endl;
    if (flag_7 > 5) {
        return 7;
    }
    _hitlayer = hitlayer;
    _hitno = 0;
    total_energy = 0;
    _big_hitno = big_hitno;
    layer_energy->clear();
    layer_hitno->clear();
    for (int i = 0; i < 40; i++) {
        layer_hitno->push_back(layerhitno[i]);
        layer_energy->push_back(lenergy[i]);
    }
    for (int i = 0; i < layer_range[particleID]; i++) {
        _hitno += layer_hitno->at(i);
        total_energy += layer_energy->at(i);
    }
    if (total_energy < 0.5 * MIP_E) {
        return 0;
    }
    if (shower_start >= shower_start_ref[particleID - 1]) {
        return 3;
    }
    if (shower_end >= shower_end_ref[particleID - 1]) {
        return 9;
    }
    // cout << particleID << "  " << beam_energy << "  " << hitnocut_down[particleID - 1][beam_energy] << "  " << _hitno << "  " << hitnocut_up[particleID - 1][beam_energy] << endl;
    if (_hitno < hitnocut_down[particleID - 1][beam_energy] || _hitno > hitnocut_up[particleID - 1][beam_energy]) {
        return 1;
    }
    if (_hitlayer < hitlayercut_down[particleID - 1][beam_energy] || _hitlayer > hitlayercut_up[particleID - 1][beam_energy]) {
        return 2;
    }
    if (Ismip()) {
        return 4;
    }
    // if (!Isstraight()) {
    //     return 5;
    // }
    if (layer_hitno->at(0) > 1 || layer_hitno->at(0) == 0) {
        return 8;
    }
    if (cont_layers < cont_layers_cut[particleID - 1][beam_energy]) {
        return 10;
    }
    return 6;
}

bool Select::Isstraight() {
    double maxx = -360, minx = 360, maxy = -360, miny = 360;
    for (int i = 0; i < 30; i++) {
        if (ishit[i] && lenergy[i] / tenergy > 0.05) {
            // if(maxenergyX[i]>maxx)maxx=maxenergyX[i];
            // if(maxenergyX[i]<minx)minx=maxenergyX[i];
            // if(maxenergyY[i]>maxx)maxy=maxenergyY[i];
            // if(maxenergyY[i]<minx)miny=maxenergyY[i];
            if (CenterX.at(i) > maxx)
                maxx = CenterX.at(i);
            if (CenterX.at(i) < minx)
                minx = CenterX.at(i);
            if (CenterY.at(i) > maxy)
                maxy = CenterY.at(i);
            if (CenterY.at(i) < miny)
                miny = CenterY.at(i);
        }
    }
    // cout<<maxx<<"  "<<minx<<"  "<<maxy<<"  "<<miny<<endl;
    if ((maxx - minx) < 200 && (maxy - miny) < 200) {
        return true;
    } else {
        return false;
    }
}

bool Select::Ismip() {
    double lastenergy = 0;
    for (int i = 30; i < 38; i++) {
        lastenergy += lenergy[i];
    }
    return (lastenergy / tenergy > 0.05);
    // int vlayer[40]={0};
    // for(int i=0;i<Hit_Z->size();i++){
    //     int layer=Hit_Z->at(i)/30;
    //     if(vlayer[layer]>2){
    //         return true;
    //     }
    //     vlayer[layer]++;
    // }
    // return false;
}
bool Select::Ishadron() {
    for (int i = 20; i < 40; i++) {
        if (layerhitno[i] >= 3) {
            return true;
        }
    }
    return false;
}
void Select::E_Hit() {
    E_hit = 0;
    for (int i = 0; i < Hit_Energy.size(); i++) {
        E_hit += Hit_Energy.at(i);
    }
    E_hit /= Hit_X.size();
}
void Select::FD() {
    // fitter->SetFormula("pol1");
    // fitter->StoreData(0);
    fd = 0;
    for (int i = 0; i < 9; i++) {
        alpha[i] = log(i + 2);
        R[i] = log(R_alpha(i + 2));
        fd += R[i] / alpha[i];
    }
    alpha[9] = log(20);
    R[9] = log(R_alpha(20));
    fd += R[9] / alpha[9];
    fd /= 10;
    // fitter->AssignData(20, 1, alpha, R);
    // fitter->Eval();
    // fd = fitter->GetParameter(1);
}
double Select::R_alpha(int _alpha) {
    double R_alpha = 0;
    m_hit.clear();
    double x = 0, y = 0;
    if (_alpha % 2 == 1) {
        x = 20, y = 20;
    }
    for (int i = 0; i < Hit_X.size(); i++) {
        int a = floor((Hit_X.at(i) - x + 20 * (_alpha - 1)) / 40 / _alpha);
        int b = floor((Hit_Y.at(i) - x + 20 * (_alpha - 1)) / 40 / _alpha);
        m_hit[int(Hit_Z.at(i) / 30) * 10000 + a * 100 + b] = 1;
    }
    R_alpha = true_hitno / double(m_hit.size());
    return R_alpha;
}
void Select::FindTrack() {
    for (int i = 2; i < 15; i++) {
        if (CenterX.at(i) != -500) {
            // grx->SetPoint(i, i, CenterX.at(i));
            // gry->SetPoint(i, i, CenterY.at(i));
            double x[1] = {double(i)};
            fitterx->AddPoint(x, CenterX.at(i), 1 / sqrt(lenergy[i]));
            fittery->AddPoint(x, CenterY.at(i), 1 / sqrt(lenergy[i]));
        }
    }
    if (fitterx->GetNpoints() > 3) {
        // grx->Fit("pol1", "q");
        // gry->Fit("pol1", "q");
        // track_chi2x = grx->GetFunction("pol1")->GetChisquare();
        // track_chi2y = gry->GetFunction("pol1")->GetChisquare();
        // track_chi2x /= grx->GetN();
        // track_chi2y /= gry->GetN();
        fitterx->Eval();
        fittery->Eval();
        track_chi2x = fitterx->GetChisquare();
        // cout << track_chi2x << endl;
        track_chi2y = fittery->GetChisquare();
        track_chi2x /= fitterx->GetNpoints() - 2;
        track_chi2y /= fittery->GetNpoints() - 2;
    }
}
void Select::dfs(int layer, int chip, int channel, vector<int> &cluster) {
    is_visited[layer][chip][channel] = true;
    cluster.push_back(is_hit[layer][chip][channel] - 1);
    // cluster.push_back(layer * 10000 + chip * 1000 + channel);
    int x = (Pos_X(channel, chip) + 360) / 40;
    int y = (Pos_Y(channel, chip) + 360) / 40;
    // cout << layer << "  " << chip << "  " << channel << "  " << endl;
    for (int i = 0; i < 28; i++) {
        // cout << x << "  " << y << "  " << layer << "  seed" << endl;

        int x_temp = x + dx[i];
        int y_temp = y + dy[i];
        int layer_temp = layer + dz[i];
        inverse(x_temp * 40 - 340, y_temp * 40 - 340, chip, channel);

        // cout << layer_temp << "  " << chip << "  " << channel << "  " << endl;
        if (layer_temp >= 0 && layer_temp < 40 && chip >= 0 && chip < 9 && channel >= 0 && channel < 36 && is_hit[layer_temp][chip][channel] && !is_visited[layer_temp][chip][channel]) {
            // cout << x << "  " << y << "  " << layer_temp << "  " << is_hit[layer_temp][chip][channel] << "  " << is_visited[layer_temp][chip][channel] << endl;
            dfs(layer_temp, chip, channel, cluster);
        }
    }
}
void Select::FindCluster() {
    for (int layer = 0; layer < 40; layer++) {
        for (int chip = 0; chip < 9; chip++) {
            for (int channel = 0; channel < 36; channel++) {
                if (is_hit[layer][chip][channel] && !is_visited[layer][chip][channel]) {
                    // cout << "new " << layer << " " << chip << "  " << channel << endl;
                    vector<int> cluster;
                    dfs(layer, chip, channel, cluster);
                    clusters.push_back(cluster);
                }
            }
        }
    }
    // cout << "******************************" << endl;
    // cout << clusters.size() << endl;
    // for (int i = 0; i < clusters.size(); i++) {
    //     cout << clusters[i].size() << endl;
    // }
    // // for (int x = 0; x < 10; x++) {
    // //     cout << x << ":\n";
    // //     for (int y = 4; y < 5; y++) {
    // //         for (int z = 0; z < 36; z++) {
    // //             cout << is_hit[x][y][z] << " ";
    // //         }
    // //         cout << "\n";// 每行 y 换行
    // //     }
    // //     cout << "\n";// 每层 x 换行
    // // }

    // cout << "******************************" << endl;
}
void Select::EvalEnergy() {
    int j = 0, max = 0;
    for (int i = 0; i < clusters.size(); i++) {
        int length = clusters[i].size();
        cluster_size.push_back(length);
        if (length > max) {
            j = i;
            max = length;
        }
    }
    hitno = 0;
    hitlayer = 0;
    fill(CenterX.begin(), CenterX.end(), 0);
    fill(CenterY.begin(), CenterY.end(), 0);
    fill(layerhitno, layerhitno + 40, 0);
    fill(lenergy, lenergy + 40, 0);
    fill(ishit, ishit + 40, 0);
    for (int i = 0; i < max; i++) {
        int shift = clusters[j][i];
        int layer = Hit_Z.at(shift) / 30;
        double e = Hit_Energy.at(shift);
        CenterX.at(layer) += Hit_X.at(shift) * e;
        CenterY.at(layer) += Hit_Y.at(shift) * e;
        hitno++;
        layerhitno[layer]++;
        lenergy[layer] += e;
    }
    for (int i = 0; i < 40; i++) {
        if (layerhitno[i]) {
            hitlayer++;
            ishit[i] = 1;
            CenterX.at(i) /= lenergy[i];
            CenterY.at(i) /= lenergy[i];
        } else {
            ishit[i] = 0;
            CenterX.at(i) = -500;
            CenterY.at(i) = -500;
        }
    }
    for (int i = 0; i < 40; i++) {
        if (layerhitno[i] >= 4) {
            shower_start = i;
            break;
        }
    }
    for (int i = 37; i >= 0; i--) {
        if (layerhitno[i] >= 4) {
            shower_end = i;
            break;
        }
    }
}