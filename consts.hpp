#ifndef _CONSTS_H_
#define _CONSTS_H_
#include <math.h>

using namespace std;

const int NAtom = 1000; //粒子数
const double NRho = 0.8; //数密度
const double LBox = sqrt(NAtom/NRho); //ボックスの一辺の長さ
const double LHalf = LBox/2.0; //その半分

const double SDB = 1.0;
const double Rcut = 4.0; //カットオフ長
const int NCycle = 200; //計算回数
const int NPrint = 10; //結果をprintする頻度
const double dT = 0.001; //時間の刻み幅
const int NInitConfig = 1+floor(sqrt(NAtom)); //粒子の初期配置における一辺あたりの粒子数

const double RList = 4.5; //VNLに入れる範囲
const double Skin = RList - Rcut;
const int MCell = 3;
const int MCell_rough = floor(LBox/RList);



#endif
