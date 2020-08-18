#include <iostream>
#include <vector>
#include <array>
#include <math.h>
#include "consts.hpp"
#include "ensemble.hpp"
#include "MT.h"

using namespace std;

//コンストラクタ
Ens::Ens(){
  int i, ix, iy;
  div_t div_i;
  double sum_ux=0, sum_uy=0, r1, r2;
  init_genrand(1);
  for(i=0; i<NAtom; i++){
    div_i = div(i, NInitConfig);
    ix = div_i.quot;
    iy = div_i.rem;
    x[i] = ((double)ix+0.5)*LBox/((double)NInitConfig);
    y[i] = ((double)iy+0.5)*LBox/((double)NInitConfig);
    x_ref[i] = x[i];
    y_ref[i] = y[i];
    r1=genrand_real3();
    r2=genrand_real3();
    ux[i]=sqrt(-2.0*log(r1))*cos(2.0*M_PI*r2);
    uy[i]=sqrt(-2.0*log(r1))*sin(2.0*M_PI*r2);
    sum_ux += ux[i];
    sum_uy += uy[i];
  }
  //系全体の運動量をゼロにする
  for(i=0; i<NAtom; i++){
    ux[i] -= sum_ux/(double)NAtom;
    uy[i] -= sum_ux/(double)NAtom;
  }
}

bool Ens::apply_pbc(){
  int i;
  for(i=0;i<NAtom;i++){
    if(x[i]<0.0){
      x[i] += LBox;
    }else if(x[i]>=LBox){
      x[i] -= LBox;
    }
    if(y[i]<0.0){
      y[i] += LBox;
    }else if(y[i]>=LBox){
      y[i] -= LBox;
    }
  }
  return true;
}

void Ens::calcu(double T1){
  int i;
  for(i=0;i<NAtom;i++){
    ux[i] += T1*Fx[i];
    uy[i] += T1*Fy[i];
  }
}

void Ens::calcr(double T1){
  int i;
  for(i=0;i<NAtom;i++){
    x[i] += T1*ux[i];
    y[i] += T1*uy[i];
  }
}

double Ens::pbc_dist(double r){
  double result;
  if(r >= LHalf){
    result = r - LBox;
  }else if(r < -LHalf){
    result = r + LBox;
  }else{
    result = r;
  }
  return result;
}

void Ens::give_label(){
  int i,im;
  for(im=0;im<MCell*MCell;im++){
    LinkedList[im].clear();
    LinkedList[im].shrink_to_fit();
  }

  for(i=0;i<NAtom;i++){
    //粒子iがどのセルに属するかを調べる →コメントでは、簡単にセル(ix,iy)と表す
    label_x[i] = floor(MCell/LBox*x[i]);
    label_y[i] = floor(MCell/LBox*y[i]);
    //セル(ix,iy)を番号ix+iy*Mで指定
    //LinkedListを(ix, iy)に入っている粒子の一覧にするため、LinkedListに番号iを追加
    LinkedList[label_x[i]+label_y[i]*MCell].push_back(i);
  }
}

void Ens::make_VNL(){
  NeighborList_i.clear();
  NeighborList_i.shrink_to_fit();
  NeighborList_j.clear();
  NeighborList_j.shrink_to_fit();
  pot=0.0;
  int i,ix,iy,i_delta,jx,jy,j_cell,k,j;
  double xij, yij, rij, Eij, Fij;
  for(i=0;i<NAtom;i++){
    Fx[i]=0.0;
    Fy[i]=0.0;
  }
  give_label();

  array<int, 5> delta_ix={0, 0, 1, 1, 1};
  array<int, 5> delta_iy={0, -1, -1, 0, 1}; //順に、同じセル、上、右上、右、右下を指定する
  for(i=0;i<NAtom;i++){
    ix=label_x[i];
    iy=label_y[i];

    for(i_delta=0;i_delta<5;i_delta++){
      jx = floored_mod((ix+delta_ix[i_delta]) , MCell);  //近接するセルを指定
      jy = floored_mod((iy+delta_iy[i_delta]) , MCell);
      j_cell=jx+jy*MCell;
      for(k=0;k<LinkedList[j_cell].size();k++){
        j=LinkedList[j_cell][k];
        if((i_delta==0) && (j<=i)){
          continue;
        }
        xij = pbc_dist(x[j] - x[i]);
        yij = pbc_dist(y[j] - y[i]);
        rij=sqrt(pow(xij,2)+pow(yij,2));
        Eij=0.0;
        Fij=0.0;

        if(rij<RList){
          NeighborList_i.push_back(i);
          NeighborList_j.push_back(j);
          if(rij<Rcut){
            Eij=4.0*(pow(rij,-12)-pow(rij,-6));
            Fij=24.0*(-2.0*pow(rij,-13)+pow(rij,-7));
            Fx[i] += Fij*xij/rij;
            Fy[i] += Fij*yij/rij;
            Fx[j] -= Fij*xij/rij;
            Fy[j] -= Fij*yij/rij;
            pot += Eij/(double)NAtom;
          }
        }
      }
    }
  }
  cout << "update VNL" <<endl;
}

/*
void Ens::make_VNL2(){
  NeighborList_i.clear();
  NeighborList_i.shrink_to_fit();
  NeighborList_j.clear();
  NeighborList_j.shrink_to_fit();
  pot=0.0;
  int i,j;
  double xij, yij, rij, Eij, Fij;
  for(i=0;i<NAtom;i++){
    Fx[i]=0.0;
    Fy[i]=0.0;
  }
  for(i=0;i<NAtom;i++){
    for(j=i+1;j<NAtom;j++){
      xij = pbc_dist(x[j] - x[i]);
      yij = pbc_dist(y[j] - y[i]);
      rij=sqrt(pow(xij,2)+pow(yij,2));
      Eij=0.0;
      Fij=0.0;

      if(rij<RList){
        NeighborList_i.push_back(i);
        NeighborList_j.push_back(j);
        if(rij<Rcut){
          Eij=4.0*(pow(rij,-12)-pow(rij,-6));
          Fij=24.0*(-2.0*pow(rij,-13)+pow(rij,-7));
          Fx[i] += Fij*xij/rij;
          Fy[i] += Fij*yij/rij;
          Fx[j] -= Fij*xij/rij;
          Fy[j] -= Fij*yij/rij;
          pot += Eij/(double)NAtom;
        }
      }
    }
  }
  cout << "update VNL" <<endl;
}*/

bool Ens::check_update_VNL(){
  double diff,diffmax; //i番目の粒子がrefと比べてどれくらいの距離離れたか
  int i;
  diffmax=0.0;
  for(i=0;i<NAtom;i++){
    diff=sqrt(pow((x[i]-x_ref[i]),2)+pow((y[i]-y_ref[i]),2));
    if(diff>diffmax){
      diffmax=diff;
    }
  }
  if(diffmax > Skin/2.0){
    return true;
  }else{
    return false;
  }
}

void Ens::use_VNL(){
  pot=0.0;
  int i,j,k,i_F;
  for(i_F=0;i_F<NAtom;i_F++){
    Fx[i_F]=0.0;
    Fy[i_F]=0.0;
  }
  double xij,yij,rij,Eij,Fij;
  for(k=0;k<NeighborList_i.size();k++){
    i=NeighborList_i[k];
    j=NeighborList_j[k];
    xij = pbc_dist(x[j] - x[i]);
    yij = pbc_dist(y[j] - y[i]);
    rij=sqrt(pow(xij,2)+pow(yij,2));
    if(rij<Rcut){
      Eij=4.0*(pow(rij,-12)-pow(rij,-6));
      Fij=24.0*(-2.0*pow(rij,-13)+pow(rij,-7));
      Fx[i] += Fij*xij/rij;
      Fy[i] += Fij*yij/rij;
      Fx[j] -= Fij*xij/rij;
      Fy[j] -= Fij*yij/rij;
      pot += Eij/(double)NAtom;
    }
  }
}

void Ens::calcF(){
  if(check_update_VNL()){
    make_VNL(); //!
    int i;
    for(i=0;i<NAtom;i++){
      x_ref[i]=(x[i]*1.0);
      y_ref[i]=(y[i]*1.0);
    }
  }else{
    use_VNL();
  }
}

void Ens::VVcycle(){
  calcu(dT/2.0);
  calcr(dT);
  apply_pbc();
  calcF();
  calcu(dT/2.0);
  temp = 0;
  int i;
  for(i=0;i<NAtom;i++){
    temp += ((pow(ux[i],2)+pow(uy[i],2))/2.0)/double(NAtom);
  }
  energy = temp + pot;
}

array<double,NAtom> Ens::getx(){
  return x;
}

array<double,NAtom> Ens::gety(){
  return y;
}

array<double,NAtom> Ens::getux(){
  return ux;
}

array<double,NAtom> Ens::getuy(){
  return uy;
}

array<double,NAtom> Ens::getFx(){
  return Fx;
}

array<double,NAtom> Ens::getFy(){
  return Fy;
}

double Ens::get_energy(){
  return energy;
}

double Ens::get_pot(){
  return pot;
}

double Ens::get_temp(){
  return temp;
}

//-------------------------------------
//切り捨て割り算（余りの符号はbの符号と一致）
int floored_mod(int a, int b){
  return a-floor((double)a/(double)b)*b;
}
