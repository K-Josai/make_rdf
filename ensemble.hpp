#ifndef _ENSEMBLE_H_
#define _ENSEMBLE_H_

#include <iostream>
#include <vector>
#include <array>
#include <math.h>
#include "consts.hpp"

using namespace std;

class Ens{
private:
  array<double, NAtom> x, y, ux, uy, Fx, Fy, x_ref, y_ref;
  array<int, NAtom> label_x, label_y;
  double pot, temp, energy;
  array< vector<int>, MCell*MCell > LinkedList;
  vector<int> NeighborList_i, NeighborList_j;

public:
  Ens(); //コンストラクタ
  int set_initial_value();
  array<double,NAtom> getx();
  array<double,NAtom> gety();
  array<double,NAtom> getux();
  array<double,NAtom> getuy();
  array<double,NAtom> getFx();
  array<double,NAtom> getFy();
  double get_energy();
  double get_pot();
  double get_temp();
  bool apply_pbc();
  void calcu(double T1);
  void calcr(double T1);
  void give_label();
  void make_VNL();
  // void make_VNL2(); //CLLを使わない簡易的なもの
  bool check_update_VNL();
  void use_VNL();
  void calcF();
  static double pbc_dist(double r);
  void VVcycle();
};

int floored_mod(int, int);

#endif
