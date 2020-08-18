#include <iostream>
#include <vector>
#include <array>
#include <math.h>
#include "consts.hpp"
#include <stdio.h>
#include "ensemble.hpp"
#include <time.h>

using namespace std;

int main(){
  std::cout << "MCellの目安: "<<MCell_rough <<std::endl;
  Ens ens;
  ens.make_VNL();
  FILE *fp, *fp2, *fp3;
  fp=fopen("config.txt","w");
  fp2=fopen("energy.txt","w");
  fp3=fopen("final_config.txt","w");
  fprintf(fp,"NATOM NRHO NCYCLE NPRINT \n");
  fprintf(fp,"%d %f %d %d \n",NAtom,NRho,NCycle,NPrint);
  int i,j;
  clock_t start = clock();
  for(i=0;i<NCycle;i++){
    cout << "\r" << i <<" "<< flush;
    ens.VVcycle();
    if(i%NPrint==0){
      for(j=0;j<NAtom;j++){
        fprintf(fp,"%f ",ens.getx()[j]);
      }
      fprintf(fp,"\n");
      for(j=0;j<NAtom;j++){
        fprintf(fp,"%f ",ens.gety()[j]);
      }
      fprintf(fp,"\n");
      fprintf(fp2,"%f %f\n",i*dT,ens.get_energy());
    }
  }

  fprintf(fp3,"NATOM NRHO NCYCLE NPRINT \n");
  fprintf(fp3,"%d %f %d %d \n",NAtom,NRho,NCycle,NPrint);
  for(j=0;j<NAtom;j++){
    fprintf(fp3,"%f ",ens.getx()[j]);
  }
  fprintf(fp3,"\n");
  for(j=0;j<NAtom;j++){
    fprintf(fp3,"%f ",ens.gety()[j]);
  }
  fprintf(fp3,"\n");

  fclose(fp);
  fclose(fp2);
  fclose(fp3);
  clock_t end = clock();
  cout << "処理時間 " << end-start <<endl;
  return 0;
}
