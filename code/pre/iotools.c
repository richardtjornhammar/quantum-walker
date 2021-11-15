#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include "define.h"

#define NIDAT   10

int load(FILE *in, float data[])
{
  int i,r=0;

  for(i=0;i<NIDAT;i++){
    if(fscanf(in,"%f",&data[i])==EOF)
      fatal("problem with input data");
  }
  return 0;
}

int anint(float r)
{
   int sign,res;
   double frac,part;
   
   sign=(r>=0)?1:-1;
   frac=modf(r,&part);

   res=(part+((frac*sign>=0.5)?1.0:0.0)*sign);

   return(res);   

}

t_matrix amat(int N, int M)
{
  int    i,j,N0;
  t_matrix mat;

  N0=((N>M)?N:M);
  mat=malloc(sizeof(matrix));
  mat->dat=malloc( (N0)*sizeof(float*) );
  if(mat->dat==NULL)
    fatal("FAILED IN ALLOCATING MEMORY :: MATRIX I");

  for(i=0;i<N0;i++){
    mat->dat[i] = malloc( (N0*N0)*sizeof(float) );
    if(mat->dat[i]==NULL)
      fatal("FAILED IN ALLOCATING MEMORY :: MATRIX J");
  }
  (mat->Ni)=N; 
  (mat->Mj)=M;

  return(mat);
}

void initMat(t_matrix psi, int N0, float R[3][3], int cmd)
{
  int   i,j;
  float sign=-1.0;

  for(i=0;i<psi->Ni;i++){
    psi->dat[i][0]=(i<N0?1.0:0.0);
    for(j=0;j<psi->Mj;j++){
      psi->dat[i][j+1]=R[cmd][j]/2*((j==(psi->Mj-2.0))?sign:1.0);
    }
    sign*=-1;
  }
}

void printMat(t_matrix psi)
{
  int i,j;

  fprintf(stdout,"\nOUTPUT OF MATRIX::\n");
  for(i=0;i<psi->Ni;i++){
    if((int)psi->dat[i][0]){
      for(j=0;j<psi->Mj-1;j++){
        fprintf(stdout,"%9.5f ",psi->dat[i][j+1]);
      }
      fprintf(stdout,"\n");
    }
  }
}

void avec(float *vec, int N)
{
  vec=malloc(sizeof(float)*N);
}

//THIS VERSION OF GASDEV USES THE MERSENNE TWISTER VIA DEFINE OF RND1
float gasdev(void)
{
//    float rnd1(void);
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;

    if (iset == 0) {
        do {
            v1=2.0*rnd1()-1.0;
            v2=2.0*rnd1()-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return (float)(v2*fac);
    } else {
        iset=0;
        return (float)gset;
    }
}


void walk(t_matrix psi, float sigma)
{
  int i,j,k,N,M;
  N=psi->Ni; M=psi->Mj;

  for(i=0;i<N;i++){
    if((int)psi->dat[i][0]) //IF IT IS ALIVE
      for(j=1;j<M;j++){
        psi->dat[i][j]+=0.0;//%sigma; //*sigma;
    }
  }
  fprintf(stdout,"\n%f %f\n",rnd1(),sigma);
}


