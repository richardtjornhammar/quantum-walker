#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include "mt19937H.h"
#include "define.h"

#define NIDAT   10

int load(FILE *in, double data[])
{
  int i,r=0;

  for(i=0;i<NIDAT;i++){
    if(fscanf(in,"%lf",&data[i])==EOF)
      fatal("problem with input data");
  }
  return 0;
}

int anint(double r)
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
  int i,j,N0;
  t_matrix mat;

  N0=((N>M)?N:M);
  mat=malloc(sizeof(matrix));
  mat->dat=malloc( (N0)*sizeof(double*) );
  if(mat->dat==NULL)
    fatal("FAILED IN ALLOCATING MEMORY :: MATRIX I");

  for(i=0;i<N0;i++){
    mat->dat[i] = malloc( N0*sizeof(double) );
    if(mat->dat[i]==NULL)
      fatal("FAILED IN ALLOCATING MEMORY :: MATRIX J");
  }
  (mat->Ni)=N;
  (mat->Mj)=M;

  return(mat);
}

void permute(int *a, int n)
{
/*      PSEUDOCODE
        for i from n downto 2
        do   di ← random element of { 0, ..., i − 1 }
            swap a[di] and a[i − 1]
*/
  int i,j,at;

  for(i=n;i>=2;i--){
     j=anint(rnd1()*(i-1));
     at=a[j];a[j]=a[i-1];a[i-1]=at;
  }
}

void initMat(t_matrix walker, int N0, double R[7][3], int cmd)
{
  int   i,j,k[3]={0,1,2};
  int   l,m,n;
  double sign=-1.0;

  switch(cmd){
    case 5:
    for(i=0;i<walker->Ni;i++){
      walker->dat[i][0]=(i<N0?1.0:0.0);
      for(j=0;j<walker->Mj-1;j++){
        walker->dat[i][j+1]= walker->dat[i][0]>0.0 ? ( R[7][j] ) :0.0;
      }
    }
    break;
    case 4:
    for(i=0;i<walker->Ni;i++){
      walker->dat[i][0]=(i<N0?1.0:0.0);
      permute(k,3);
      for(j=XX;j<=ZZ;j++){ // H
        walker->dat[i][j+1]= walker->dat[i][0]>0.0 ? ( R[5][k[j]]* pow(sign,(double) anint(rnd1()+1.0) ) ) :0.0; 
      }
      for(l=1;l<=9;l++){   // O-
        permute(k,3);
        for(j=XX;j<=ZZ;j++){
          walker->dat[i][j+1+(ZZ+1)*l] = walker->dat[i][0]>0.0 ? ( R[6][k[j]]* pow(sign,(double) anint(rnd1()+1.0) ) ) :0.0; 
        }
      }
    }
    break;
    case 3:
    for(i=0;i<walker->Ni;i++){
      walker->dat[i][0]=(i<N0?1.0:0.0);
      permute(k,3);
      for(j=XX;j<=ZZ;j++){ 
        walker->dat[i][j+1]= walker->dat[i][0]>0.0 ? ( R[cmd][k[j]]* pow(sign,(double) anint(rnd1()+1.0) ) ) :0.0; 
      }
      for(j=XX;j<=ZZ;j++){ 
        walker->dat[i][j+1+ZZ+1]= walker->dat[i][0]>0.0 ? ( R[cmd+1][k[j]]* pow(sign,(double) anint(rnd1()+1.0) ) ) :0.0; 
      }      
    }
    break;
    case 2:
    for(i=0;i<walker->Ni;i++){
      walker->dat[i][0]=(i<N0?1.0:0.0);
      permute(k,3);
      for(j=0;j<walker->Mj-1;j++){ 
        walker->dat[i][j+1]= walker->dat[i][0]>0.0 ? ( R[cmd][k[j]]* pow(sign,(double) anint(rnd1()+1.0) ) ) :0.0; 
      }
    }
    break;
    default:
    for(i=0;i<walker->Ni;i++){
      walker->dat[i][0]=(i<N0?1.0:0.0);
      for(j=0;j<walker->Mj-1;j++){ 
        walker->dat[i][j+1]= walker->dat[i][0]>0.0 ? ( R[cmd][j] ) :0.0;
      }
    }
    break;
  }
}

void printMat(t_matrix walker, FILE *out, int DIM)
{ //FOR PRINTING THE GSWF
  int i,j;
  int M;

  M=DIM+2;//( ((walker->Mj)<DIM )?(walker->Mj):(DIM+1));

  for(i=0;i<walker->Ni;i++){
    for(j=0;j<M;j++){
      fprintf(out,"% 12.5f ",walker->dat[i][j]);
    }
    fprintf(out,"\n");  
  }
}

void avec(double *vec, int N)
{
  vec=malloc(sizeof(double)*N);
}

//THIS GAUSSIAN RANDOM NUMBER GENERATOR USES THE MERSENNE TWISTER VIA DEFINE OF RND1
double gaussrnd(void)
{
  double fac=0.0,rsq=100.0,v1=0.0,v2=0.0,gset=0.0;

  while (rsq >= 1.0 || rsq == 0.0){
    v1=2.0*rnd1()-1.0;
    v2=2.0*rnd1()-1.0;
    rsq=v1*v1+v2*v2;
  }
  fac=sqrt(-2.0*log(rsq)/rsq);
  gset=v1*fac;
  return (double)(v2*fac);
}

int birth_death(double W)
{
  int i;
  double dnr,part,frac;

  dnr=(W+rnd1());
  frac=modf(dnr,&part);
  //i=(part>=3)?3:part; //CAP REPLICA CREATION AT 2 FOR EACH STEP
  i=part;
  return(i);
}


double Potential(double X[], double R1[], double R2[], int type, int d, int s, double dtau)
{
    double V=0.0,omega,K=1.0,a=1.0,L=0.0;
    int i,j,k,l,m,n;

    switch(type){
       case 0:
         K=1.0;
         for(i=s;i<=d;i++){
           V+=0.5*K*X[i]*X[i];
         }
         break;
       case 1:
         K=0.5; a=1.0;
         for(i=s;i<=d;i++){
           V+=K*( exp(-2.0*X[i]*a)-2.0*exp(-1.0*X[i]*a));
         }
         break;
       case 2:
         K=0.0;
         for(i=s;i<=d;i++){
           omega=pow((X[i]-R1[i-1]),2.0);
           K+=omega;
         }         
         V-=1.0/(double)sqrt(K);
         break;
       case 3:
         K=0.0;
         for(i=s;i<=d;i++){
           omega=pow((X[i]-R1[i-1]),2.0);
           K+=omega;
         }         
         V-=1.0/(double)sqrt(K);
         K=0.0;
         for(i=s;i<=d;i++){
           omega=pow((X[i]-R2[i-1]),2.0);
           K+=omega;
         }         
         V-=1.0/(double)sqrt(K);
	 K=0.0;
         for(i=s;i<=d;i++){
           omega=pow((X[i+d]-R1[i-1]),2.0);
           K+=omega;
         }         
         V-=1.0/(double)sqrt(K);
	 K=0.0;
         for(i=s;i<=d;i++){
           omega=pow((X[i+d]-R2[i-1]),2.0);
           K+=omega;
         }         
         V-=1.0/(double)sqrt(K);
	 K=0.0;
         for(i=s;i<=d;i++){
           omega=pow((X[i]-X[i+d]),2.0);
           K+=omega;
         }         
         V+=1.0/(double)sqrt(K);
         break;
         case 4: //R1 is H R2 is O
           n=9;
           for(j=0;j<n;j++){  //REPULSION TERMS
             for(l=j+1;l<=n;l++){
               K=0.0;
               for(i=XX+1;i<=ZZ+1;i++){
                 omega=pow((X[i+(ZZ+1)*l]-X[i+(ZZ+1)*j]),2.0);
                 K+=omega;
               }          
               V+=1.0/(double)sqrt(K);
             }
           }
           for(j=0;j<=n;j++){ // ion-electron attraction
             K=0.0;
             for(i=XX+1;i<=ZZ+1;i++){
               omega=pow((R1[i-1]-X[i+(ZZ+1)*j]),2.0);
               K+=omega;
             }          
             V-=1.0/(double)sqrt(K); // H
             K=0.0;
             for(i=XX+1;i<=ZZ+1;i++){
               omega=pow((R2[i-1]-X[i+(ZZ+1)*j]),2.0);
               K+=omega;
             }          
             V-=8.0/(double)sqrt(K); // O
           }
         break;
         case 5://WOOD SAXON POTENTIAL
           a=1.0;
           V=-1/(1+exp( (X[1]-R1[0])/a ) );
         break;
    }

    return(V);
}

double Prob(double X1[], double X0[],double dtau, int d, int s) //s=1
{
   //X1 is the present and X0 the preceding step
   int i;
   double r;

   r=1/sqrt(PI2*dtau);

   for(i=s;i<=d;i++)
     r*=exp(-pow(X1[i]-X0[i],2.0)/2/dtau);

   return(r);
}

double Weight(double V, double ER, double dtau)
{
  double W=0.0;

  W = exp( -(V-ER)*dtau );

  return( W );
 
}

void walk(t_matrix walker, t_matrix owalker, t_matrix E, double R[][3] , double dtau, int DIM, int type)
{
  int i,j,k,N,M,F,RR;
  double V,W,sigma,PROD;

  sigma=sqrt(dtau);
  N=walker->Ni; 
  M=walker->Mj;

  for(i=0;i<N;i++){
    owalker->dat[i][0]=walker->dat[i][0];
    if((int)walker->dat[i][0]){ //IF IT IS ALIVE
      //RR=XX-1;
      for(j=1;j<M;j++){
        owalker->dat[i][j]=walker->dat[i][j];
        walker->dat[i][j]+=gaussrnd()*sigma;
        //V = Potential(walker->dat[i],R[type],R[type+1],type,d,1,scale);
        PROD=Prob(walker->dat[i],owalker->dat[i],dtau,j,j);
        F=rnd1()>PROD;
        if(F){ 
          walker->dat[i][j]=owalker->dat[i][j];
        }
        else {
          owalker->dat[i][j]=walker->dat[i][j];
        }
      }
    }
  }

}

int branch(t_matrix walker, t_matrix owalker, double R[][3], double scale, t_matrix E, int type, int d, int N[])
{
  int i,j,k,l,Nt,M,mn,flag=100;
  double V,V0,W,P,mV1=0.0,mV0=0.0,mW=0.0,alpha=0.0,trash=0.0;

  Nt=walker->Ni; M=walker->Mj;
  N[1]=0;
  alpha=1.0;//E->dat[0][5];
  
  for(i=0;i<Nt;i++){
    if( (int)walker->dat[i][0]>=1){ //IF IT IS ALIVE
      V = Potential(walker->dat[i],R[type],R[type+1],type,d,1,scale);
      W = Weight(V,E->dat[0][0],scale);

      //REPLICATE
      mn=birth_death(W);
      N[1]+=mn;

      if(N[1]>Nt){
        fprintf(stderr,"\nN[1](%10d) = %d\n",i,N[1]);
        fatal("PROCESS DIED [ OUT OF BOUNDS ]");
      }
      switch(mn){
        case 0:
          walker->dat[i][0]=0;
          break;
        case 1:
          break;
        default: //SPAWN NEW PROCESSES
          k=1;
          for(j=0;j<Nt;j++){     
            if( !((int)walker->dat[j][0]) && k < mn ){
              walker->dat[j][0]=1;
              for(l=1;l<M;l++){
                walker->dat[j][l]=walker->dat[i][l];
              }
              k++;
            }
          }
          break;
      }  
      mV1+=V*mn;
      mW+=W*mn;
    }
  }
  mV1/=N[1]>0?(double)N[1]:1.0;
  mW/=N[1]>0?(double)N[1]:1.0;

  E->dat[0][3]+=mV1;
  E->dat[0][4]+=1.0;
  E->dat[0][1]=E->dat[0][3]/E->dat[0][4];
  mn=E->dat[0][4]<1000?1:0;
  E->dat[0][0]=E->dat[0][mn]+alpha*(1-(double)N[1]/(double)N[0]);
  E->dat[0][6]=mW;

  if(!N[1]){
    fprintf(stderr,"%10.6f%10.6f%10.6f%10d\n",V,mV1,E->dat[0][0],N[1]);
    fatal("PROCESS DIED [ NO WALKERS LEFT ]");
  }
  N[0]=N[1];

  return(N[1]);
}

void setAxis(t_matrix GSWF, int nb,double min,double max,int typ,int DIM)
{                                                    //  CURRENTLY ONLY 1D MAPPED PROBLEMS
  int i,j,k,l,N,M,tmp,D;
  double temp,sqSUM=0.0;
  
  N=GSWF->Ni; M=GSWF->Mj;
  D=(typ==2)?1:DIM;                                  //  WE HAVE ALREADY COLLECTED THIS SPATIAL DISTRIBUTION ON A RADIAL AXIS

  for(i=0;i<N;i++){
    GSWF->dat[i][0]=(double)i*(max-min)/(double)nb+min;
    for(j=1;j<=D;j++){
      if(typ==2){
        GSWF->dat[i][j]/=GSWF->dat[i][0];
      }
      sqSUM+=GSWF->dat[i][j]*GSWF->dat[i][j];        //  NORMALISING SQUARE OF WAVE FUNCTION
    }
  }
  temp=sqrt( pow( ((max-min)/(double)nb),D) * sqSUM); //  SIZE OF ONE ELEMENT
  switch(typ){
    case 0:
      sqSUM=1.0/sqrt(sqrt(PI));
      break;
    case 1:
      sqSUM=sqrt(2.0);
      break;
    case 2:
      sqSUM=2.0;
      break;
    default:
      sqSUM=1.0;
      break;
  }
  for(i=0;i<N;i++){
    GSWF->dat[i][1]/=temp;
    switch(typ){
      case 0:
        GSWF->dat[i][2]=sqSUM*exp(-0.5*GSWF->dat[i][0]*GSWF->dat[i][0]);
        break;
      case 1:
        GSWF->dat[i][2]=sqSUM*exp(-1.0*exp(-1.0*GSWF->dat[i][0])-GSWF->dat[i][0]/2.0);
        break;
      case 2:
        GSWF->dat[i][2]=((GSWF->dat[i][0]>=0)?sqSUM:0.0)*exp(-1.0*GSWF->dat[i][0]);//*GSWF->dat[i][0]
        GSWF->dat[i][3]=((GSWF->dat[i][0]>=0)?sqSUM:0.0)*exp(-1.0*GSWF->dat[i][0])*GSWF->dat[i][0];
        GSWF->dat[i][4]=GSWF->dat[i][1]/GSWF->dat[i][0];
        break;
    }
  }
}

void count(t_matrix walker, t_matrix GSWF, int nb, double min, double max, int typ, int DIM)
{
  int i,j,k,l,N,M,tmp;
  double temp,r=0.0,fact;

  N=walker->Ni; M=walker->Mj;
  temp=1.0;  fact=2.0;

  for(i=0;i<N;i++){
    if(i<=nb)
      GSWF->dat[i][0]=1.0;
    if((int)walker->dat[i][0]){ //IF IT IS ALIVE
      for(l=1;l<2;l++){
        r=0.0;
        if(typ==2){ // radial binning for the H atom
          for(k=1;k<=DIM;k++){
            r+=pow(walker->dat[i][k],2.0);
          }
          r=sqrt(r);
        }
        else{
          r=walker->dat[i][l];
        }
        tmp=floor( (r-min)*temp / ( (max-min)/((double)nb) ) );
        tmp=tmp>nb?nb:tmp;
        tmp=tmp<0?0:tmp;
        GSWF->dat[tmp][l]+=1.0;
      }
    }
  } 
}
