//	gcc *.c -o qmc -lm
#define NM        0.000000001
#define KB        1.3806503*NM*NM*NM*10000
#define EQ        1.60217646*NM*NM*0.1
#define EPS       78.000000
#define ELS       -2.837297
#define C0        299792458.00
#define PI        3.14159265358979
#define PI2       6.28318530717958
#define PI4       12.5663706143592
#define MU0       (PI4)*0.0000001
#define EP0       1/(C0*C0*MU0)
#define NA        6.02214179/NM/NM*100000
#define ME	  9.10938215*NM*NM*NM*NM*100000
#define H         6.62606896*NM*NM*NM*NM*100
#define HBAR      1.054571628*NM*NM*NM*NM*100
#define SQ(X)     (X)*(X)
#define CB(X)     (X)*SQ(X)
#define SGN(X)    ((X)>=0?1.0:-1.0)
#define rnd1	  genrand_real1
#define XX        0
#define YY        1
#define ZZ        2


typedef struct vect{
  double* dat ;
  int i;
} vector;

typedef vector *t_vector;

typedef struct mat{
  double** dat;
  int     Ni,Mj;
} matrix;

typedef matrix *t_matrix;

enum {
  etGRO=0, etPDB=1, etXYZ=2, etTOP=3, etITP=4, etDAT=5, etNFO=6
};

enum {
  cmdHARM=0, cmdMORSE=1, cmdH=2, cmdH2=3
};
void setAxis(t_matrix, int,double,double,int,int);
void fatal(char*);

/*
int load(FILE*,double[]);
int anint(double);
double gaussrnd(void);
void walk(t_matrix, t_matrix, t_matrix, double[][3], double, int, int);
int branch(t_matrix, t_matrix,double[][3], double, t_matrix, int, int, int[]);
double Weight(double, double, double);
double Prob(double[], double[],double, int,int);
double Potential(double[],double[],double[],int, int,int,double);
int birth_death(double);
void count(t_matrix,t_matrix,int,double,double,int,int);
*/
