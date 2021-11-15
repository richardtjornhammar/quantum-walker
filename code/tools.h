#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include "mt19937ar.h"
//#include "define.h"

#define NIDAT   10

int       load(FILE *, double[] );
int       anint(double);
t_matrix  amat(int , int );
void      permute(int *, int );
void      initMat(t_matrix , int, double (*)[3], int);
void      printMat(t_matrix walker, FILE *out, int DIM);
void      avec(double *vec, int N);
double    gaussrnd(void);
int       birth_death(double W);
double    Potential(double X[], double R1[], double R2[], int type, int d, int s, double dtau);
double    Prob(double X1[], double X0[],double dtau, int d, int s);
double    Weight(double V, double ER, double dtau);
void      walk(t_matrix walker, t_matrix owalker, t_matrix E, double R[][3] , double dtau, int DIM, int type);
int       branch(t_matrix walker, t_matrix owalker, double R[][3], double scale, t_matrix E, int type, int d, int N[]);
void      setAxis( t_matrix , int , double , double , int , int );
void      count(t_matrix walker, t_matrix GSWF, int nb, double min, double max, int typ, int DIM);
