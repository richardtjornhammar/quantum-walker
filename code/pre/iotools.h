#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include "define.h"

#define NIDAT   10

int       load(FILE *in, float data[]);
int       anint(float r);
t_matrix  amat(int N, int M);
void      initMat(t_matrix psi, int N0, float R[3][3], int cmd);
void      printMat(t_matrix psi);
void      avec(float *vec, int N);
float     gasdev(void);
void      walk(t_matrix psi, float sigma);

