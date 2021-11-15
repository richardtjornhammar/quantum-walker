#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include <math.h>
#include "define.h"
#include "tools.h"
#include "mt19937H.h"

const char *argp_program_version =
"QMC v1.0 by Richard Tjornhammar";

const char *argp_program_bug_address =
"<richard.tjornhammar@gmail.com>";

/* This structure is used by main to communicate with parse_opt. */
struct arguments
{
  char *args[1];            /* ARG1 and ARG2 */
  int verbose;              /* The -v flag */
  char *outfile,*infile;
  float RAD;
};

void fatal(char errstr[])
{
  fprintf(stderr,"ERROR::  %s\n",errstr);
  exit(0);
}

int file_ext(char filenm[], char filext[])
{
  char* ptr;

  ptr=strpbrk(filenm,".");
  ptr;
  if(ptr==NULL)
    return(-1);
  else{
    ptr++;
    return(strcmp(ptr,filext));
  }
}

/*
   OPTIONS.  Field 1 in ARGP.
   Order of fields: {NAME, KEY, ARG, FLAGS, DOC}.
*/
static struct argp_option options[] =
{
  {"verbose", 'v', 0, 0, "Produce verbose output"},
  {"input",  'i', "INFILE", 0,
   "Input from INFILE instead"},
  {"output",  'o', "OUTFILE", 0,
   "Output to OUTFILE instead of to standard output"},
  {"radius",  'R', "RADIUS", 0,
   "Use this radial separation distance"},
  {0}
};
/*
   PARSER. Field 2 in ARGP.
   Order of parameters: KEY, ARG, STATE.
*/
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'R':
      arguments->RAD = atof(arg);
      break;
    case 'v':
      arguments->verbose = 1;
      break;
    case 'o':
      arguments->outfile = arg;
      break;
    case 'i':
      arguments->infile = arg;
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num >= 1)
	{
	  argp_usage(state);
	}
      arguments->args[state->arg_num] = arg;
      break;
    case ARGP_KEY_END:
      if (state->arg_num < 1)
	{
	  argp_usage (state);
	}
      break;
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/*
   ARGS_DOC. Field 3 in ARGP.
   A description of the non-option command-line arguments
     that we accept.
*/
static char args_doc[] = "CMD";//"ARG1 ARG2";

/*
  DOC.  Field 4 in ARGP.
  Program documentation.
*/
static char doc[] =
"QMC a small path integral Quantum Monte Carlo code";
/*
   The ARGP structure itself.
*/
static struct argp argp = {options, parse_opt, args_doc, doc};

int main (int argc, char **argv)
{
  //GENERAL
  int i,j,k,nr;

  //IO
  struct arguments arguments;
  FILE *out,*in,*outWF,*outRE;
  int  typ,cmdtyp;
  char *commands[]={"HARMONIC","MORSE","H","H2","OH","WOODS"};
  char *cmd;
  double data[11];

  /* COORDINATES */
  double CELL[6],a,b,c;
  int   LA,LB;
  double R0[8][3]={{0,0,0},{0,0,0},{0,0,0},{0,0,0.7},{0,0,-0.7},{0,0,0.625},{0,0,-0.625},{0.2,0.0,0.0}}; //ION POSITION
  double RE[8][3]={{0,0,0},{0,0,0},{0,0,1.1},{0,0,1.0},{0,0,-1.0},{0,0,-0.50},{0,0,0.50},{0.5,0.0,0.5}}; //RELATIVE ELECTRON POSITION

  /* QMC ALGO */
  int      N0,Nmax,nb,DIM,seed,NSTEPS,iT,N[2],P=1;
  double   tau0,dtau,rmin,rmax,V,ER,nmean;
  double   sigm=0.0,time;
  t_matrix walker,owalker,GSWF,Energy;
  double   frac,fract,part,parti;

  /* DEFAULTS */
  arguments.outfile = NULL;
  arguments.infile = NULL;
  arguments.verbose = 0;
  arguments.RAD=0.0;

  /* PARSER */
  argp_parse (&argp, argc, argv, 0, 0, &arguments);
  cmd=arguments.args[0];

  i=0;
  while(!(cmd[i]=='\0')){
    cmd[i]=toupper(cmd[i]);
    i++;
  }
  for(i=0;i<6;i++){
    if(!strcmp(commands[i],cmd))
      cmdtyp=i;
  }
  if(!(cmdtyp==0||cmdtyp==1 || cmdtyp==2 || cmdtyp==3 || cmdtyp==4 || cmdtyp==5))
    fatal("ONLY THE HARMONIC, MORSE, H, H2 , OH AND WOODS COMMANDS ARE IMPLEMENTED"); // OH WILL NOT CONVERGE
  DIM=((cmdtyp>=2 && cmdtyp<5 )?3:1);
  P=((cmdtyp==3)?(2):(((cmdtyp==4)?(10):(1))));

  if( (cmdtyp>=3) && arguments.RAD){
    fprintf(stderr,"\nCAUGHT R=%f\n",arguments.RAD);
    R0[(cmdtyp==4?5:cmdtyp)][ZZ]=arguments.RAD*0.5;
    R0[(cmdtyp==4?6:cmdtyp+1)][ZZ]=-arguments.RAD*0.5;
  }

  for(i=XX;i<=ZZ;i++){
    RE[(cmdtyp==4?5:cmdtyp==5?7:cmdtyp)][i]+=R0[(cmdtyp==4?5:cmdtyp==5?7:cmdtyp)][i];
    if(cmdtyp==3||cmdtyp==4){
      RE[(cmdtyp==4?6:cmdtyp+1)][i]+=R0[(cmdtyp==4?6:cmdtyp+1)][i];
    }
  }

  /* OUTPUT */
  if (arguments.outfile){
    out = fopen (arguments.outfile, "w");
    if(out==NULL)
      fatal("COULD NOT OPEN OUTPUT FILE");
  }
  else
    out = stdout;

  outWF = fopen ("GSWF.dat", "w");
  outRE = fopen ("REPL.dat", "w");

  /* VERBOSE */
  if (arguments.verbose)
    fprintf(stderr,"%s\n","VERBOSE");

  if(arguments.infile){
    if(!file_ext(arguments.infile,"nfo"))
      typ=etNFO;
    if(file_ext(arguments.infile,"nfo"))
      fatal("INPUT FILE IS NOT A RECONGIZED FILETYPE");
    in=fopen(arguments.infile,"r");
    if(in==NULL)
      fatal("FAILED TO OPEN INPUT FILE");
  }
  else
    fatal("NO INPUT FILE GIVEN");

  load(in,data);
  Nmax=anint(data[0]); nb=anint(data[5]); NSTEPS=data[1];
  dtau=data[2]; rmin=data[3]; rmax=data[4]; tau0=NSTEPS*dtau;
  V=anint(data[6]); N0=anint(data[8]); seed=anint(data[9]);
  NSTEPS=anint(tau0/dtau); sigm=sqrt(dtau);

  if(arguments.verbose){
    fprintf(stderr,"COLLECTED THE FOLLOWING DATA:\n");
    fprintf(stderr,"N0    = %d\n",N0);
    fprintf(stderr,"Nmax  = %d\n",Nmax);
    fprintf(stderr,"Tau0  = %f\n",tau0);
    fprintf(stderr,"dTau  = %f\n",dtau);
    fprintf(stderr,"rmin  = %f\n",rmin);
    fprintf(stderr,"rmax  = %f\n",rmax);
    fprintf(stderr,"nb    = %d\n",nb);
    fprintf(stderr,"V     = %f\n",V);
    fprintf(stderr,"Sigma = %f\n",sigm);
    fprintf(stderr,"d     = %d\n",DIM);
    fprintf(stderr,"seed  = %d\n",seed);
    fprintf(stderr,"Steps = %d\n",NSTEPS);
  }
  init_genrand((long)seed);

  /* ALLOC */
  //first walker[i][0] is the existance flag, the rest are coordinates
  walker=amat(Nmax,DIM*P+1);
  owalker=amat(Nmax,DIM*P+1);
  iT=10;
  Energy=amat(1,iT);
  for(i=0;i<iT;i++)
    Energy->dat[0][i]=0.0;
  Energy->dat[0][5]=dtau;
  Energy->dat[0][9]=V;

  GSWF=amat(nb+1,(DIM+1)*(2+P-1));

  if(arguments.verbose)
    fprintf(stderr,"\nALLOCATED:: (%d,%d)-MATRIX\n",walker->Ni,walker->Mj);

  initMat(walker,N0,RE,cmdtyp);
  initMat(owalker,N0,RE,cmdtyp);
  N[0]=N0; N[1]=N0;
  N0=0; nmean=1000.0;

  time=0.0;
  for(iT=0;iT<1000*(P*cmdtyp+1)*data[7];iT++){ //LET THE SYSTEM DIFFUSE
    walk(walker,owalker,Energy,R0,0.1,DIM,cmdtyp);
    time+=1.0;
    if(time>200.0 && arguments.verbose){
       fprintf(stderr,".");
       time=0.0;
    }
  }
  if(arguments.verbose)
    fprintf(stderr,"\n");

  time=0.0;
  for(iT=0;iT<NSTEPS;iT++){
    time+=dtau;
    walk(walker,owalker,Energy,R0,dtau,DIM,cmdtyp);
    N[1]=branch(walker,owalker,R0,dtau,Energy,cmdtyp,DIM,N);
    fprintf(outRE,"%7.4f %10d %10d\n",time,N[0],N[1]);
    if(time>tau0*0.5){
      count(walker,GSWF,nb,rmin,rmax,cmdtyp,DIM);
    }
    part=Energy->dat[0][3]/Energy->dat[0][4];
    fprintf(out,"%7.4f % 20.10f % 20.10f % 20.10f  % 20.10f \n",time,part,Energy->dat[0][0],Energy->dat[0][6],Energy->dat[0][5]);
  }
  setAxis(GSWF,nb,rmin,rmax,cmdtyp,DIM);     //THIS ALSO NORMALISES THE HISTOGRAM
  printMat(((cmdtyp==3||cmdtyp==4)?walker:GSWF),outWF,DIM*P);

  if(arguments.verbose)
    fprintf(stderr,"\nCLOSING STREAMS\n\n");

  if(cmdtyp==3)
    fprintf(stderr,"\n%s\n","load GSWF.dat\nplot3(GSWF(:,2),GSWF(:,3),GSWF(:,4),'bo'),hold on,plot3(GSWF(:,5),GSWF(:,6),GSWF(:,7),'rd')\n" );

  fclose(in);
  fclose(out);
  fclose(outWF);
  fclose(outRE);
  return 0;
}
