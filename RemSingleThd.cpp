//============================================================================
// Name        : RemSingleThd.cpp
// Author      : Tim
// Version     :
/* ********************************************
    program to do remodelling simulations
      using finite element analysis

This version implements a remodeling rate equation connected with

minimizing a weighted sum of the mass and the strain energy in a structure


           by Timothy P. Harrigan

    DISCLAIMER: I guarantee nothing.  The responsibility for
         any use or modification of this code lies entirely with
         the user. Check everything.

 *********************************************  */
/* modified to time step with the euler backward method  free */
/* trusses added 8 dec-1992 */

/* main program for finite element analysis  */
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <malloc.h>
#include <errno.h>
//#include <ulocks.h>
//#include <mpi.h>
#include <stdlib.h>
void diagtrans(int *, int);
void colheight(int, int *, int *);
double elstend(int *,double *,double *,int,double,double *);
void eldisdot(int *,double *,double *,int,double *);
void elcalc(int,int,int,int,int *);
void elrate(int,int,int,int,int *);

void elforce(int eltyp,int numel,int maxnodes,int elsize,int eldat[]);

int threedsolid(int nel,int nnodes,int nint,int matprp,double xx[][3],
		double b[][6],double s[],FILE* output,double volume[],int phindex);
int twodsolid(int nel,int nnodes,int itype,int nint,int matprp,double xx[][2],
		double b[][4],double s[],FILE* output,double volume[],int phindex);
int trussel(int nel,int matprp,double xx[][3],double s[],FILE* output);

int addel(int dof[],double elstf[],double matrix[],int diags[],int neq,int nel);
int blockalloc(int diags[],int ndof,FILE* output);
int plinsolv(double matrix[],double force[],int diags[],int ndof,double work[],int flag,FILE* output);
int fgetc(FILE *);
//overall simulation memory pointers

int maxcolht,ndoftot,nlcase;
double *coorptr, *stifptr, *loadptr, *workptr, *matptr, *nsptr;
double modexpnt,transexpnt,nulstate,dtime,gamnorm,ratenorm,dgamnorm,tol,itol;
double oblnorm;
double *gamptr,*dgamptr,*oldgamptr,*rateptr,*relvol,*ebmatrow;

int *colptr;
double *elnodf,*eldiag;
double *alpha;
FILE *output,*ebfile;
FILE *input,*disdens,*loadfile;

int eldone[4], netmaxdel,blalcdone,numrel;
int *dofptr,*eldat,*eldof,*remdiags;
int *active,netact;

int main()
{
char title[80];
int c,c2,c3,c4,c5,n,dof,doftot,dmy,eldmy,df,dofx,dofy,dofz;
double x,y,z,lmag,e,nu,ux,uy,uz,thic,volume,area,ddmy;
int elnod, mpnum,elsize,offset,memsize,nodnum;
int numel, eltyp, maxnodes, memelem,maxmelem;
int lnode, ldof,matlsize,mtype,maxdel,itype;
int ncase, netdof,nloads, fstatus, rstatus, wstatus, solstat;
int remflag,rind,r,ind,iter,itermax, doneflag,ebstart,nupper;
int ntnodes,masterdof,nelgps,matsize,nmatp,error;

int c3a,c4a;

int prnum,eldstart,tst,ct;

int **elarry;

double gfrac1,gfrac2,modinv,phi,phinorm;


/* **************************************************************************


     INPUT PHASE:

   get input data:
       number of nodes,
       DOF per node,
       number of load cases,
       number of element groups

   **********************************************************************   */

maxmelem=0;
netmaxdel=0;
itermax=20;
iter=0;
blalcdone=0;
/* open files */

input=fopen("/home/tim/test/FEINFILE", "r");
output=fopen("FEOUTFILE", "w");
disdens=fopen("FEDISDENS", "w");

ebfile=fopen("EBMAT", "w+");
loadfile=fopen("FELOADS", "w+");

/* read in overal control data */

c=0;
while((title[c++]=fgetc(input))!= '\n');

fscanf(input,"%d %d %d %d",&ntnodes,&masterdof,&nlcase,&nelgps);

// elarry = calloc(nelgps*2,sizeof(int *));
	elarry=new int*[nelgps*2];
/* print up a nice looking header for the output file */

fprintf(output," ********************************************************\n");
fprintf(output,"\n\n      A LINEAR ELASTIC FINITE ELEMENT PROGRAM\n");
fprintf(output,"  FOR BONE REMODELLING STUDIES AND ASSORTED OTHER PRANKS\n");
fprintf(output,"     AT TRUMAN MEDICAL CENTER/UMKC MED SCHOOL, KC MO\n");
fprintf(output,"     WRITTEN BY TIM HARRIGAN - FIRST REVISION 6/6/91\n\n");
fprintf(output,"     THIS VERSION USES MULTIPLE LOADS AND GAMMA VARIABLES\n");
fprintf(output,"     TIM HARRIGAN, 6 NOVEMBER 1992\n");

/* *************************************************************************
   get memory for nodal point coordinates
              and degrees of freedom
   **********************************************************************    */

//coorptr=calloc(ntnodes*3, sizeof(double));
//if(coorptr == NULL)printf("coordinate allocation failed (main)\n");
	coorptr= new double[ntnodes*3];

//dofptr =calloc(ntnodes*3, sizeof(int));
//if(dofptr == NULL)printf("degree of freedom allocation failed (main)\n");
	dofptr=new int[ntnodes*3];

	/* **********************************************************************
   read in nodal points
   **********************************************************************    */

fprintf(output," Nodal Point Degrees of Freedom and Coordinates\n\n");
fprintf(output,
   "  n  xdof  ydof  zdof         x             y            z\n\n");

ndoftot=1;
for(c=0;c<ntnodes;c++){
  fscanf(input, " %d %d %d %d %lf %lf %lf",&n,&dofx,&dofy,&dofz,&x,&y,&z);
    dmy=3*(n-1);
  if(n<=ntnodes){
    coorptr[dmy]=x;
    coorptr[dmy+1]=y;
    coorptr[dmy+2]=z;
    dofptr[dmy]=dofx;
    dofptr[dmy+1]=dofy;
    dofptr[dmy+2]=dofz;
    if(dofx>0)ndoftot++;
    if(dofy>0)ndoftot++;
    if(dofz>0)ndoftot++;
  } // if
fprintf(output," %5d %5d %5d %5d %12.9lf %12.9lf,%12.9lf\n",
   n,dofptr[dmy],dofptr[dmy+1],dofptr[dmy+2],x,y,z);
  }  //for(c=0

/* ***************************************************************
   get memory for column heights
   ***************************************************************     */

//colptr= calloc(ndoftot+2, sizeof(int));
colptr= new int[ndoftot+2];
if(colptr == NULL)printf("column pointer allocation failed (main)\n");

/* ***************************************************************

   LOOP THROUGH ELEMENT GROUPS

   read in element information
      to get the matrix profile
      and remodelling element indices
      then store the information on disk

   ***************************************************************    */

numrel = 0;
rind=0;

for(c=0;c<nelgps;c++){

  fscanf(input," %d %d %d %d", &numel,&eltyp,&maxnodes,&remflag);

  if(remflag == 1)numrel+= numel;

 //   elarry[c*2]=calloc(4,sizeof(int));
  	  elarry[c*2]= new int[4];
    eldat=elarry[c*2];

    if(eldat == NULL)printf("element group data allocation failed (main)\n");

  eldat[0]=numel;
  eldat[1]=eltyp;
  eldat[2]=maxnodes;
  eldat[3]=remflag;


/* define netmaxdel = maximum element degrees of freedom
   for later use in remodelling iterations */

  if(eltyp == 1)maxdel=6;
  if(eltyp/10 == 2){
      maxdel=maxnodes*2;}
  if(eltyp == 3){
      maxdel=maxnodes*3;}

  netmaxdel = (maxdel >netmaxdel) ? maxdel : netmaxdel;


/* *****************************************************************
   get memory for element info
   ******************************************************************   */

  if(eltyp == 1)elsize= 4 + maxnodes*4;
  if(eltyp == 3)elsize= 4 + maxnodes*4;
  if(eltyp/10 == 2)elsize= 4 + maxnodes*3;

  memelem = numel*elsize;

//  elarry[c*2 +1]=calloc(memelem, sizeof(int));
    elarry[c*2+1]= new int[memelem];
   eldat=elarry[c*2 + 1];

  if(eldat == NULL)printf("element data allocation failed (main)\n");

 /* ********************************************************************
     read in elements,
       send info to column height routine,
       store elements for later
    *******************************************************************  */

  for(c2=0;c2<numel;c2++){

  offset=c2*elsize;

  fscanf(input," %d %d %d", &eldat[offset], &eldat[offset+1],
          &eldat[offset+2]);

  eldat[offset+3] = (remflag == 1)? rind++ : -1;


/* ********************************************************************
   for each element
   eldat[offset]=element number
   eldat[offset+1]=material property set number
   eldat[offset+2]=number of nodes
   eldat[offset+3]=element density index (-1 if no remodelling)
   ********************************************************************   */

  fprintf(output,"\nElement Number %d degrees of freedom\n", c2+1);
  fprintf(output,"local node, global node, xdof, ydof, zdof\n\n");

  eldof= eldat + offset + (maxnodes+4);

  eldmy = 0;
  netdof=0;
    for(c3=0;c3<eldat[offset+2];c3++){

      fscanf(input," %d", &eldat[offset+c3+4]);

      dmy=3*(eldat[offset+c3+4]-1); /* index to nodal point dof */

      if(eltyp/10 != 2)eldof[eldmy++]=dofptr[dmy];
      eldof[eldmy++]=dofptr[dmy+1];
      eldof[eldmy++]=dofptr[dmy+2];

    fprintf(output," %5d %5d %5d %5d %5d\n",c3, eldat[offset+c3+4],
              dofptr[dmy],dofptr[dmy+1],dofptr[dmy+2]);

    }  //node


/* adjust matrix profile to account for this element */

    colheight(eldmy,eldof,colptr);


  } /* each element */


} /* loop through element groups */


/* *******************************************
   Read in number of material property sets.
      If the element is to be remodelled, the elastic modulus
      in its property set must be the modulus for volumetric
      density equal to one.

   *******************************************    */

fscanf(input," %d", &nmatp);

/* *********************************************
   get memory for material property information
   *********************************************    */

fprintf(output,"\n\n Material Property Sets:\n\n");

matlsize=3*nmatp;

//matptr= calloc(matlsize, sizeof(double));
matptr= new double[matlsize];
if(matptr == NULL)printf("material data allocation failed (main)\n");

/* **************************************
   read in material property information
   **************************************    */

for(c=0;c<nmatp;c++){
fscanf(input," %d", &mtype);
if(mtype == 2)fscanf(input," %d %lf %lf %lf",&c2, &e, &nu, &thic);
if(mtype == 2)fprintf(output,
"\n Set number= %d (%dD) E = %lf Nu =%lf Thickness=%lf\n",c2,mtype,e,nu,thic);

if(mtype == 3)fscanf(input," %d %lf %lf",&c2, &e, &nu);
if(mtype == 3)fprintf(output,
"\n Set number= %d (%dD) E = %lf Nu =%lf\n",c2,mtype,e,nu);

if(mtype == 1)fscanf(input," %d %lf %lf %lf",&c2, &e, &nu, &area);
if(mtype == 1)fprintf(output,
"\n Set number= %d (%dD) E = %lf Nu =%lf Area=%lf\n",c2,mtype,e,nu,area);



dmy=3*(c2-1);
matptr[dmy]=e;
matptr[dmy+1]=nu;
if(mtype == 2)matptr[dmy+2]=thic;
if(mtype == 1)matptr[dmy+2]=area;
} //nmatp

/* ***************************************
   get memory for volumetric density and rate info
   *************************************** */

 //  gamptr=calloc(numrel, sizeof(double)); /* actual gamma in iterations */
 //  if(gamptr == NULL)printf("gamma allocation failed (main)\n");

	gamptr= new double[numrel];

//   rateptr=calloc(numrel, sizeof(double)); /* rate of change */
//   if(rateptr == NULL)printf("rate allocation failed (main)\n");
     rateptr= new double[numrel];

//     dgamptr=calloc(numrel, sizeof(double));  /* increment in gamma */
//   if(dgamptr == NULL)printf("gamma increment allocation failed (main)\n");
     dgamptr= new double[numrel];

//   oldgamptr=calloc(numrel, sizeof(double));
		/* converged gamma from previous iteration */
//   if(oldgamptr == NULL)printf("old gamma allocation failed (main)\n");
     oldgamptr = new double[numrel];


//   relvol=calloc(numrel,sizeof(double)); /* element volumes */
//   if(relvol == NULL)printf("element volume allocation failed (main)\n");
      relvol = new double[numrel];

//   remdiags=calloc(numrel+1,sizeof(int));
//               /* addresses of remodelling matrix diagonals */
//   if(remdiags == NULL)
//             printf("remodelling diagonal index allocation failed (main)\n");
      remdiags= new int[numrel+1];

//      active = calloc(numrel,sizeof(int));
//   if(active == NULL)
//             printf("remodelling active index allocation failed (main)\n");
      active = new int[numrel];
      /* index for active elements  */
      netact=numrel;
      for(c=0;c<netact;c++)active[c]=c;


/* set up remodelling matrix as full */

for(c=0;c<=numrel;c++)remdiags[c]=c*(c+1)/2;

/* print out number of remodelling elements */

fprintf(output,"number of remodelling elements = %d\n",numrel);

/* ************************************
   read in initial density distribution
   ************************************  */

   for(c=0;c<numrel;c++)
        fscanf(input," %lf",gamptr+c);

   for(c=0;c<numrel;c++)oldgamptr[c]=gamptr[c];


/* *****************************************
   read in loading information:
       for each case,
       create the load vector and store it.
   ***************************************** */

//alpha = calloc(nlcase,sizeof(double));
alpha= new double[nlcase];
if(alpha == NULL)printf("alpha allocation failed\n");


fprintf(output,"\n\nLoad Cases:\n");

for(c=0;c<nlcase;c++){

/* ****************************
   get memory for loading case
   ****************************     */


//loadptr= calloc(ndoftot, sizeof(double));
loadptr = new double[ndoftot]();
if(loadptr == NULL)printf("load allocation failed (main)\n");


fscanf(input," %d %d %lf", &ncase, &nloads, &alpha[c]);
fprintf(output,"\n load case %d has %d nodal loads\n", ncase, nloads);
fprintf(output, "with a remodeling weight of %lf\n", alpha[c]);
fprintf(output,"\n node    direction    value\n");

for(c2=0;c2<nloads;c2++){
    fscanf(input," %d %d %lf", &lnode,&ldof,&lmag);
    netdof=dofptr[3*(lnode-1)+ldof];
    loadptr[netdof]= lmag;
    fprintf(output," %5d %5d %lf\n", lnode,ldof,lmag);
  } //each nodal load
    wstatus=fwrite(loadptr,sizeof(double),ndoftot,loadfile);

    //free(loadptr);
    delete[] loadptr;
    loadptr = nullptr;
} /* for c ... each load case */

/* ************************************
   read in remodelling simulation data
   ************************************ */
fscanf(input," %lf %lf %lf %lf %lf %lf"
   ,&modexpnt,&transexpnt,&nulstate,&dtime,&tol,&itol);

/* ****************************************************
   calculate addresses of diagonal elements in finite element matrix
   ************************************************  */

  diagtrans(colptr,ndoftot);

fprintf(output, "\nDiagonal Addresses:\n");
fprintf(output, "\nDOF:  Address: Minimum Row: \n");

for(c=0;c<=ndoftot;c++){
          r = colptr[c]+c+1-colptr[c+1];
          fprintf(output, "%d %d %d\n",c,colptr[c],r);  }

/* *****************************************
   allocate memory for the stiffness matrix
   *****************************************     */

  matsize= colptr[ndoftot];

//  loadptr= calloc((ndoftot+1)*nlcase, sizeof(double));
  loadptr=new double[(ndoftot+1)*nlcase];
  if(loadptr == NULL)printf("load allocation failed (main)\n");

  //workptr= calloc(ndoftot+1, sizeof(double));
  workptr= new double[ndoftot+1];

if(workptr == NULL)printf("diagonal allocation failed (main)\n");

/* print out overall matrix data */

fprintf(output,"\ntotal degrees of freedom = %d\n",ndoftot);
fprintf(output,"\ntotal matrix size = %d\n",matsize);

prnum=4;
//m_set_procs(prnum);

/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   main loop for remodelling iterations

   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
doneflag=0;

for(;;){

iter=0;

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   loop for euler backward iterations

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  */

do{

iter += 1;

//  stifptr = calloc(matsize, sizeof(double));
  stifptr = new double[matsize]();
if(stifptr == NULL)printf("stiffness matrix allocation failed (main)\n");

/* *****************************
   loop through element data
     to assemble overall matrix
   *****************************   */

eldstart=0;

printf(" starting assembly\n");

for(c=0;c<nelgps;c++){

/* **********************************
   get memory for element group data
   **********************************   */

eldone[0]=eldone[1]=eldone[2]=eldone[3] = 0;

  eldat=elarry[c*2];

  numel=eldat[0];
  eltyp=eldat[1];
  maxnodes=eldat[2];
  remflag=eldat[3];

   if(eltyp == 1)elsize= 4 + maxnodes*4;
   if(eltyp == 3)elsize= 4 + maxnodes*4;
   if(eltyp/10 == 2)elsize= 4 + maxnodes*3;

   memelem = numel*elsize;

  eldat= elarry[c*2 + 1];

	fflush (stdout);
	fflush (stderr);

 //m_fork(elcalc,eltyp,numel,maxnodes,elsize,eldat);
  elcalc(eltyp,numel,maxnodes,elsize,eldat);
//do{
//tst=1;
//for(ct=0;ct<prnum;ct++)tst *= ((eldone[ct] == 0)? 0:1);
//}
//while(tst == 0);

 } /* for (c=... element group */

/* *******************************

   solve for the displacements

   *******************************    */

printf(" starting solution\n");
fflush(stdout);

/* ******************************
   get matrix blocking info
   ******************************** */
   blockalloc(colptr,ndoftot,output);


/* triangularize the stiffness matrix */

solstat=plinsolv(stifptr,loadptr,colptr,ndoftot,workptr,1,output);
printf("solstat = %d\n",solstat);

rewind(loadfile);

/* solve for displacements for the selected load case */

for(c=0;c<nlcase;c++){

/* read in the force vector */

rstatus=fread(loadptr+c*(ndoftot+1),sizeof(double),ndoftot,loadfile);
if(rstatus != ndoftot){
	printf("force vector read failed (main)\n");
	fflush (stdout);}

solstat=plinsolv(stifptr,loadptr+c*(ndoftot+1),colptr,ndoftot,workptr,3,output);
} // for c... each load case

/* **********************************************************************
calculate the average strain energy density
summed overthe load cases
for each remodelling element
and the time stepping matrix
 ************************************************************************/

/* get memory for element nodal point forces */

// elnodf = calloc(numrel*netmaxdel*nlcase,sizeof(double));
elnodf= new double[numrel*netmaxdel*nlcase];
if(elnodf == NULL)printf("element nodal force allocation failed (main)\n");

for(c=0;c<nelgps;c++){

/* calculate element stiffness matrices and incremental matrix */
/* **********************************
   get memory for element group data
   **********************************   */

eldat=elarry[c*2];

  numel=eldat[0];
  eltyp=eldat[1];
  maxnodes=eldat[2];
  remflag=eldat[3];

   if(eltyp == 1)elsize= 4 + maxnodes*4;
   if(eltyp == 3)elsize= 4 + maxnodes*4;
   if(eltyp/10 == 2)elsize= 4 + maxnodes*3;

   memelem = numel*elsize;

eldat=elarry[c*2 + 1];

  if(remflag == 1){

/***********************************************************
   calculate element stiffness matrices, strain energy density,
   and incremental matrix rows in parallel
 ***********************************************************/

  printf("calculating element forces\n");
  fflush(stdout);

// eldone[0]=eldone[1]=eldone[2]=eldone[3] = 0;

 //m_fork(elforce,eltyp,numel,maxnodes,elsize,eldat);
elforce(eltyp,numel,maxnodes,elsize,eldat);
//do{
//tst=1;
//for(ct=0;ct<prnum;ct++)tst *= ((eldone[ct] == 0)? 0:1);
//}
//while(tst == 0);

 //m_kill_procs();

printf("calculating Euler backward matrix components \n");
fflush(stdout);

//   eldone[0]=eldone[1]=eldone[2]=eldone[3] = 0;

   //m_fork(elrate,eltyp,numel,maxnodes,elsize,eldat);
elrate(eltyp,numel,maxnodes,elsize,eldat);
//do{
//tst=1;
//for(ct=0;ct<prnum;ct++)tst *= ((eldone[ct] == 0)? 0:1);
//}
//while(tst == 0);

 //m_kill_procs();

  } /* if (remflag == 1...  */

} // for (c=... element groups

/************* c loops through element groups ***************/

  //free(elnodf);
  delete[] elnodf;
  elnodf=nullptr;

/* calculate new values for the densities */
/* first, let go of the memory for the overall stiffness matrix */

  //free(stifptr);
  delete[] stifptr;
  stifptr=nullptr;

/* then allocate memory for the (symmetric) time stepping matrix */

  nupper = numrel*(numrel+1)/2;

//stifptr = calloc(nupper,sizeof(double));
stifptr = new double[nupper]();
if(stifptr == NULL)printf("incremental density matrix allocation failed (main)\n");

//ebmatrow = calloc(numrel,sizeof(double));
ebmatrow = new double[numrel]();
if(ebmatrow == NULL)printf("Euler Backward row allocation failed (main)\n");

printf("reading back in Euler backward matrix ");
fflush(stdout);

for(c3=0;c3<netact;c3++){

  ebstart = c3*(netact)*sizeof(double);

  fstatus=fseek(ebfile,ebstart,SEEK_SET);

  if(fstatus != 0)printf("euler backward file seek failed\n");

  rstatus=fread(ebmatrow,sizeof(double), netact,ebfile);

  if(rstatus != netact){
	printf("Euler Backward matrix read failed (main)\n");
	fflush (stdout);}

    for(c4=c3;c4<netact;c4++){

    ind=c4*(c4+1)/2 + (c4-c3);
    c4a=active[c4];
    c3a=active[c3];

    stifptr[ind]= modexpnt * ebmatrow[c4]/
               (gamptr[c4a]*gamptr[c3a]*transexpnt);

    } /* for(c4... */

  } /* for(c3... */

//free(ebmatrow);
delete[] ebmatrow;
ebmatrow=nullptr;

fprintf(output,"\n\n    volumetric densities \n\n");
gamnorm=0.0;
dgamnorm=0.0;
ratenorm=0.0;
oblnorm=0.0;

gfrac1 = 1.0/transexpnt - 1.0;
gfrac2 = 1.0/transexpnt - 2.0;

/* adding on diagonal terms to remodeling matrix
  and defining the remodeling load vector */

for(c2=0;c2<netact;c2++){

  ind=c2*(c2+1)/2;

 c= active[c2];

ddmy = (transexpnt*relvol[c])/(modexpnt*dtime);
ddmy +=
 (transexpnt-modexpnt)*rateptr[c]*relvol[c]
  /(transexpnt*gamptr[c]*gamptr[c]);

/* minimizes total mass (not used currently)
ddmy +=
 (1.0-transexpnt)*nulstate*relvol[c]*pow(gamptr[c],gfrac2)
          /(transexpnt*modexpnt); */

stifptr[ind] += ddmy;

dgamptr[c2] = rateptr[c]*relvol[c]/gamptr[c];

/* minimizes total mass
dgamptr[c2] -= nulstate*relvol[c]*pow(gamptr[c],gfrac1)/transexpnt;
*/
/* minimizes integrated gamma */
dgamptr[c2] -= transexpnt*nulstate*relvol[c]/modexpnt;


dgamptr[c2] +=
   (oldgamptr[c]- gamptr[c])*relvol[c]*transexpnt/(modexpnt*dtime);

   oblnorm += dgamptr[c2]*dgamptr[c2];

} /* for (c2... */

if(netact< numrel){
for(c=netact;c<numrel;c++)dgamptr[c]=0.0; }


if(iter > itermax){

printf(" maximum iteration count exceeded \n");

/* check for positive definiteness of the density system */

/*  for(c=0;c<nupper;c++)stifptr[c] += nsptr[c];  */


  for(c=0;c<netact;c++){
     ind=c*(c+1)/2;
     c2=active[c];
     stifptr[ind] -= (transexpnt*relvol[c2])/(modexpnt*dtime);}

  blockalloc(remdiags,netact,output);

printf("starting solution to test stability\n\n");
fflush(stdout);

solstat=plinsolv(stifptr,dgamptr-1,remdiags-1,netact+1,workptr,1,output);

  printf("density evolution system is unstable if negative pivots are found\n");

  printf("Otherwise,the system is apparently stable - try smaller time steps\n");


  fclose(input);
  fclose(output);

  fclose(loadfile);
  fclose(ebfile);
  exit(0);

  } /* if(iter... */

if(doneflag == 1){

/* check for positive definiteness of the density system */

  for(c=0;c<netact;c++){
     c2=active[c];
     ind=c*(c+1)/2;
     stifptr[ind] -= (transexpnt*relvol[c2])/(modexpnt*dtime);}

  blockalloc(remdiags,netact,output);

printf("starting solution to test stability\n\n");
fflush(stdout);

solstat=plinsolv(stifptr,dgamptr-1,remdiags-1,netact+1,workptr,1,output);

  printf("Done: density evolution system is unstable if negative pivots are found\n");

  printf("Otherwise,the system is apparently stable - try smaller time steps\n");

  fclose(input);
  fclose(output);
  fclose(loadfile);
  fclose(ebfile);
  exit(0);

  } /* if(doneflag... */


  blockalloc(remdiags-1,netact+1,output);

solstat=plinsolv(stifptr,dgamptr-1,remdiags-1,netact+1,workptr,1,output);
solstat=plinsolv(stifptr,dgamptr-1,remdiags-1,netact+1,workptr,3,output);

printf("density status = %d\n",solstat);

/* * * * * * * * * * * * * * * * *

debugging code for the iteration loop  */

for(c2=0;c2<ntnodes;c2++){
   ux=uy=uz=0.;
   dmy=c2*3;
   if((df=dofptr[dmy]) != 0)ux=loadptr[df];
   if((df=dofptr[dmy+1]) != 0)uy=loadptr[df];
   if((df=dofptr[dmy+2]) != 0)uz=loadptr[df];
   fprintf(disdens,"%8d, %20.15lf, %20.15lf, %20.15lf\n",c2+1,ux,uy,uz);
}


for(c=0;c<numrel;c++){
fprintf(disdens,"   %8d   %15.12lf   %15.12lf %8d\n",
                       c,gamptr[c],dgamptr[c],active[c]);
}

/*  end of debugging code

* * * * * * * * * * * * * * * * * */


for(c=0;c<netact;c++){
   c2=active[c];
   if(gamptr[c2]+dgamptr[c] >1.0)dgamptr[c]=1.0-gamptr[c2];
   if(gamptr[c2]+dgamptr[c] < 0.000001)dgamptr[c]=0.000001-gamptr[c2];
   gamptr[c2] += dgamptr[c];
}

/* * * * * * * * * * * * * * * * * * * * *

 debugging code for the iteration loop

for(c2=0;c2<ntnodes;c2++){
   ux=uy=uz=0.;
   dmy=c2*3;
   if((df=dofptr[dmy]) != 0)ux=loadptr[df];
   if((df=dofptr[dmy+1]) != 0)uy=loadptr[df];
   if((df=dofptr[dmy+2]) != 0)uz=loadptr[df];
   fprintf(disdens,"%8d, %20.15lf, %20.15lf, %20.15lf\n",c2+1,ux,uy,uz);
}

for(c=0;c<numrel;c++){

fprintf(disdens,"   %8d   %15.12lf   %15.12lf\n",c,gamptr[c],dgamptr[c]);
}
fflush(disdens);

 debugging code for the iteration loop

* * * * * * * * * * * * * * * * * * * * * * */


for(c=0;c<numrel;c++){
         gamnorm += gamptr[c]*gamptr[c];
         dgamnorm += dgamptr[c]*dgamptr[c];
          }
gamnorm = sqrt(gamnorm/numrel);
dgamnorm = sqrt(dgamnorm/numrel);
oblnorm = sqrt(oblnorm/numrel);

printf("dgamma norm= %15.12f, gamma norm = %15.12f\n",dgamnorm, gamnorm);
fprintf(output,"dgamma norm= %15.12f, gamma norm = %15.12f\n",dgamnorm, gamnorm);

 //m_kill_procs();

//free(stifptr);
delete[] stifptr;
stifptr=nullptr;

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end of timestep iteration loop

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   */

} while ((dgamnorm/gamnorm) > itol);

modinv = 1.0/modexpnt;

phinorm=0.0;
ratenorm=0.0;

netact=0;

for(c=0;c<numrel;c++){

  rateptr[c] *= relvol[c]/gamptr[c];

/* minimized integrated gamma  */

  rateptr[c] -= nulstate*relvol[c]*transexpnt/modexpnt;

/* minimized total mass
  rateptr[c] -= nulstate*relvol[c]*pow(gamptr[c],gfrac1)/modexpnt;
*/

   if(gamptr[c]>0.000001 && gamptr[c]<1.0){
         phi = pow(gamptr[c],modinv);
         phinorm += phi*phi;
         ratenorm += rateptr[c]*rateptr[c]/(relvol[c]*relvol[c]);
         active[netact++]=c; }

   if(gamptr[c]<=0.000001){
       gamptr[c]=0.000001;
	if(rateptr[c]> 0.0)active[netact++]=c;}

   if(gamptr[c]>=1.0){
	gamptr[c]=1.0;
	if(rateptr[c]< 0.0)active[netact++]=c;}

 } /* for(c... numrel */

phinorm = sqrt(phinorm/netact);
ratenorm = sqrt(ratenorm/netact);

printf("rate norm = %15.12f, phinorm = %15.12f\n", ratenorm, phinorm);


 /* write out the displacement solution */

fprintf(output,"\n\n DISPLACEMENTS ");

fprintf(output,
     "\n             ux                   uy                   uz\n\n");
for(c2=0;c2<ntnodes;c2++){
   ux=uy=uz=0.;
   dmy=c2*3;
   if((df=dofptr[dmy]) != 0)ux=loadptr[df];
   if((df=dofptr[dmy+1]) != 0)uy=loadptr[df];
   if((df=dofptr[dmy+2]) != 0)uz=loadptr[df];
   fprintf(output,"%8d, %20.15lf, %20.15lf, %20.15lf\n",c2+1,ux,uy,uz);

   fprintf(disdens,"%8d, %20.15lf, %20.15lf, %20.15lf\n",c2+1,ux,uy,uz);
} /* for (c2... */

/* write our the densities as converged for this time step */

fprintf(output,"\n\n    volumetric densities \n\n");

for(c=0;c<numrel;c++){
oldgamptr[c]=gamptr[c];
dgamptr[c]=0.0;

         phi = pow(gamptr[c],modinv);
fprintf(output,"   %8d   %15.12lf   %15.12lf\n",c,phi,rateptr[c]);

fprintf(disdens,"   %8d   %15.12lf   %15.12lf\n",c,phi,rateptr[c]);
} //for c=0  numrel

fprintf(output, "  density norm = %15.12lf, rate norm = %15.12lf\n"
         ,phinorm,ratenorm);
fflush(disdens);



if(ratenorm/phinorm < tol)doneflag=1;


/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
   end of overall remodelling simulation loop
   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

} /* for loop */


/* end of program   */

}

/****************************************************************************
*****************************************************************************/

void diagtrans (int mcolptr[],int ndoftot)
{
int c,h,hlast;
maxcolht=0;
mcolptr[1]=0;
hlast=mcolptr[2];
mcolptr[2]=1;

for(c=3;c<ndoftot+2;c++){
  maxcolht = (maxcolht >hlast) ? maxcolht: hlast;
  h=mcolptr[c];
  mcolptr[c]=mcolptr[c-1]+hlast+1;
  hlast=h;
  }
}
/****************************************************************************
*****************************************************************************/

void colheight (int edof,int eldof[],int gbldof[])
{
int c,min,colht,dmy;
/* find smallest dof */
  min=10000000;
for(c=0;c<edof;c++){
  if(eldof[c]>0){ min= (eldof[c]<min) ? eldof[c] : min;
  }
}
/* update column data */
for(c=0;c<edof;c++){
  if((dmy=eldof[c])>0){
    colht=dmy-min;
    gbldof[dmy] = (gbldof[dmy] > colht) ? gbldof[dmy] : colht;
  }

}
}
/****************************************************************************
*****************************************************************************/

/* **************************************************
   routine to calculate element strain energy density
   ************************************************** */

double elstend(int eldof[],double elstf[],double disp[],int neldof,double volume,double work[])
{
int c, c2,dind,elind;
double net;

net=0.;

for(c=0;c<neldof;c++){

work[c]=0.;

   if(eldof[c] != 0){

      for(c2=0;c2<neldof;c2++){

/* ********************
   element matrix index
   ******************** */

          elind= (c>c2)? (c*(c+1)/2 + c2): (c2*(c2+1)/2 + c);
          dind= eldof[c2];  /*  displacement index */
          if(eldof[c2] !=0)work[c] += elstf[elind]*disp[dind];
        } /* for(c2 ... */

      }  /* if(eldof[c]... */

    }  /* for(c...  */


for(c=0;c<neldof;c++){

   if(eldof[c] != 0){

       dind=eldof[c];

       net += work[c]*disp[dind];

       }  /* if(eldof[c]... */

    }  /* for(c...  */

   net /= 2.0*volume;

return(net);

}

/****************************************************************************
*****************************************************************************/

/* **************************************************
   routine to calculate element matrix - total displacement dot products
   ************************************************** */

void eldisdot(int eldof[],double elstf[],double disp[],int neldof,double work[])
{
int c, c2,dind,elind;
double net;

for(c=0;c<neldof;c++){

work[c]=0.;

   if(eldof[c] != 0){

      for(c2=0;c2<neldof;c2++){

/* ********************
   element matrix index
   ******************** */

          elind= (c>c2)? (c*(c+1)/2 + c2): (c2*(c2+1)/2 + c);
          dind= eldof[c2];  /*  displacement index */
          if(eldof[c2] !=0)work[c] += elstf[elind]*disp[dind];

        } /* for(c2 ... */

      }  /* if(eldof[c]... */

    }  /* for(c...  */

return;

} /* end of routine */

/****************************************************************************
*****************************************************************************/
/****************************************************************************
*****************************************************************************/
/* routine to do parallel backsubstitutions to get Euler Backward Matrix
   components  */

void elrate(int eltyp,int numel,int maxnodes,int elsize,int eldat[])
//int eltyp,numel,maxnodes,elsize,*eldat;
{
int maxdel,memsize,*eldof,offset;
double *elstf,*elstdm,*elnpc,volume,netdot;
int id,c2,c3,c4,itype,nodnum,solstat,rind,offorc,wstatus,fstatus,ebstart;
int c5,c2a,c3a;
extern int nlcase;
double *elload,*ebmatrow;
extern FILE *output;

//   elload=calloc(ndoftot,sizeof(double));
   elload = new double[ndoftot];
if(elload == NULL){
           //m_lock();
           printf("element load vector allocation failed (elrate)\n");
           //m_unlock();
           } //elload

  if(eltyp/10 == 2){
      maxdel=maxnodes*2;}
  if(eltyp == 3){
      maxdel=maxnodes*3;}

/* *************************
   for each element
   eldat[offset]=element number
   eldat[offset+1]=material property set number
   eldat[offset+2]=number of nodes
   eldat[offset+3]=remodelling index
   *************************   */

/* get memory for rows of (nonsymmetric) Euler Backward matrix */

//ebmatrow = calloc(numrel,sizeof(double));
ebmatrow = new double[numrel]();

if(ebmatrow == NULL){

         //m_lock();
         printf("Euler Backward row allocation failed (elrate)\n");
         //m_unlock();
         } // if(ebmatrow

for(c2=0;c2<netact;c2++){

//c2=m_next();
/* for now, all the remodeling elements must be in one group */

//if(c2 >= netact){
//id = m_get_myid();
//eldone[id]=1;
//m_lock();
//printf("elrate proc = %d\n",id);
//fflush(stdout);
//m_unlock();
//free(ebmatrow);
//free(elload);
//  return;
//}

	c2a=active[c2];

	for(c5=0;c5<nlcase;c5++){

		offset=c2a *elsize;
		eldof= eldat + offset + (maxnodes+4);
		offorc= c2a *netmaxdel + netmaxdel*numrel*c5;;

/* take dot products with calculated displacements
   to get strain energy density and rate of change */

		for(c3=0;c3<ndoftot;c3++)elload[c3]=0.0;

		for(c3=0;c3<maxdel;c3++){
			rind = eldof[c3];
			elload[rind]=elnodf[offorc+c3];} //for(c3=0
			solstat=plinsolv(stifptr,elload,colptr,ndoftot,workptr,3,output);

		for(c3=c2;c3<netact;c3++){
			c3a = active[c3];

			netdot = 0.0;
			offset=c3a *elsize;
			eldof= eldat + offset + (maxnodes+4);
			offorc= c3a *netmaxdel + netmaxdel*numrel*c5;

			for(c4=0;c4<maxdel;c4++){
				rind = eldof[c4];
				netdot +=  elload[rind]*elnodf[offorc+c4]; }  //for(c4...
			ebmatrow[c3]+=alpha[c5]*netdot;

		} // for(c3...
	} /* for(c5: load cases */

	ebstart = c2*netact*sizeof(double);

//m_lock();

	fseek(ebfile,ebstart,SEEK_SET);

	wstatus=fwrite(ebmatrow,sizeof(double),netact,ebfile);

//m_unlock();

/* printf("rind= %d,rateptr[rind]= %lf,eldiag[rind]= %lf\n",
        rind,rateptr[rind],eldiag[rind]);  */

	for(c3=0;c3<numrel;c3++)ebmatrow[c3]=0.0;

  	}  /* for(;;... c2 */

//free(ebmatrow);
delete[] ebmatrow;
ebmatrow = nullptr;
delete[] elload;
elload = nullptr;
//free(elload);
return;
} // end of routine

