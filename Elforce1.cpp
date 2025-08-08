/*
 * Elforce1.cpp
 *
 *  Created on: Aug 7, 2025
 *      Author: tim
 */
/****************************************************************************
****************************************************************************/
#include <stdio.h>
#include <math.h>
extern FILE *output;
//int nmatp,mtype,matlsize;
//int *eldat;
//double *volume;
int twodsolid(int nel,int nnodes,int itype,int nint,double matprops[4],double xx[][2],double s[],FILE *output,
		double *volume,double gammar,double gexpnt);
int threedsolid(int nel,int nnodes,int nint,double matprops[3],double xps[][3],double s[],FILE *output,double *volume,double gammar,double gexpnt);
int trussel(int nel,int matprp,double xx[][3],double s[],FILE *output);
int addel(int dof[],double elstf[],double matrix[],int diags[],int neq,int nel);
void eldisdot(int eldof[],double elstf[],double disp[],int neldof,double work[]);
double elstend(int eldof[],double elstf[],double disp[],int neldof,double volume,double work[]);

/* routine to do parallel dot products to get strain energy density
   and element nodal point force info  */

void elforce(int eltyp,int numel,int maxnodes,int elsize,int eldat[])
//int eltyp,numel,maxnodes,elsize,*eldat;
{
int maxdel,memsize,*eldof,offset;
double *elstf,*elstdm,*elnpc,elnpc2[9][2],elnpc3[21][3],volume;
int c2,c3,itype,nodnum,solstat,rind,elfoff,id,index,dfind;
extern int nlcase;
double elstdm2[9][4],elstdm3[21][6];
double *s;
int eNodDF[63];

extern double *matptr,*gamptr,*coorptr,*stifptr,*relvol,*rateptr,*alpha,*loadptr;
extern int *colptr,ndoftot,*dofptr,netmaxdel;
extern double modexpnt,transexpnt,*elnodf;
extern FILE *output;
extern int blalcdone,numrel;

/* ******************************************
   get memory for element stiffness matrices
   ******************************************   */

  if(eltyp/10 == 2){
      maxdel=maxnodes*2;
      memsize=maxdel*(maxdel+1)/2 + 6*maxdel;
      elstf= new double [memsize];
      if(elstf == NULL)printf("element stiffness matrix allocation failed (elforce)\n"); //if(elstf
      elstdm=elstf + maxdel*(maxdel+1)/2;
      elnpc=elstdm+4*maxdel; } // if(eltyp...

  if(eltyp == 3){
      maxdel=maxnodes*3;
      memsize=maxdel*(maxdel+1)/2 + 9*maxdel;

      elstf = new double[memsize];
      if(elstf == NULL)printf("element stiffness matrix allocation failed (elforce)\n");
      elstdm=elstf + maxdel*(maxdel+1)/2;
      elnpc=elstdm+6*maxdel; } //if(eltyp

/* *************************
   for each element
   eldat[offset]=element number
   eldat[offset+1]=material property set number
   eldat[offset+2]=number of nodes
   eldat[offset+3]=remodelling index
   *************************   */
	//set up material data to pass to the element threads
    //offset is an offset within eldat - eldat points to an array of individual element data

  for(c2=0;c2<numel;c2++){ //replaces for(;;

  offset=c2*elsize;
  eldof= eldat + offset + (maxnodes+4);

  if(eldat[offset+3]>=0){

	  int matnumber=eldat[offset+1];
	  double matpass[4];
	  int elemdof=0;
	  if (eltyp/10==2){
      elemdof=2*eldat[offset+2];
	  s=new double[elemdof*(elemdof+1)/2];
	  matpass[0]=matptr[(matnumber-1)*4];
	  matpass[1]=matptr[(matnumber-1)*4+1];
	  matpass[2]=matptr[(matnumber-1)*4+2];
	  matpass[3]=matptr[(matnumber-1)*4+3];
	  }
	  if(eltyp==3){
	   elemdof=3*eldat[offset+2];
	   s=new double[elemdof*(elemdof+1)/2];
	   matpass[0]=matptr[(matnumber-1)*4];
	   matpass[1]=matptr[(matnumber-1)*4+1];
	   matpass[2]=matptr[(matnumber-1)*4+2];
	  }
	  int denind=eldat[offset+3];

	  double gammar= denind>=0 ? gamptr[denind] : -1;


	  if(eltyp == 3){
		  dfind=0;
		  for(c3=0;c3<eldat[offset+2];c3++){
			  nodnum=eldat[offset+c3+4];
		      elnpc3[c3][0]=coorptr[3*(nodnum-1)];
		      elnpc3[c3][1]=coorptr[3*(nodnum-1)+1];
		      elnpc3[c3][2]=coorptr[3*(nodnum-1)+2];
		  	  eNodDF[dfind++]=dofptr[(nodnum-1)*3+0];
		  	  eNodDF[dfind++]=dofptr[(nodnum-1)*3+1];
		  	  eNodDF[dfind++]=dofptr[(nodnum-1)*3+2];


		  } // for(c3

		  threedsolid(eldat[offset],eldat[offset+2],3,
				 matpass,elnpc3,s,output,&volume,gammar,transexpnt);

	      }  /* if(eltyp == 3... */

	  if(eltyp/10 == 2){
		  dfind=0;
		  for(c3=0;c3<eldat[offset+2];c3++){
			  nodnum=eldat[offset+c3+4];
	      	  elnpc2[c3][0]=coorptr[3*(nodnum-1)+1];
	      	  elnpc2[c3][1]=coorptr[3*(nodnum-1)+2];
			  eNodDF[dfind++]=dofptr[(nodnum-1)*3+1];
			  eNodDF[dfind++]=dofptr[(nodnum-1)*3+2];

		  	  	  } // for(c3
		  itype=eltyp-20;

		  twodsolid(eldat[offset],eldat[offset+2],itype,3,
				  matpass,elnpc2,s,output,&volume,gammar,transexpnt);
	  }  /* if(eltyp/10 == 2... */

/* take dot products with calculated displacements
   to get strain energy density and rate of change */

	  int rind=eldat[offset+3];

	  rateptr[rind]=0.0;

	  for(c3=0;c3<nlcase;c3++){

		  index=c3*(ndoftot+1);
		  rateptr[rind] += alpha[c3]*
				  elstend(eldof,elstf,&loadptr[index],maxdel,volume,elnpc);

		  elfoff = rind*netmaxdel + c3*netmaxdel*numrel;

		  eldisdot(eldof,elstf,&loadptr[index],maxdel,elnodf+elfoff);

/* printf("rind= %d,rateptr[rind]= %lf,eldiag[rind]= %lf\n",
        rind,rateptr[rind],eldiag[rind]); */
	  } //for (c3=0 .. nlcase

  	} /* if(eldat[offset+3]>=0 ...  */
  } //for(c2=0...

//free(elstf);
delete[] elstf;
elstf=nullptr;

} /* end of routine */



