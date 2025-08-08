#include <stdio.h>
#include <math.h>
extern FILE *output;
int nmatp,mtype,matlsize;
extern int *eldat;
double *volume;
int twodsolid(int nel,int nnodes,int itype,int nint,double matprops[4],double xx[][2],double s[],FILE *output,
		double *volume,double gammar,double gexpnt);
int threedsolid(int nel,int nnodes,int nint,double matprops[3],double xps[][3],
        double s[],FILE *output,double *volume,double gammar,double gexpnt);
int trussel(int nel,int matprp,double xx[][3],double s[],FILE *output);
int addel(int dof[],double elstf[],double matrix[],int diags[],int neq,int nel);

void elcalc(int eltyp,int numel,int maxnodes,int elsize,int eldat[])
{

extern double *matptr,*gamptr,*coorptr,*stifptr,*relvol;
extern int *colptr,ndoftot,*dofptr;
extern double modexpnt,transexpnt;
extern FILE *output;


int maxdel,memsize,*eldof;
double *elstf,*elstdm,elnpc,elnpc2[9][2],elnpc3[21][3],volume;
double elstdm2[9][4],elstdm3[21][6];
int c2,c3,itype,nodnum,offset;
int remind,id,estat,dfind;
int eNodDF[63];
double *s;
/* subroutine to do parallel calculation of element matrices */

/* ******************************************
   get memory for element stiffness matrices
   ******************************************   */
/* truss */

  if(eltyp == 1){
      maxdel=maxnodes*3;
      memsize=maxdel*(maxdel+1)/2 + 3*maxdel;
      elstf=new double[memsize];

      	  if(elstf == NULL){
          printf("element stiffness matrix allocation failed (elcalc)\n");
          } //if(elstf == Null
          
           }  //if(eltyp == 1

/* plane */

  if(eltyp/10 == 2){
      maxdel=maxnodes*2;
      memsize=maxdel*(maxdel+1)/2 + 6*maxdel;
      elstf= new double [memsize];

      	  if(elstf == NULL){
      		  printf("element stiffness matrix allocation failed (elcalc)\n");
      		  }  //if(elstf==Null
      		  
           } //if(eltyp/10

/* volume */

  if(eltyp == 3){
      maxdel=maxnodes*3;
      memsize=maxdel*(maxdel+1)/2 + 9*maxdel;
 
      elstf= new double[memsize];

      if(elstf == NULL){ 
          	  printf("element stiffness matrix allocation failed (elcalc)\n");
          	  } //if(elstf == Null

  	  	  }  //if(eltyp == 3


/* *************************
   for each element
   eldat[offset]=element number
   eldat[offset+1]=material property set number
   eldat[offset+2]=number of nodes
   eldat[offset+3]=remodelling index
   *************************   */
   //offset is an offset within eldat - eldat points to an array of individual element data
	  /* ********************************************************************
	     for each element, setting up the data for 2D solids
	     eldat[offset]=element number
	     eldat[offset+1]=material property set number
	     eldat[offset+2]=number of nodes
	     eldat[offset+3]=element density index (-1 if the element is not remodeling)
	     ********************************************************************   */

  for(c2=0;c2<numel;c2++){ //each element loop
	  offset=c2*elsize;
	  eldof= eldat + offset + (maxnodes+4);
	  remind= eldat[offset+3];

	//set up material data to pass to the element threads
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

	  if(eltyp == 1){
		 dfind=0;
	  	 for(c3=0;c3<eldat[offset+2];c3++){
			  nodnum=eldat[offset+c3+4];
	      		elnpc3[c3][0]=coorptr[3*(nodnum-1)];
	     		elnpc3[c3][1]=coorptr[3*(nodnum-1)+1];
	     		elnpc3[c3][2]=coorptr[3*(nodnum-1)+2];
		  		eNodDF[dfind++]=dofptr[(nodnum-1)*3+0];
		  		eNodDF[dfind++]=dofptr[(nodnum-1)*3+1];
		  		eNodDF[dfind++]=dofptr[(nodnum-1)*3+2];
	  	  	  	  } // for(c3..
	  	  estat=trussel(eldat[offset],eldat[offset+1],elnpc3,elstf,output);
	  	  addel(eNodDF,elstf,stifptr,colptr,ndoftot,eldat[offset+2]*3);
  	  }  // if(eltyp == 1...

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
             } //for(c3..
	       
	  	  estat=threedsolid(eldat[offset],eldat[offset+2],3,matpass,elnpc3,s,output,&volume,gammar,transexpnt);
	  	  if(estat != 0)printf("negative jacobian in element # %d\n",eldat[offset]);
	  	  addel(eNodDF,s,stifptr,colptr,ndoftot,eldat[offset+2]*3);
  	  }  // if(eltyp == 3... 

  if(eltyp/10 == 2){   //2D element
     dfind=0;
	 for(c3=0;c3<eldat[offset+2];c3++){
		  nodnum=eldat[offset+c3+4];
     	    elnpc2[c3][0]=coorptr[3*(nodnum-1)+1];
      	    elnpc2[c3][1]=coorptr[3*(nodnum-1)+2];
		    eNodDF[dfind++]=dofptr[(nodnum-1)*3+1];
		    eNodDF[dfind++]=dofptr[(nodnum-1)*3+2];
		    }  // for(c3...

	  	  int itype=eltyp-20;
	  	  estat=twodsolid(eldat[offset],eldat[offset+2],itype,3,matpass,elnpc2,s,output,&volume,gammar, transexpnt);
	  	  if(estat != 0)printf("negative jacobian in element # %d\n",eldat[offset]);
	  	  addel(eNodDF,s,stifptr,colptr,ndoftot,eldat[offset+2]*2);
  	  }  // if(eltyp/10 == 2... 


   	   if (remind>=0)relvol[remind]=volume;

  }  /* (c2) for loop  */

} /* routine end */

