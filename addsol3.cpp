/* this file contains routines to assemble and solve a finite element matrix */

/* this is a single threaded implamantation of the parallel routine */

#include <stdio.h>
#include <malloc.h>

//int plinsolv(double *,double *,int *,int,double *,int,FILE *);
int plinsolv(double matrix[],double force[],int diags[],int ndof,double work[],int flag,FILE *output);

//int addel(int *,double *,double *,int *,int,int);
int addel(int dof[],double elstf[],double matrix[],int diags[],int neq,int nel);

//void blockreduce(double *,int *,int,double *, int *);
void blockreduce(double matrix[],int diags[],int ndof,double work[],int dflag[]);

//int blockalloc(int *,int,FILE *);
int blockalloc(int diags[],int ndof,FILE *output);

int *doneflag,*mincpl,*bstartrow;

int nblocks,solverr;

int blstart;

/* *****************
   addel
   ***************** */

/* routine to add element stiffness to global stiffness matrix */

/* assumes an element matrix num,bering scheme which is:

 elstf[0]
 elstf[1] elstf[2]
 elstf[3] elstf[4] elstf[5]
 elstf[6] elstf[7] elstf[8] elstf[9]
 .
 .
 .

and a global matrix numbering scheme which is

 matrix[0]
 matrix[2] matrix[1] 
           matrix[4] matrix[3]
 matrix[8] matrix[7] matrix[6] matrix[5]
 .
 .
 .
   with diags holding the addresses of the diagonals in the matrix 
*/

int addel(int dof[],double elstf[],double matrix[],int diags[],int neq,int nel)
//int *dof,neq,nel,*diags;
//double *matrix,*elstf;
{
int i1,i2,ind2,elind,df1dmy, df2dmy; /* indices */
for(i1=0;i1<nel;i1++){
if(dof[i1] !=0 ){
   for(i2=0;i2<nel;i2++){
      if(dof[i2]!=0){
        if(dof[i2]<=dof[i1]){
        df1dmy=dof[i1];
        df2dmy=dof[i2];
        ind2=diags[df1dmy]+(df1dmy-df2dmy);
         elind=(i1>=i2)? i1*(i1+1)/2 + i2 : i2*(i2+1)/2 + i1 ;
         matrix[ind2]+=elstf[elind];
         } /* if (dof[i2 ... */
      } /* if (dof[i2]... */
    } /* for(i2=0... */
} /* if(dof[i1]... */
} /* for (i1...) */
return(0);
} /* routine end */

/* routine to divide up the matrix into blocks which have roughly the
     same amount of work in each  */

int blockalloc(int diags[],int ndof,FILE* output)
//int *diags,ndof;
//FILE *output;
{
int c,c2,c3,r,ncalc,ncb,npblock,h,i1,havg,hmax,blocksize;
extern int *mincpl,*bstartrow,*doneflag;
extern int blalcdone;

nblocks=1;


hmax=0;
havg=diags[ndof]/ndof;

ncalc=0;

/* calculate total number of multiplies  */

for(c=0;c<ndof;c++){

   h=diags[c+1]-1-(diags[c]+1);
   ncalc += h*(h+1)/2 + h;  

   hmax = (h>hmax) ? h : hmax;

    } /* for(c= ... */

/* blocksize = havg/4; */
blocksize=ndof;

nblocks = ndof/blocksize;

nblocks += (ndof%blocksize == 0) ? 0:1;

fprintf(output,"Matrix Blocking Information:\n\n");
fprintf(output," %d blocks, hmax = %d\n\n",nblocks,hmax);


if(blalcdone != 0){
//free(mincpl);
delete[] mincpl;
mincpl=nullptr;
//free(bstartrow);
delete[] bstartrow;
bstartrow=nullptr;
//free(doneflag);
delete[] doneflag;
doneflag=nullptr;
}

//mincpl =calloc(nblocks+1, sizeof((int)0));
mincpl = new int[nblocks+1];
//bstartrow =calloc(nblocks+1, sizeof((int)0));
bstartrow=new int[nblocks+1];
//doneflag =calloc(nblocks, sizeof((int)0));
doneflag=new int[nblocks];

blalcdone = 1;


/* multiplies per block  */

npblock=ncalc/nblocks;

ncb=0;
c2=1;

/* divide up the matrix */


for(c=0;c<nblocks;c++)bstartrow[c]=c*(ndof/nblocks);

bstartrow[0]=1;
bstartrow[nblocks]=ndof;


/* establish block coupling  */


for(c=0;c<nblocks;c++){

mincpl[c]=c;

   for(c2=bstartrow[c];c2<bstartrow[c+1];c2++){

      r = diags[c2]+c2+1-diags[c2+1];

      for(c3=0;c3<c;c3++){
       if (r<bstartrow[c3+1])mincpl[c]=(mincpl[c]<c3)? mincpl[c]:c3;
       } /* for(c3... (each previous block)*/
      }  /* for(c2... (each column) */
   }   /* for(c...  (each block)*/

fprintf(output,"block  Minimum equation   Minimum Coupling Block\n");
for(c=0;c<nblocks;c++)
fprintf(output,"%5d   %5d     %5d\n",c, bstartrow[c],mincpl[c]);
fflush(output);

return(0);
}


/* **********************************************
    c routine to solve a set of linear equations
   ********************************************** */


int plinsolv(double matrix[],double force[],int diags[],int ndof,double work[],int flag,FILE* output)
//double *matrix,*force,*work;
//int *diags,ndof,flag;
//FILE *output;
{
int i,i1,i2,i3,sdot1,sdot2,edot,ldot,h,h2,dd,r,istart,n,m,cc,bstat;
double b,c;
int nprocs;
extern int *doneflag,nblocks,solverr;

int ibdmy;
double blasdmy;

solverr = 0;

if(flag == 1 || flag ==2){

/* factorization */
/* parallel calls */

for(i=0;i<nblocks;i++)doneflag[i]=0;

/*
m_fork(blockreduce,matrix,diags,ndof,work,doneflag);
  */

blockreduce(matrix,diags,ndof,work,doneflag);


/* while ((doneflag[nblocks-1]==0) && (solverr == 0)); */

/*
m_kill_procs();
   */

if (solverr != 0)return(-solverr);


} /* if (flag == 1 ...  */

/***** reducing right-hand side vector *****/

if(flag == 2 || flag == 3){

for(i1=0;i1<ndof;i1++){
        sdot1=diags[i1]+1;
        sdot2=i1-1;
        h=diags[i1+1]-diags[i1]-1;
        c=0;
        if(h>0){

            
                for(i2=0;i2<h;i2++)c+=matrix[sdot1+i2]*force[sdot2-i2];
             

               /* i2=1;
                ibdmy=-1;
                blstart=sdot2-h+1;
                c = sdot(&h,&matrix[sdot1],&i2,&force[blstart],&ibdmy); */

        } /* if(h>0)... */

        force[i1]-=c;
} /* for(i1... */

/* back-substitution */

for(i1=0;i1<ndof;i1++)force[i1]*=work[i1];

if(ndof == 1)return(1);

for(i1=ndof-1;i1>0;i1--){

   sdot1=diags[i1]+1;
   sdot2=i1-1;

   h=diags[i1+1]-diags[i1]-1;

   if(h>0){


     for(i2=0;i2<h;i2++)force[sdot2-i2]-=matrix[sdot1+i2]*force[i1];


    /* i2=1;
     ibdmy=-1;
     blasdmy= -force[i1];
     blstart= sdot2-h+1;

     saxpy(&h,&blasdmy,&matrix[sdot1],&i2,&force[blstart],&ibdmy);  */
 }

  } /* for(i1... */

} /* if (flag == 2... */

return(0);

}/* routine end */




void blockreduce(double matrix[],int diags[],int ndof,double work[],int dflag[])
//double *matrix,*work;
//int *diags,ndof,*dflag;
{
int i,i1,i2,i3,sdot1,sdot2,edot,ldot,h,h2,dd,r,istart,m,n,cc,id;
double b,c;
extern int solverr;


/*
id = m_get_myid();

m_lock();
printf("solve id = %d\n",id);
fflush(stdout);
m_unlock();
*/

for(n=0;n<nblocks;n++){

/*
n= m_next();


if((n>=nblocks) || (solverr != 0))return;
*/

   if(mincpl[n]<n){

/*  ********************************************
    factorize current block with previous blocks
    ******************************************** */ 

      for(m=mincpl[n];m<n;m++){


/*          while(dflag[m] == 0);  */


  /* wait until the needed block is done */

          for(i1=bstartrow[n];i1<bstartrow[n+1];i1++){

             r = diags[i1]+i1+1-diags[i1+1];

             if(r < bstartrow[m+1]){

               istart = (r<bstartrow[m]) ? bstartrow[m]-r-1 : 0;

               h=diags[i1+1]-1-(diags[i1]+1);

               if(h>0 && istart<h){
                   for(i2=istart;i2<(bstartrow[m+1]-r-1);i2++){

 /* set up dot products for g(ij) factors */

                     sdot1=diags[i1-h+i2]+1;  /* other column */
                     sdot2=diags[i1+1]-i2-1;  /* current column */

                     h2=diags[i1-h+i2+1]-sdot1; /* height of other column */

                       if(h2>0){

                         ldot=(i2+1>h2)?h2:i2+1;

/* dot product length is the smaller of the two quantities */

                         c=0.;

/* serial */                      
                         for(i3=0;i3<ldot;i3++)
                             c+=matrix[sdot1+i3]*matrix[sdot2+i3];

 /* substitute BLAS call
                          i3=1;
                          c=sdot(&ldot,&matrix[sdot1],&i3,&matrix[sdot2],&i3); */

/* calculate g(ij) factor  */

                         matrix[sdot2-1]-=c;

                       } /* if(h2>0)... */

                  } /* for(i2=bstartrow[... */

                }  /* if(h>0... */

             } /* if(r<bstartrow[...   */


          }  /* for(i1... */

      } /* for(m=mincpl... */

   }  /* if(mincpl...  */

/* ************************************
   loop to reduce the block by itself 
   ************************************ */

for(i1=bstartrow[n];i1<bstartrow[n+1];i1++){

/* generate g(ij) factors from bathe's book */

/* h=column height -1 */

   r = diags[i1]+i1+1-diags[i1+1];   /* minimum coupled equation */

   istart = (r<bstartrow[n]) ? bstartrow[n]-r-1 : 0;

   h=diags[i1+1]-1-(diags[i1]+1);

   if( h>0 && istart<h){
      for(i2=istart;i2<h;i2++){

        /* set up dot products for g(ij) factors */

        sdot1=diags[i1-h+i2]+1;  /* other column */
        sdot2=diags[i1+1]-i2-1;  /* current column */

        h2=diags[i1-h+i2+1]-sdot1; /* height of other column */

        if(h2>0){

                ldot=(i2+1>h2)?h2:i2+1;

                /* dot product length is the smaller of the two quantities */

                c=0.;

/* serial code */ 

                for(i3=0;i3<ldot;i3++)c+=matrix[sdot1+i3]*matrix[sdot2+i3]; 
  
 /* substitute with BLAS call 
                i3=1;
                c=sdot(&ldot,&matrix[sdot1],&i3,&matrix[sdot2],&i3); */

                /* calculate g factor  */

                matrix[sdot2-1]-=c;

        } /* if(h2>0)... */

        } /* for(i2=istart... */

    } /* if(h>0... */

/* generate l(ij) factors from the g(ij) factors as in bathe's book */

/* we only need to do this for the last block */


   if(h>=0){

        b=0.;

        dd=diags[i1]+1;

        for(i3=0;i3<=h;i3++){

                c=matrix[dd+i3]*work[i1-i3-1];

                b += c*matrix[dd+i3];

                matrix[dd+i3]=c;

        } /* for(i3=... */

        dd=diags[i1];

        matrix[dd]-=b;

   } /* if (h>=0... */

dd=diags[i1];

work[i1]=1.0/matrix[dd];

if(work[i1] <=0.0){

/* m_lock();  */
fprintf(stdout," error: negative pivot in equation %d\n",i1);
fflush(stdout);
/* m_unlock(); */

solverr = i1;
return;
}

} /* for(i1... */


dflag[n]=1;

}  /* endless for loop */

} /* end of routine  */
