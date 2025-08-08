/* this routine calculates the stiffness matrix for a 2-noded
 straight truss element */

#include <stdio.h>
#include <math.h>

//int trussel(int,int,double [][3],double *,FILE *);
int trussel(int nel,int matprp,double xx[][3],double s[],FILE* output);

int trussel(int nel,int matprp,double xx[][3],double s[],FILE* output)
//int nel,matprp;
//double xx[][3],s[];
//FILE *output;
{

double l,k,ym,a,cs[6];
int c,c2,index;

extern double *matptr;

c=(matprp-1)*3;
ym=matptr[c];
a=matptr[c+2];

l= sqrt((xx[1][0] - xx[0][0])*(xx[1][0] - xx[0][0])+
        (xx[1][1] - xx[0][1])*(xx[1][1] - xx[0][1])+
        (xx[1][2] - xx[0][2])*(xx[1][2] - xx[0][2]) );

k=ym*a/l;

cs[0]=(xx[1][0] - xx[0][0])/l;
cs[1]=(xx[1][1] - xx[0][1])/l;
cs[2]=(xx[1][2] - xx[0][2])/l;
cs[3]= -cs[0];
cs[4]= -cs[1];
cs[5]= -cs[2];

for(c=0;c<6;c++){

for(c2=0;c2<c;c2++){

index=c*(c+1)/2 + c2;

s[index]= k*cs[c]*cs[c2];}}
return(0);
}
