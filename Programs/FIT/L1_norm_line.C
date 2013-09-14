#include<stdio.h>
#include <stdlib.h>
#include<iostream>
#include<math.h>
#include <string>
#include <unistd.h>
using namespace std;

#define NMAX 10

void arg(double *a, double *b, int *n,int x,int y) {
   int k,l,i,j;
   for(i=0,k=0;i<*n;i++,k++) {
      for(j=0,l=0;j<*n;j++,l++) {
         if(i==x)
            i++;
         if(j==y)
            j++;
         *(b+NMAX*k+l)=*(a+NMAX*i+j);
      }
   }
   *n=*n-1;
}

double det(double *p,int *n) {
   int i,j,m;
   double d[NMAX][NMAX], sum=0;
   m=*n;
   if(*n==2)
      return(*p**(p+NMAX+1)-*(p+1)**(p+NMAX));
   for(i=0,j=0;j<m;j++) {
      *n=m;
      arg(p,&d[0][0],n,i,j);
      sum=sum+*(p+NMAX*i+j)*pow(-1,(i+j))*det(&d[0][0],n);
   }

   return(sum);
}

int Inverse( double **datin, int n, double **datout )
{
   //void arg(int *,int *, int *,int ,int );
   //int det(int *,int *);
   int i,j,m;
   double a[NMAX][NMAX],b[NMAX][NMAX],c[NMAX][NMAX],d;
   //clrscr();
   for(i=0;i<n;i++) for(j=0;j<n;j++) a[i][j] = datin[i][j];
   if(n==2) {
      c[0][0]=a[1][1];
      c[1][1]=a[0][0];
      c[0][1]=-a[0][1];
      c[1][0]=-a[1][0];
      d=a[0][0]*a[1][1]-a[0][1]*a[1][0];
      //printf("Determinant: %lf\n",d);
      if(d==0) {
         //getch();
         //std::cin.get();
         return 0; //exit((int)d-'0');
      }

      for(i=0;i<n;i++) {
         for(j=0;j<n;j++)
            datout[i][j] = c[i][j]/d;
      }
   }
   else {
      m=n;
      for(i=0;i<m;i++) {
         for(j=0;j<m;j++) {
            n=m;
            arg(&a[0][0],&b[0][0],&n,i,j);
            c[j][i]=pow(-1,(i+j))*det(&b[0][0],&n);
         }
      }
      n=m;
      d=det(&a[0][0],&n);
      //printf("Determinant is :%d\n",d);
      if(d==0) {
         printf("INVERSE DOES NOT EXIST");
         //getch();
         //std::cin.get();
         return 0; //exit((int)d-'0');
      }
      for(i=0;i<m;i++) {
         for(j=0;j<m;j++)
            datout[i][j] = c[i][j]/d; //printf(" %f",c[i][j]/d);
      }
   } //std::cin.get();
   return 1;
}

#define NDAT 5000

int L1_fit_l_s1_s2 (double *datx, double *daty, double *weit, int ndat, double *a, double *b, double *std)
{
   int i, j, ii, iter;
   int dim=2;
   
//-------------initialize data matrices-----------//
   double A[dim][ndat], AC[dim], coef[dim], **AA, **AAI;
   AA = (double **) malloc ( dim * sizeof(double *) );
   AAI = (double **) malloc ( dim * sizeof(double *) );
   for(i=0;i<dim;i++) AA[i] = (double *) malloc ( dim * sizeof(double) );
   for(i=0;i<dim;i++) AAI[i] = (double *) malloc ( dim * sizeof(double) );
//-------------define G matrix-------------------//
   for(i=0;i<ndat;i++) {
      A[0][i] = datx[i]*weit[i];
      A[1][i] = weit[i];
   }
//-------------IRLS loop starts------------------//
   double datf[ndat], R[ndat], coef_o[dim], dm2, m2, tolr=1e-20, tolm=1e-20;
   for(i=0;i<ndat;i++) R[i] = 1.;
   for(i=0;i<dim;i++) coef[i] = 0;
   for(iter=0;;iter++) {
      //----------G transpose * R * G--------------//
      for(i=0;i<dim;i++)
         for(j=i;j<dim;j++) {
            AA[i][j] = 0.;
            for(ii=0;ii<ndat;ii++) AA[i][j] += A[i][ii]*R[ii]*A[j][ii];
         }
      for(i=1;i<dim;i++) for(j=0;j<i;j++) AA[i][j] = AA[j][i];
      //----------G'RG inverse-------------------//
      Inverse( AA, dim, AAI );
      //----------G transpose * R * d---------------//
      for(i=0;i<dim;i++) {
         AC[i]=0;
         for(ii=0;ii<ndat;ii++) AC[i] += A[i][ii]*R[ii]*daty[ii]*weit[ii];
      }
      //----------(G'RG) inverse * G'Rd-----------//
      for(i=0;i<dim;i++) {
         coef_o[i] = coef[i];
         coef[i] = 0;
         for(j=0;j<dim;j++) {
             coef[i] += AAI[i][j]*AC[j];
         }
      }
      dm2=0; m2=0;
      //----------residual vector converge?-------//
      for(i=0;i<dim;i++) {
         dm2 += pow(coef[i]-coef_o[i],2);
         m2 += pow(coef_o[i],2);
      }
      if(dm2/(1.+m2) < tolm) break;
      //----------compute new R------------------//
      for(i=0;i<ndat;i++) {
         datf[i] = 0;
         for(ii=0;ii<dim;ii++) datf[i]+=A[ii][i]*coef[ii];
         datf[i]/=weit[i];
         R[i] = fabs(daty[i] - datf[i]);
         if(R[i]<tolr) R[i]=fabs(1./tolr);
         else R[i] = fabs(1./R[i]);
      }
   }
   cout<<iter<<" iterations performed for tolr="<<tolr<<" and tolm="<<tolm<<endl;

/*/------------------compute variance reduction and reduced chi-squre value--------------------------//
   double dd;
   for(i=0;i<ndat;i++) {
      datf[i] = 0;
      for(ii=0;ii<dim;ii++) datf[i]+=A[ii][i]*coef[ii];
      datf[i]/=weit[i];
   }
   dd=0; *vr=0; *rc=0;
   for(i=0;i<ndat;i++) {
      dd += pow(daty[i],2);
      *vr += pow(daty[i]-datf[i],2);
      *rc += pow((daty[i]-datf[i])*weit[i],2);
   }
   *vr = *vr/dd;

*/
//---------------------------compute fitting coefficients---------------------------------//
   *a=coef[0];
   *b=coef[1];
   *std = 0;
   for(i=0;i<ndat;i++) *std += pow((daty[i] - (*a * datx[i] + *b)),2);
   *std = sqrt(*std/(ndat-1));
   return 1;
}


int main (int argc, char *argv[])
{
   if(argc != 2){
      printf("Usage: L1_norm_line.C [input file]\n");
      exit(-1);
   }

   FILE *ff;
   char buff[300], ctmp[300];
   int i, j, ndat;
   double datx[NDAT], daty[NDAT], weit[NDAT], ftmp, a, b, std;

   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file: "<<argv[1]<<endl;
      exit(0);
   }
   for(i=0;;i++) {
      if(fgets(buff,300,ff)==NULL) break;
      sscanf(buff,"%lf %lf %lf", &datx[i], &daty[i], &weit[i]);
      weit[i] = 1./weit[i];
   }
   fclose(ff);
   ndat=i;

   //cout.precision(15);

   L1_fit_l_s1_s2 (&datx[0], &daty[0], &weit[0], ndat, &a, &b, &std);
   cout<<a<<"*x + "<<b<<"  std: "<<std<<endl;
   return 1;
}
