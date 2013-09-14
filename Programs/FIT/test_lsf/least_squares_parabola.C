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
         printf("\n");
         for(j=0;j<n;j++)
            printf(" %f",c[i][j]/(float)d);
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

#define NDAT 500
int least_fit_parabola (double *datx, double *daty, double *weit, int ndat, double *a, double *b, double *c, double *std, double *stdw)
{
   int i, j, ii;
   double pi=3.14159265359;


   int dim=3;
   double A[dim][ndat], AC[dim], coef[dim], **AA, **AAI;
   AA = (double **) malloc ( dim * sizeof(double *) );
   AAI = (double **) malloc ( dim * sizeof(double *) );
   for(i=0;i<dim;i++) AA[i] = (double *) malloc ( dim * sizeof(double) );
   for(i=0;i<dim;i++) AAI[i] = (double *) malloc ( dim * sizeof(double) );

   for(i=0;i<ndat;i++) {
      A[0][i] = pow(datx[i],2)*weit[i];
      A[1][i] = datx[i]*weit[i];
      A[2][i] = weit[i];
   }
   for(i=0;i<dim;i++) 
      for(j=i;j<dim;j++) {
         AA[i][j] = 0.;
         for(ii=0;ii<ndat;ii++) AA[i][j] += A[i][ii]*A[j][ii];
      }
   for(i=1;i<dim;i++) for(j=0;j<i;j++) AA[i][j] = AA[j][i];
   Inverse( AA, dim, AAI );
   for(i=0;i<dim;i++) {
      AC[i]=0;
      for(ii=0;ii<ndat;ii++) AC[i] += A[i][ii]*daty[ii]*weit[ii];
   }

   for(i=0;i<dim;i++) {
      coef[i] = 0;
      for(j=0;j<dim;j++) {
          coef[i] += AAI[i][j]*AC[j];
      }
      //cout<<"coef["<<i<<"]: "<<coef[i]<<endl;
   }

   *a=coef[0]; 
   *b=coef[1]; 
   *c=coef[2];
   *std = 0; *stdw = 0;
   double v1=0, v2=0;
   for(i=0;i<ndat;i++) {
      *std = *std + pow(daty[i] - (*a * pow(datx[i],2) + *b * datx[i] + *c),2);
      *stdw = *std + pow(daty[i] - (*a * pow(datx[i],2) + *b * datx[i] + *c),2) * weit[i];
      v1 += weit[i];
      v2 += pow(weit[i],2);
   }
   *std = sqrt(*std/(ndat-3.));
   *stdw = sqrt(*stdw*v1/(v1*v1-3.*v2));
   //cout<<*A0<<" + "<<*A1<<" * sin( theta + "<<*phi<<" )"<<endl;

   return 1;
}

int main (int argc, char *argv[])
{
   if(argc != 2){
      printf("Usage: least_squares_parabola.C [input file (x y optional_weight)]\n");
      exit(-1);
   }

   FILE *ff;
   char buff[300];
   int i, ndat, weitflag=1;
   double datx[NDAT], daty[NDAT], weit[NDAT], a, b, c, std, stdw;

   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file: "<<argv[1]<<endl;
      exit(0);
   }
   for(i=0;;i++) {
      if(fgets(buff,300,ff)==NULL) break;
      if( sscanf(buff,"%lf %lf %lf",&datx[i],&daty[i],&weit[i])==2 ) weitflag=0;
      else weit[i] = 1./weit[i];
   }
   fclose(ff);
   ndat=i;
   if( !weitflag ) for(i=0;i<ndat;i++) weit[i] = 1.;

   least_fit_parabola (&datx[0], &daty[0], &weit[0], ndat, &a, &b, &c, &std, &stdw);
   cout<<"1) "<<a<<" * x^2 + "<<b<<" * x + "<<c<<endl<<endl;
   cout<<"2) data variance (with denominator n-3): "<<std*std<<endl<<endl<<endl;

   return 1;
}
