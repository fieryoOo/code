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
   double A[dim][ndat], AC[dim+1], coef[dim], **AA, **AAI;
   AA = (double **) malloc ( (dim+1) * sizeof(double *) );
   AAI = (double **) malloc ( (dim+1) * sizeof(double *) );
   for(i=0;i<dim+1;i++) AA[i] = (double *) malloc ( (dim+1) * sizeof(double) );
   for(i=0;i<dim+1;i++) AAI[i] = (double *) malloc ( (dim+1) * sizeof(double) );

   for(i=0;i<ndat;i++) {
      A[2][i] = pow(datx[i],2)*weit[i];
      A[1][i] = datx[i]*weit[i];
      A[0][i] = weit[i];
   }
   for(i=0;i<dim;i++) 
      for(j=i;j<dim;j++) {
         AA[i][j] = 0.;
         for(ii=0;ii<ndat;ii++) AA[i][j] += A[i][ii]*A[j][ii];
      }
   for(i=1;i<dim;i++) for(j=0;j<i;j++) AA[i][j] = AA[j][i];

   int nc=2;
   for(i=0;i<dim;i++) {
      AA[i][dim] = A[i][nc];
      AA[dim][i] = A[i][nc];
   }
   AA[dim][dim] = 0;
/*
   for(i=0;i<dim+1;i++) {
      for(j=0;j<dim+1;j++) cout<<AA[i][j]<<" ";
      cout<<endl;
   }
*/
   Inverse( AA, dim+1, AAI );
/*
   for(i=0;i<dim+1;i++) {
      for(j=0;j<dim+1;j++) cout<<AAI[i][j]<<" ";
      cout<<endl;
   }
*/
   for(i=0;i<dim;i++) {
      AC[i]=0;
      for(ii=0;ii<ndat;ii++) AC[i] += A[i][ii]*daty[ii]*weit[ii];
   }
   AC[dim] = daty[nc];

   for(i=0;i<dim;i++) {
      coef[i] = 0;
      for(j=0;j<dim+1;j++) {
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
      *std = *std + pow(daty[i] - (*c * pow(datx[i],2) + *b * datx[i] + *a),2);
      *stdw = *std + pow(daty[i] - (*c * pow(datx[i],2) + *b * datx[i] + *a),2) * weit[i];
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
      printf("Usage: least_squares_parabola.C [input file (x y weight)]\n");
      exit(-1);
   }

   FILE *ff;
   char buff[300];
   int i, ndat;
   double datx[NDAT], daty[NDAT], weit[NDAT], a, b, c, std, stdw;

   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file: "<<argv[1]<<endl;
      exit(0);
   }
   for(i=0;;i++) {
      if(fgets(buff,300,ff)==NULL) break;
      sscanf(buff,"%lf %lf %lf",&datx[i],&daty[i],&weit[i]);
      weit[i] = 1./weit[i];
   }
   fclose(ff);
   ndat=i;

   double weit0[ndat];
   for(i=0;i<ndat;i++) weit0[i] = 1;
   least_fit_parabola (&datx[0], &daty[0], &weit0[0], ndat, &a, &b, &c, &std, &stdw);
   cout<<"(2) By constrain the inversion to go through the 3rd point, we get:"<<endl;
   cout<<"   "<<c<<" * x^2 + "<<b<<" * x + "<<a<<"  (See page 2 for plots)";
   cout<<"    with a data variance (with denominator n-3) of: "<<std*std<<endl;
   cout<<"    In HW3 we see that the weighted fit went higher than the unweighted fit in the middle part. That was because the 3rd measurements was weighted heavier thus dragged the mid part up. Here the solution is even higher as the parabola is constrained to go through the 3rd point, which functions like an infinite weighting."<<endl<<endl;

   cout<<"(3) Linear programming is a technique for optimizing a linear objective function under a set of linear equality and inequality constraints. This technique can be applied to improve the accuracy and to reduce the computational cost in computing earthquake mechanisms when both the first-motion polarities and amplitude rations were taken into consideration."<<endl;
   cout<<"    A relevant paper: <Earthquake Mechanisms from Linear-Programming Inversion of Seismic-Wave Amplitude Ratios> by Bruce R. Julian and G. R. Foulger. Published on August 1996. This is the paper where the linear-programming method was extended to work on amplitude ratios in earthquake mechanism computation method."<<endl;

   return 1;
}
