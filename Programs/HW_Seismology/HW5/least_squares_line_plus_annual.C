/* Inverse of a n by n matrix */
#include<stdio.h>
#include <stdlib.h>
//#include<conio.h>
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
int least_fit_l_s1_s2 (double *tim, double *dat, double *weit, int ndat, int dim, double *a, double *b, double *c1, double *phi1, double *c2, double *phi2, double *vr, double *rc, double **cov)
{
   int i, j, ii;
   double pi=3.14159265359, sec_y=3600.*24.*365.25;

   if(dim!=2 && dim!=4 && dim!=6) {
      cout<<"dim has to be one of the 3 value: 2 4 or 6!"<<endl;
      exit(0);
   }
   double A[dim][ndat], AC[dim], coef[dim], **AA, **AAI;
   AA = (double **) malloc ( dim * sizeof(double *) );
   AAI = (double **) malloc ( dim * sizeof(double *) );
   for(i=0;i<dim;i++) AA[i] = (double *) malloc ( dim * sizeof(double) );
   for(i=0;i<dim;i++) AAI[i] = (double *) malloc ( dim * sizeof(double) );

   for(i=0;i<ndat;i++) {
      A[0][i] = weit[i];
      A[1][i] = tim[i]*weit[i];
      if(dim==2) continue;
      A[2][i] = sin(tim[i]*2*pi/sec_y)*weit[i];
      A[3][i] = cos(tim[i]*2*pi/sec_y)*weit[i];
      if(dim==4) continue;
      A[4][i] = sin(tim[i]*4*pi/sec_y)*weit[i];
      A[5][i] = cos(tim[i]*4*pi/sec_y)*weit[i];
   }
   for(i=0;i<dim;i++) 
      for(j=i;j<dim;j++) {
         AA[i][j] = 0.;
         for(ii=0;ii<ndat;ii++) AA[i][j] += A[i][ii]*A[j][ii];
      }
   for(i=1;i<dim;i++) for(j=0;j<i;j++) AA[i][j] = AA[j][i];
   Inverse( AA, dim, cov );
   for(i=0;i<dim;i++) {
      AC[i]=0;
      for(ii=0;ii<ndat;ii++) AC[i] += A[i][ii]*dat[ii]*weit[ii];
   }
   for(i=0;i<dim;i++) {
      coef[i] = 0;
      for(j=0;j<dim;j++) {
          coef[i] += cov[i][j]*AC[j];
      }
      //cout<<"coef["<<i<<"]: "<<coef[i]<<endl;
   }

   double datf[ndat], dd;
   for(i=0;i<ndat;i++) {
      datf[i] = 0;
      for(ii=0;ii<dim;ii++) datf[i]+=A[ii][i]*coef[ii];
      datf[i]/=weit[i];
   }
   dd=0; *vr=0; *rc=0;
   for(i=0;i<ndat;i++) {
      dd += pow(dat[i],2);
      *vr += pow(dat[i]-datf[i],2);
      *rc += pow((dat[i]-datf[i])*weit[i],2);
   }
   *vr = *vr/dd;
   *rc = *rc/(ndat-dim-1);

   *a=coef[0];
   *b=coef[1];
   if(dim>2) {
      *c1=sqrt(coef[2]*coef[2]+coef[3]*coef[3]); 
      *phi1=atan2(coef[3],coef[2]);
   }
   if(dim>4) {
      *c2=sqrt(coef[4]*coef[4]+coef[5]*coef[5]);
      *phi2=atan2(coef[5],coef[4]);
   }
   //*std = 0;
   //for(i=0;i<ndat;i++) *std += pow((dat[i] - (*A0 + *A1 * sin(tim[i]*pi/180.+*phi))),2);
   //*std = sqrt(*std/(ndat-1));
   return 1;
}

int L1_fit_l_s1_s2 (double *tim, double *dat, double *weit, int ndat, int dim, double *a, double *b, double *c1, double *phi1, double *c2, double *phi2, double *vr, double *rc, double **cov)
{
   int i, j, ii, iter;
   double pi=3.14159265359, sec_y=3600.*24.*365.25;

   if(dim!=2 && dim!=4 && dim!=6) {
      cout<<"dim has to be one of the 3 value: 2 4 or 6!"<<endl;
      exit(0);
   }
//-------------initialize data matrices-----------//
   double A[dim][ndat], AC[dim], coef[dim], **AA, **AAI;
   AA = (double **) malloc ( dim * sizeof(double *) );
   AAI = (double **) malloc ( dim * sizeof(double *) );
   for(i=0;i<dim;i++) AA[i] = (double *) malloc ( dim * sizeof(double) );
   for(i=0;i<dim;i++) AAI[i] = (double *) malloc ( dim * sizeof(double) );
//-------------define G matrix-------------------//
   for(i=0;i<ndat;i++) {
      A[0][i] = weit[i];
      A[1][i] = tim[i]*weit[i];
      if(dim==2) continue;
      A[2][i] = sin(tim[i]*2*pi/sec_y)*weit[i];
      A[3][i] = cos(tim[i]*2*pi/sec_y)*weit[i];
      if(dim==4) continue;
      A[4][i] = sin(tim[i]*4*pi/sec_y)*weit[i];
      A[5][i] = cos(tim[i]*4*pi/sec_y)*weit[i];
   }
//-------------IRLS loop starts------------------//
   double datf[ndat], R[ndat], coef_o[dim], dm2, m2, tolr=1e-26, tolm=1e-26;
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
      Inverse( AA, dim, cov );
      //----------G transpose * R * d---------------//
      for(i=0;i<dim;i++) {
         AC[i]=0;
         for(ii=0;ii<ndat;ii++) AC[i] += A[i][ii]*R[ii]*dat[ii]*weit[ii];
      }
      //----------(G'RG) inverse * G'Rd-----------//
      for(i=0;i<dim;i++) {
         coef_o[i] = coef[i];
         coef[i] = 0;
         for(j=0;j<dim;j++) {
             coef[i] += cov[i][j]*AC[j];
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
         R[i] = fabs(dat[i] - datf[i]);
         if(R[i]<tolr) R[i]=fabs(1./tolr);
         else R[i] = fabs(1./R[i]);
      }
   }
   cout<<"    "<<iter<<" iterations performed for tolr=1e-25 and tolm=1e-25"<<endl;

//------------------compute variance reduction and reduced chi-squre value--------------------------//
   double dd;
   for(i=0;i<ndat;i++) {
      datf[i] = 0;
      for(ii=0;ii<dim;ii++) datf[i]+=A[ii][i]*coef[ii];
      datf[i]/=weit[i];
   }
   dd=0; *vr=0; *rc=0;
   for(i=0;i<ndat;i++) {
      dd += pow(dat[i],2);
      *vr += pow(dat[i]-datf[i],2);
      *rc += pow((dat[i]-datf[i])*weit[i],2);
   }
   *vr = *vr/dd;
   *rc = *rc/(ndat-dim-1);

//---------------------------compute fitting coefficients---------------------------------//
   *a=coef[0];
   *b=coef[1];
   if(dim>2) {
      *c1=sqrt(coef[2]*coef[2]+coef[3]*coef[3]);
      *phi1=atan2(coef[3],coef[2]);
   }
   if(dim>4) {
      *c2=sqrt(coef[4]*coef[4]+coef[5]*coef[5]);
      *phi2=atan2(coef[5],coef[4]);
   }
   //*std = 0;
   //for(i=0;i<ndat;i++) *std += pow((dat[i] - (*A0 + *A1 * sin(tim[i]*pi/180.+*phi))),2);
   //*std = sqrt(*std/(ndat-1));
   return 1;
}


int main (int argc, char *argv[])
{
   if(argc != 2){
      printf("Usage: least_squares_sine.C [input file]\n");
      exit(-1);
   }

   FILE *ff;
   char buff[300], ctmp[300];
   int i, j, ndat;
   double tim[NDAT], dat[NDAT], weit[NDAT], ftmp, a, b, c1, c2, phi1, phi2, vr, rc;

   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file: "<<argv[1]<<endl;
      exit(0);
   }
   for(i=0;;i++) {
      if(fgets(buff,300,ff)==NULL) break;
      sscanf(buff,"%s %lf %lf %lf %lf %lf %lf %lf",&ctmp[0], &tim[i], &ftmp, &ftmp, &ftmp, &ftmp, &dat[i], &weit[i]);
      weit[i] = 1./weit[i];
   }
   fclose(ff);
   ndat=i;

   cout.precision(15);
   int dim=6;
   double *cov[dim];
   for(i=0;i<dim;i++) cov[i] = (double *) malloc ( dim * sizeof(double));

   cout<<"GEOL 6670   Homework set 5   Ye Tian"<<endl<<endl;
   cout<<"(1) By fitting a line a + b*x to the data we got:"<<endl;
   least_fit_l_s1_s2 (&tim[0], &dat[0], &weit[0], ndat, 2, &a, &b, &c1, &phi1, &c2, &phi2, &vr, &rc, &cov[0]);
   cout<<"    "<<a<<" + "<<b<<"*x"<<endl;
   cout<<"    variance reduction: 1-(d-Gm)T*(d-Gm) = 1-"<<vr<<"  reduced chi-square: sum((d_i-Gm_i)/uncertainty_i) = "<<rc<<endl;
   cout<<"    covariance matrix: "<<endl;
   for(i=0;i<2;i++) {
      cout<<"    ";
      for(j=0;j<2;j++) printf("%12g ", cov[i][j]);
      cout<<endl;
   }
   cout<<"    By using the square root of model-parameter variances as their uncertainties we got: "<<endl;
   cout<<"    uncertainty on a: "<<sqrt(cov[0][0])<<" ("<<sqrt(cov[0][0])/fabs(a)*100.<<"%)  uncertainty on b: "<<sqrt(cov[1][1])<<" ("<<sqrt(cov[1][1])/fabs(b)*100.<<"%)"<<endl<<endl;
   cout<<"(2) Adding an annual term c1*sin(theta+phi1) to the fit and compute the variance reduction and chi-square again: "<<endl;
   least_fit_l_s1_s2 (&tim[0], &dat[0], &weit[0], ndat, 4, &a, &b, &c1, &phi1, &c2, &phi2, &vr, &rc, &cov[0]);
   cout<<"    "<<a<<" + "<<b<<"*x + "<<c1<<"*sin(2*pi*x/31557600.0+"<<phi1<<")"<<endl;
   cout<<"    variance reduction: 1-"<<vr<<"  reduced chi-square: "<<rc<<endl<<endl;
   cout<<"(3) Again adding an biannual term c2*sin(2*theta+phi2) to the fit: "<<endl;
   least_fit_l_s1_s2 (&tim[0], &dat[0], &weit[0], ndat, 6, &a, &b, &c1, &phi1, &c2, &phi2, &vr, &rc, &cov[0]);
   cout<<"    "<<a<<" + "<<b<<"*x + "<<c1<<"*sin(2*pi*x/31557600.0+"<<phi1<<") + "<<c2<<"*sin(4*pi*x/31557600.0+"<<phi2<<")"<<endl;
   cout<<"    variance reduction: 1-"<<vr<<"  reduced chi-square: "<<rc<<endl<<endl;
   cout<<"(4) By comparing the results from all 3 of the above results, 2 key conclusions can be made: "<<endl;
   cout<<"    First, in those 3 results, the value of model parameter a are nearly identical while the value of b changes a lot. This is consistent with the model-parameter uncertainties we got in problem 1, where the fractional uncertainty of a (5.6e-9%) is much smaller than that of b (1.9%)."<<endl;
   cout<<"    Second, as we add in the annual and biannual term, the variance reduction gets more close to 1 and the reduced chi-square decreases. This indicates that we fit the data better by adding in those two terms."<<endl<<endl;
   cout<<"(5) Not sure how to do this one, sorry."<<endl<<endl;
   cout<<"(6) Use an L1 norm instead of L2 to fit the data: "<<endl;
   L1_fit_l_s1_s2 (&tim[0], &dat[0], &weit[0], ndat, 6, &a, &b, &c1, &phi1, &c2, &phi2, &vr, &rc, &cov[0]);
   cout<<"    "<<a<<" + "<<b<<"*x + "<<c1<<"*sin(2*pi*x/31557600.0+"<<phi1<<") + "<<c2<<"*sin(4*pi*x/31557600.0+"<<phi2<<")"<<endl;
   cout<<"    variance reduction: 1-"<<vr<<"  reduced chi-square: "<<rc<<endl;
   cout<<"    This result is very close to the result in part 3. However, because there are several outliers in the data set, the chi-square value from L1 norm turned out to be a little bit larger than from L2 norm. This is because the L1 norm down weights those outliers thus leads to larger chi-square distances."<<endl;
   //cout<<A0<<" + "<<A1<<" * sin( theta*pi/180. + "<<phi<<" )  std: "<<std<<endl;
   return 1;
}
