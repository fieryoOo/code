#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
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
         return 0; //exit((int)d-'0');
      }
      for(i=0;i<m;i++) {
         for(j=0;j<m;j++)
            datout[i][j] = c[i][j]/d; //printf(" %f",c[i][j]/d);
      }
   } 
   return 1;
}

#define NDAT 500
int least_fit_parabola (double *datx, double *daty, double *weit, int ndat, double *a, double *b, double *c, double *STD, double *STDw)
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
      //A[0][i] = pow(datx[i],2)*weit[i];
      A[0][i] = datx[i]*weit[i];
      A[1][i] = sqrt(datx[i])*weit[i];
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
   *STD = 0; *STDw = 0;
   double v1=0, v2=0;
   for(i=0;i<ndat;i++) {
      *STD = *STD + pow(daty[i] - (*a * pow(datx[i],2) + *b * datx[i] + *c),2);
      *STDw = *STD + pow(daty[i] - (*a * pow(datx[i],2) + *b * datx[i] + *c),2) * weit[i];
      v1 += weit[i];
      v2 += pow(weit[i],2);
   }
   *STD = sqrt(*STD/(ndat-3.));
   *STDw = sqrt(*STDw*v1/(v1*v1-3.*v2));
   //cout<<*A0<<" + "<<*A1<<" * sin( theta + "<<*phi<<" )"<<endl;

   return 1;
}

int readdata (char *infile, double **data, int iage, int ndep, float depstep) {

   FILE *ff;
   char buff[300];
   int i, idep;
   double a, b, c, STD, STDw;
   double dep, vel;

   if((ff=fopen(infile,"r"))==NULL) {
      std::cout<<"Can't open file: "<<infile<<endl;
      return 0;
   }
   for(i=0;;i++) {
      if(fgets(buff,300,ff)==NULL) break;
      if( sscanf(buff,"%lf %lf",&dep,&vel)!=2 ) { cout<<"formmat error!"<<endl; return 0; }
      idep = int(dep/depstep);
      data[idep][iage] = vel;
   }
   fclose(ff);

   return 1;
}

#define NFILE 10
int main (int argc, char *argv[]) {
   //10 files max
   int i, nfile;
   char ** filelst = new char*[NFILE];

   //define model files
   char const dir[50] = "./SV_age_model";
   float age, magemin=0.5, magemax=3.5, magestep = 0.5;
   nfile = 0;
   for(age=magemin; age<=magemax; age+=magestep) {
      filelst[nfile] = new char[100];
      sprintf(filelst[nfile++], "%s/%.1f.mod.out_sm", dir, age);
   }

   //initialize data matrix
   float agestep = 0.1, agemax = 10.1, depstep = 0.3, depmax = 80;
   int ndep = int(depmax/depstep)+1, nage = int(agemax/agestep)+1;
   double **dataall = (double **) calloc (ndep , sizeof(double *));
   for(i=0;i<ndep;i++) dataall[i] = (double *) calloc (nage , sizeof(double));

   int iage, idep;
   for(age=magemin,i=0;i<nfile;i++) {
      iage = int(age/agestep);
      if( ! readdata(filelst[i], dataall, iage, ndep, depstep) ) {
	 cout<<"readdata failed!"<<endl;
	 exit(-1);
      }
      age+=magestep;
   }

   //compute/update one layer at a time
   int ivalid;
   double datx[nfile], daty[nfile], weit[nfile];
   double a, b, c, STD, STDw, vel;
   for(i=0;i<nfile;i++) weit[i] = 1.;
   for(idep=0;idep<ndep;idep++) {
      for(ivalid=0,age=magemin,i=0;i<nfile;i++) {
         iage = int(age/agestep);
	 vel = dataall[idep][iage];
	 if(vel<0.01) continue;
	 datx[ivalid] = age; daty[ivalid++] = vel;
	 age+=magestep;
      }
      if(ivalid<6) continue;
      least_fit_parabola (datx, daty, weit, nfile, &a, &b, &c, &STD, &STDw);
      //cerr<<idep*depstep<<":  vel = "<<a<<"*age + "<<b<<"*sqrt(age) + "<<c<<"   std: "<<STD<<" "<<STDw<<endl;
      for(iage=0;iage<nage;iage++) {
         age = iage*agestep;
         dataall[idep][iage] = a*age+b*sqrt(age)+c;
      }
   }
//for(idep=0;idep<ndep;idep++) {
//   cout.precision(4);
//   for(iage=0;iage<nage;iage++) cout<< setw(1) <<dataall[idep][iage]<<" ";
//   cout<<endl;
//}

   //output results for each age
   FILE *fout = NULL;
   char outname[50];
   char const outdir[40] = "./SV_allage_model";
   int started;
   for(age=0;age<agemax;age+=agestep) {
      iage = int(age/agestep);
      sprintf(outname, "%s/%.1f.mod.out", outdir, age);
      fout = fopen(outname, "w");
      started=0;
      for(idep=0;idep<ndep;idep++) {
	 vel = dataall[idep][iage];
	 if( vel>=0.01 ) started=1;
	 else if( started && vel<0.01 ) continue;
	 fprintf(fout, "%.1lf %lf\n", idep*depstep, vel);
      }
      fclose(fout); fout=NULL;
   }
   
   //cleanup
   for(i=0;i<ndep;i++) free(dataall[i]);
   free(dataall);
   for(i=0;i<nfile;i++) delete [] filelst[i];
   delete [] filelst;
   return 1;
}
