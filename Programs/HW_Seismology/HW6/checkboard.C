#include<stdio.h>
#include <stdlib.h>
#include<iostream>
#include<iomanip>
#include<math.h>
#include <string>
#include <unistd.h>
using namespace std;

#define NMAX 10
#define M_define(A, ni, nj) \
   A = (double **) malloc ( ni * sizeof(double *) ); \
   for(i=0;i<ni;i++) A[i] = (double *) malloc ( nj * sizeof(double) );


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

int Transpose(double **A, int ni, int nj, double **A_T)
{
   int i,j;
   for(i=0;i<ni;i++)
      for(j=0;j<nj;j++)
         A_T[j][i] = A[i][j];
   return 1;
}

int Multiply(double **A, int ni, int nj, int nk, double **B, double **out)
{
   int i, j, k;
   for(i=0;i<ni;i++)
      for(k=0;k<nk;k++) {
         out[i][k] = 0.;
         for(j=0;j<nj;j++) out[i][k] += A[i][j]*B[j][k];
      }
   return 1;
}

int Print_Matrix (double **A, int ni, int nj, int nw)
{
   int i, j;
   cout.precision(4);
   for(i=0;i<ni;i++) {
      for(j=0;j<nj;j++) {
         if(fabs(A[i][j])<1e-5) A[i][j] = 0;
         cout<<setw(nw)<<A[i][j]<<" ";
      }
      cout<<endl;
   }
   cout<<endl;
   return 1;
}

int Print_Matrix_33 (double **A, int nw)
{
   int i, j;
   for(i=0;i<9;i++) {
      if(fabs(A[i][0])<1e-5) A[i][0] = 0;
      cout<<setw(nw)<<A[i][0]<<" ";
      if((i+1)%3==0)cout<<endl;
   }
   cout<<endl;
   return 1;
}


int main (int argc, char *argv[])
{
   if(argc != 5){
      printf("Usage: checkboard.C [G] [Up] [Sp] [Vp_T]\n");
      exit(-1);
   }

   FILE *ff;
   char buff[300];
   double dtmp;
/*/-------Definition of matrices-----------
G: forward problem matrix
Up: left orthogonal matrix from SVD decomposition in compact form
Sp: diagonal matrix from SVD decomposition in compact form
Vp: right orthogonal matrix from SVD decomposition in compact form
Up_T: transpose of Up
Sp_T: transpose pf Sp
Vp_I: inverse of Vp
G_I: Pseudo inverse of G
G_I_T: transpose of G_I
M_t: true model
M_r: recovered model
T: travel time measurements 
----------------------------------------*/
   double **G, **Up, **Up_T, **Sp, **Sp_I, **Vp, **Vp_T, **G_I, **M_t, **M_r, **T, **G_I_T, **Mtmp;
//--------------------------//
   int i, j, k;
   M_define(G, 8, 9)
   M_define(Up, 8, 7)
   M_define(Up_T, 7, 8)
   M_define(Sp, 7, 7)
   M_define(Sp_I, 7, 7)
   M_define(Vp, 9, 7)
   M_define(Vp_T, 7, 9)
   M_define(G_I, 9, 8)
   M_define(M_t, 9, 1)
   M_define(M_r, 9, 1)
   M_define(T, 8, 1)
   M_define(G_I_T, 8, 9)
   M_define(Mtmp, 9, 7)

//--------------------read in G-----------------------//
   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file "<<argv[2]<<endl;
      return 0;
   }
   for(i=0;i<8;i++) {
      if(fgets(buff,300,ff)==NULL) return 0;
      if((sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf", &G[i][0], &G[i][1], &G[i][2], &G[i][3], &G[i][4], &G[i][5], &G[i][6], &G[i][7], &G[i][8]))!=9) return 0;
   }
   fclose(ff);
   cout<<"G matrix: "<<endl;
   Print_Matrix (G, 8, 9, 7);

//-------------------read in Up----------------------//
   if((ff=fopen(argv[2],"r"))==NULL) {
      cout<<"Can't open file "<<argv[2]<<endl;
      return 0;
   }
   for(i=0;i<8;i++) {
      if(fgets(buff,300,ff)==NULL) return 0;
      if((sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf", &Up[i][0], &Up[i][1], &Up[i][2], &Up[i][3], &Up[i][4], &Up[i][5], &Up[i][6]))!=7) return 0;
   }
   fclose(ff);
   //cout<<"Up matrix: "<<endl;
   //Print_Matrix (Up, 8, 7, 12);

//-------------------read in Sp----------------------//
   if((ff=fopen(argv[3],"r"))==NULL) {
      cout<<"Can't open file "<<argv[3]<<endl;
      return 0;
   }
   for(i=0;i<7;i++) {
      if(fgets(buff,300,ff)==NULL) return 0;
      if((sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf", &Sp[i][0], &Sp[i][1], &Sp[i][2], &Sp[i][3], &Sp[i][4], &Sp[i][5], &Sp[i][6]))!=7) return 0;
   }
   fclose(ff);
   //cout<<"Sp matrix: "<<endl;
   //Print_Matrix (Sp, 7, 7, 8);

//------------------read in Vp_T---------------------//
   if((ff=fopen(argv[4],"r"))==NULL) {
      cout<<"Can't open file "<<argv[4]<<endl;
      return 0;
   }
   for(i=0;i<7;i++) {
      if(fgets(buff,300,ff)==NULL) return 0;
      if((sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf", &Vp_T[i][0], &Vp_T[i][1], &Vp_T[i][2], &Vp_T[i][3], &Vp_T[i][4], &Vp_T[i][5], &Vp_T[i][6], &Vp_T[i][7], &Vp_T[i][8]))!=9) return 0;
   }
   fclose(ff);
   //cout<<"Vp_T matrix: "<<endl;
   //Print_Matrix (Vp_T, 7, 7, 12);

//------define true model-----//
   for(i=0;i<9;i++) M_t[i][0] = pow(-1.,i+1);
   cout<<"True model: "<<endl;
   Print_Matrix_33 (M_t, 5);

//-----compute travel time from G-----//
   Multiply(G, 8, 9, 1, M_t, T);  
   cout<<"Travel times computed from G and the true model: "<<endl;
   Print_Matrix (T, 8, 1, 8);

//-----compute generalized inverse of G-----//
   Transpose(Vp_T, 7, 9, Vp);
   Transpose(Up, 8, 7, Up_T);
   Inverse(Sp, 7, Sp_I);
   Multiply(Vp, 9, 7, 7, Sp_I, Mtmp);
   Multiply(Mtmp, 9, 7, 8, Up_T, G_I);
   cout<<"Generalized inverse of G:"<<endl;
   Print_Matrix (G_I, 9, 8, 8);

//--------------compute model from G_I---------------//
   Multiply(G_I, 9, 8, 1, T, M_r);
   cout<<"recovered model: "<<endl;
   Print_Matrix_33 (M_r, 8);

//---------compute difference between true and recovered model---------//
   double **M_diff;
   M_define(M_diff, 9, 1)
   for(i=0;i<9;i++) M_diff[i][0] = M_r[i][0]-M_t[i][0];
   cout<<"Difference between true and recovered model: "<<endl;
   Print_Matrix_33 (M_diff, 8);

//------compute covariance matrix-------//
   double **Cov;
   M_define(Cov, 9, 9)
   Transpose(G_I, 9, 8, G_I_T);
   Multiply(G_I, 9, 8, 9, G_I_T, Cov);
   cout<<"Covariance: "<<endl;
   Print_Matrix(Cov, 9, 9, 7);  

   return 1;
}
