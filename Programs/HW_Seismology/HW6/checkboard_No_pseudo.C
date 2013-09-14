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
   for(i=0;i<ni;i++) {
      for(j=0;j<nj;j++) cout<<setw(nw)<<A[i][j]<<" ";
      cout<<endl;
   }
   return 1;
}


int main (int argc, char *argv[])
{
   if(argc != 2){
      printf("Usage: checkboard.C [input G file]\n");
      exit(-1);
   }

   FILE *ff;
   char buff[300];
   double dtmp;
   double **G, **G_T, **S_t, **S_p, **GG, **GGI, **GT, **T, **Itmp;
   int i, j, k;
   M_define(G, 8, 9)
   M_define(G_T, 9, 8)
   M_define(S_t, 9, 1)
   M_define(S_p, 9, 1)
   M_define(GG, 9, 9)
   M_define(GGI, 9, 9)
   M_define(GT, 9, 1)
   M_define(T, 8, 1)
   M_define(Itmp, 9, 9)

   if((ff=fopen(argv[1],"r"))==NULL) {
      cout<<"Can't open file "<<argv[2]<<endl;
      return 0;
   }
   for(i=0;i<8;i++) {
      if(fgets(buff,300,ff)==NULL) return 0;
      if((sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf", &G[i][0], &G[i][1], &G[i][2], &G[i][3], &G[i][4], &G[i][5], &G[i][6], &G[i][7], &G[i][8], &G[i][9]))!=9) return 0;
   }
   fclose(ff);
   cout<<"G matrix: "<<endl;
   Print_Matrix (G, 8, 9, 8);

   for(i=0;i<9;i++) S_t[i][0] = pow(-1.,i+1);
   cout<<"True model: "<<endl;
   Print_Matrix (S_t, 9, 1, 5);

   Multiply(G, 8, 9, 1, S_t, T);  
   cout<<"Travel times: "<<endl;
   Print_Matrix (T, 8, 1, 5);

   Transpose(G, 8, 9, G_T);
   Multiply(G_T, 9, 8, 9, G, GG);
   cout<<"G'G: "<<endl;
   Print_Matrix (GG, 9, 9, 5);

   Inverse( GG, 9, GGI);
   if((ff=fopen("GG_inv.txt","r"))==NULL) {
      cout<<"Can't open file GG_inv.txt"<<endl;
      return 0;
   }
   for(i=0;i<9;i++) {
      if(fgets(buff,300,ff)==NULL) return 0;
      if((sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf", &GGI[i][0], &GGI[i][1], &GGI[i][2], &GGI[i][3], &GGI[i][4], &GGI[i][5], &GGI[i][6], &GGI[i][7], &GGI[i][8], &GGI[i][9]))!=9) return 0;
   }
   fclose(ff);

   cout<<"Inverse(G'G): "<<endl;
   Print_Matrix (GGI, 9, 9, 12);

   Multiply(GG, 9, 9, 9, GGI, Itmp);
   cout<<"Identity?: "<<endl;
   Print_Matrix (Itmp, 9, 9, 12);

   Multiply(G_T, 9, 8, 1, T, GT);
   cout<<"G'T: "<<endl;
   Print_Matrix (GT, 9, 1, 5);

   Multiply(GGI, 9, 9, 1, GT, S_p);
   cout<<"recovered model: "<<endl;
   Print_Matrix (S_p, 9, 1, 5);

   return 1;
}
