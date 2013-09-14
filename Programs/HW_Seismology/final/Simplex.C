#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
using namespace std;

/*/-------------------------variable definition-------------------------//
nv: number of variables in objective function
nc: number of constraints
Ni: number of rows of the tableau
Nj: number of columns of the tableau
F_opt: variable indicating whether the optimal solution has been found
P1: row number of the entering basic variable
P2: column number of the entering basic variable
XERR: variable indicating the error type of the problem
PH: phase number
Tb1: tableau for the 1st phase
Tb: tableau for the 2nd phase
//----------------------------------------------------------------------//
*/
int Ni, Nj, nc, nv, F_opt, P1, P2, XERR, PH;
double **Tb;
//-------subroutine to output simplex tableau--------//
void PrintData(double **A, int ni, int nj) {
   int i,j;
   for(i=0;i<ni;i++) {
      for(j=0;j<nj;j++) cout<<setw(10)<<A[i][j]<<" ";
      cout<<endl;
   }
   cout<<endl;
}
//----------read in objective function & constraints and set up simplex tableau-------//
int Tableau(char *fname) {
   FILE *ff;
   double R2;
   char R, buff[300], *pbuff, ctmp[10];
   int i,j,R1,flag;
//----------open input data file------------//
   if((ff=fopen(fname,"r"))==NULL) {
      cout<<"Cannot open dat file: "<<fname<<endl;
      exit(0);
   }
//-----------read in max/min choice, # of variables, and # of constraints------------//
   fgets(buff, 300, ff);
   sscanf(buff, "%d %d %d", &R1, &nv, &nc);
   if(R1!=-1 && R1!=1) {
      cout<<"Wrong value for min/max choice. Expecting 1 for max or -1 for min."<<endl;
      exit(0);
   }
//----------assign memory to the Tableau-----------//
   Ni=1+nc; Nj=2+nv+nc;
   Tb = (double **) malloc ( Ni * sizeof(double *) );
   for(i=0;i<Ni;i++) Tb[i] = (double *) malloc ( Nj * sizeof(double) );
//----------read in coefs and set values for tableau--------//
   Tb[0][0]=1; Tb[0][Nj-1]=0;
   for(i=1;i<Ni;i++) Tb[i][0] = 0.;
  //---read in objective function--//
   fgets(buff, 300, ff);
   pbuff = buff;
   for(j=1;j<=nv;j++) {
      if((sscanf(pbuff,"%s", &ctmp)) != 1) {
         cout<<"Less then "<<nv<<" variable-coefs found in data file. Stopped."<<endl;
         exit(0);
      }
      Tb[0][j] = atof(ctmp)*R1;
      pbuff += (strlen(ctmp)+1);
   }
  //---read in constraints---//
   for(i=1;i<=nc;i++) {
      if((fgets(buff, 300, ff))==NULL) {
         cout<<"Less then "<<nc<<" constraints found in data file. Stopped."<<endl;
         exit(0);
      }
      pbuff = buff;
      flag=0;
      for(j=1;j<=nv;j++) {
         if((sscanf(pbuff,"%s",&ctmp)) != 1) {
            flag=1;
            break;
         }
         Tb[i][j] = atof(ctmp);
         pbuff += (strlen(ctmp)+1);
      }
      if((sscanf(pbuff,"%s",&ctmp)) != 1) flag=1;
      Tb[i][Nj-1] = atof(ctmp);
      if(flag) {
         cout<<"No enough coefs found for the "<<j<<"th constraint. Stopped."<<endl;
         exit(0);
      }
   }
   fclose(ff);
  //---define the slack and error parts---//
   for(i=0;i<Ni;i++)
      for(j=nv+1;j<Nj-1;j++) {
         if(i+nv==j) Tb[i][j]=1;
         //else if(i+nv+nc==j && i>0) Tb[i][j]=-1;
         else Tb[i][j]=0;
      }
//-------------print initial tableau------------------//
   cout<<"Initial tableau:"<<endl;
   PrintData(Tb,Ni,Nj);
   return 1;
}

void Select(double **A, int ni, int nj);
void Transform(double **A, int ni, int nj);
//-------perform pivot operations to a initialized (starting vertex found) tableau------//
int Pivot(double **A, int ni, int nj) {
   for(;;) {
   //-------Select entering and leaving basic variable-----//
      Select(A, ni, nj);
      if(XERR==-2) {
         cout<<"unbounded. Stopped."<<endl;
         exit(0);
      }
   //-------transform to move from one vertex to an adjacent one------//
      Transform(A, ni, nj);
      //PrintData(A ,ni, nj);
   //-----check if the optimal solution has been found-----//
      if (F_opt) break;
   }
   return 1;
}
//-----subroutine of selecting entering and leaving basic variable-----//
void Select(double **A, int ni, int nj) {
   double ratio, vcmax;
   int i,j;
   //---------searching for the maximum coefficient in the objective function to decide the entering variable------//
   vcmax = 0.;
   for(j=1; j<nj-1; j++) {
      if (A[0][j] > vcmax) {
         vcmax = A[0][j];
         P2 = j;
      }
   }
   //---------searching for the minimum ratio to decide the leaving variable-----//
   ratio = 1e10; P1=0;
   for (i=3-PH; i<ni; i++) {
      if (A[i][P2] <= 0.0) continue;
      vcmax = A[i][nj-1] / A[i][P2];
      if (vcmax < ratio) {
         ratio = vcmax;
         P1 = i;
      }
    }
    if(P1==0) XERR = -2;
    return;
}

void Transform(double **A, int ni, int nj) {
   int i,j;
   double dtmp;
   //--------normalize entering row------//
   dtmp=A[P1][P2];
   for (j=1; j<nj; j++) {
      A[P1][j] /= dtmp;
   }
   //--------make substractions to other rows-------------//
   for (i=0; i<ni; i++) {
      if (i == P1) continue;
      dtmp=A[i][P2];
      for (j=1; j<nj; j++) {
         A[i][j] -= (A[P1][j] * dtmp);
         if(fabs(A[i][j])<1e-10) A[i][j]=0.;
      }
   }
   //for (i=1; I<Ni; i++)
     // if (Tb[i][Nj-1] < 0.0)  XERR = 1;
   F_opt = 1;
   //if (XERR == 1)  return;
   for (j=1; j<nj-1; j++)
      if (A[0][j] > 0.0)  F_opt = 0;
   return;
}
//-----subroutine to extract the current value of variables from the tableau-----//
void Pick_Var(double **A, int ni, int nj) {
   int i, j, ii, flag;
   double Var[nv];
   for(i=1;i<ni;i++)
      for(j=1;j<=nv;j++) {
         if(A[i][j]!=1) continue;
         flag=0;
         for(ii=1;ii<ni;ii++) {
            if(ii==i) continue;
            if(A[ii][j]!=0) {flag=1; break;}
         }
         if(flag==0) { 
            Var[j-1]=A[i][nj-1]; 
            break; 
         }
      }
   cout<<"Variables: ";
   for(j=0;j<nv;j++) printf("%.3f ",Var[j]);
   cout<<endl;
   cout<<"Moment tensor: ";
   for(j=0;j<6;j++) printf("%.3g ",(Var[j]-Var[j+6])*1.46e+24/(Var[0]-Var[6]));
   return;
}
//---------Phase one: searching for the starting vertex and checking feasibility--------//
int Phase1() {
   int i, j, jj, flag=0;
   PH=1;
//----check----//
   for(i=1;i<Ni;i++) if(Tb[i][Nj-1]<0) flag++;
   if(flag==0) return 1;
//------construct new tableau for phase 1-------//
   int ni=Ni+1, nj=Nj+flag+1;
   double **Tb1;
   Tb1 = (double **) malloc ( ni * sizeof(double *) );
   for(i=0;i<ni;i++) Tb1[i] = (double *) malloc ( nj * sizeof(double) );
   for(i=0;i<ni;i++) for(j=0;j<nj;j++) Tb1[i][j] = 0;
   Tb1[0][0]=1;
   for(i=0;i<Ni;i++) {
      for(j=0;j<Nj-1;j++) Tb1[i+1][j+1] = Tb[i][j];
      Tb1[i+1][nj-1] = Tb[i][Nj-1];
   }
//------define new objective function and error colums------//
   flag=0;
   for(i=2;i<ni;i++) if(Tb1[i][nj-1]<0) {
      Tb1[0][Nj+flag] = -1;
      for(j=2;j<nj;j++) Tb1[i][j] = -Tb1[i][j];
      Tb1[i][Nj+flag] = 1;
      flag++;
   }
//------pricing out--------//
   cout<<"Phase1 initial:"<<endl;
   PrintData(Tb1,ni,nj);
   for(i=2;i<ni;i++)for(j=Nj;j<Nj+flag;j++)
      if(Tb1[i][j]==1)for(jj=1;jj<nj;jj++) Tb1[0][jj] += Tb1[i][jj];
   cout<<"Phase1 mid:"<<endl;
   //PrintData(Tb1,ni,nj);
//------perform pivot operations------//
   Pivot(Tb1,ni,nj);
   cout<<"Phase1 optimal:"<<endl;
   PrintData(Tb1,ni,nj);
   double error = Tb1[0][nj-1];
   for(i=0;i<Ni;i++) {
      for(j=0;j<Nj-1;j++) Tb[i][j] = Tb1[i+1][j+1];
      Tb[i][Nj-1] = Tb1[i+1][nj-1];
   }
//-------output Phase one results---------//
   Pick_Var(Tb,Ni,Nj);
   if(Tb1[0][nj-1]>0) {
      cout<<"infeasible with error "<<error<<endl;
      XERR = -1;
   }
   return 1;
}
//---------Phase two: find the optimal solution if the problem is feasible------//
int Phase2() {
   if(XERR==-1) return 0;
   PH=2;
   Pivot(Tb,Ni,Nj);
   cout<<"Phase2 optimal:"<<endl;
   PrintData(Tb,Ni,Nj);
   Pick_Var(Tb,Ni,Nj);
   return 1;
}

int main(int argc, char *argv[])  
{
  if(argc!=2) {
     cout<<"Usage: Simplex [input_dat_file]"<<endl;
     return 0;
  }
  Tableau(argv[1]);
  Phase1();
  Phase2();
  //Results();

  return 1;

}
