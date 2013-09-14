#include <iostream>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <time.h>
#include <sys/time.h>
using namespace std;

//--------Clock, swap & median----------//
uint64_t ClockGetTime()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
}

void swap(double *a, double *b){
     double t=*a; *a=*b; *b=t;
}

int median(double *arr, int a, int b, int c){
     int sm, lg;
     if(arr[a]<arr[b]) {sm=a; lg=b;}
     else {sm=b; lg=a;}
     if(arr[c]<arr[sm]) return sm;
     else if(arr[c]>arr[lg]) return lg;
     else return c;
}

//----------optimized median-of-3 quicksort---------------//
int Opt_QuickSort(double *arr, int elements) {
  int max_l=elements;
  //if(max_l>1000) max_l=1000;
  int  beg[max_l], end[max_l], i=0, j, L, R;
  double piv;
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (L<R) {
      if (i==max_l-1) {cout<<"Increase max layer!"<<endl; return 0;}
      j=median(arr,L,(int)((L+R)/2),R);
      if(j!=L) swap(&arr[L],&arr[j]);
      piv=arr[L]; 
      while (L<R) {
        while (arr[R]>=piv && L<R) R--; if (L<R) arr[L++]=arr[R];
        while (arr[L]<=piv && L<R) L++; if (L<R) arr[R--]=arr[L]; 
      }
      arr[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L; 
    }
    else i--;
  }
  return 1; 
}

//-------------optimized median-of-3 quicksort combined with insertion sort---------------//
void InserSort(double *arr, int p, int r);
int Opt_ins_QuickSort(double *arr, int elements) {
  int max_l=elements;
  //if(max_l>1000) max_l=1000;
  int  beg[max_l], end[max_l], i=0, j, L, R, thres=20;
  double piv;
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; R=end[i]-1;
    if (R-L>thres) {
      if (i==max_l-1) {cout<<"Increase max layer!"<<endl; return 0;}
      j=median(arr,L,(int)((L+R)/2),R);
      if(j!=L) swap(&arr[L],&arr[j]);
      piv=arr[L];
      while (L<R) {
        while (arr[R]>=piv && L<R) R--; if (L<R) arr[L++]=arr[R];
        while (arr[L]<=piv && L<R) L++; if (L<R) arr[R--]=arr[L];
      }
      arr[L]=piv; beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
    }
    else i--;
  }
  InserSort(arr, 0, elements-1);
  return 1;
}

//-------------optimized 3-way median-of-3 quicksort combined with insertion sort---------------//
void InserSort(double *arr, int p, int r);
int Opt_ins_3w_QuickSort(double *arr, int elements) {
  int max_l=elements;
  //if(max_l>1000) max_l=1000;
  int  beg[max_l], end[max_l], i=0, j, L, R, l, r, thres=20;
  double piv;
  beg[0]=0; end[0]=elements;
  while (i>=0) {
    L=beg[i]; l=L; R=end[i]-1; r=R;
    if (R-L>thres) {
      if (i==max_l-1) {cout<<"Increase max layer!"<<endl; return 0;}
      j=median(arr,L,(int)((L+R)/2),R);
      if(j!=L) swap(&arr[L],&arr[j]);
      piv=arr[L];
      while (L<R) {
        for(;R>L;R--) {if(arr[R]<piv){arr[L++]=arr[R]; break;} if(arr[R]==piv)arr[R]=arr[r--];}
        for(;L<R;L++) {if(arr[L]>piv){arr[R--]=arr[L]; break;} if(arr[L]==piv)arr[L]=arr[l++];}
      }
      arr[L]=piv;
      for(j=beg[i],L=L-1;j<l;j++,L--) {if(L<l){L=j-1;break;}arr[j]=arr[L];arr[L]=piv;}
      for(j=end[i]-1,R=R+1;j>r;j--,R++) {if(R>r){R=j+1;break;}arr[j]=arr[R];arr[R]=piv;}
      beg[i+1]=R; end[i+1]=end[i]; end[i++]=L+1;
    }
    else i--;
  }
  InserSort(arr, 0, elements-1);
  return 1;
}

//----------------median-of-3 quicksort (combined with insertion sort)-------------------//
int M_Partition(double *arr, int p, int r) {
   int j=median(arr, p,(int)((p+r)/2),r);
   if(j!=r) swap(&arr[r],&arr[j]);
   double x=arr[r];
   int i=p-1;
   for(j=p;j<r;j++)
      if(arr[j]<=x){
         i++;
         swap(&arr[i],&arr[j]);
      }
   swap(&arr[i+1],&arr[r]);
   return i+1;
}
void QuickSort( double *arr, int p, int r) {
    int q;
    if(p<r) {
       q=M_Partition(arr,p,r);
       QuickSort(arr,p,q-1);
       QuickSort(arr,q+1,r);
    }
}
void ins_QuickSort( double *arr, int p, int r) {
    int q;
    if(p<r) {
       if(r-p>20) {
          q=M_Partition(arr,p,r);
          ins_QuickSort(arr,p,q-1);
          ins_QuickSort(arr,q+1,r);
       }
       else InserSort(arr,p,r);
    }
}

//--------------merge sort------------------//
/*
void Merge(double *arr, int p, int q, int r) {
   int i, j, k, n1 = q-p+1, n2 = r-q;
   double L[n1+1], R[n2+1];
   for(i=0; i<n1; i++) L[i] = arr[p+i];
   for(j=0; j<n2; j++) R[j] = arr[q+j+1];
   L[i] = 1e30; R[j] = 1e30;
   for(i=0,j=0,k=p;k<=r;k++)
      if (L[i]<=R[j]) arr[k] = L[i++];
      else arr[k] = R[j++];
}
*/
void Merge(double *arr, int p, int q, int r) {
   int i, j, k, n1 = q-p+1, n2 = r-q;
   double L[n1], R[n2];
   for(i=0; i<n1; i++) L[i] = arr[p+i];
   for(j=0; j<n2; j++) R[j] = arr[q+j+1];
   for(i=0,j=0,k=p;;k++)
      if(L[i]<=R[j]) {
         arr[k] = L[i++];
         if(i==n1) break;
      }
      else {
         arr[k] = R[j++];
         if(j==n2) break;
      }
   if(i==n1) for(k++;k<=r;k++,j++) arr[k] = R[j];
   else for(k++;k<=r;k++,i++) arr[k] = L[i];
}
void MergeSort(double *arr, int p, int r) {
   if(p<r) {
      if(r-p>20) {
         int q=(p+r)/2;
         MergeSort(arr, p, q);
         MergeSort(arr, q+1, r);
         Merge(arr,p,q,r);
      }
      else InserSort(arr,p,r);
   }
}

//--------------insertion sort-----------------//
void InserSort(double *arr, int p, int r) {
   int i,j;
   double dtmp;
   for(i=p+1;i<=r;i++) {
      dtmp=arr[i];
      for(j=i;j>p && dtmp<arr[j-1]; j--) arr[j]=arr[j-1];
      arr[j]=dtmp;
   }
}

//---------------selection sort------------------//
int selectsort( double *datax, int ndat ) {
  int i, j, ii;
  double datamin;

  for(i=0;i<ndat;i++) {
     datamin=datax[i]; ii=i;
     for(j=i+1;j<ndat;j++)
        if(datamin>datax[j]){
           datamin=datax[j];
           ii=j;
        }
     if(ii==i) continue;
     swap(&datax[i],&datax[ii]);
     //dataxx=datax[i];
     //datax[i]=datax[ii];
     //datax[ii]=dataxx;
  }

  return 1;
}

int compare (const void *a, const void *b) {
   double aa = *(double *)a, bb = *(double *)b;
   if(aa > bb) return 1;
   if(aa == bb) return -1;
   return 0;
}

int main(int argc, char *argv[]) {
   if(argc!=2){
      cout<<"usage: "<<argv[0]<<" [array size]"<<endl;
      return 0;
   }
   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::normal_distribution<> normal(0,100);
   boost::variate_generator< RNGType, boost::normal_distribution<> > Gauss(gener, normal);
   boost::uniform_int<> uniform(0,atoi(argv[1])/200);
   boost::variate_generator< RNGType, boost::uniform_int<> > Uniint(gener, uniform);

   int i, j, ndat=atoi(argv[1]);// thres=atoi(argv[2]);
   double tb, ts, te;
   double datG[ndat], datI[ndat], dat1[ndat];
   for(i=0;i<ndat;i++) { datG[i] = Gauss(); datI[i] = Uniint(); }
   cout<<"Array size: "<<ndat<<endl;

   for(i=0;i<ndat;i++) dat1[i]=datG[i];
   ts=ClockGetTime();
   qsort(dat1, ndat, sizeof(double), compare);
   te=ClockGetTime();
   cout<<"Internal qsort (Gaussian/Sorted/FewUnique): "<<(te-ts)/1000.<<"msec ";
   ts=ClockGetTime();
   qsort(dat1, ndat, sizeof(double), compare);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec ";
   for(i=0;i<ndat;i++) dat1[i]=datI[i];
   ts=ClockGetTime();
   qsort(dat1, ndat, sizeof(double), compare);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec"<<endl;

   for(i=0;i<ndat;i++) dat1[i]=datG[i];
   ts=ClockGetTime();
   Opt_ins_3w_QuickSort(dat1, ndat);
   te=ClockGetTime();
   cout<<"Opt inser 3w-quicksort (Gaussian/Sorted/FewUnique): "<<(te-ts)/1000.<<"msec ";
   ts=ClockGetTime();
   Opt_ins_3w_QuickSort(dat1, ndat);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec ";
   for(i=0;i<ndat;i++) dat1[i]=datI[i];
   ts=ClockGetTime();
   Opt_ins_3w_QuickSort(dat1, ndat);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec"<<endl;

   for(i=0;i<ndat;i++) dat1[i]=datG[i];
   ts=ClockGetTime();
   Opt_ins_QuickSort(dat1, ndat);
   te=ClockGetTime();
   cout<<"Opt inser quicksort (Gaussian/Sorted/FewUnique): "<<(te-ts)/1000.<<"msec ";
   ts=ClockGetTime();
   Opt_ins_QuickSort(dat1, ndat);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec ";
   for(i=0;i<ndat;i++) dat1[i]=datI[i];
   ts=ClockGetTime();
   Opt_ins_QuickSort(dat1, ndat);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec"<<endl;

   for(i=0;i<ndat;i++) dat1[i]=datG[i];
   ts=ClockGetTime();
   Opt_QuickSort(dat1, ndat);
   te=ClockGetTime();
   cout<<"Opt quicksort (Gaussian/Sorted/FeqUnique): "<<(te-ts)/1000.<<"msec ";
   ts=ClockGetTime();
   Opt_QuickSort(dat1, ndat);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec ";
   for(i=0;i<ndat;i++) dat1[i]=datI[i];
   ts=ClockGetTime();
   Opt_QuickSort(dat1, ndat);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec"<<endl;

   for(i=0;i<ndat;i++) dat1[i]=datG[i];
   ts=ClockGetTime();
   ins_QuickSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<"inser Quicksort (Gaussian/Sorted/FewUnique): "<<(te-ts)/1000.<<"msec ";
   ts=ClockGetTime();
   ins_QuickSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec ";
   for(i=0;i<ndat;i++) dat1[i]=datI[i];
   ts=ClockGetTime();
   ins_QuickSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec"<<endl;

   for(i=0;i<ndat;i++) dat1[i]=datG[i];
   ts=ClockGetTime();
   QuickSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<"Quicksort (Gaussian/Sorted/FeqUnique): "<<(te-ts)/1000.<<"msec ";
   ts=ClockGetTime();
   QuickSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec ";
   for(i=0;i<ndat;i++) dat1[i]=datI[i];
   ts=ClockGetTime();
   QuickSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec"<<endl;

   for(i=0;i<ndat;i++) dat1[i]=datG[i];
   ts=ClockGetTime();
   MergeSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<"inser Mergesort (Gaussian/Sorted/FeqUnique): "<<(te-ts)/1000.<<"msec ";
   ts=ClockGetTime();
   MergeSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec ";
   for(i=0;i<ndat;i++) dat1[i]=datI[i];
   ts=ClockGetTime();
   MergeSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec"<<endl;

   for(i=0;i<ndat;i++) dat1[i]=datG[i];
   ts=ClockGetTime();
   InserSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<"Insertion sort (Gaussian/Sorted/FewUnique): "<<(te-ts)/1000.<<"msec ";
   ts=ClockGetTime();
   InserSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec ";
   for(i=0;i<ndat;i++) dat1[i]=datI[i];
   ts=ClockGetTime();
   InserSort(dat1, 0, ndat-1);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec"<<endl;

   for(i=0;i<ndat;i++) dat1[i]=datG[i];
   ts=ClockGetTime();
   selectsort(dat1, ndat);
   te=ClockGetTime();
   cout<<"selection sort (Gaussian/Sorted/FewUnique):"<<(te-ts)/1000.<<"msec ";
   ts=ClockGetTime();
   selectsort(dat1, ndat);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec ";
   for(i=0;i<ndat;i++) dat1[i]=datI[i];
   ts=ClockGetTime();
   selectsort(dat1, ndat);
   te=ClockGetTime();
   cout<<(te-ts)/1000.<<"msec"<<endl;

/*
   for(j=0;j<=100;j=j+10) {
      for(i=0;i<ndat;i++) dat1[i]=dat0[i];
      ts=ClockGetTime();
      ins_QuickSort(dat1, 0, ndat-1, j);
      te=ClockGetTime();
      cout<<"thres "<<j<<": "<<(te-ts)/1000.<<"msec"<<endl;   
    }
*/
   return 1;
}

