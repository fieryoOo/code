#include <iostream>
using namespace std;

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
