#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <stdlib.h>
#include <float.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <sys/time.h>
using namespace std;

#define LMAX 15

//Function to get the current clock time in nanosecond (to provide input for the random number generator)
uint64_t ClockGetTime()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
}

//Define a node which consist of a 'int' component that stores level#, a 'double' component that stores the key, and a pointer to an array of pointers that point to the next node at each level
struct nd {
   int nl;
   double key;
   struct nd ** pnext;
};
typedef struct nd node ;

//Function that returns a radom level number decided by flip a coin at each level.
int RandLevel() {
   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::uniform_int<> uni(0,1);
   boost::variate_generator< RNGType, boost::uniform_int<> > coin(gener, uni);
   int i;
   for(i=1;i<LMAX;i++) if(coin()) break;
   return i;
}

//Creat a node that has nl level.
node *CreatNode (int nl, double key, int *Ta) {
   node *p = (node *) malloc (sizeof(node)); *Ta = *Ta+3;
   int i;
   p->nl = nl; *Ta=*Ta+1;
   p->key = key; *Ta=*Ta+1;
   p->pnext = (node **) malloc (nl * sizeof(node *)); *Ta=*Ta+4;
   for(i=0;i<nl;i++) p->pnext[i] = NULL; *Ta=*Ta+nl*2;
   return p;
}

//Search for a key value in the SkipList. Return a pointer that points to an array of pointers which are associated with nodes prior to the key value on each level.
node **FindNode (node *head, double key, int *Ta) {
   int i = head->nl-1, T=0;
   node **prev, *p = head; T++;
   prev = (node **) malloc ( (i+1) * sizeof(node *)); T += 5;
   for(;i>=0;i--) {
      T++;
      while(p->pnext[i]!=NULL && p->pnext[i]->key < key) { p = p->pnext[i]; T += 3; }
      prev[i] = p; T++;
   }
   *Ta = *Ta + T;
 if(prev[0]==NULL) {cout<<head->nl<<endl;}
   return prev;
}

node *Find (node *head, double key, int *Ta) {
   if(head==NULL || head->nl==0) { *Ta = 1; return NULL;} //cout<<"Find: empty list"<<endl;
   int T = 1;
   node **prev = FindNode(head, key, &T); T++;
   T++;
   *Ta = T;
   if(prev[0]->pnext[0]==NULL || prev[0]->pnext[0]->key!=key) {cout<<key<<" not exists in List"<<endl; return NULL;}
   else return prev[0]->pnext[0];
}

//Insert a node into the SkipList and return a pointer to the inserted node. Creat a head node if the SkipList is empty.
node *Insert (node **head, double key, int *Ta) {
   int i, nl=RandLevel(), T=1;
   node *pinser = CreatNode(nl, key, &T); T++;
   T++;
   if(*head == NULL) {
      *head = CreatNode(LMAX, -DBL_MAX, &T); T++;
      (*head)->nl = nl; T++;
      for(i=nl-1;i>=0;i--) { pinser->pnext[i] = NULL; (*head)->pnext[i] = pinser; T+=3; }
   }
   else {
      i = (*head)->nl; T+=2; if(nl > i) { (*head)->nl = nl; T++; }
      if(i==0) for(i=nl-1;i>=0;i--) { pinser->pnext[i] = NULL; (*head)->pnext[i] = pinser; T+=3; }
      else { 
         node **pprev = FindNode(*head, key, &T); T++;
         T+=2;
         if(pprev[0]->pnext[0]!=NULL && pprev[0]->pnext[0]->key==key) { 
            //cout<<key<<" already in list, skipped."<<endl; 
            (*head)->nl = i; T++; *Ta = T; return NULL; 
         }
         for(i=nl-1;i>=0;i--) { pinser->pnext[i] = pprev[i]->pnext[i]; pprev[i]->pnext[i] = pinser; T+=3; }
      }
   }
   *Ta = T;
   return pinser;
}

//Delet a node from the SkipList and free the memory.
void Delete (node *head, double key, int *Ta) {
   int T=1;
   if(head==NULL || head->nl==0) {*Ta = T; return;} //cout<<"Delete: empty list"<<endl;
   node **pprev = FindNode(head, key, &T); T++;
   T++;
   if(pprev[0]->pnext[0]==NULL || (pprev[0]->pnext[0])->key!=key) { 
      //cout<<key<<" not exists in list, skipped."<<endl; 
      *Ta = T; return; 
   }
   int i, nl=pprev[0]->pnext[0]->nl; T++;
   for(i=nl-1;i>=0;i--) {
      pprev[i]->pnext[i] = pprev[i]->pnext[i]->pnext[i]; 
      T += 2;
   }
   free(pprev);
   T+=2;
   if(nl == head->nl) {
      for(i=nl-1;i>=0;i--) { T+=2; if(head->pnext[i]!=NULL) break; }
      head->nl = i+1; T++;
   }
   *Ta = T;
}

//Function that prints a whole skip list.
void PrintList(node *head) {
   node *p = head;
   for(int i=head->nl-1;i>=0;i--) {
      while(p->pnext[0] != NULL) {
         p = p->pnext[0];
         if(p->nl<=i) cout<<setw(10)<<"NaN"<<" ";
         else cout<<setw(10)<<p->key<<" "; 
      }
      cout<<endl;
      p = head;
   }
}

int LLength(node *head) {
   int L=0;
   node *p = head;
   if(p!=NULL)while(p->pnext[0] != NULL) { p = p->pnext[0]; L++; }
   return L;
}

void PlotKeys(node *head) {
   node *p = head;
   if(p!=NULL) {
      while(p->pnext[0] != NULL) { p = p->pnext[0]; cout<<p->key<<" "; }
      cout<<endl;
   }
}

#define nkey (int)pow(2.,5.)
int main (int argc, char *argv[]) {
   //initialize uniform integer generator
   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::uniform_int<> uniform(0,nkey);
   boost::variate_generator< RNGType, boost::uniform_int<> > Uniint(gener, uniform);

   int iiter, i, L, Ta;
   int NIkey[nkey]={ 0 }, NFkey[nkey]={ 0 }, NDkey[nkey]={ 0 }, TaI[nkey]={ 0 }, TaF[nkey]={ 0 }, TaD[nkey]={ 0 };
   double key[nkey];
   node *SL = NULL, *p = NULL;
   //run the experiment 100 times
   for(iiter=0;iiter<1;iiter++) {
      //produce input data randomly
      for(i=0;i<nkey;i++) { key[i] = Uniint(); cout<<key[i]<<" "; } cout<<endl;
      //perform Insert, Find and Delete operation and record the atomic operation number
      for(i=0;i<nkey;i++) { 
         L = LLength(SL); NIkey[L]++;
         Insert(&SL, key[i], &Ta); TaI[L] += Ta;
         PlotKeys(SL);
         L = LLength(SL); NFkey[L]++;
         p = Find(SL, key[i], &Ta); TaF[L] += Ta;
         cout<<p->key<<endl;
      }
      for(i=0;i<nkey;i++) { 
         L = LLength(SL); NDkey[L]++; 
         Delete(SL, key[i], &Ta); TaD[L] += Ta; 
         PlotKeys(SL);
      }
   }
   //compute avrage atomic operation number of all the 100 trials
   for(L=0;L<nkey;L++) {
      if(NIkey[L]!=0) TaI[L] /= NIkey[L]; 
      if(NFkey[L]!=0) TaF[L] /= NFkey[L]; 
      if(NDkey[L]!=0) TaD[L] /= NDkey[L];
   }

   FILE *ff;
   //output results
   ff=fopen("Insert_Atomic.txt","w");
   for(i=0;i<nkey;i++) if(NIkey[i]>100) fprintf(ff,"%d %d\n", i, TaI[i]);
   fclose(ff);
   ff=fopen("Find_Atomic.txt","w");
   for(i=0;i<nkey;i++) if(NFkey[i]>100) fprintf(ff,"%d %d\n", i, TaF[i]);
   fclose(ff);
   ff=fopen("Delete_Atomic.txt","w");
   for(i=0;i<nkey;i++) if(NDkey[i]>100) fprintf(ff,"%d %d\n", i, TaD[i]);
   fclose(ff);

   return 1;
}
