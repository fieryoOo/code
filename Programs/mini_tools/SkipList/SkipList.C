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
   for(i=1;i<=LMAX;i++) if(coin()) break;
   return i;
}

//Creat a node that has nl level.
node *CreatNode (int nl, double key) {
   node *p = (node *) malloc (sizeof(node));
   int i;
   p->nl = nl;
   p->key = key;
   p->pnext = (node **) malloc (nl * sizeof(node *));
   for(i=0;i<nl;i++) p->pnext[i] = NULL;
   return p;
}

//Search for a key value in the SkipList. Return a pointer that points to an array of pointers which are associated with nodes prior to the key value on each level.
node **FindPreNode (node *head, double key) {
   int i = head->nl-1;
   node **prev, *p = head;
   prev = (node **) malloc ( (i+1) * sizeof(node *));
   for(;i>=0;i--) {
      while(p->pnext[i]!=NULL && p->pnext[i]->key < key) p = p->pnext[i];
      prev[i] = p;
   }
   return prev;
}

node *Find (node *head, double key) {
   if(head==NULL || head->nl==0) {cout<<"Find: empty list"<<endl; return NULL;}
   node **prev = FindPreNode(head, key);
   if(prev[0]->pnext[0]==NULL || prev[0]->pnext[0]->key!=key) {cout<<key<<" not exists in List"<<endl; return NULL;}
   else return prev[0]->pnext[0];
}

//Insert a node into the SkipList and return a pointer to the inserted node. Creat a head node if the SkipList is empty.
node *Insert (node **head, double key) {
   int i, nl=RandLevel();
   node *pinser = CreatNode(nl, key);
   if(*head == NULL) {
      *head = CreatNode(LMAX, -DBL_MAX);
      (*head)->nl = nl;
      for(i=nl-1;i>=0;i--) { pinser->pnext[i] = NULL; (*head)->pnext[i] = pinser; }
   }
   else {
      i = (*head)->nl; if(nl > i) (*head)->nl = nl;
      if(i==0) for(i=nl-1;i>=0;i--) { pinser->pnext[i] = NULL; (*head)->pnext[i] = pinser; }
      else {
         node **pprev = FindPreNode(*head, key);
         if(pprev[0]->pnext[0]!=NULL && pprev[0]->pnext[0]->key==key) {
            (*head)->nl = i; return NULL;
         }
         for(i=nl-1;i>=0;i--) { pinser->pnext[i] = pprev[i]->pnext[i]; pprev[i]->pnext[i] = pinser; }
      }
   }
   return pinser;
}

//Delet a node from the SkipList and free the memory.
void Delete (node *head, double key) {
   if(head==NULL || head->nl==0) {cout<<"Delete: empty list"<<endl; return;}
   node **pprev = FindPreNode(head, key);
   if(pprev[0]->pnext[0]==NULL || pprev[0]->pnext[0]->key!=key) { 
      cout<<key<<" not exists in list, skipped."<<endl; return; 
   }
   int i, nl=pprev[0]->pnext[0]->nl;
   for(i=nl-1;i>=0;i--) {
      //if(pprev[i]->pnext[i]==NULL || pprev[i]->pnext[i]->key!=key) break;
      pprev[i]->pnext[i] = pprev[i]->pnext[i]->pnext[i]; 
   }
   if(nl == head->nl) {
      for(i=nl-1;i>=0;i--) if(head->pnext[i]!=NULL) break;
      head->nl = i+1;
   }
   //free()
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

int main (int argc, char *argv[]) {
   int i;
   node *SL = NULL;
   double key[10] = {38., 59., 22., -78.23897, 1099.88, -2538999, 59., -78.23897, 0., 10};
   for(i=0;i<10;i++) Insert(&SL, key[i]);
   cout<<"After Inser "<<endl; PrintList(SL);
   Delete(SL, 22);
   Delete(SL, 57.9); 
   cout<<"After Deletion:"<<endl; PrintList(SL);

   return 1;
}
