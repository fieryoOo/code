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

//Define a node which consist of a 'int' component that stores level#, a 'double' component that stores the key, and a pointer to an array of pointers that point to the next node at each level

struct nd {
   int nl;
   double key;
   void *p;
   struct nd ** pnext;
};
typedef struct nd node ;


class SkipList {
public:
   SkipList();
   int ONum;
   node *head;
   node **FindNode (double key);
   node *Find (double key);
   node *Insert (double key);
   void Delete (double key);  
   void PrintList();
   void PrintKeys();
   int LLength();
   uint64_t ClockGetTime();

private:
   int RandLevel();
   node *CreatNode (int nl, double key);
};

SkipList::SkipList() {
   head = NULL;
   ONum = 0;
}
//Function to get the current clock time in nanosecond (to provide input for the random number generator)
uint64_t SkipList::ClockGetTime()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
}

//Function that returns a radom level number decided by flip a coin at each level.
int SkipList::RandLevel() {
   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::uniform_int<> uni(0,1);
   boost::variate_generator< RNGType, boost::uniform_int<> > coin(gener, uni);
   int i;
   for(i=1;i<LMAX;i++) if(coin()) break;
   return i;
}

//Creat a node that has nl level.
node* SkipList::CreatNode (int nl, double key) {
   node *p = (node *) malloc (sizeof(node));
   int i, T=3;
   p->nl = nl; T++;
   p->key = key; T++;
   p->pnext = (node **) malloc (nl * sizeof(node *)); T += 4;
   for(i=0;i<nl;i++) p->pnext[i] = NULL; T += nl*2;
   ONum += T;
   return p;
}

//Search for a key value in the SkipList. Return a pointer that points to an array of pointers which are associated with nodes prior to the key value on each level.
node** SkipList::FindNode (double key) {
   int i = head->nl-1, T=0;
   node **prev, *p = head; T++;
   prev = (node **) malloc ( (i+1) * sizeof(node *)); T += 5;
   for(;i>=0;i--) {
      T++;
      while(p->pnext[i]!=NULL && p->pnext[i]->key < key) { p = p->pnext[i]; T += 3; }
      prev[i] = p; T++;
   }
   ONum += T;
 if(prev[0]==NULL) {cout<<head->nl<<endl;}
   return prev;
}

node* SkipList::Find (double key) {
   if(head==NULL || head->nl==0) { ONum += 1; return NULL;} //cout<<"Find: empty list"<<endl;
   int T = 1;
   node **prev = FindNode(key); T++;
   T++;
   ONum += T;
   if(prev[0]->pnext[0]==NULL || prev[0]->pnext[0]->key!=key) {cout<<key<<" not exists in List"<<endl; return NULL;}
   else return prev[0]->pnext[0];
}

//Insert a node into the SkipList and return a pointer to the inserted node. Creat a head node if the SkipList is empty.
node* SkipList::Insert (double key) {
   int i, nl=RandLevel(), T=1;
   node *pinser = CreatNode(nl, key); T++;
   T++;
   if(head == NULL) {
      head = CreatNode(LMAX, -DBL_MAX); T++;
      head->nl = nl; T++;
      for(i=nl-1;i>=0;i--) { pinser->pnext[i] = NULL; head->pnext[i] = pinser; T+=3; }
   }
   else {
      i = head->nl; T+=2; if(nl > i) { head->nl = nl; T++; }
      if(i==0) for(i=nl-1;i>=0;i--) { pinser->pnext[i] = NULL; head->pnext[i] = pinser; T+=3; }
      else { 
         node **pprev = FindNode(key); T++;
         T+=2;
         if(pprev[0]->pnext[0]!=NULL && pprev[0]->pnext[0]->key==key) { 
            cout<<key<<" already in list, skipped."<<endl; 
            head->nl = i; T++; ONum += T; return NULL; 
         }
         for(i=nl-1;i>=0;i--) { pinser->pnext[i] = pprev[i]->pnext[i]; pprev[i]->pnext[i] = pinser; T+=3; }
      }
   }
   ONum += T;
   return pinser;
}

//Delet a node from the SkipList and free the memory.
void SkipList::Delete (double key) {
   if(head==NULL || head->nl==0) {ONum += 1; return;} //cout<<"Delete: empty list"<<endl;
   int T=1;
   node **pprev = FindNode(key); T++;
   T++;
   if(pprev[0]->pnext[0]==NULL || (pprev[0]->pnext[0])->key!=key) { 
      //cout<<key<<" not exists in list, skipped."<<endl; 
      ONum += T; return; 
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
   ONum += T;
}

//Function that prints a whole skip list.
void SkipList::PrintList() {
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

//Function that print the bottom-most linked list
void SkipList::PrintKeys() {
   node *p = head;
   if(p!=NULL) {
      while(p->pnext[0] != NULL) { p = p->pnext[0]; cout<<p->key<<" "; }
      cout<<endl;
   }
}

//Function that returns the length of the bottom-level linked list
int SkipList::LLength() {
   int L=0;
   node *p = head;
   if(p!=NULL)while(p->pnext[0] != NULL) { p = p->pnext[0]; L++; }
   return L;
}


#define nkey (int)pow(2.,5.)
int main (int argc, char *argv[]) {
   SkipList sl;
   node *p = NULL;

   double key[10] = {38., 59., 22., -78.23897, 1099.88, -2538999, 59., -78.23897, 0., 10};
   int i, iter;
   for(iter=0;iter<1;iter++) {
      sl.ONum = 0;
      for(i=0;i<10;i++) sl.Insert(key[i]);
   cout<<"After Inser "<<endl; sl.PrintList(); 
      cout<<sl.ONum<<endl;
      sl.Delete(22);
      sl.Delete(57.9);
      sl.Delete(-78.23897);
   cout<<"After Deletion:"<<endl; sl.PrintList(); 
//      cout<<(iter+1)*12<<" "<<sl.ONum<<endl;
   }

   return 1;
}
