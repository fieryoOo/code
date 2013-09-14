#include "VoronoiFortunes.h"

SkipList::SkipList() {
   head = NULL;
   Onum = 0;
   Msize = 0;
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
node* SkipList::CreatNode (int nl, struct Bisector *bs, double key1, double key2) {
cout<<"CreatNode"<<endl;
   node *p = (node *) malloc (sizeof(node));
cout<<"!!!created at "<<p<<endl;
   Msize += sizeof(node);
   int i, T=3;
   p->nl = nl; T++;
   p->key1 = key1; 
   p->key2 = key2;
   p->p = bs;
   T += 3;
   p->pnext = (node **) malloc (nl * sizeof(node *)); T += 4;
   p->plast = (node **) malloc (nl * sizeof(node *));
   Msize += nl * sizeof(node *);
   for(i=0;i<nl;i++) { p->pnext[i] = NULL; p->plast[i] = NULL;} T += nl*2;
   Onum += T;
cout<<"Creatend"<<endl;
   return p;
}

//Search for a key value in the SkipList. Return a pointer that points to an array of pointers which are associated with nodes prior to the key value on each level.
node** SkipList::FindPreNode (double key1, double key2) {
   int i = head->nl-1, T=0;
   node **prev, *p = head; T++;
   prev = (node **) malloc ( (i+1) * sizeof(node *)); T += 5;
   Msize += (i+1) * sizeof(node *);
   for(;i>=0;i--) {
      T++;
      while(p->pnext[i]!=NULL && ( p->pnext[i]->key1 < key1 || p->pnext[i]->key1 == key1 && p->pnext[i]->key2 < key2 ) ) 
         { p = p->pnext[i]; T += 3; }
      prev[i] = p; T++;
   }
   Onum += T;
   //if(prev[0]==NULL) {cout<<head->nl<<endl;}
   return prev;
}

//Search for a key value, where the key of each node in the list is calculated on the fly.
double SkipList::Breakpoint(struct Bisector *el,double ysl) {
cout<<"in Breakpoint"<<endl;
   struct Edge *e = el -> ELedge;
   double x1 = e->site[0]->coord.x, y1 = e->site[0]->coord.y;
   double x2 = e->site[1]->coord.x, y2 = e->site[1]->coord.y;
//cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<" "<<ysl<<endl;
   double a = y2-y1;
   double b = 2*(x1*(ysl-y2)+x2*(y1-ysl));
   double c = -(ysl-y1)*(ysl-y2)*(y2-y1)+x1*x1*(y2-ysl)-x2*x2*(y1-ysl);
   double k = b*b-4*a*c;
   Onum += 31;
   if( k<0 ) {
      printf("Error! need debug!\n");
      return -1;
   }
   Onum += 4;
cout<<"Breakpend"<<endl;
   if(el->Bort == le) return (-b-sqrt(k))/2./a;
   else return (-b+sqrt(k))/2./a;
}
node** SkipList::FindPreNodeOTF (struct Point pnt) {
   int i;
   node **prev = (node **) malloc ( LMAX * sizeof(node *));
   Msize += LMAX * sizeof(node *);
   for(i=0;i<LMAX;i++) prev[i] = head;
   if( head==NULL || head->nl==0 ) return prev;
   i = head->nl-1;
   int T=6;
   double key = pnt.x;
   node *p=head;
   for(;i>=0;i--) {
      T++;
      while(p->pnext[i]!=NULL && Breakpoint( (struct Bisector *)(p->pnext[i]->p), pnt.y ) < key )
         { p = p->pnext[i]; T += 3; }
      prev[i] = p; T++;
   }
   Onum += T;
 if(prev[0]==NULL) {cout<<head->nl<<endl;}
//   struct Edge *tmp = ((struct Bisector *)(prev[0]->pnext[0]->p))->ELedge;
//cout<<"cc"<<endl;
//   if(tmp==NULL) cout<<"NULL"<<endl;
//   else cout<<"Found: "<<tmp->a<<" "<<tmp->b<<" "<<tmp->c<<"  bis(site): "<< Breakpoint( (struct Bisector *)(prev[0]->p), pnt.y )<<" (vertex): "<<Breakpoint( (struct Bisector *)(prev[0]->pnext[0]->p), pnt.y )<<endl;
   return prev;
}

node* SkipList::Find (double key1, double key2) {
   if(head==NULL || head->nl==0) { Onum += 1; return NULL;} //cout<<"Find: empty list"<<endl;
   int T = 1;
   node **prev = FindPreNode(key1, key2); T++;
   T++;
   Onum += T;
   if(prev[0]->pnext[0]==NULL || !((prev[0]->pnext[0])->key1==key1 && (prev[0]->pnext[0])->key2==key2)) 
      {cout<<key1<<" - "<<key2<<" not exists in List"<<endl; return NULL;}
   else return prev[0]->pnext[0];
}

//Insert a node into the SkipList and return a pointer to the inserted node. Create a head node if the SkipList is empty.
node* SkipList::Insert (struct Bisector *he, double key1, double key2) {
cout<<"Insert"<<endl;
   int i, nl=RandLevel(), T=1;
   node *pinser = CreatNode(nl, he, key1, key2); T++;
  // he->parent = pinser;
   T++;
   if(head == NULL) {
      head = CreatNode(LMAX, (struct Bisector *) NULL, -DBL_MAX, -DBL_MAX); T++;
      head->nl = nl; T++;
      for(i=nl-1;i>=0;i--) { pinser->pnext[i] = NULL; pinser->plast[i] = head; head->pnext[i] = pinser; T+=3; }
   }
   else {
      i = head->nl; T+=2; 
      if(nl > i) { head->nl = nl; T++; }
      if(i==0) for(i=nl-1;i>=0;i--) { pinser->pnext[i] = NULL; pinser->plast[i] = head; head->pnext[i] = pinser; T+=3; }
      else { 
         node **pprev = FindPreNode(key1, key2); T++;
         T+=2;
         if(pprev[0]->pnext[0]!=NULL && pprev[0]->pnext[0]->key1==key1 && pprev[0]->pnext[0]->key2==key2) { 
            cout<<key1<<" - "<<key2<<" already in list, skipped."<<endl; 
            head->nl = i; T++; Onum += T; return NULL; 
         }
         for(i=nl-1;i>=0;i--) { 
            pinser->pnext[i] = pprev[i]->pnext[i]; pprev[i]->pnext[i] = pinser; T+=3; 
            pinser->plast[i] = pprev[i];
            if(pinser->pnext[i]!=NULL) pinser->pnext[i]->plast[i] = pinser;
         }
      }
   }
   Onum += T;
cout<<"Ins end"<<endl;
   return pinser;
}

//Insert a node by by FindPredNodeOTF
node* SkipList::InsertOTF (struct Bisector *bs, struct Point pnt) {
cout<<"InserOTF"<<endl;
   int i, nl=RandLevel(), T=1;
   node *pinser = CreatNode(nl, bs, 0, 0), **pprev;
   bs->parent = pinser;
cout<<"I 11"<<endl;
   T+=2;
   if(head == NULL) {
cout<<"null head"<<endl;
      head = CreatNode(LMAX, (struct Bisector *) NULL, -DBL_MAX, -DBL_MAX); T++;
      head->nl = nl; T++;
      for(i=nl-1;i>=0;i--) { pinser->pnext[i] = NULL; pinser->plast[i] = head; head->pnext[i] = pinser; T+=3; }
   }
   else {
cout<<"not null"<<endl;
      i = head->nl; T+=2;
      if(nl > i) { head->nl = nl; T++; }
      if(i==0) for(i=nl-1;i>=0;i--) { 
         pinser->pnext[i] = NULL; pinser->plast[i] = head;
         head->pnext[i] = pinser; T+=4; 
      }
      else {
cout<<"I 12"<<endl;
         pprev = FindPreNodeOTF(pnt);
cout<<"I 13"<<endl;
         T+=2;
         for(i=nl-1;i>=0;i--) { 
            pinser->pnext[i] = pprev[i]->pnext[i]; pprev[i]->pnext[i] = pinser; T+=3; 
            pinser->plast[i] = pprev[i];
            if(pinser->pnext[i]!=NULL) pinser->pnext[i]->plast[i] = pinser;
         }
      }
   }
   Onum += T;
cout<<"InsOTF end"<<endl;
   return pprev[0];
}

//Insert a new node right next to a given node. Return an array of pointers to the previous position of the next node.
node** SkipList::Insertnext (node **pprev, struct Bisector *bs) {
cout<<"Insertnext"<<endl;
   int i, nl=RandLevel(), T=1;
   node *pinser = CreatNode(nl, bs, 0, 0); T++;
   bs->parent = pinser;
   T++;
   if(pprev==NULL || (pprev[0]==NULL && head!=NULL) ) {
      cout<<"Skiplist.Insertnext error: NULL reference node"<<endl;
      return (node **)NULL;
   }
   if(head==NULL) {
      if(pprev[0]!=NULL) { cout<<"Insertnext to NULL head!"<<endl; return NULL; }
      head = CreatNode(LMAX, (struct Bisector *) NULL, -DBL_MAX, -DBL_MAX); T++;
      head->nl = nl; T++;
      for(i=nl-1;i>=0;i--) { 
         pinser->pnext[i] = NULL; pinser->plast[i] = head;
         head->pnext[i] = pinser; pprev[i] = pinser; T+=4; 
      }
      for(i=nl;i<LMAX;i++) pprev[i] = head;
   }
   else {
      i = head->nl; T+=2;
      if(nl > i) { head->nl = nl; T++; }
      for(i=nl-1;i>=0;i--) {
         pinser->pnext[i] = pprev[i]->pnext[i]; pprev[i]->pnext[i] = pinser;
         pinser->plast[i] = pprev[i]; pprev[i] = pinser;
         if(pinser->pnext[i]!=NULL) (pinser->pnext[i])->plast[i] = pinser; T+=5;
      }
   }
   Onum += T;
cout<<"Insnext end"<<endl;
   return pprev;
}

//Delete a node from the SkipList and free the memory.
void SkipList::Delete (double key1, double key2) {
   if(head==NULL || head->nl==0) {Onum += 1; return;} //cout<<"Delete: empty list"<<endl;
   int T=1;
   node **pprev = FindPreNode(key1, key2); T++;
   T++;
   if(pprev[0]->pnext[0]==NULL || !((pprev[0]->pnext[0])->key1==key1 && (pprev[0]->pnext[0])->key2==key2)) { 
      //cout<<key<<" not exists in list, skipped."<<endl; 
      Onum += T; return; 
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
   Onum += T;
}

//Delete a node on the fly at Point pnt, which point to Bisector he
struct Bisector* SkipList::DeleteOTF (struct Bisector *he, Point pnt) {
cout<<"DeleteOTF"<<endl;
   if(head==NULL || head->nl==0) {Onum += 1; return NULL;} //cout<<"Delete: empty list"<<endl;
   int T=1;
   node **pprev = FindPreNodeOTF(pnt);
   node *ptmp = pprev[0]->pnext[0];
   T += 2;
   int i, nl=pprev[0]->pnext[0]->nl; T++;
   struct Bisector *pbs = (struct Bisector *)((pprev[0]->pnext[0])->p);
   for(i=nl-1;i>=0;i--) {
      pprev[i]->pnext[i] = pprev[i]->pnext[i]->pnext[i];
      T += 2;
   }
   free(pprev);
   free(ptmp);
   T+=2;
   if(nl == head->nl) {
      for(i=nl-1;i>=0;i--) { T+=2; if(head->pnext[i]!=NULL) break; }
      head->nl = i+1; T++;
   }
   Onum += T;
   return pbs;
}

//Delete a node give the previous node of it
void SkipList::DeleteNode (struct Bisector *bs) {
cout<<"DeleteNode"<<endl;
   node *pdelete = bs->parent;
   int i, nl = pdelete->nl, T=0;
   if(nl == head->nl) {
      for(i=nl-1;i>=0;i--) { T+=2; if(head->pnext[i]!=NULL) break; }
      head->nl = i+1; T++;
   }
   for(i=nl-1;i>=0;i--) {
      if(pdelete->plast[i]==NULL) {cout<<"!!!shit is at"<<pdelete<<endl;if(pdelete==head)cout<<"is head"<<endl;}
      if(pdelete->plast[i]!=NULL) (pdelete->plast[i])->pnext[i] = pdelete->pnext[i];
      if(pdelete->pnext[i]!=NULL) (pdelete->pnext[i])->plast[i] = pdelete->plast[i];
      T += 5;
   }
   free(pdelete);
   Onum += T;
cout<<"end"<<endl;
}

//return the first node in the SkipList.
struct Bisector* SkipList::GetMin () {
   if(head==NULL || head->nl==0) {Onum += 1; return (struct Bisector*) NULL;}
   node* hmin = head->pnext[0];
   if(hmin==NULL) {Onum += 2; return (struct Bisector*) NULL;}
   return (Bisector*)(hmin->p);
}

//Pop out the first node in the SkipList.
struct Bisector* SkipList::PopMin () {
   if(head==NULL || head->nl==0) {Onum += 1; return (struct Bisector*) NULL;}
   node* hmin = head->pnext[0];
   if(hmin==NULL) {Onum += 2; return (struct Bisector*) NULL;}
   for(int i=hmin->nl-1;i>=0;i--) {
      head->pnext[i] = hmin->pnext[i];
      Onum += 2;
   }
   return (Bisector*)(hmin->p);
}

//Function that prints a whole skip list.
void SkipList::PrintList() {
   node *p = head;
   for(int i=head->nl-1;i>=0;i--) {
      while(p->pnext[0] != NULL) {
         p = p->pnext[0];
         if(p->nl<=i) cout<<setw(10)<<"NaN"<<" ";
         else cout<<setw(10)<<p->key1<<"-"<<p->key2<<" "; 
      }
      cout<<endl;
      p = head;
   }
}

//Function that print the bottom-most linked list
void SkipList::PrintKeys() {
   node *p = head;
   if(p!=NULL) {
      while(p->pnext[0] != NULL) { p = p->pnext[0]; cout<<p->key1<<"-"<<p->key2<<" "; }
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
