#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <limits.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <sys/time.h>
using namespace std;

struct edge {
   struct vertex *v; //neighboring vertex
   struct edge *next; //next edge in the adj list
};

struct vertex {
   int num;
   int vflag; //visit flag
   int dist; //path distance
   struct vertex *pre; //previous vertex on the path
   struct edge *e1; //first neighboring edge in the adj list
};

struct Queue {
   struct vertex *v;
   struct Queue *next;
};

// Enqueue a vertex
void EnQueue (struct vertex *u, struct Queue **Q, struct Queue **Qtail, int *T) {
   struct Queue *Qpn;
   Qpn = (struct Queue*) malloc (sizeof(struct Queue));
   Qpn->v = u;
   Qpn->next = NULL;
   *T = *T + 4;
   if(*Qtail==NULL) {*Qtail = Qpn; *Q = Qpn; return; }
   (*Qtail)->next = Qpn;
   *Qtail = Qpn;
   *T = *T + 2;
   return;
}

// Dequeue and return a vertex
struct vertex* DeQueue (struct Queue **Q, struct Queue **Qtail, int *T) {
   if(*Q==NULL) return NULL;
   struct Queue *head;
   struct vertex *v;
   head = *Q;
   v = head->v;
   *Q = (*Q)->next;
   if(*Q==NULL) *Qtail = NULL;
   free(head);
   *T = *T + 4;
   return v;
}

// Insert an edge to a vertex in the adj list
void InsertEdge (struct vertex *u, struct vertex *v) {
   struct edge *ep, *epn;
   epn = (struct edge*) malloc (sizeof(struct edge));
   epn->v = v;
   epn->next = NULL;
   if(u->e1==NULL) { u->e1 = epn; return; }
   for(ep=u->e1;ep->next!=NULL;ep=ep->next) {}
   ep->next = epn;
}

// Build the adj list from an adj matrix with n vertices
struct vertex* AdjList (int **A, int n, int *SP) {
   int i, j;
   struct vertex *u;
   u = (struct vertex *) malloc (n * sizeof(struct vertex));
   for(i=0;i<n;i++) { u[i].num = i; u[i].e1 = NULL; *SP = *SP + 2;}
   for(i=0;i<n;i++)
      for(j=i;j<n;j++)
         if(A[i][j]>0) {
            InsertEdge(&u[i], &u[j]);
            InsertEdge(&u[j], &u[i]);
            *SP = *SP + 1;
         }
   return u;
}

// Breadth-first search given a source vertex. Return the maximum distance
int BFS(struct vertex *s, struct vertex *v, int n, int *T, int *SP) {
   struct Queue *Q = NULL, *Qtail = NULL;
   struct vertex u;
   struct edge *e;
   int i, dm=0;
   //float inf = std::numeric_limits<float>::infinity();
   int inf = std::numeric_limits<int>::max();
   *T = *T + 2;
   //initialize vertices
   for(i=0;i<n;i++) {
      v[i].vflag=0;
      v[i].dist=inf;
      v[i].pre=NULL;
      *SP = *SP + 3;
      *T = *T + 4;
   }
   s->vflag = 1;
   s->dist = 0;
   *T = *T + 2;
   EnQueue(s, &Q, &Qtail, T);
   //cout<<Q->v->num<<" "<<Qtail->v->num<<endl;
   *SP = *SP + 1;
   while(Q) {
      u = *DeQueue(&Q, &Qtail, T); //visit a vertex u
      for(e=u.e1;e;e=e->next) { //EnQueue a neighboring vertex if it is not visited
         *T = *T + 2;
         if(e->v->vflag==0) {
            e->v->vflag = 1;
            e->v->dist = u.dist + 1;
            e->v->pre = &u;
            *T = *T + 3;
            EnQueue(e->v, &Q, &Qtail, T);
            *SP = *SP + 1;
         }
      }
      if(u.dist > dm) dm = u.dist;
      *T = *T + 2;
   }
//   for(i=0;i<n;i++) if(v[i].dist > dm) dm = v[i].dist;
   return dm;
}

// read the AdjMatrix from an input file given the number of vertices
int ** ReadMatrix(char *fname, int nv) {
   FILE *fin;
   char buff[300], *ptk;
   int i, j, **A;
   A = (int **) malloc (nv * sizeof(int *));
   for(i=0;i<nv;i++) A[i] = (int *) malloc (nv * sizeof(int));
   if((fin=fopen(fname,"r"))==NULL) {
      cout<<"Cannot open file "<<fname<<endl;
      exit(0);
   }
   cout<<"AdjMatrix read in: "<<endl;
   for(i=0;i<nv;i++) {
      if(fgets(buff, 300, fin)==NULL) { 
         cout<<"Only "<<i<<" rows in file "<<fname<<". Stopped!"<<endl;
         exit(0);
      }
      ptk = strtok(buff," \n");
      A[i][0] = atoi(ptk);
      cout<<A[i][0]<<" ";
      for(j=1;j<nv;j++) {
         ptk = strtok(NULL," \n");
         if(ptk==NULL) {
            cout<<"Only "<<j<<" columns in file "<<fname<<". Stopped!"<<endl;
            exit(0);
         }
         A[i][j] = atoi(ptk);
         cout<<A[i][j]<<" ";
      }
      cout<<endl;
   }
   fclose(fin);
   return A;
}

//Function to get the current clock time in nanosecond (to provide input for the random number generator)
uint64_t ClockGetTime()
{
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
}

//Produce a symmetric matrix with its component set to 1 or 0 randomly with a probability computed from the input expectation
int ** RandomMatrix (int c, int nv) {
   typedef boost::mt19937 RNGType;
   RNGType gener( ClockGetTime() );
   boost::uniform_real<> uni(0,1);
   boost::variate_generator< RNGType, boost::uniform_real<> > dcoin(gener, uni);
   int i, j, **A;
   double p = (double)c/(nv-1);
   A = (int **) malloc (nv * sizeof(int *));
   for(i=0;i<nv;i++) A[i] = (int *) malloc (nv * sizeof(int));
   for(i=0;i<nv;i++) {
      A[i][i] = 0;
      for(j=i+1;j<nv;j++) {
         if(dcoin()<p) A[i][j] = 1;
         else A[i][j] = 0;
         A[j][i] = A[i][j];
      }
   }
//   cout<<"Matrix produced: "<<endl;
//   for(i=0;i<nv;i++) {
//      for(j=0;j<nv;j++) cout<<A[i][j]<<" ";
//      cout<<endl;
//   }
   return A;
}

//Free the memory of AdjMatrix and AdjList
void ReleaseMemory(int **A, struct vertex *G, int n) {
   int i;
   struct edge *ep, *eptmp;
   struct vertex *vp;
   for(i=0;i<n;i++) free(A[i]);
   free(A);
   for(i=0;i<n;i++)
      for(ep=G[i].e1;ep!=NULL;) {
         eptmp = ep;
         ep=ep->next;
         free(eptmp);
      }
   free(G);
}


int main (int argc, char *arg[]) {
//   if(argc!=3) {
//      cout<<"usage: "<<arg[0]<<" [Input_AdjMatrix_file] [number of vertices]"<<endl;
//      return -1;
//   }
//   int nv = atoi(arg[2]), n;
   //int **A = ReadMatrix(arg[1], nv); // read in the adj matrix from file
   int i, j, dist, dmax=0, c=5, n, niter;
   int **B;
   int Ttmp, SPmax, SPtmp;
   double davg, Tavg, SPavg;
   struct vertex* G;
   FILE *fout = fopen("vertex.num_diameter_T_SP.txt","w");
   for(n=6;n<=900;n++) {
      cout<<"Working on vertex# "<<n<<"..."<<endl;
      davg = 0.; Tavg = 0.; SPavg = 0.;
      niter = 300-n;
      if(niter<20) niter = 20;
      for(j=0;j<niter;j++) {
         B = RandomMatrix(c, n); // produce an Erdos-Renyi random graph with n vertices and expectation c
         G=AdjList (B, n, &SPtmp); // convert the matrix into an adj list
         SPavg += SPtmp;
         dmax = 0; SPmax = 0;
         for(i=0;i<n;i++) {
            Ttmp = 0; SPtmp = 0;
            dist = BFS(&G[i], G, n, &Ttmp, &SPtmp); // run a BFS on each vertex to find the maximum distance
            Tavg += Ttmp;
            if(SPmax<SPtmp) SPmax = SPtmp;
            if(dmax<dist) dmax = dist;
         }
         ReleaseMemory(B,G,n);
         SPavg += SPmax;
         davg += dmax;
         //cout<<dmax<<endl;
      }
      davg /= niter;
      SPavg /= niter;
      Tavg /= niter;
      fprintf(fout, "%d %lf %lf %lf\n", n, davg, Tavg, SPavg);
   }
   fclose(fout);
   return 1;
}
