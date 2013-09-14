#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "SkipList.h"

#define DELETED -2
#define le 0
#define re 1



struct	Memnode	{
   struct Memnode *Nnext;
};

struct MemBlock {
   struct Memnode *memory;
   struct MemBlock *next;
};

struct	Memlst {
   struct Memnode *head;
   int nodesize;
};

struct Point {
   double x;
   double y;
};

//node for site or vertex
struct Site {
   struct Point coord;
   int id;
   int ntask;
};

//node for the geometry of an bisector
struct Edge {
   double a,b,c;
   struct Site *endp[2];
   struct Site *site[2];
};

//node to store finalized waiting for output bisectors
struct FBisector {
   double x1,y1,x2,y2;
   struct FBisector* next;
};

//node to store the bisectors from both site and vertex events
struct Bisector {
   struct Edge *ELedge;
   struct Bisector *Blast, *Bnext, *Vnext;
   node * parent;
   char Bort;
   struct Site *vertex;
   double yvev;
};




class VoronoiFortunes {
 public:
   VoronoiFortunes();
   ~VoronoiFortunes();

   int Onum, OnumVL, OnumBL; //Recorders of memory sizes and
   int Msize, MsizeVL, MsizeBL; // atomic operation numbers
   void VoronoiGenerator(double *xValues, double *yValues, int numPoints, double minX, double maxX, double minY, double maxY);
   bool GetBisector(double *x1, double *y1, double *x2, double *y2);
   void Rewind();


 private:
   //Variables
   MemBlock* MemBlockVF;
   MemBlock* MemBlockcur;
   struct Memlst mBL;
   struct Memlst mST;
   struct Memlst mED;
   int MemSize;

   struct FBisector* BisectorList;
   struct FBisector* Bcurrent;
   struct Bisector *BLhead, *BLtail;
   int sitenum;
   int isite;
   struct Site *sites;
   struct Site *sbot;

   SkipList VL;
   SkipList BL;

   double RegionXl, RegionXh, RegionYl, RegionYh;

   //Subroutines
   bool RunFortune();
   //system & memories
   void MemRefresh();
   void BSFree();
   char *BuildMemLink(struct Memlst *fl);	
   void InitialMemList(struct Memlst *fl,int size);
   void NodeLink(struct Memnode *curr,struct Memlst *fl);
   void SiteInitialize(double* xdat, double* ydat, int datnum);
   void EventInitialize();
   int Pcmp(Site* a, Site* b);
   void Merge(Site *arr, int p, int q, int r);
   void MergeSort(Site *arr, int p, int r);

   //bisector list
   struct Bisector *NewBisector(struct Edge *e, int ort);
   void BLDelete(struct Bisector *bs);//(Point p);
   bool BLInitialize();
   node** BLInsertnext(node **prevs, struct Bisector *Bnew);
   void BLInsert(struct Bisector *Bnew, Point vcoord);
	
   //Bisector geometries
   double dist(struct Site *s,struct Site *t);
   struct Site* GetSite(struct Bisector *he, int flag);
   struct Edge *CreateEdge(struct Site *s1,struct Site *s2);

   //Site events
   struct Site *SiteNext();

   //Vertex events
   int VLempty();
   struct Point VLmin();
   void VLInsert(struct Bisector *he,struct Site * v, double offset);
   void VLDelete(struct Bisector *he);
   struct Site *NewVertex(struct Bisector *el1, struct Bisector *el2, struct Point *p=0);
   void release(struct Site *v);
   void SetEndp(struct Edge *e,int lr,struct Site * s);

   //Finalize and output
   void FinalizeBisector(struct Edge *e);
   void OutputBisector(double x1, double y1, double x2, double y2);

};
