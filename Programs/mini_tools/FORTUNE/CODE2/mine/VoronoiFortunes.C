#include "VoronoiFortunes.h"

// Constructor, initialize memories
VoronoiFortunes::VoronoiFortunes() {
   isite = 0;
   sites = NULL;

   //initialize memory block
   BisectorList = NULL;
   Bcurrent = NULL;
   MemBlockVF = new MemBlock; 
   MemBlockVF->memory = NULL;
   MemBlockVF->next = NULL;
   MemBlockcur = MemBlockVF;
}

//Deconstructor, cleanup memories.
VoronoiFortunes::~VoronoiFortunes() {
   MemRefresh();
   BSFree();
   if(MemBlockVF != 0) delete MemBlockVF;
}

//The public Voronoi generator function. Takes 2-D points and region to be computed as inputs
void VoronoiFortunes::VoronoiGenerator(double *xdat, double *ydat, int datnum, double minX, double maxX, double minY, double maxY) {
   if( minX > maxX || minY > maxY ) {
      cout<<"Incorrect region inputs!"<<endl;
      exit(0);
   }
   RegionXl = minX; RegionXh = maxX;
   RegionYl = minY; RegionYh = maxY;

   Onum = 15; OnumVL = 0; OnumBL = 0;
   Msize = sizeof(VoronoiFortunes) + sizeof(MemBlock);
   MsizeVL = 0; MsizeBL = 0;
   

   MemRefresh();
   BSFree();

   SiteInitialize(xdat, ydat, datnum);
   MemSize = 4*(int)sqrt((double)sitenum);
   InitialMemList(&mED, sizeof(Edge));

   RunFortune();
   Bcurrent = BisectorList;
}

//The public function to access bisectors after 'VoronoiGenerator' is done. Returns the next bisector.
//set rewind = 1 to go back to the first bisector.
void VoronoiFortunes::Rewind() {
   Bcurrent = BisectorList;
}
bool VoronoiFortunes::GetBisector(double *x1, double *y1, double *x2, double *y2) {
   if( Bcurrent == NULL ) return false;

   *x1 = Bcurrent->x1;
   *x2 = Bcurrent->x2;
   *y1 = Bcurrent->y1;
   *y2 = Bcurrent->y2;

   Bcurrent = Bcurrent->next;

   Onum += 5;
   return true;
}

//Initialize 'sites' and site region
void VoronoiFortunes::SiteInitialize(double* xdat, double* ydat, int datnum) {
   sitenum = datnum;
   InitialMemList(&mST, sizeof(Site));
   sites = new Site[sitenum];
   Msize += sitenum * sizeof(Site);
   Onum += 5;
   if(sites == NULL) {
      cout<<"SiteInitialize: Cannot allocate memory for sites"<<endl;
      exit(-1);
   }

   isite = 0;
   int i, j;
   for(i = 0; i< sitenum; i++) {
      sites[i].coord.x = xdat[i];
      sites[i].coord.y = ydat[i];
   }

   MergeSort(sites, 0, sitenum-1);

   sites[0].ntask = 0; sites[0].id = 0;
   Onum += 5 + 3*sitenum;
   for(i=1; i< sitenum; ) {
      if( sites[i].coord.x==sites[i-1].coord.x && sites[i].coord.y==sites[i-1].coord.y ) {
         for(j=i; j<sitenum--; j++) sites[j] = sites[j+1];
         Onum += 3 + 2*(sitenum-i);
      }
      else { 
         sites[i].ntask = 0; sites[i].id = i; i++; 
         Onum += 4;
      }
   }
   cout<<"number of effecient sites: "<<sitenum<<endl;
}

//MergeSort
int VoronoiFortunes::Pcmp(Site* a, Site* b) {
   struct Point pa = a->coord, pb = b->coord;
   if(pa.y < pb.y) { Onum += 3; return -1; }
   if(pa.y > pb.y) { Onum += 4; return 1; }
   if(pa.x < pb.x) { Onum += 5; return -1; }
   if(pa.x > pa.x) { Onum += 6; return 1; }
   Onum += 7; return 0;
}
void VoronoiFortunes::Merge(Site *arr, int p, int q, int r) {
   int i, j, k, n1 = q-p+1, n2 = r-q;
   Site L[n1], R[n2];
   Onum += 5;
   for(i=0; i<n1; i++) L[i] = arr[p+i];
   for(j=0; j<n2; j++) R[j] = arr[q+j+1];
   Onum += n1+n2;
   for(i=0,j=0,k=p;;k++)
      if(Pcmp(&L[i], &R[j])<=0) {
         arr[k] = L[i++];
         if(i==n1) break;
      }
      else {
         arr[k] = R[j++];
         if(j==n2) break;
      }
   if(i==n1) for(k++;k<=r;k++,j++) arr[k] = R[j];
   else for(k++;k<=r;k++,i++) arr[k] = L[i];
   Onum += 2 * (n1+n2);
}
void VoronoiFortunes::MergeSort(Site *arr, int p, int r) {
   if(p<r) {
      int q=(p+r)/2;
      MergeSort(arr, p, q);
      MergeSort(arr, q+1, r);
      Merge(arr,p,q,r);
      Onum += 6;
   }
   Onum ++;
}

bool VoronoiFortunes::BLInitialize() {
   InitialMemList(&mBL, sizeof *BLhead);
   BLhead = NewBisector( (struct Edge *)NULL, 0);
   BLtail = NewBisector( (struct Edge *)NULL, 0);
   BLhead -> Blast = (struct Bisector *)NULL;
   BLhead -> Bnext = BLtail;
   BLtail -> Blast = BLhead;
   BLtail -> Bnext = (struct Bisector *)NULL;
   Onum += 7;
   return true;
}

//get the next site in the list
struct Site * VoronoiFortunes::SiteNext() {
   struct Site *s;
   Onum ++;
   if(isite < sitenum) {
      s = &sites[isite];
      isite += 1;
      Onum += 2;
      return(s);
   }
   else return( (struct Site *)NULL);
}

struct Bisector* VoronoiFortunes::NewBisector(struct Edge *e,int ort) {
   struct Bisector *Bnew;
   Bnew = (struct Bisector *) BuildMemLink(&mBL);
   Bnew -> ELedge = e;
   Bnew -> Bort = ort;
   Bnew -> Vnext = (struct Bisector *) NULL;
   Bnew -> vertex = (struct Site *) NULL;
   Onum += 5;
   return(Bnew);
}


node** VoronoiFortunes::BLInsertnext(node **prevs, struct Bisector *Bnew) {
   struct Bisector *Bpre;
   if(prevs[0]==NULL || prevs[0]->p==NULL) Bpre = BLhead;
   else Bpre  = (struct Bisector *)(prevs[0]->p);
   Bnew -> Blast = Bpre;
   Bnew -> Bnext = Bpre->Bnext;
   Bpre->Bnext->Blast = Bnew;
   Bpre->Bnext = Bnew;
   Onum += 7;
   return BL.Insertnext(prevs, Bnew);
}

void VoronoiFortunes::BLInsert(struct Bisector *Bnew, Point vcoord) {
   struct Bisector *Bpre = (struct Bisector *)(BL.InsertOTF(Bnew, vcoord)->p);
   if(Bpre==NULL) Bpre = BLhead;
   Bnew->Blast = Bpre;
   Bnew->Bnext = Bpre->Bnext;
   Bpre->Bnext->Blast = Bnew;
   Bpre->Bnext = Bnew;
   Onum += 7;
}

void VoronoiFortunes::BLDelete(struct Bisector *bs) {
   //vcoord.x -= 1.e-8;
//   struct Bisector *bs = 
   BL.DeleteNode(bs);//vcoord);
   (bs -> Blast) -> Bnext = bs -> Bnext;
   (bs -> Bnext) -> Blast = bs -> Blast;
   bs -> ELedge = (struct Edge *)DELETED;
   Onum += 10;
}

//Return site of the input bisector: flag=0 for the left one and flag=1 for the right one
struct Site * VoronoiFortunes::GetSite(struct Bisector *bs, int ort) {
   Onum += 2;
   if(bs->ELedge == (struct Edge *)NULL) return(sbot);
   Onum += 3;
   return( bs->Bort == ort ? bs -> ELedge -> site[le] : bs -> ELedge -> site[re]);
}

//calculate the line formula of the bisector in the form ax+by-c=0
struct Edge * VoronoiFortunes::CreateEdge(struct Site *s1, struct Site *s2) {
   double dx,dy;
   struct Edge *enew = (struct Edge *) BuildMemLink(&mED);

   enew -> site[0] = s1; //set site[0] to be the lower site
   enew -> site[1] = s2; //site[1] to be the higher
   (s1->ntask)++;
   (s2->ntask)++;
   enew -> endp[0] = (struct Site *) NULL;
   enew -> endp[1] = (struct Site *) NULL;
   dx = s2->coord.x - s1->coord.x;
   dy = s2->coord.y - s1->coord.y;
   enew->a = dx/dy;
   enew->b = 1.;
   double x1 = s1->coord.x, x2 = s2->coord.x, y1 = s1->coord.y, y2 = s2->coord.y;
   enew -> c = 0.5 * (x2*x2 - x1*x1 + y2*y2 - y1*y1) / dy;
   Onum += 28;
	
    return(enew);
}

//Calculate the intersection of two bisectors and produce a new vertex event if necessary
struct Site * VoronoiFortunes::NewVertex(struct Bisector *el1, struct Bisector *el2, struct Point *p) {
   struct Edge *e1,*e2, *e;
   struct Bisector *el;
   double d, xint, yint;
   struct Site *v;
	
   e1 = el1 -> ELedge;
   e2 = el2 -> ELedge;
   Onum += 5;
   if(e1 == (struct Edge*)NULL || e2 == (struct Edge*)NULL) return ((struct Site *) NULL);
   if (e1->site[1] == e2->site[1]) return ((struct Site *) NULL); //from same site

   d = e1->a * e2->b - e1->b * e2->a;
   Onum == 5;
   if ( fabs(d)<1.0e-10 ) return ((struct Site *) NULL); //parallel
	
   xint = (e1->c*e2->b - e2->c*e1->b)/d;
   yint = (e2->c*e1->a - e1->c*e2->a)/d;
	
   if( (e1->site[1]->coord.y < e2->site[1]->coord.y) || (e1->site[1]->coord.y == e2->site[1]->coord.y && e1->site[1]->coord.x < e2->site[1]->coord.x) ) {	
      el = el1; 
      e = e1;
   }
   else {	
      el = el2; 
      e = e2;
   }
	
   int fpos = xint >= e -> site[1] -> coord.x;
   Onum += 21;
   if ((fpos && el -> Bort == le) || (!fpos && el -> Bort == re)) return ((struct Site *) NULL);
	
   //create a new vertex event.
   v = (struct Site *) BuildMemLink(&mST);
   v -> ntask = 0;
   v -> coord.x = xint;
   v -> coord.y = yint;
   Onum == 4;
   return(v);
}

void VoronoiFortunes::SetEndp(struct Edge *e,int lr,struct Site * s) {
   e -> endp[lr] = s;
   (s->ntask)++;
   Onum += 3;
   if(e -> endp[re-lr]== (struct Site *) NULL) return;

   FinalizeBisector(e);
   release(e->site[le]);
   release(e->site[re]);
   NodeLink((Memnode*)e, &mED);
   Onum += 4;
}

//compute distance between two sites
double VoronoiFortunes::dist(struct Site *s1,struct Site *s2) {
   double dx,dy;
   dx = s1->coord.x - s2->coord.x;
   dy = s1->coord.y - s2->coord.y;
   Onum += 8;
   return sqrt(dx*dx + dy*dy);
}

//reduce task number of input site by 1 and delete if all tasks are done.
void VoronoiFortunes::release(struct Site *v) {
   if( (v->ntask)-- == 0 ) NodeLink((Memnode*)v, &mST);
   Onum += 2;
}

//Store the vertex in the input bisector and insert into the vertex-events skiplist
void VoronoiFortunes::VLInsert(struct Bisector *bs, struct Site * v, double Rc) {
   bs->vertex = v;
   (v->ntask)++;
   bs->yvev = v->coord.y + Rc;
   VL.Insert(bs, bs->yvev, v->coord.x);
   Onum += 5;
}

//delete the bisector from the vertex-event skip list
void VoronoiFortunes::VLDelete(struct Bisector *bs) {
   if(bs -> vertex == (struct Site *) NULL) return;
   VL.Delete(bs->yvev, bs->vertex->coord.x);
   release(bs -> vertex);
   Onum += 3;
   bs->vertex = (struct Site *) NULL;
}

//Check if the vertex list is empty
int VoronoiFortunes::VLempty() {
   Onum += 2;
   return ( VL.head == NULL || VL.head->pnext[0] == NULL );
}

//return the coordinates of the lower-most vertex event
struct Point VoronoiFortunes::VLmin() {
   struct Bisector* min = VL.GetMin();
   //if( min==NULL ) return (struct Point)NULL;
   struct Point pt;
   pt.x = min->vertex->coord.x;
   pt.y = min->yvev;
   Onum == 5;
   return pt;
}

void VoronoiFortunes::InitialMemList(struct Memlst *mlst,int size) {
   mlst -> head = (struct Memnode *) NULL;
   mlst -> nodesize = size;
   Onum += 2;
}

//allocate memory block, divide it into linked nodes according to the nodesize of the input memory list head
char * VoronoiFortunes::BuildMemLink(struct Memlst *mlst) {
   int i; 
   struct Memnode *mNew;

   if(mlst->head == (struct Memnode *) NULL) {	
      mNew =  (struct Memnode *) malloc(MemSize * mlst->nodesize);
      Onum += 2;
      if(mNew == NULL) return 0;
      Msize += MemSize * mlst->nodesize;
		
      MemBlockcur->next = new MemBlock;
      Msize += sizeof(MemBlock);
      MemBlockcur = MemBlockcur->next;
      MemBlockcur->memory = mNew;
      MemBlockcur->next = NULL;
      for(i=0; i<MemSize; i++) NodeLink((struct Memnode *)((char *)mNew+i*mlst->nodesize), mlst);
      Onum += 8+MemSize*4;
   }
   mNew = mlst->head;
   mlst->head = (mlst->head)->Nnext;
   Onum += 3;
   return((char *)mNew);
}


//link a memory node at the front of the input memory list
void VoronoiFortunes::NodeLink(struct Memnode *mNew, struct Memlst *mlst) {
   mNew -> Nnext = mlst -> head;
   mlst -> head = mNew;
   Onum += 2;
}

void VoronoiFortunes::MemRefresh() {
   if(sites != NULL) { delete[] sites; sites = NULL; }

   MemBlock* MBcur = NULL, *Mpre = NULL;

   Onum += 4;
   for(MBcur = MemBlockVF; MBcur->next!=NULL;) {
      Mpre = MBcur;
      MBcur = MBcur->next;
      free(Mpre->memory);
      delete Mpre;
      Mpre = NULL;
      Onum += 6;
   }
   if(MBcur != 0 && MBcur->memory != 0) {
      free(MBcur->memory);
      delete MBcur;
      MBcur = NULL;
      Onum += 3;
   }

   MemBlockVF = new MemBlock;
   Msize += sizeof(MemBlock);
   MemBlockVF->next = NULL;
   MemBlockVF->memory = NULL;
   MemBlockcur = MemBlockVF;
   Onum += 8;
}

void VoronoiFortunes::BSFree() {
   FBisector* BScur = NULL, *BSpre = NULL;
   BScur = BSpre = BisectorList;

   while(BScur != NULL && BScur->next != NULL) {
      BSpre = BScur;
      BScur = BScur->next;
      delete BSpre;
      BSpre = NULL;
      Onum += 5;
      //Msize -= sizeof(FBisector *);
   }

   BisectorList = NULL;
   Onum += 5;
}

void VoronoiFortunes::OutputBisector(double x1, double y1, double x2, double y2) {
   FBisector* Bnew = new FBisector;
   Msize += sizeof(FBisector);
   Bnew->next = BisectorList;
   BisectorList = Bnew;
   Bnew->x1 = x1;
   Bnew->y1 = y1;
   Bnew->x2 = x2;
   Bnew->y2 = y2;
   Onum += 8;
}

//Finalize the bisectors by cutting the non-bounded side at the boundart of the input region
void VoronoiFortunes::FinalizeBisector(struct Edge *e) {
   struct Site *sl, *sr;
   double xl, yl, xr, yr;
   double a = e->a, b = e->b, c = e->c;

   sl = e->endp[0]; sr = e->endp[1];
   if( sl == (struct Site *)NULL ) {
      xl = RegionXl;
      yl = (-a*xl+c)/b;
   }
   else {
      xl = sl->coord.x;
      yl = sl->coord.y;
   }
   if( sr == (struct Site *)NULL ) {
      xr = RegionXh;
      yr = (-a*xr+c)/b;
   }
   else {
      xr = sr->coord.x;
      yr = sr->coord.y;
   }
   Onum += 20;
   //skip if the edge is out of the interested region
   if( xl>=RegionXh || xr<=RegionXl ) return;
   if( yl<=yr ) { if( yr<=RegionYl || yl>=RegionYh ) return; }
   else if( yl<=RegionYl || yr>=RegionYh) return;

   OutputBisector(xl,yl,xr,yr);
   Onum += 9;
}


bool VoronoiFortunes::RunFortune() {
   struct Site *sitecur, *v, *bot, *top, *stmp, *p;
   struct Point Csite, Cvertex;
   struct Bisector *bsl, *bsr, *bsll, *bsrr, *bsnew;
   int ort, eflag;
   struct Edge *e;
   node **prevs;

   sbot = SiteNext();
   BLInitialize();

   sitecur = SiteNext();
   Onum += 3;
   while( !(sitecur==NULL && VLempty()) ) {
      //Check to see if the next event is a site or vertex event.
      if(sitecur==NULL) { eflag = 1; Cvertex = VLmin(); }
      else if(VLempty()) eflag = 0;
      else {
         Csite = sitecur->coord;
         Cvertex = VLmin();
         if( Csite.y<Cvertex.y || (Csite.y==Cvertex.y && Csite.x<Cvertex.x) ) eflag = 0;
         else eflag = 1;
      }
      Onum += 10;
                
      if(eflag) { //vertex event
	 bsl = VL.PopMin();  //Pop out the Bisector to the left of the vertex event (which stores the vertex event)
	 bsll = bsl -> Blast;
	 bsr = bsl -> Bnext;
	 bsrr = bsr -> Bnext;

	 bot = GetSite(bsl, le);  //Get the site to the left of the left bisector
	 top = GetSite(bsr, re);  //Get the site to the right of the right bisector
	 v = bsl->vertex; //Get the vertex location of the event
	 SetEndp(bsl->ELedge,bsl->Bort,v);    //set the endpoint of the left Bisector to be this vector
	 SetEndp(bsr->ELedge,bsr->Bort,v);    //set the endpoint of the right Bisector to be this vector
	 VLDelete(bsr);
         BLDelete(bsl); //Cvertex//Delete the 2 bisectors that intersect at this vertex. Potential error may occur
         BLDelete(bsr); // when 3 bisectors intersect together. Can be fixed by using a doublly linked skip list
         ort = le;
         if (bot->coord.y > top->coord.y) { stmp = bot; bot = top; top = stmp; ort = re; } //decide which site is higher
         e = CreateEdge(bot, top); //create a new bisector between the left and right site
         bsnew = NewBisector(e, ort); // that heading to the direction decided by ort
         BLInsert(bsnew, Cvertex); //Insert it into the Bisector list
         SetEndp(e, re-ort, v); //set the starting point of the new Bisector to v
         release(v); //delete vertex v from the Vertex list
         //if left and new bisector donnot intersect, reset the left bisector
         if((p = NewVertex(bsll, bsnew)) != (struct Site *) NULL) {
            VLDelete(bsll);
            VLInsert(bsll, p, dist(p,bot));
         }
         //if right and new bisector donnot intersect, reset the new bisector
         if ((p = NewVertex(bsnew, bsrr)) != (struct Site *) NULL) 
	    VLInsert(bsnew, p, dist(p,bot));
         Onum += 30;
      }
      else { //site
	 prevs = BL.FindPreNodeOTF(sitecur->coord); //Get the first bisector to the left of the site event
	 if( prevs[0]==NULL ) bsl = BLhead; // and a list of the previous node in the SkipList
	 else if( prevs[0]->p==NULL ) bsl = BLhead;
	 else bsl = (struct Bisector *)(prevs[0] -> p);
	 bsr = bsl -> Bnext; //Get the first bisector to the right of the site event
	 bot = GetSite(bsl, re); //Get the higher site of the left bisector
	 e = CreateEdge(bot, sitecur); //create a new Bisector between the new site and
	 bsnew = NewBisector(e, le); // that point to the left
	 prevs = BLInsertnext(prevs, bsnew); //Insert the new bisector and update the 'previous' list

	 if ((p = NewVertex(bsl, bsnew)) != (struct Site *) NULL) {//update the vertex event recorded on the left edge
	    VLDelete(bsl);
	    VLInsert(bsl, p, dist(p,sitecur));
	 }
	 bsnew = NewBisector(e, re); //The node that stores the same bisector but point to the right.
	 BLInsertnext(prevs, bsnew); //Insert it
	 if ((p = NewVertex(bsnew, bsr)) != (struct Site *) NULL)	//update the vertex event record on the right edge
	    VLInsert(bsnew, p, dist(p,sitecur));
	 sitecur = SiteNext();
         Onum += 18;
      }
   }

   //Finalilze the remaining bisectors in the list
   for(bsl=BLhead -> Bnext; bsl != BLtail; bsl=bsl -> Bnext) {	
      e = bsl -> ELedge;
      FinalizeBisector(e);
      Onum += 3;
   }

   //Update memory-size and operation-number recorders;
   OnumVL += VL.Onum; OnumBL += BL.Onum;
   MsizeVL += VL.Msize; MsizeBL += BL.Msize;
cout<<Msize<<" "<<VL.Msize<<" "<<BL.Msize<<endl;

   MemRefresh();
   return true;
}
