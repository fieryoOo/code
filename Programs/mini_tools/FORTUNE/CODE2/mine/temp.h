#ifndef VORONOI_DIAGRAM_GENERATOR
#define VORONOI_DIAGRAM_GENERATOR

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "SkipList.h"


#ifndef NULL
#define NULL 0
#endif
#define DELETED -2

#define le 0
#define re 1



struct	Freenode	
{
	struct	Freenode *nextfree;
};

struct FreeNodeArrayList
{
	struct	Freenode* memory;
	struct	FreeNodeArrayList* next;

};

struct	Freelist	
{
	struct	Freenode	*head;
	int		nodesize;
};

struct Point	
{
	float x,y;
};

// structure used both for sites and for vertices 
struct Site	
{
	struct	Point	coord;
	int		sitenbr;
	int		refcnt;
};



struct Edge	
{
	float   a,b,c;
	struct	Site 	*ep[2];
	struct	Site	*reg[2];
	int		edgenbr;

};

struct GraphEdge
{
	float x1,y1,x2,y2;
	struct GraphEdge* next;
};




struct Halfedge 
{
	struct	Halfedge	*ELleft, *ELright;
	struct	Edge	*ELedge;
	int		ELrefcnt;
	char	ELpm;
	struct	Site	*vertex;
	float	ystar;
	struct	Halfedge *Vnext;
};




class VoronoiDiagramGenerator
{
public:
	VoronoiDiagramGenerator();
	~VoronoiDiagramGenerator();

	bool generateVoronoi(float *xValues, float *yValues, int numPoints, float minX, float maxX, float minY, float maxY);

	void resetIterator()
	{
		iteratorEdges = allEdges;
	}

	bool getNext(float& x1, float& y1, float& x2, float& y2)
	{
		if(iteratorEdges == 0)
			return false;
		
		x1 = iteratorEdges->x1;
		x2 = iteratorEdges->x2;
		y1 = iteratorEdges->y1;
		y2 = iteratorEdges->y2;

		iteratorEdges = iteratorEdges->next;

		return true;
	}


private:
	bool voronoi();
	//system & memories
	void cleanup();
	void cleanupEdges();
	char *myalloc(unsigned n);
	char *getfree(struct Freelist *fl);	
	void freeinit(struct Freelist *fl,int size);
	void makefree(struct Freenode *curr,struct Freelist *fl);
	void geominit();

	//bisector list
	struct	Halfedge **ELhash;
	struct Halfedge *ELleftbnd(struct Point *p);
	struct	Halfedge *HEcreate(struct Edge *e,int pm);
	void ELdelete(struct Halfedge *he);
	bool ELinitialize();
        void ELinsert(struct    Halfedge *lb, struct Halfedge *newHe);
	struct Halfedge *ELgethash(int b);

	//Bisector geometries
        float dist(struct Site *s,struct Site *t);
        struct Site *leftreg(struct Halfedge *he);
        struct Site *rightreg(struct Halfedge *he);
        struct Edge *bisect(struct Site *s1,struct Site *s2);
        float Breakpoint(struct Halfedge *el,float ysl);

	//Site events
	struct Site *nextone();

	//Vertex events
	int VLempty();
	struct Point VLmin();
	void VLinsert(struct Halfedge *he,struct Site * v, float offset);
        void VLdelete(struct Halfedge *he);
	void endpoint(struct Edge *e,int lr,struct Site * s);
	void clip_line(struct Edge *e);
	struct Site *NewVertex(struct Halfedge *el1, struct Halfedge *el2, struct Point *p=0);

	void ref(struct Site *v);
	void deref(struct Site *v);


	void pushGraphEdge(float x1, float y1, float x2, float y2);


	struct  Freelist	hfl;
	struct	Halfedge *ELleftend, *ELrightend;
	int 	ELhashsize;

	float	xmin, xmax, ymin, ymax, deltax;

	struct	Site	*sites;
	int		nsites;
	int		siteidx;
	int		sqrt_nsites;
	int		vertnum;
	struct 	Freelist sfl;
	struct	Site	*bottomsite;

	int		edgenum;
	struct	Freelist efl;

	struct SkipList VL;

	float	pxmin, pxmax, pymin, pymax;

	float borderMinX, borderMaxX, borderMinY, borderMaxY;

	FreeNodeArrayList* allMemoryList;
	FreeNodeArrayList* currentMemoryBlock;

	GraphEdge* allEdges;
	GraphEdge* iteratorEdges;

};

int scomp(const void *p1,const void *p2);


#endif


