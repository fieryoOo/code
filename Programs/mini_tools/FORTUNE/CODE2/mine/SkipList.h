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

#define LMAX 20

//Define a node which consist of a 'int' component that stores level#, a 'double' component that stores the key, and a pointer to an array of pointers that point to the next node at each level

struct nd {
   int nl;
   double key1, key2;
   void *p;
   struct nd ** pnext;
   struct nd ** plast;
};
typedef struct nd node;


class SkipList {
public:
   SkipList();
   int Onum, Msize, Msizec;
   node *head;
   node **FindPreNode (double key1, double key2);
   node **FindPreNodeOTF (struct Point pnt);
   node *Find (double key1, double key2);
   node *Insert (struct Bisector *he, double key1, double key2);
   node *InsertOTF (struct Bisector *he, struct Point pnt);
   node **Insertnext (node **pprev, struct Bisector *he);
   void Delete (double key1, double key2);
   struct Bisector* DeleteOTF (struct Bisector *bs, Point pnt);
   void DeleteNode (struct Bisector *bs);
   struct Bisector* GetMin ();
   struct Bisector* PopMin ();
   void PrintList();
   void PrintKeys();
   int LLength();
   uint64_t ClockGetTime();

private:
   int RandLevel();
   node *CreatNode (int nl, struct Bisector *he, double key1, double key2);
   double Breakpoint(struct Bisector *el,double ysl);
};
