#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <stdlib.h>
#include <limits.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <sys/time.h>
using namespace std;

#define Parent(i) ( (i-1)>>1 )
#define Right(i) ( (i<<1)+2 )
#define Left(i) ( (i<<1)+1 )

//---Tree Node structure---//
struct Node {
   unsigned f;
   char c;
   Node *left, *right;
};

//---Min Heap structure---//
struct MinHeap {
   unsigned length;
   unsigned size;
   Node **node;
};

//---swap function---//
void swap(Node **a, Node **b) {
   Node *Ntmp = *a;
   *a = *b; *b = Ntmp;
}

//---Create single Tree Node---//
Node *CreateNode(char c, unsigned f) {
   Node *np = (Node *) malloc(sizeof(Node));
   np->c = c; np->f = f;
   np->left = NULL; np->right = NULL;
   return np;
}

//---Heapify function for min heap---//
void MinHeapify(MinHeap *mh, unsigned i) {
   unsigned il = Left(i), ir = Right(i), is = i;
   if(il < mh->size && mh->node[il]->f < mh->node[is]->f) is = il;
   if(ir < mh->size && mh->node[ir]->f < mh->node[is]->f) is = ir;
   if(is != i) {
      swap(&(mh->node[i]), &(mh->node[is]));
      MinHeapify(mh, is);
   }
}

//---Extract min node---//
Node *ExtractMin(MinHeap *mh) {
   if(mh->size < 1) { cout<<"Extraction failed: heap underflow!"<<endl; return NULL; }
   Node *np = mh->node[0];
   //swap(&mh->node[0], &mh->node[--(mh->size)]);
   mh->node[0] = mh->node[--(mh->size)];
   MinHeapify(mh, 0);
   return np;
}

//---Decrease key---//
void Decrease(MinHeap *mh, unsigned i, unsigned f) {
   if(f > mh->node[i]->f) { cout<<"Decrease failed: "<<mh->node[i]->f<<" -> "<<f<<endl; return; }
   mh->node[i]->f = f;
   unsigned ip = Parent(i);
   Node *Ntmp = mh->node[i];
   while(i>0 && mh->node[ip]->f > f) {
      mh->node[i] = mh->node[ip];
      i = ip; ip = Parent(i);
   }
   mh->node[i] = Ntmp;
}

//---Insert a new node---//
void InsertMin(MinHeap *mh, Node *nd) {
   if(mh->size == mh->length) {cout<<"Insertion failed: heap overflow!"<<endl; return; }
   mh->node[mh->size] = nd;
   Decrease(mh, (mh->size)++, nd->f);
}

//---Create min heap---//
MinHeap *CreateMinHeap(char c[], unsigned f[], unsigned nc) {
   int i;
   MinHeap *mh = (MinHeap *) malloc (sizeof(MinHeap));
   mh->size = nc; mh->length = nc;
   mh->node = (Node **) malloc (nc * sizeof(Node *));
   for(i=0; i<nc; i++) mh->node[i] = CreateNode(c[i], f[i]);
   for(i=Parent(nc-1); i>=0; i--) MinHeapify(mh, i);
   return mh;
}

//---Build the Tree of Huffman code. Store the depth of each non-leaf node in node.c---//
Node *HuffmanTree(char c[], unsigned f[], unsigned nc) {
   Node *ndl, *ndr, *ndp;
   MinHeap *mh = CreateMinHeap(c, f, nc);
   unsigned depl, depr;
int j;
   for(int i=0; i<nc-1; i++) {
      ndl = ExtractMin(mh);
      if(ndl->left==NULL && ndl->right==NULL) depl = 0;
      else depl = ndl->c;
      ndr = ExtractMin(mh);
      if(ndr->left==NULL && ndr->right==NULL) depr = 0;
      else depr = ndr->c;
      ndp = CreateNode((depl>depr ? depl:depr) + 1, ndl->f + ndr->f);
      ndp->left = ndl; ndp->right = ndr;
      InsertMin(mh, ndp);
    //cout<<"After iteration "<<i<<": "; for(j=0;j<mh->size;j++) cout<<mh->node[j]->c<<" "<<mh->node[j]->f<<" "<<mh->node[j]<<" "; cout<<endl;
   }
   return ExtractMin(mh);
}

//---Traverse the HuffmanTree and produce an array of Node pointers for all the 128 possible chars---//
void Traverse(Node *mh, int code[], int iroot, int **cp) {
   int i;
   if(mh->left) { code[iroot] = 0; Traverse(mh->left, code, iroot+1, cp); }
   if(mh->right) { code[iroot] = 1; Traverse(mh->right, code, iroot+1, cp); }
   if((mh->left)==NULL && (mh->right)==NULL) {
      cp[mh->c] = (int *) malloc ((iroot+1) * sizeof(int));
      for(i=0; i<iroot; i++) cp[mh->c][i] = code[i];
      cp[mh->c][i] = -1;
      return;
   }
}
int **CharIndex(Node *mh) {
   int **cp = (int **) malloc (128 * sizeof(int *));
   for(int i=0; i<128; i++) cp[i] = NULL;
   int code[mh->c];
   Traverse(mh, code, 0, cp);
   return cp;
}

//---Print out a Tree on the screen
void gotoxy(int x, int y)
{
    printf("%c[%d;%df",0x1B,y,x);
}
void PrintT(Node *mh, int x, int y, int flag) {
   char ctmp = mh->c;
   if(mh->left || mh->right) ctmp=' ';
   gotoxy(x-1,y);
   if(flag==0) cout<<"[]";
   else if(flag==-1) cout<<"["<<ctmp;
   else cout<<ctmp<<"]";
   if(mh->left) PrintT(mh->left, x-(int)(50./pow(2.,y-1.)+0.5), y+1, -1);
   if(mh->right) PrintT(mh->right, x+(int)(50/pow(2.,y-1.)+0.5), y+1, 1);
}
void PrintTree(Node *mh) {
   system("clear");
   PrintT(mh, 100, 1, 0);
   gotoxy(1,(mh->c)+2);
   cout<<endl;
}

//---count through the input string to get all input characters and their frequencies
void count_freq(char *str, char *c, unsigned *f, unsigned *nc) {
   int i, j;
   for(i=0;i<128;i++) f[i] = 0;
   for(i=0;;i++) {
      if(str[i]=='\0') break;
      c[str[i]] = str[i];
      f[str[i]]++;
   }
   *nc = 128;
   for(i=0;i<*nc;) {
      if(f[i]==0) {
         *nc = *nc-1;
         for(j=i;j<*nc;j++) {
            f[j] = f[j+1];
            c[j] = c[j+1];
         }
      }
      else i++;
   }
}


int main(int argc, char *arg[]) {
   if(argc!=2) {
      cout<<"usage: "<<arg[0]<<" [input_file]"<<endl;
      return 0;
   }
   FILE *fi;
   char str[1000];
   int i;
   if((fi = fopen(arg[1],"r"))==NULL) {
      cout<<"Cannot open file "<<arg[1]<<endl;
      return -1;
   }
   fgets(str, 1000, fi);
   int ns = strlen(str)-1;
   str[ns]='\0';
   fclose(fi);

   char *c = (char *) malloc (128*sizeof(char));
   unsigned *f = (unsigned *) malloc (128*sizeof(unsigned)), nc;
   cout<<"Input string: "<<str<<endl;
   cout<<"Length of string: "<<ns<<endl;
   count_freq(str, c, f, &nc); // count for the symbols and their frequencies
   cout<<nc<<" symbols in total."<<endl;
   Node *mh = HuffmanTree(c, f, nc); // Constructing the Huffman encoding tree.
   int is, ii, **cp = CharIndex(mh); // extract the Huffman code of each symbol from the constructed Huffman tree.
   cout<<"The input symbols and their frequencies: ";
   for(i=0;i<nc;i++) cout<<c[i]<<" "<<f[i]<<"   "; cout<<endl;
   double H=0.;
   for(i=0;i<nc;i++) H = H - ((double)f[i]/ns)*log((double)f[i]/ns)/log(2.); // computing the theoretical lower limit
   cout<<"The entropy of the frequency of the symbols in the string: "<<H<<endl;
   cout<<"Thus the theoretical lower limit for the number of bits required is given by H * ns =  "<<H*ns<<endl;
   PrintTree(mh);
   for(ii=0,is=0;is<ns;is++) for(i=0;;i++) { if(cp[str[is]][i]<0) break; cout<<cp[str[is]][i]; ii++;} cout<<endl;
   cout<<"number of bits in the encoded string: "<<ii<<endl;

   return 1;
}
