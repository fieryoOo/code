#include<stdio.h>
#include<stdlib.h>
#include<iostream>
using namespace std;

#define SIZE 10
#define T 1
#define F 0
struct edge
{
 int terminal;
 struct edge *next;
};
struct vertex
{
 int visit;
 int vertex_no;
 char info;
 int path_lenth;
 struct edge *edge_ptr;
};
struct q
{
 int info;
 struct q *next;
};
void table(int,int matrix[SIZE][SIZE],struct vertex vert[SIZE]);
struct edge *insert_vertex(int,struct edge *);
void bfs(int ,int *dist,struct vertex vert[SIZE]);
void output(int,int a[SIZE][SIZE]);
void input(int,int a[SIZE][SIZE]);
struct q *insert_queue(int vertex_no,struct q *first);
struct q *delet_queue(int *vertex_no,struct q *first);
struct edge *insert_vertex(int vertex_no,struct edge *first)
{
 struct edge *new1,*current ;
 new1=(struct edge*)malloc(sizeof(struct edge));
 new1->terminal=vertex_no;
 new1->next=NULL;
 if(!first)
   return new1;
 for(current=first;current->next;current=current->next);
 current->next=new1;
 return first;
}
struct q *insert_queue(int vertex_no,struct q *first)
{
 struct q *new1,*current;
 new1=(struct q*)malloc(sizeof(struct q));
 new1->info=vertex_no;
 new1->next=NULL;
 if(!first)
   return new1;
 for(current=first;current->next;current=current->next);
 current->next=new1;
 return first;
}
struct q *delet_queue(int *vertex_no,struct q *first)
{
 struct q *previous;
 if(!first)
   return NULL;
 *vertex_no=first->info;
 first=first->next;
 free(previous);
 return first;
}
void table(int vertex_num,int matrix[SIZE][SIZE],struct vertex vert[SIZE])
{
 int i,j;
 for(i=0;i<vertex_num;i++)
 {
  vert[i].visit=F;
  vert[i].vertex_no=i+1;
  vert[i].info='A'+i;
  vert[i].path_lenth=0;
  vert[i].edge_ptr=NULL;
 }
 for(i=0;i<vertex_num;i++)
   for(j=0;j<vertex_num;j++)
    if(matrix[i][j]>0)
     vert[i].edge_ptr=insert_vertex(j,vert[i].edge_ptr);
}
void bfs(int index,int *dist,struct vertex vert[SIZE])
{
 struct q *queue=NULL;
 struct edge *link;
 vert[index].visit=T;
 queue=insert_queue(index,queue);
 while(queue)
 {
  queue=delet_queue(&index,queue);
  for(link=vert[index].edge_ptr;link;link=link->next)
   if(vert[link->terminal].visit==F)
   {
    vert[link->terminal].visit=T;
    vert[link->terminal].path_lenth=vert[index].path_lenth+1;
    queue=insert_queue(link->terminal,queue);
   }
  }
}
void input(int number,int a[SIZE][SIZE])
{
 int i,j;
 printf("input the adjacency matrix is:\n");
 for(i=0;i<number;i++)
 {
  for(j=0;j<number;j++)
   scanf("%d",&a[i][j]);
  printf("\n");
 }
}
void output(int number,int a[SIZE][SIZE])
{
 int i,j;
 printf("\n adjacency matrix is:\n");
 for(i=0;i<number;i++)
 {
  for(j=0;j<number;j++)
   printf("%d\t",a[i][j]);
  printf("\n");
 }
}
int main()
{
 int i,number,index,dist,a[SIZE][SIZE];
 struct vertex vert[SIZE];
 struct edge *list;
 system("clear");
 printf("\ninput the number of vertices in graph:");
 scanf("%d",&number);
 input(number,a);
 output(number,a);
 table(number,a,vert);
 printf("\ninput the strating vertex 0-%d:",number-1);
 scanf("%d",&index);
 dist=0;
 bfs(index,&dist,vert);
 printf("\n path lenth of the vertex %c",vert[index].info);
 printf("\n vertex length vertex complexity\n");
 for(i=0;i<number;i++)
 {
  printf("\n%c\t\t%d",vert[i].info,vert[i].path_lenth);
  for(list=vert[i].edge_ptr;list;list=list->next)
  {
    printf(" ");
    cout<<list->terminal+'A';
  }
 }
 //getch();
 //cout<<endl;
 return 1;
}

