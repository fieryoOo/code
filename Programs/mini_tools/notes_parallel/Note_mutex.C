#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <pthread.h>
using namespace std;

#define NT 4
#define NDAT 10000000

struct DAT {
   float sig[NDAT];
   pthread_mutex_t lock[NDAT];
};
struct DAT data;

void *calc (void *x) {
   int i;
   //update signal array with locking and unlocking mutex(lock) array
   for(i=0;i<NDAT;i++) {
      pthread_mutex_lock(&(data.lock[i]));
      data.sig[i]++;
      pthread_mutex_unlock(&(data.lock[i]));
   }
   int tid = *((int *) x);
   printf("Calculation completed from thread %d!\n", tid);
   //return NULL;
   pthread_exit(x);
}


/* stack management
pthread_attr_getstacksize (attr, stacksize)
pthread_attr_setstacksize (attr, stacksize)
pthread_attr_getstackaddr (attr, stackaddr)
pthread_attr_setstackaddr (attr, stackaddr) 
*/
int main ()
{
   int i;

   //Initialize signal and mutex(lock) array
   for(i=0; i<NDAT; i++) {
      data.sig[i] = 0.;
      pthread_mutex_init(&(data.lock[i]), NULL);
   }
   
   //thread ids
   pthread_t threads[NT];
   //thread attributes
   pthread_attr_t attr_j, attr_d;
   //status attribute
   void *status;

   //attributes initialize and state set
   pthread_attr_init(&attr_j);
   pthread_attr_init(&attr_d);
   pthread_attr_setdetachstate(&attr_j, PTHREAD_CREATE_JOINABLE);
   pthread_attr_setdetachstate(&attr_d, PTHREAD_CREATE_DETACHED);
   int thread_args[NT];
   for(i=0; i<NT; i++) {
      thread_args[i] = i;
      printf("spawning thread %d\n", i);
      //create
      if(pthread_create(&threads[i], &attr_j, calc, (void *) &thread_args[i])) cout<<"Create failed!"<<endl;
   }
   //free attributes
   pthread_attr_destroy(&attr_j);
   pthread_attr_destroy(&attr_d);
   //free mutex array
   for(i=0; i<NT; i++) pthread_mutex_destroy(&data.lock[i]);

   for (i=0; i<NT; ++i) {
      if(pthread_join(threads[i], &status)) cout<<"Join failed!"<<endl;
      else printf("Main: completed join with thread %d having a status of %d\n", i, *((int *)status));
   }

   for(i=0;i<NDAT;i++) if(data.sig[i]!=NT) cout<<i<<" "<<data.sig[i]<<endl;
   //block main and keep it alive to support created threads until they are done
   pthread_exit(NULL);
}
