#include "Param.h"

int NTHRDS;
long MemAvail;

void EstimateMemAvail (long *MemAvail);

void InitialPthread() {
   size_t stacksize = 2097152;
   NTHRDS = sysconf( _SC_NPROCESSORS_ONLN );
   NTHRDS += NTHRDS/10;
   pthread_attr_init(&attr_j);
   pthread_attr_setdetachstate(&attr_j, PTHREAD_CREATE_JOINABLE);
   pthread_attr_setstacksize (&attr_j, stacksize);
   pthread_mutex_init(&cevlock, NULL);
   pthread_mutex_init(&fftlock, NULL);
   pthread_mutex_init(&fiolock, NULL);
   EstimateMemAvail (&MemAvail);
}

void CleanupPthread() {
   pthread_attr_destroy(&attr_j);
   pthread_mutex_destroy(&cevlock);
   pthread_mutex_destroy(&fftlock);
   pthread_mutex_destroy(&fiolock);
}
