#include <pthread.h>

extern   pthread_attr_t attr_j;
extern   pthread_mutex_t cevlock, fftlock, fiolock;

void InitialPthread() {
   size_t stacksize = 67108864;
   pthread_attr_init(&attr_j);
   pthread_attr_setdetachstate(&attr_j, PTHREAD_CREATE_JOINABLE);
   pthread_attr_setstacksize (&attr_j, stacksize);
   pthread_mutex_init(&cevlock, NULL);
   pthread_mutex_init(&fftlock, NULL);
   pthread_mutex_init(&fiolock, NULL);
}

void CleanupPthread() {
   pthread_attr_destroy(&attr_j);
   pthread_mutex_destroy(&cevlock);
   pthread_mutex_destroy(&fftlock);
   pthread_mutex_destroy(&fiolock);
}
