#define NFCMAX 100000
char *buffdel, *buffend, *taildel;
pthread_mutex_t dellock, dexlock;

void DeleteInitial () {
   //lock(dellock);
   buffdel = new char[NFCMAX];
   buffend = buffdel + (NFCMAX-300)*sizeof(char);
   taildel = buffdel;
   //unlock(dellock);
}

void Delete (char *fname) {
   //lock(dellock);
   taildel += sprintf(taildel, "%s\n", fname);
   if( taildel > buffend ) DeleteExcute(0);
   else //unlock(dellock);
}

void DeleteExcute (int flag) { //flag=1 if called from outside
   if(flag) //lock(dellock);
   FILE *fdel = fopen("fdel_lst.tmp", "w");
   fprintf(fdel, "%s fdel_lst.tmp\n", buffdel);
   fclose(fdel);
   taildel = buffdel;
   //unlock(dellock);
   //lock(dexlock);
   system("more fdel_lst.tmp | xargs rm -f");
   //unlock(dexlock);
}

void DeleteCleanup () {
   DeleteExcute(1);
   //lock(dellock);
   delete [] buffdel; buffdel = NULL;
   taildel = NULL;
   //unlock(dellock);
}
