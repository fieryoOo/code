#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <ftw.h>
#include <unistd.h>
//#include <string.h>
//#include <iostream>
#include <errno.h>
#include <sys/time.h>
//using namespace std;

void TimedContinue (int time) {
   fd_set fdset;
   struct timeval timeout;
   int rc;
   timeout.tv_sec = time;
   timeout.tv_usec = 0;
   FD_ZERO(&fdset);
   FD_SET(0, &fdset);
   //cerr<<"### Sure to continue (continue in "<<time<<" sec if no input received)?  "<<endl;
   rc = select(1, &fdset, NULL, NULL, &timeout);
   if(rc == -1) {
      //cerr<<"### Error: select failed! ###"<<endl;
      exit(0);
   }
   else if (rc && FD_ISSET(0, &fdset)) {
      char cin, buff[300];
      fgets(buff, 300, stdin);
      sscanf(buff, "%c", &cin);
      if( cin!='Y' && cin!='y' ) exit(0);
   }
}

void fMove (char *oldname, char *newname) {
   if( rename(oldname, newname) == 0 ) return; //succeed
   int errsv = errno;
perror(newname); printf("errno: %d\n", errsv);

   if( errsv == 2 ) return;
   if( errsv == 21 || errsv == 39 ) {
      char crctname[150];
      sprintf(crctname, "%s/%s", newname, oldname);
      fMove(oldname, crctname); return;
   }
   perror("### Warning: Moving failed");
exit(0);
//   cerr<<"### Warning: failed moving file "<<oldname<<" ###"<<endl;
   TimedContinue(10);
}

int main() {
   int numCPU = sysconf( _SC_NPROCESSORS_ONLN );
   printf("numCPU: %d\n", numCPU);
   fMove("temp1", "temp2");

   fMove("temp2", "temp3");

   return 1;
}
