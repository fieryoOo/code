#include <stdio.h>
#include <stdlib.h>
#include <ftw.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <iostream>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <math.h>
using namespace std;

/* ---------------------- Move file or directory ----------------------------- */
#define PLENMAX 150
void Move (char *oldname, char *newname) {
   if( rename(oldname, newname) == 0 ) return; //succeed
   int errsv = errno;
   if( errsv == 2 ) return; // old file not exist
   if( errsv == 21 || errsv == 39 ) { // newfile is a directory
      perror(newname);
      char crctname[PLENMAX];
      sprintf(crctname, "%s/%s", newname, oldname);
      Move(oldname, crctname); return;
   }
   perror("### Warning: Moving failed");
   //TimedContinue(10);
}

/* ------------------------ Listing (wildcards matching) ----------------------------- */
#include <sys/types.h>
#include <fts.h>
#include <fnmatch.h>

#define BLKSIZE 1024
//int namecmp(const FTSENT **f1, const FTSENT **f2) { return strcmp((*f1)->fts_name, (*f2)->fts_name); }
char * List(char *dir, const char *pattern, int type, int *nfile) {
   /* type value decides how sub-directories are handdled
   0: list files int the root dir only
   1: list files and dir names in the root dir
   2: list all files
   3: list all files and directories */
   if( type>3 || type<0 ) {
      cerr<<"Unknow list type: "<<type<<endl;
      return NULL;
   }
   FTS *tree;
   FTSENT *file;
   char *dirlist[] = { dir, NULL }; //may send in multiple dirs
   //get handle of the file hierarchy; FTS_LOGICAL follows symbolic links and detects cycles.
   //replace '0' with 'namecmp' to sort files by name
   tree = fts_open(dirlist, FTS_LOGICAL | FTS_NOSTAT | FTS_NOCHDIR, 0);
   if (tree == NULL) perror("fts_open");

   char *sblk = NULL;
   int sleng = 0, bsize = 0;
   int outflag = 1; // if listing within current directory
   if( type<2 && strcmp(dir, ".")==0 ) outflag=0; // path will not be printed
   *nfile = 0;
   //ignores '.' and '..' as FTS_SEEDOT is not set
   while ((file = fts_read(tree))) {
      switch (file->fts_info) { //current node
         case FTS_DNR: // is a non-readable dir
         case FTS_ERR: // has common errors
         case FTS_NS: // has no stat info
         case FTS_DC: // causes cycle
            perror(file->fts_path);
         case FTS_DP: // is a post-order dir
            continue; //skip all above cases

         case FTS_D: // is a directory
            if(file->fts_level>0) switch(type) {
                case 0:
                   fts_set(tree, file, FTS_SKIP); //no descend
                   continue; // and skip
                case 1:
                   fts_set(tree, file, FTS_SKIP); //no descend
                   break; // and stop switch
                case 2:
                   continue; //skip directories
                case 3:;
            }
      }

      if (fnmatch(pattern, file->fts_name, FNM_PERIOD) == 0) {
         if( sleng > bsize-PLENMAX ) {
            bsize += BLKSIZE;
            sblk = (char *) realloc (sblk, bsize * sizeof(char));
         }
         if(outflag) sleng += sprintf(&sblk[sleng], "%s\n", file->fts_path);
         else sleng += sprintf(&sblk[sleng], "%s\n", file->fts_name);
	 *nfile = *nfile+1;
      }
   }

   if (errno != 0) perror("fts_read");
   if (fts_close(tree) < 0) perror("fts_close");
   return sblk;
}
// Wildcards Moving
char * wMove (char *pattern, char *tdir, int retlst, int *nfile) {
  // return a list of moved files when retlst==1
   //list files matching pattern
   char *list = List((char *)".", pattern, 0, nfile);
   if(list==NULL) return NULL;
   //move and rename;
   char *filelst = NULL;
   if(retlst) filelst = new char[strlen(list) + *nfile*(strlen(tdir)+3)];
   char list_name[PLENMAX], tname[PLENMAX];
   int offset, curp = 0, sleng = 0;
   while( (sscanf(&list[curp], "%s%n", list_name, &offset)) == 1 ) {
      sprintf(tname, "%s/%s", tdir, list_name);
      Move(list_name, tname);
      if(retlst) sleng += sprintf(&filelst[sleng], "%s\n", tname);
      curp += offset;
   }
   free(list);
   return filelst;
}

int main () {
   int nfile;
   wMove("*.tmp", "testwMove", 0, &nfile);
   char *list = List("testwMove","*.tmp", 3, &nfile);
   cout<<list<<endl;
   wMove("testwMove/*.tmp", ".", 0, &nfile);
   list = List(".","*.tmp", 3, &nfile);
   cout<<list<<endl;
   free(list);
   return 1;
}
