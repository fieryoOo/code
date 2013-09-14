#include <stdio.h>
#include <sys/types.h>
#include <fts.h>
#include <fnmatch.h>
#include <errno.h>
#include <iostream>
using namespace std;

//int namecmp(const FTSENT **f1, const FTSENT **f2) { return strcmp((*f1)->fts_name, (*f2)->fts_name); }

#define BLKSIZE 1024
#define PLENMAX 150
char * List(char *dir, const char *pattern, int type) {
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
      }
   }
 
   if (errno != 0) perror("fts_read");
   if (fts_close(tree) < 0) perror("fts_close");
   return sblk;
}
 
int main(int argc, char *argv[])
{
   if( argc != 3 ) {
      cerr<<"Usage: "<<argv[0]<<" [list type (int)] [pattern]"<<endl;
      exit(-1);
   }
   char *list = List(".", argv[2], atoi(argv[1]));
   int offset, curp = 0;
   char buff[PLENMAX];
   while( (sscanf(&list[curp], "%s%n", buff, &offset)) == 1 ) {
      curp += offset;
      cout<<buff<<endl;
   }
   free(list);
   return 0;
}
