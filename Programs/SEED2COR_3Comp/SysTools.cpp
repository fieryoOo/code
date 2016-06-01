#include "SysTools.h"
//#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <stdlib.h>
#include <ftw.h>
#include <unistd.h>
#include <errno.h>
#include <cstring>
#include <iostream>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <math.h>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

/* -------------- Read memory info from /proc/meminfo -------------------- */
struct SYSINFO {
   long MemTotal;
   long MemFree;
   long SwapFree;
   long Cached;
};

int GetSysinfo (long &MemAvail) {
   struct sysinfo SysInfoSTD;
   sysinfo(&SysInfoSTD);
   MemAvail = SysInfoSTD.freeram;

   FILE *fmem;
   if( ! (fmem=fopen("/proc/meminfo", "r")) ) return 0;

   char buff[300], *phead;
   struct SYSINFO SysInfo;
   SysInfo.MemTotal = -1;
   SysInfo.MemFree = -1;
   SysInfo.SwapFree = -1;
   SysInfo.Cached = -1;
   while( fgets(buff, 300, fmem) ) {
      if( SysInfo.MemTotal==-1 && (phead=strstr(buff, "MemTotal")) ) {
	 phead += strlen("MemTotal ");
	 sscanf(phead, "%ld", &(SysInfo.MemTotal));
	 SysInfo.MemTotal *= 1024;
	 if( fabs(SysInfoSTD.totalram-SysInfo.MemTotal)/double(SysInfoSTD.totalram) > 0.01 ) return 0;
      }
      else if ( SysInfo.MemFree==-1 && (phead=strstr(buff, "MemFree")) ) {
	 phead += strlen("MemFree ");
	 sscanf(phead, "%ld", &(SysInfo.MemFree));
	 SysInfo.MemFree *= 1024;
      }
      else if ( SysInfo.SwapFree==-1 && (phead=strstr(buff, "SwapFree")) ) {
         phead += strlen("SwapFree ");
         sscanf(phead, "%ld", &(SysInfo.SwapFree));
         SysInfo.SwapFree *= 1024;
      }
      else if ( SysInfo.Cached==-1 && (phead=strstr(buff, "Cached")) ) {
         phead += strlen("Cached ");
         sscanf(phead, "%ld", &(SysInfo.Cached));
         SysInfo.Cached *= 1024;
      }
   }
   fclose(fmem);

   if( SysInfo.MemFree==-1 || SysInfo.Cached==-1 ) return 0;
   MemAvail = SysInfo.MemFree + SysInfo.Cached;
   return 1;
}

void TimedContinue (int time) {
   fd_set fdset;
   struct timeval timeout;
   int rc;
   timeout.tv_sec = time;
   timeout.tv_usec = 0;
   FD_ZERO(&fdset);
   FD_SET(0, &fdset);
   cout<<"### Sure to continue (continue in "<<time<<" sec if no input received)?  "<<endl;
   rc = select(1, &fdset, NULL, NULL, &timeout);
   if(rc == -1) {
      cerr<<"Error(TimedContinue): select failed!"<<endl;
      exit(0);
   }
   else if (rc && FD_ISSET(0, &fdset)) {
      char cin, buff[300];
      fgets(buff, 300, stdin);
      sscanf(buff, "%c", &cin);
      if( cin!='Y' && cin!='y' ) exit(0);
   }
}

void EstimateMemAvail (long &MemAvail) {
   //Estimating succeed
   if( GetSysinfo(MemAvail) ) return;
   //else
   cout<<"### Warning: Failed to get cached memory size from /proc/meminfo!! ###"<<endl;
   cout<<"### The current freeRAM ("<<MemAvail/1024./1024.<<" Mb) will be used ###"<<endl;
   TimedContinue(10); 
}

/* ------------------------- make dir ----------------------------------------------- */
#include <sys/stat.h>
bool MKDir(const char *dirname) {
  //create dir if not exists
  //with read/write/search permissions for owner and group, and with read/search permissions for others if not already exists
	if( mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0 ) return true;
	switch(errno) {
		case EEXIST:
			return false;
		default:
			perror("### Error: MKDir failed"); //failed. prompt to continue
			exit(0);
   }
}

bool MKDirs(const char *dirname) {
	std::stringstream sin(dirname);
	std::string pathname = "";
	bool succeed = false;
	for( std::string dircur; std::getline(sin, dircur, '/'); ) {
		pathname += dircur + "/";
		bool suc = MKDir( pathname.c_str() );
		succeed = succeed || suc;
	}
	return succeed;
}

/* --------------------- Delete file or directory ---------------------------- */
//#define _XOPEN_SOURCE 500
//#include <ftw.h>
//#include <unistd.h>
int fdprompt = -10;
void fRemove (const char *fname) {
   //cerr<<"Removing "<<fname<<endl;
   if( remove( fname ) == 0 ) return; //succeed
   int ersv = errno;
   if( ersv == ENOENT ) return; //file not exists
   if( fdprompt>0 ) return; fdprompt++;
   perror("### Warning: Deleting failed"); //failed. prompt to continue
   TimedContinue(10);
}

int Unlink(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf) {
   fRemove((char *)fpath);
   return 0;
}
int dRemove(const char *dirname) {
   return nftw(dirname, Unlink, 64, FTW_DEPTH | FTW_PHYS);
}

/* ---------------------- Move file or directory ----------------------------- */
#define PLENMAX 200
void Move (const char *oldname, const char *newname) {
   if( rename(oldname, newname) == 0 ) return; //succeed
   int errsv = errno;
   if( errsv == 2 ) return; // old file not exist
   if( errsv == 21 || errsv == 39 ) { // newfile is a directory
		std::string strtmp = "SysTools::Move: "+std::string(newname)+" is a directory";
      perror( strtmp.c_str() );
      char crctname[PLENMAX];
      sprintf(crctname, "%s/%s", newname, oldname);
      Move(oldname, crctname); return;
   }
   perror("### Warning: Moving failed");
   TimedContinue(10);
}

// Wildcards Moving
#include <libgen.h>
//char * List(const char *dir, const char *pattern, int type, int *nfile);
bool List(const char *dir, const char *pattern, int type, std::vector<std::string> &filelist);
bool wMove (const char *odir, const char *pattern, const char *tdir, std::vector<std::string> &outlist) {
  // return a list of moved files when retlst==1
   //list files matching pattern
   //char *list = List(odir, pattern, 0, nfile);
   //if(list==NULL) return NULL;
   std::vector<std::string> list;
	if( ! List(odir, pattern, 0, list) ) return false;
   //if( ! List(odir, pattern, 0, list) ) return false;
   //move and rename;
   char *filelst = NULL;
   //if( retlst == 1 ) filelst = new char[strlen(list) + *nfile*(strlen(tdir)+3)];
   char list_name[PLENMAX], tname[PLENMAX], *path, *bname;
   int offset, curp = 0, sleng = 0;
   //while( (sscanf(&list[curp], "%s%n", list_name, &offset)) == 1 ) {
   for(int i=0; i<list.size(); i++) {
      path = strdup(list.at(i).c_str()); bname = basename(path);
      sprintf(tname, "%s/%s", tdir, bname);
      Move(list.at(i).c_str(), tname);
      //if( retlst == 1 ) sleng += sprintf(&filelst[sleng], "%s\n", tname);
      outlist.push_back(tname);
      curp += offset;
   }
   //free(list);
   return true;
}

/* ------------------------ copy a file ------------------------------- */
#define BSZ 8192
void Copy(const char *oldname, const char *newname) {
   FILE *fin, *fou;
   if( (fin = fopen(oldname, "r")) == NULL ) {
      perror("### Warning from Copy (fopen)");
      TimedContinue(10);
      return;
   }
   if( (fou = fopen(newname, "w")) == NULL ) {
      perror("### Warning from Copy (fopen)");
      fclose(fin);
      TimedContinue(10);
      return;
   }
   int result;
   char buff[BSZ];
   while ( (result = fread(buff, 1, BSZ, fin)) ) {
       fwrite(buff, 1, result, fou);
   }
   fclose(fin); fclose(fou);
}

/* ------------------------ Listing (wildcards matching) ----------------------------- */
#include <sys/types.h>
#include <fts.h>
#include <fnmatch.h>
//#define BLKSIZE 1024
//int namecmp(const FTSENT **f1, const FTSENT **f2) { return strcmp((*f1)->fts_name, (*f2)->fts_name); }
bool List(const char *dir, const char *pattern, int type, std::vector<std::string> &filelist) {
   /* type value decides how sub-directories are handdled
   0: list files in the root dir only
   1: list files in the root dir with dir paths
   2: list all file names
   3: list all files with dir paths */
   //*nfile = 0;
   if( type>3 || type<0 ) {
      cerr<<"ERROR(List): Unknow list type: "<<type<<endl;
      return false;
   }
   FTS *tree;
   FTSENT *file;
   char *dirlist[] = { (char *)dir, NULL }; //may send in multiple dirs
   //get handle of the file hierarchy; FTS_LOGICAL follows symbolic links and detects cycles.
   //replace '0' with 'namecmp' to sort files by name
   tree = fts_open(dirlist, FTS_LOGICAL | FTS_NOSTAT | FTS_NOCHDIR, 0);
   if (tree == NULL) perror("fts_open");

   //char *sblk = NULL;
   //int sleng = 0, bsize = 0;
   //empty the input filelist
   filelist.clear();
   int outflag = 1; // if listing within current directory
   if( type<2 && (strcmp(dir, ".")==0||strcmp(dir, "./")==0) ) outflag=0; // path will not be printed
   //ignores '.' and '..' as FTS_SEEDOT is not set
	while ((file = fts_read(tree))) {
		std::string strtmp;
		switch (file->fts_info) { //current node
			case FTS_DNR: // is a non-readable dir
			case FTS_ERR: // has common errors
			case FTS_NS: // has no stat info
			case FTS_DC: // causes cycle
				//strtmp = "SysTools::List: "+std::string(file->fts_path)+" causes cycle!";
				//perror( strtmp.c_str() );
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

		if( strcmp(file->fts_name,"") == 0 ) continue;	// skip empty name
		if (fnmatch(pattern, file->fts_name, FNM_PERIOD) == 0) {
			/*
				if( sleng > bsize-PLENMAX ) {
            bsize += BLKSIZE;
            sblk = (char *) realloc (sblk, bsize * sizeof(char));
         }
         if(outflag) sleng += sprintf(&sblk[sleng], "%s\n", file->fts_path);
         else sleng += sprintf(&sblk[sleng], "%s\n", file->fts_name);
			*/
			if(outflag) {
				filelist.push_back(file->fts_path);
			} else if( strcmp(file->fts_name, file->fts_path) != 0 ) {
				filelist.push_back(file->fts_name);
			}
			//*nfile = *nfile+1;
      }
   }

   if (errno != 0) perror("fts_read");
   if (fts_close(tree) < 0) perror("fts_close");
	return filelist.size()!=0;
}

