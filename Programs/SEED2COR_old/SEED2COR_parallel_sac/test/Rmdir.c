#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <ftw.h>
#include <unistd.h>
#include <errno.h>

int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
{
printf("%s\n", fpath);
    int rv = remove(fpath);
    if (rv)
        perror(fpath);

    return rv;
}

int rmrf(char *path)
{
    return nftw(path, unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
}

int main() {
   printf("%d\n", rmrf("temp_dir"));
   int rs = remove(NULL), ensv = errno;
   printf("NULLrt er: %d %d\n", rs, ensv);
   printf("EFAULT: %d\n", EFAULT);
   perror("Nuller");
   rs = remove("temp1"); ensv = errno;
   printf("temp1rt: %d %d\n", rs, ensv);
   perror("temp1er");
   printf("ENOENT: %d\n", ENOENT);
   return 1;
}
