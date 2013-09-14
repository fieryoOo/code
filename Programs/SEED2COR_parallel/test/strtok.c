/* strtok example */
#include <stdio.h>
#include <string.h>
#include <libgen.h>

int main ()
{
  char str[] =". aab ,aa a-ccc/";
  char * pch;
  printf ("Splitting string \"%s\" into tokens:\n",str);
  pch = strtok (str," ,.-");
  while (pch != NULL)
  {
    printf ("%s\t\t%s\n",pch,str);
    pch = strtok (NULL, " ,.-");
  }

  char *dirc, *basec, *bname, *dname;
  char *path = "/etc/passwd/";
  dirc = strdup(path);
  basec = strdup(path);
  dname = dirname(dirc);
  bname = basename(basec);
  printf("path=%s, dirname=%s, basename=%s\n", path, dname, bname);
  return 0;
}
