#include <iostream>

using namespace std;

int main(int na, char *arg[])
{ 
  if(na!=2)
    {
      cout<<"usage:touch_sac filelst"<<endl;
      return 0;
    }
  FILE *f1,*f2;
  char file_name[20000];
  if((f1=fopen(arg[1],"r"))==NULL)
    {
      cout<<"no filelst exist!!!"<<endl;
      return 1;
    }
  for(;;)
    {
      if(fscanf(f1,"%s",file_name)==EOF) break;
      printf("Touching file: %s...\n",file_name);
      f2=fopen("runsac.csh","w");
      fprintf(f2,"sac <<END\n");
      fprintf(f2,"r %s\n",file_name);
      fprintf(f2,"w %s\n",file_name);
      fprintf(f2,"quit\n");
      fprintf(f2,"END\n\n");
      fclose(f2);
      system("csh runsac.csh");
    }
  fclose(f1);
  system("rm -f runsac.csh");
  return 0;
}
