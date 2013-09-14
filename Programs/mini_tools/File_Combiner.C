#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

main(int na, char *arg[])
{
  if (na!=4)
    {
      cout<<"usage: File_Combiner [in_file_lst] [column # (0 for entire row; -1 for non-comp join)] [out_file]"<<endl;
      return 0;
    }
  FILE *f1, *f2, *fout;
  char *filename, *tmp, *out, temp[100], str[1000000][100], buff[1000000][100], num_buff[33];
  int i, j, ifrow, irow=0, flag[1000000];

  if(atoi(arg[2]) != atof(arg[2]) || atoi(arg[2])<-1){
     cout<<"Integer(>=-1) expected for column #"<<endl;
     return 0;
    } 
  if(!(f1=fopen(arg[1],"r"))){
     cout<<"Can't open file "<<arg[1]<<" to read!"<<endl;
     return 0;
    }
  for(;;){
     if(fgets(temp,100,f1) == NULL)break;
     //itoa (ifl, num_buff, 10);
     strtok(temp, "\n");
     tmp = strdup(temp);
     filename = strsep(&tmp," ");
     cout<<"Reading "<<filename<<endl;
     if(!(f2=fopen(filename,"r"))){
        cout<<"Can't open file "<<filename<<" to read. Skipped!"<<endl;
        continue;
       }
     ifrow=0;
     for(;;irow++){
        if(fgets(buff[irow],100,f2) == NULL)break;
        strtok(buff[irow],"\n");
        if(atoi(arg[2])==0){strcpy(str[irow],buff[irow]);cout<<str[irow]<<endl;}
        else if(atoi(arg[2])==-1) {sprintf(str[irow],"%d",ifrow);}
        else {
           tmp = strdup(buff[irow]);
           out = strtok(tmp, " ");
           if(atoi(arg[2])>1)
              for (i=2;out=strtok(NULL, " ");i++)
                 if(i == atoi(arg[2]))break;
           if(out != NULL)strcpy(str[irow],out);
           else {cout<<"Column # out of range!"<<endl; return 0;}
          }
          ifrow++;
       }
     fclose(f2);
    }
  fclose(f1);
  
  if(!(fout=fopen(arg[3],"w"))){
     cout<<"Can't open file "<<arg[3]<<" to write!"<<endl;
     return 0;
    }
  for(i=0;i<irow;i++)flag[i]=0;
  for(i=0;i<irow;i++){
     if(flag[i])continue;
     fprintf(fout,"%s",buff[i]);
     for(j=i+1;j<irow;j++)
        if(strcmp(str[i],str[j])==0){
           fprintf(fout,"\t%s",buff[j]);
           flag[j]=1;
          }
     fprintf(fout,"\n");
    }
  fclose(fout);
}
