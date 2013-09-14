#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;

main(int na, char *arg[])
{
  if (na!=7)
    {
      cout<<"usage: Requesting_ASN [sta_list (as 'M12A TA')] [month_list (as '2008 6')] [label] [email] [name] [channel list]"<<endl;
      return 0;
    }
  FILE *f1,*f2,*fout;
  int i, nch, year, year2, month, month2, day, day2, month_dn, day_num[12];
  char month_name[12][4], intemp[100], station[6], net[3];
  char buff[300], tp, ctmp, ch[10][5];

  sprintf(month_name[0], "JAN\0"); day_num[0] = 31;
  sprintf(month_name[1], "FEB\0"); day_num[1] = 28;
  sprintf(month_name[2], "MAR\0"); day_num[2] = 31;
  sprintf(month_name[3], "APR\0"); day_num[3] = 30;
  sprintf(month_name[4], "MAY\0"); day_num[4] = 31;
  sprintf(month_name[5], "JUN\0"); day_num[5] = 30;
  sprintf(month_name[6], "JUL\0"); day_num[6] = 31;
  sprintf(month_name[7], "AUG\0"); day_num[7] = 31;
  sprintf(month_name[8], "SEP\0"); day_num[8] = 30;
  sprintf(month_name[9], "OCT\0"); day_num[9] = 31;
  sprintf(month_name[10], "NOV\0"); day_num[10] = 30;
  sprintf(month_name[11], "DEC\0"); day_num[11] = 31;

  if(!(fout=fopen(arg[6],"r"))){
     cout<<"Can't access channel.lst file "<<arg[6]<<endl;
     return 0;
  }
  cout<<"Channels: "<<endl;
  for(nch=0; fgets(buff, 300, fout)!=NULL; nch++) {
     sscanf(buff, "%s", ch[nch]);
     cout<<ch[nch]<<endl;
  }
  fclose(fout);
  cout<<"Coninue? ";
  scanf("%c", &ctmp);
  if(ctmp!='y' && ctmp!='Y') exit(0);
  
  if(!(f1=fopen(arg[2],"r"))){
     cout<<"Can't open month.lst file "<<arg[2]<<" to read."<<endl;
     return 0;
    }
  if(!(f2=fopen(arg[1],"r"))){
     cout<<"Can't open sta.lst file "<<arg[1]<<" to read."<<endl;
     return 0;
    }
  for(;;){
     if(fscanf(f1,"%d %d",&year,&month)!=2)break;
     month = month-1;
     month_dn = day_num[month];
     if( month == 1)if((year%4==0 && year%100!=0) || year%400==0)month_dn=29;
     if ( month >= 0 && month < 11 ){
        year2=year;
        month2=month+1;
       }
     else if ( month == 11 ){
        year2=year+1;
        month2=0;
       }
     else {
        cout<<"Wrong month info in the month.lst: "<<month+1<<"th month!"<<endl;
        continue;
       }
     for(day=1;day<=month_dn;day++){
        //if(day!=31) continue;
        if(!(fout=fopen("requesting_email","w"))){
           cout<<"Can't open month.lst file "<<arg[2]<<" to read."<<endl;
           return 0;
          }
        fprintf(fout,".NAME %s\n.INST CU\n.MAIL University of Colorado at Boulder\n.EMAIL %s\n.PHONE\n.FAX\n.MEDIA: Electronic (FTP)\n.ALTERNATE MEDIA: Electronic (FTP)\n.LABEL %s_%d.%s.%d\n.END\n", arg[5], arg[4], arg[3], year, month_name[month],day);
        if ( day < month_dn ){
           day2=day+1;
           sprintf(intemp,"%d %d %d 0 0 0 %d %d %d 0 0 0 1\0", year, month+1, day, year, month+1, day2);
          }
        else
           sprintf(intemp,"%d %d %d 0 0 0 %d %d 1 0 0 0 1\0", year, month+1, day, year2, month2+1);
        rewind(f2);
        for(;;){
           if(fscanf(f2,"%s %s", &station, &net)!=2)break;
           /*
           sscanf(arg[6], "%c%c%c", &tp, &ctmp, &ch);
           if(tp!='L'&&tp!='B'||ctmp!='H'||ch!='Z'&&ch!='E'&&ch!='N'&&ch!='A') {
              cout<<"Channel input error: "<<arg[6]<<endl;
              exit(0);
           }
           if(ch == 'A') {
              fprintf(fout,"%s %s %s %cHZ\n", station, net, intemp, tp);
              fprintf(fout,"%s %s %s %cHE\n", station, net, intemp, tp);
              fprintf(fout,"%s %s %s %cHN\n", station, net, intemp, tp);
           }
           else fprintf(fout,"%s %s %s %s\n", station, net, intemp, arg[6]);
           */
           for(i=0;i<nch;i++) fprintf(fout,"%s %s %s %s\n", station, net, intemp, ch[i]);
        }
        fclose(fout);
        system("cat requesting_email | mail -s 'Requesting Data' breq_fast@iris.washington.edu");
        cout<<"Request of "<<month_name[month]<<" "<<day<<"th "<<year<<" sent..."<<endl;
       }
    }
  fclose(f1);
  fclose(f2);
}
