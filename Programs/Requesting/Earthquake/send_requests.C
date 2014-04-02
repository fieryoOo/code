#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void find_end (int addh,int addmin,int y, int m, int d,int h,int min,int *y2, int *m2, int *d2,int *h2,int *min2);
/*////////////////////////////////////////////////////////////////////////////////////////*/
/*----------------------------------------------------------------------------------------*/
int jday (int y, int m, int d)
/*----------------------------------------------------------------------------------------*/
{
  int jd = 0;
  int i;

  for ( i = 1; i < m; i++ )
    {
      if ( (i==1) || (i==3) || (i==5) || (i==7) || (i==8) || (i==10) ) jd += 31;
      else if (i==2)
	{
	  if ( y == 4*(y/4) ) jd += 29;
	  else jd += 28;
	}
      else jd += 30;
    }

  jd += d;

  return jd;
}

/*////////////////////////////////////////////////////////////////////////////////////////*/
/*----------------------------------------------------------------------------------------*/
void dm_fr_jday (int jd, int *y, int *m, int *d)
/*----------------------------------------------------------------------------------------*/
{
  int jdmax = 365, dy = 0, dy0, dmmax, i;

  if ( *y == 4*((*y)/4) ) jdmax = 366;

  if ( jd > jdmax)
    {
      (*y)++;
      jd -= jdmax;
    }

  for ( i = 1; i <= 12; i++ )
    {
      dy0 = dy;

        if ( (i==1) || (i==3) || (i==5) || (i==7) || (i==8) || (i==10)|| (i==12) ) dmmax = 31;
	else if (i==2)
	  {
	    if ( *y == 4*((*y)/4) ) dmmax = 29;
	    else dmmax = 28;
	  }
	else dmmax = 30;

	dy += dmmax;

	if ( jd <= dy )
	  {
	    *m = i;
	    *d = jd - dy0;
	    break;
	  }
    }
}

/*////////////////////////////////////////////////////////////////////////////////////////*/
char infline[20][200];
/*----------------------------------------------------------------------------------------*/
void one_event (int y, int m, int d, int h, int min,int sec,int ndays )
/*----------------------------------------------------------------------------------------*/
{
  int y2, m2, d2, h2,min2,sec2, jd;
  int y3, m3, d3, h3,min3,sec3;
  int nil, len, i;
  FILE *fst, *fre, *finf;
  static char sta[200];
  static char net[200];
  char buff1[200];
  finf = fopen("request.info","r");
  fst = fopen("station.loc","r");
  fre = fopen("request_mail","w");

  jd = jday( y, m, d) + ndays;


  find_end(0,-15, y,m,d,h,min, &y2,&m2,&d2,&h2,&min2);
  find_end(2, 30, y,m,d,h,min, &y3,&m3,&d3,&h3,&min3);

  printf(" this month %d\n", m);
  for ( nil = 0;; nil++ )
    {
      if ( !fgets(&(infline[nil][0]),200,finf) ) break;
    }

  if ( nil != 10 ) 
    {
      fprintf(stderr,"improper number of lines in request.info file\n");
      exit(0);
    }

  len = strlen(&(infline[nil-1][0]));
  sprintf(&(infline[nil-1][len-1]),"_%d_%d_%d_%d_%d\n", y, m, d, h,min );

  for ( i = 0; i < nil; i++ ) fputs(&(infline[i][0]),fre);
  fprintf(fre,".END\n");

  i = 0;
  for (;;)
    {
      i = i+1;
    
      //      if ( !fgets(sta,20,fst) ) break;
      // if ( strlen(sta) < 8 ) break;
      if(fscanf(fst,"%s %s",sta,net)==EOF) 
	break;
      sprintf(buff1,"%s %s %d %d %d %d %d %d.0 %d %d %d %d %d %d.0 1 LHE\n\0",sta, net,y2, m2, d2, h2,min2,sec, y3, m3, d3, h3,min3,sec );
      fputs(buff1,fre);
      sprintf(buff1,"%s %s %d %d %d %d %d %d.0 %d %d %d %d %d %d.0 1 LHN\n\0",sta, net,y2, m2, d2, h2,min2,sec, y3, m3, d3, h3,min3,sec );
      fputs(buff1,fre);
      sprintf(buff1,"%s %s %d %d %d %d %d %d.0 %d %d %d %d %d %d.0 1 LHZ\n\0",sta, net,y2, m2, d2, h2,min2,sec, y3, m3, d3, h3,min3,sec );
      fputs(buff1,fre);

      //      sprintf(&(sta[15])," %d %d %d %d %d %d.0 %d %d %d %d %d %d.0 1 LHN\n\0", y2, m2, d2, h2,min2,sec, y3, m3, d3, h3,min3,sec );
      //fputs(sta,fre);

      //      sprintf(&(sta[7])," %d %d %d %d %d %d.0 %d %d %d %d %d %d.0 1 BHN\n\0", y2, m2, d2, h2,min2,sec, y3, m3, d3, h3,min3,sec );
      //      fputs(sta,fre);

      //      sprintf(&(sta[7])," %d %d %d %d %d %d.0 %d %d %d %d %d %d.0 1 BHZ\n\0", y2, m2, d2, h2,min2,sec, y3, m3, d3, h3,min3,sec );
      //      fputs(sta,fre);

    }

  fclose(finf);
  fclose(fst);
  fclose(fre);
}



char buff[300];

/*////////////////////////////////////////////////////////////////////////////////////////*/
/*----------------------------------------------------------------------------------------*/
int main (int na, char *arg[])
/*----------------------------------------------------------------------------------------*/
{
  int y, m, d, h,min,sec, ndays;
  FILE *ff;

  if ( na < 2 ) ndays = 1;
  else sscanf(arg[1],"%d", &ndays );

  ff = fopen("events.loc","r");
  for(;;)
    {
      if ( !fgets(buff,300,ff) ) break;
      if ( strlen(buff) < 5 ) break;

      sscanf(&(buff[0]), "%4d", &y );
      sscanf(&(buff[4]), "%2d", &m );
      sscanf(&(buff[6]), "%2d", &d );
      sscanf(&(buff[8]), "%2d", &h );
      sscanf(&(buff[10]), "%2d", &min );
      sscanf(&(buff[12]), "%2d", &sec );
      one_event ( y, m, d, h, min,sec, ndays);
      system("mail  BREQ_FAST@iris.washington.edu < request_mail");  
    }

  fclose(ff);
}






/*////////////////////////////////////////////////////////////////////////////////////////*/
/*----------------------------------------------------------------------------------------*/
void find_end (int addh,int addmin,int y, int m, int d,int h,int min,int *y2, int *m2, int *d2,int *h2,int *min2)
/*----------------------------------------------------------------------------------------*/
{
  
  *min2 = min+addmin;
  *h2 = h + addh;
  *d2 = d;
  *m2 = m;
  *y2 = y;
  printf(" initial %d\n",m);

  if( addmin > 0 )
    {
       if( *min2 >= 60 ) 
       {
         *min2 = *min2 - 60;
        *h2 = *h2 + 1;
       }

       if( *h2 >= 24 ) 
       {
	 *h2 = *h2 - 24;
	 *d2 = d + 1;
       }
    

       if (  m == 2 )
	 { 
	   if (  y == 4*((y)/4) )
	     {
	        if( *d2 > 29 ) 
		 { *d2 = *d2-29;
		 *m2 = m+1;
		 }
	     }
	   else 
	     {
	       if( *d2 > 28 ) 
		 { *d2 = *d2-28;
		   *m2 = m+1;
		 }
	     }
	 }

 
       if ( (m == 1) || (m == 3) || (m == 5) || (m ==7) || (m == 8) || (m == 10) || (m == 12))
	 { printf(" here %d  %d\n",m,d);
         if ( *d2 > 31) 
         {
          *d2 = *d2 -31;
          *m2 = m + 1;
         }
       }   

      if ( m == 4 || m == 6 || m == 9 || m ==11 )
       {
        if ( *d2 > 30) 
        {
          *d2 = *d2 -30;
          *m2 = m + 1;
        }
       }   

      if ( *m2 == 13)
	{ *y2 = y +1;
	*m2 = *m2 - 12;
	}

    }

 





  if( addmin <  0 )
 {
   if( *min2 < 0  ) 
     {
       *min2 = *min2 + 60;
       *h2 = *h2 - 1;
     }

   if( *h2 < 0 ) 
     {
       *h2 = *h2 + 24;
       *d2 = d - 1;
     }
    



   if( *d2 == 0 )
     {
          if (  m == 3 )
          { 
               if (  y == 4*((y)/4) )
               {
                *d2 = *d2 + 29;
                *m2 = m-1;
               }
              else 
              {
    	       *d2 = *d2+28;
               *m2 = m-1;
              }
          }

  
          if ( (m == 5) || (m == 7) || (m == 10) || (m == 12))
          {
               *d2 = *d2 +30;
               *m2 = m - 1;
          }   

         if ( (m ==1) || (m == 2) || (m ==4) || (m == 6) || (m == 8) || (m == 9) || (m ==11) )
         {
            *d2 = *d2 +31;
            *m2 = m - 1;
         }
    }

     if( *m2 == 0 )
     {
	   *m2 = 12;
	   *y2 = y-1;
     }
  }
  
 
}
