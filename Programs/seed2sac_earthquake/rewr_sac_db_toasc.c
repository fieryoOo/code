#define MAIN

#include <stdio.h>
#include <math.h>
#include "mysac.h"
#include "sac_db.h"


/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
	void sac_db_write_to_asc ( SAC_DB *sdb, char *fname )
/*--------------------------------------------------------------------------*/
{
  int ie, is;
  FILE *fi, *ff;
  static SAC_HD shd;

  ff = fopen(fname,"w");

  for ( ie = 0; ie < sdb->nev; ie++ ) for ( is = 0; is < sdb->nst; is++ )
    {
      fprintf(ff,"%s  %s  ", sdb->ev[ie].name, sdb->st[is].name );

      if ( sdb->rec[ie][is].n <= 0 ) fprintf(ff,"NO DATA\n");

      else 
	{
	  fi = fopen(sdb->rec[ie][is].fname,"rb");
	  fread(&shd, sizeof(SAC_HD), 1, fi );
	  fclose(fi);

	  fprintf(ff,"%s  t0: %d/%d:%d:%d:%g  %g s of record\n", sdb->rec[ie][is].fname, 
		   shd.nzyear, shd.nzjday, shd.nzhour, shd.nzmin, 
		   (shd.nzsec + 0.001*shd.nzmsec), shd.delta*shd.npts );
	}
    } 

  fclose(ff);
}




SAC_DB sdb;

/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
 int main (int na, char *arg[])
/*--------------------------------------------------------------------------*/
{
  FILE *ff;

  ff = fopen(arg[1],"rb");
  fread(&sdb, sizeof(SAC_DB), 1, ff );
  fclose(ff);

  sac_db_write_to_asc ( &sdb, arg[2] ); 
}
