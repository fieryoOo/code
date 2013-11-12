void SetName(int ne, int ns) {
   sprintf(sdb->rec[ne][ns].fname,"%s/%s/%s.%s.%s.SAC", sdb->mo[imonth].name, sdb->ev[ne].name, sdb->ev[ne].name, sdb->st[ns].name, ch);
   sprintf(sdb->rec[ne][ns].ft_fname,"%s/%s/ft_%s.%s.%s.SAC", sdb->mo[imonth].name, sdb->ev[ne].name, sdb->ev[ne].name, sdb->st[ns].name, ch);
   sprintf(sdb->rec[ne][ns].chan,"%s", ch );
}


int  main () 
{
   /* store old SAC and RESP files if there's any */
   bool oldfiles = false;
   MKDir("old_sac_files");
   if( wMove(".", "*.SAC", "old_sac_files", 0, &itmp) ) oldfiles = true;
   if( wMove(".", "RESP.*", "old_sac_files", 0, &itmp) ) oldfiles = true;
   int ithread;
 
   /* create the DailyRec object */
   DailyRec dailyrec;
   dailyrec.Load();
   dailyrec.ExtractSac();
  
   /* fetch back old files and remove temporary dir */
   if( oldfiles ) {
      wMove("old_sac_files", "*.SAC", ".", 0, &itmp);
      wMove("old_sac_files", "RESP.*", ".", 0, &itmp);
   }
   fRemove("old_sac_files");

   return 0;
}
