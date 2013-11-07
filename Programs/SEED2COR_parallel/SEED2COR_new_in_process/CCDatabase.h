class CCDatabase {
protected:
   // input parameters
   using namespace Seed2CorParam;
   //constructor (read in parameters and station/seed list. Create DailyRecs for all stations all days)
   CCDatabase();
   int GetParameters(char *fname);
   int isTermi(int deno);
   void InitialPthread();
   void FillStations();
   void FillMonths();
   void FillEvents();

   // operations
   class DailyRec * PrepareNextRec(int *nmonth); // assign file names and make directory if necessary

   //destructor (deallocate memories)
   ~CCDatabase();
   void CleanupPthread();
}

class DailyRec {
   // list of record objects with information of event(fseed), station(lon & lat) and singlestation-files(f_original_sac, f_ft_sac, f_am, f_ph)
   void ExtractSac();
   void RmRESP();
   void TempSpecNorm ();
}

typedef struct month
{
   char seedf[31][200];
   char name[9];
   int y, m;
}
MONTH;

typedef struct event
{
  float lat, lon;
  int yy, mm, dd, h, m, s, ms, jday;
  double t0;
  char name[40];
}
EVENT;

typedef struct station
{
  int flag;
  float lat, lon;
  char name[10];
}
STATION;

typedef struct record
{
  char fname[150];
  char ft_fname[150];
  char *resp_fname;
  char chan[7];
  double t0;
  float dt;
  int n;
}
RECORD;

#define NSTATION 2000
#define NEVENTS 31
#define NMONTH 100

typedef struct sac_dbase
{
  MONTH mo[NMONTH];
  EVENT ev[NEVENTS];
  STATION st[NSTATION];
  RECORD rec[NEVENTS][NSTATION];
  int nev, nst, nmo;
}
SAC_DB;

typedef struct sac_dbase3
{
  EVENT ev[NEVENTS];
  STATION st[NSTATION];
  RECORD rz[NEVENTS][NSTATION];
  RECORD rn[NEVENTS][NSTATION];
  RECORD re[NEVENTS][NSTATION];
  int nev, nst;
}
SAC_DB3;
