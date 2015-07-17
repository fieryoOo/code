#ifndef MYSAC
#define MYSAC
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cctype>

struct SAC_HD {
	float	delta,     depmin,    depmax,    scale,     odelta;    //0 - 4
	float	b,         e,         o,         a,         internal1; //5 - 9
	float	t0,        t1,        t2,        t3,        t4;        //10 - 14
	float	t5,        t6,        t7,        t8,        t9;        //15 - 19
	float	f,         resp0,     resp1,     resp2,     resp3;     //20 - 24
	float	resp4,     resp5,     resp6,     resp7,     resp8;     //25 - 29
	float	resp9,     stla,      stlo,      stel,      stdp;      //30 - 34
	float	evla,      evlo,      evel,      evdp,      unused1;   //35 - 39
	float	user0,     user1,     user2,     user3,     user4;     //40 - 44
	float	user5,     user6,     user7,     user8,     user9;     //45 - 49
	float	dist,      az,        baz,       gcarc,     internal2; 
	float	internal3, depmen,    cmpaz,     cmpinc,    unused2;   
	float	unused3,   unused4,   unused5,   unused6,   unused7;   
	float	unused8,   unused9,   unused10,  unused11,  unused12;  
	int	nzyear,    nzjday,    nzhour,    nzmin,     nzsec;     
	int	nzmsec,    internal4, internal5, internal6, npts;      
	int	internal7, internal8, unused13,  unused14,  unused15;  
	int	iftype,    idep,      iztype,    unused16,  iinst;     
	int	istreg,    ievreg,    ievtyp,    iqual,     isynth;    
	int	unused17,  unused18,  unused19,  unused20,  unused21;  
	int	unused22,  unused23,  unused24,  unused25,  unused26;  
	int	leven,     lpspol,    lovrok,    lcalda,    unused27;  
	char	kstnm[8],  kevnm[16];           
	char	khole[8],  ko[8],     ka[8];               
	char	kt0[8],    kt1[8],    kt2[8];              
	char	kt3[8],    kt4[8],    kt5[8];              
	char	kt6[8],    kt7[8],    kt8[8];              
	char	kt9[8],    kf[8],     kuser0[8];           
	char	kuser1[8], kuser2[8], kcmpnm[8];           
	char	knetwk[8], kdatrd[8], kinst[8];            

	/* ------------------------------- header streams ------------------------------- */
	bool StreamTo( std::string field, std::ostream& o ) const {
		std::transform(field.begin(), field.end(), field.begin(), ::tolower);

		std::string prefix, suffix;
		auto dumpfield = [&](const std::string& fdump) -> bool {
			bool dump = (field=="all");
			if( dump ) {
				prefix = fdump + " ";
				suffix = "\n";
			}
			return (dump || field==fdump);
		};

		if( dumpfield("dist") ) o << prefix << dist << suffix;
		if( dumpfield("az") ) o << prefix << az << suffix;
		if( dumpfield("baz") ) o << prefix << baz << suffix;
		if( dumpfield("gcarc") ) o << prefix << gcarc << suffix;
		if( dumpfield("b") ) o << prefix << b << suffix;
		if( dumpfield("e") ) o << prefix << e << suffix;

		if( dumpfield("knetwk") ) o << prefix << std::string(knetwk, 8) << suffix;
		if( dumpfield("kstnm") ) o << prefix << std::string(kstnm, 8) << suffix;
		if( dumpfield("stlo") ) o << prefix << stlo << suffix;
		if( dumpfield("stla") ) o << prefix << stla << suffix;
		if( dumpfield("stel") ) o << prefix << stel << suffix;
		if( dumpfield("stdp") ) o << prefix << stdp << suffix;

		if( dumpfield("kevnm") ) o << prefix << std::string(kevnm, 16) << suffix;
		if( dumpfield("evlo") ) o << prefix << evlo << suffix;
		if( dumpfield("evla") ) o << prefix << evla << suffix;
		if( dumpfield("evel") ) o << prefix << evel << suffix;
		if( dumpfield("evdp") ) o << prefix << evdp << suffix;

		if( dumpfield("nzyear") ) o << prefix << nzyear << suffix;
		if( dumpfield("nzjday") ) o << prefix << nzjday << suffix;
		if( dumpfield("nzhour") ) o << prefix << nzhour << suffix;
		if( dumpfield("nzmin") ) o << prefix << nzmin << suffix;
		if( dumpfield("nzsec") ) o << prefix << nzsec << suffix;
		if( dumpfield("nzmsec") ) o << prefix << nzmsec << suffix;

		if( dumpfield("kcmpnm") ) o << prefix << std::string(kcmpnm, 8) << suffix;
		if( dumpfield("cmpaz") ) o << prefix << cmpaz << suffix;
		if( dumpfield("cmpinc") ) o << prefix << cmpinc << suffix;

		if( dumpfield("o") ) o << prefix << o << suffix;
		if( dumpfield("ko") ) o << prefix << std::string(ko, 8) << suffix;
		if( dumpfield("a") ) o << prefix << a << suffix;
		if( dumpfield("ka") ) o << prefix << std::string(ka, 8) << suffix;
		if( dumpfield("f") ) o << prefix << f << suffix;
		if( dumpfield("kf") ) o << prefix << std::string(kf, 8) << suffix;

		if( dumpfield("user0") ) o << prefix << user0 << suffix;
		if( dumpfield("user1") ) o << prefix << user1 << suffix;
		if( dumpfield("user2") ) o << prefix << user2 << suffix;
		if( dumpfield("user3") ) o << prefix << user3 << suffix;
		if( dumpfield("user4") ) o << prefix << user4 << suffix;
		if( dumpfield("user5") ) o << prefix << user5 << suffix;
		if( dumpfield("user6") ) o << prefix << user6 << suffix;
		if( dumpfield("user7") ) o << prefix << user7 << suffix;
		if( dumpfield("user8") ) o << prefix << user8 << suffix;
		if( dumpfield("user9") ) o << prefix << user9 << suffix;

		else return false;

		return true;
	}

	friend std::ostream& operator<<( std::ostream& o, const SAC_HD& shd ) {
		shd.StreamTo("all", o);
		return o;
	}

	friend std::istream& operator>>( std::istream& i, SAC_HD& shd ) {
		for( std::string line; std::getline(i, line); ) {
			std::istringstream sin(line);
			std::string field; sin >> field;
			std::transform(field.begin(), field.end(), field.begin(), ::tolower);
			bool succeed = false;

			if( field == "dist" ) succeed = sin >> shd.dist;
			else if( field == "az" ) succeed = sin >> shd.az;
			else if( field == "baz" ) succeed = sin >> shd.baz;
			else if( field == "gcarc" ) succeed = sin >> shd.gcarc;
			else if( field == "b" ) succeed = sin >> shd.b;
			else if( field == "e" ) succeed = sin >> shd.e;

			else if( field == "knetwk" ) succeed = sin >> shd.knetwk;
			else if( field == "kstnm" ) succeed = sin >> shd.kstnm;
			else if( field == "stlo" ) succeed = sin >> shd.stlo;
			else if( field == "stla" ) succeed = sin >> shd.stla;
			else if( field == "stel" ) succeed = sin >> shd.stel;
			else if( field == "stdp" ) succeed = sin >> shd.stdp;

			else if( field == "kevnm" ) succeed = sin >> shd.kevnm;
			else if( field == "evlo" ) succeed = sin >> shd.evlo;
			else if( field == "evla" ) succeed = sin >> shd.evla;
			else if( field == "evel" ) succeed = sin >> shd.evel;
			else if( field == "evdp" ) succeed = sin >> shd.evdp;

			else if( field == "nzyear" ) succeed = sin >> shd.nzyear;
			else if( field == "nzjday" ) succeed = sin >> shd.nzjday;
			else if( field == "nzhour" ) succeed = sin >> shd.nzhour;
			else if( field == "nzmin" ) succeed = sin >> shd.nzmin;
			else if( field == "nzsec" ) succeed = sin >> shd.nzsec;
			else if( field == "nzmsec" ) succeed = sin >> shd.nzmsec;

			else if( field == "kcmpnm" ) succeed = sin >> shd.kcmpnm;
			else if( field == "cmpaz" ) succeed = sin >> shd.cmpaz;
			else if( field == "cmpinc" ) succeed = sin >> shd.cmpinc;

			else if( field == "o" ) succeed = sin >> shd.o;
			else if( field == "ko" ) succeed = sin >> shd.ko;
			else if( field == "a" ) succeed = sin >> shd.a;
			else if( field == "ka" ) succeed = sin >> shd.ka;
			else if( field == "f" ) succeed = sin >> shd.f;
			else if( field == "kf" ) succeed = sin >> shd.kf;

			else if( field == "user0" ) succeed = sin >> shd.user0;
			else if( field == "user1" ) succeed = sin >> shd.user1;
			else if( field == "user2" ) succeed = sin >> shd.user2;
			else if( field == "user3" ) succeed = sin >> shd.user3;
			else if( field == "user4" ) succeed = sin >> shd.user4;
			else if( field == "user5" ) succeed = sin >> shd.user5;
			else if( field == "user6" ) succeed = sin >> shd.user6;
			else if( field == "user7" ) succeed = sin >> shd.user7;
			else if( field == "user8" ) succeed = sin >> shd.user8;
			else if( field == "user9" ) succeed = sin >> shd.user9;

			if( ! succeed )
				throw std::runtime_error(std::string("Error(") + __FUNCTION__ + "): Invalid header field or value (" + line + ")");
		}
		return i;
	}

};

#ifndef MAIN
extern SAC_HD SAC_HEADER;
#else
SAC_HD SAC_HEADER;
#endif

static SAC_HD sac_null = {
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345., -12345., -12345., -12345., -12345.,
-12345, -12345, -12345, -12345, -12345,
-12345, -12345, -12345, -12345, -12345,
-12345, -12345, -12345, -12345, -12345,
-12345, -12345, -12345, -12345, -12345,
-12345, -12345, -12345, -12345, -12345,
-12345, -12345, -12345, -12345, -12345,
-12345, -12345, -12345, -12345, -12345,
-12345, -12345, -12345, -12345, -12345,
{ '-','1','2','3','4','5','\0','\0' },
{ '-','1','2','3','4','5','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0' },
{ '-','1','2','3','4','5','\0','\0' }, { '-','1','2','3','4','5','\0','\0' },
{ '-','1','2','3','4','5','\0','\0' }, { '-','1','2','3','4','5','\0','\0' },
{ '-','1','2','3','4','5','\0','\0' }, { '-','1','2','3','4','5','\0','\0' },
{ '-','1','2','3','4','5','\0','\0' }, { '-','1','2','3','4','5','\0','\0' },
{ '-','1','2','3','4','5','\0','\0' }, { '-','1','2','3','4','5','\0','\0' },
{ '-','1','2','3','4','5','\0','\0' }, { '-','1','2','3','4','5','\0','\0' },
{ '-','1','2','3','4','5','\0','\0' }, { '-','1','2','3','4','5','\0','\0' },
{ '-','1','2','3','4','5','\0','\0' }, { '-','1','2','3','4','5','\0','\0' },
{ '-','1','2','3','4','5','\0','\0' }, { '-','1','2','3','4','5','\0','\0' },
{ '-','1','2','3','4','5','\0','\0' }, { '-','1','2','3','4','5','\0','\0' },
{ '-','1','2','3','4','5','\0','\0' }
};

/* defines for logical data types */
#define TRUE    1
#define FALSE   0

/* defines for enumerated data types */
#define IREAL   0 
#define ITIME   1 
#define IRLIM   2 
#define IAMPH   3 
#define IXY     4 
#define IUNKN   5 
#define IDISP   6 
#define IVEL    7 
#define IACC    8 
#define IB      9 
#define IDAY   10 
#define IO     11 
#define IA     12 
#define IT0    13 
#define IT1    14 
#define IT2    15 
#define IT3    16 
#define IT4    17 
#define IT5    18 
#define IT6    19 
#define IT7    20 
#define IT8    21 
#define IT9    22 
#define IRADNV 23 
#define ITANNV 24 
#define IRADEV 25 
#define ITANEV 26 
#define INORTH 27 
#define IEAST  28 
#define IHORZA 29 
#define IDOWN  30 
#define IUP    31 
#define ILLLBB 32 
#define IWWSN1 33 
#define IWWSN2 34 
#define IHGLP  35 
#define ISRO   36 
#define INUCL  37 
#define IPREN  38 
#define IPOSTN 39 
#define IQUAKE 40 
#define IPREQ  41 
#define IPOSTQ 42 
#define ICHEM  43 
#define IOTHER 44 
#define IGOOD  45 
#define IGLCH  46 
#define IDROP  47 
#define ILOWSN 48 
#define IRLDTA 49 
#define IVOLTS 50 
#define INIV51 51 
#define INIV52 52 
#define INIV53 53 
#define INIV54 54 
#define INIV55 55 
#define INIV56 56 
#define INIV57 57 
#define INIV58 58 
#define INIV59 59 
#define INIV60 60 

#define FCS "%15.7f%15.7f%15.7f%15.7f%15.7f\n"
#define ICS "%10d%10d%10d%10d%10d\n"
#define CCS1 "%-8.8s%-8.8s%-8.8s\n"
#define CCS2 "%-8.8s%-16.16s\n"


#endif
