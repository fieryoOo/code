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


	/* ------------------------------- fix broken strings ------------------------------- */
	const std::string ntname() const {
		std::stringstream ss(knetwk);
		std::string ntname; ss >> ntname;
		return ntname;
	}
	const std::string evname() const {
		std::stringstream ss(kevnm);
		std::string evname; ss >> evname;
		return evname;
	}
	const std::string stname() const {
		std::stringstream ss(kstnm);
		std::string stname; ss >> stname;
		return stname;
	}

	const std::string chname() const {
		std::stringstream ss(kcmpnm);
		std::string chname; ss >> chname;
		return chname;
	}

	/* ------------------------------- header streams ------------------------------- */
	bool StreamTo( std::string field, std::ostream& os ) const {
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

		if( dumpfield("dist") ) os << prefix << dist << suffix;
		if( dumpfield("az") ) os << prefix << az << suffix;
		if( dumpfield("baz") ) os << prefix << baz << suffix;
		if( dumpfield("gcarc") ) os << prefix << gcarc << suffix;
		if( dumpfield("b") ) os << prefix << b << suffix;
		if( dumpfield("e") ) os << prefix << e << suffix;
		if( dumpfield("delta") ) os << prefix << delta << suffix;
		if( dumpfield("npts") ) os << prefix << npts << suffix;
		if( dumpfield("depmin") ) os << prefix << depmin << suffix;
		if( dumpfield("depmax") ) os << prefix << depmax << suffix;

		if( dumpfield("knetwk") ) os << prefix << ntname() << suffix;
		if( dumpfield("kstnm") ) os << prefix << stname() << suffix;
		if( dumpfield("stlo") ) os << prefix << stlo << suffix;
		if( dumpfield("stla") ) os << prefix << stla << suffix;
		if( dumpfield("stel") ) os << prefix << stel << suffix;
		if( dumpfield("stdp") ) os << prefix << stdp << suffix;

		if( dumpfield("kevnm") ) os << prefix << evname() << suffix;
		if( dumpfield("evlo") ) os << prefix << evlo << suffix;
		if( dumpfield("evla") ) os << prefix << evla << suffix;
		if( dumpfield("evel") ) os << prefix << evel << suffix;
		if( dumpfield("evdp") ) os << prefix << evdp << suffix;

		if( dumpfield("nzyear") ) os << prefix << nzyear << suffix;
		if( dumpfield("nzjday") ) os << prefix << nzjday << suffix;
		if( dumpfield("nzhour") ) os << prefix << nzhour << suffix;
		if( dumpfield("nzmin") ) os << prefix << nzmin << suffix;
		if( dumpfield("nzsec") ) os << prefix << nzsec << suffix;
		if( dumpfield("nzmsec") ) os << prefix << nzmsec << suffix;

		if( dumpfield("kcmpnm") ) os << prefix << chname() << suffix;
		if( dumpfield("cmpaz") ) os << prefix << cmpaz << suffix;
		if( dumpfield("cmpinc") ) os << prefix << cmpinc << suffix;

		if( dumpfield("o") ) os << prefix << o << suffix;
		if( dumpfield("ko") ) os << prefix << std::string(ko, 8) << suffix;
		if( dumpfield("a") ) os << prefix << a << suffix;
		if( dumpfield("ka") ) os << prefix << std::string(ka, 8) << suffix;
		if( dumpfield("f") ) os << prefix << f << suffix;
		if( dumpfield("kf") ) os << prefix << std::string(kf, 8) << suffix;

		if( dumpfield("user0") ) os << prefix << user0 << suffix;
		if( dumpfield("user1") ) os << prefix << user1 << suffix;
		if( dumpfield("user2") ) os << prefix << user2 << suffix;
		if( dumpfield("user3") ) os << prefix << user3 << suffix;
		if( dumpfield("user4") ) os << prefix << user4 << suffix;
		if( dumpfield("user5") ) os << prefix << user5 << suffix;
		if( dumpfield("user6") ) os << prefix << user6 << suffix;
		if( dumpfield("user7") ) os << prefix << user7 << suffix;
		if( dumpfield("user8") ) os << prefix << user8 << suffix;
		if( dumpfield("user9") ) os << prefix << user9 << suffix;

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
			//bool succeed = false;

			if( field == "dist" ) sin >> shd.dist;
			else if( field == "az" ) sin >> shd.az;
			else if( field == "baz" ) sin >> shd.baz;
			else if( field == "gcarc" ) sin >> shd.gcarc;
			else if( field == "b" ) sin >> shd.b;
			else if( field == "e" ) sin >> shd.e;
			else if( field == "delta" ) sin >> shd.delta;
			else if( field == "npts" ) sin >> shd.npts;

			else if( field == "knetwk" ) sin >> shd.knetwk;
			else if( field == "kstnm" ) sin >> shd.kstnm;
			else if( field == "stlo" ) sin >> shd.stlo;
			else if( field == "stla" ) sin >> shd.stla;
			else if( field == "stel" ) sin >> shd.stel;
			else if( field == "stdp" ) sin >> shd.stdp;

			else if( field == "kevnm" ) sin >> shd.kevnm;
			else if( field == "evlo" ) sin >> shd.evlo;
			else if( field == "evla" ) sin >> shd.evla;
			else if( field == "evel" ) sin >> shd.evel;
			else if( field == "evdp" ) sin >> shd.evdp;

			else if( field == "nzyear" ) sin >> shd.nzyear;
			else if( field == "nzjday" ) sin >> shd.nzjday;
			else if( field == "nzhour" ) sin >> shd.nzhour;
			else if( field == "nzmin" ) sin >> shd.nzmin;
			else if( field == "nzsec" ) sin >> shd.nzsec;
			else if( field == "nzmsec" ) sin >> shd.nzmsec;

			else if( field == "kcmpnm" ) sin >> shd.kcmpnm;
			else if( field == "cmpaz" ) sin >> shd.cmpaz;
			else if( field == "cmpinc" ) sin >> shd.cmpinc;

			else if( field == "o" ) sin >> shd.o;
			else if( field == "ko" ) sin >> shd.ko;
			else if( field == "a" ) sin >> shd.a;
			else if( field == "ka" ) sin >> shd.ka;
			else if( field == "f" ) sin >> shd.f;
			else if( field == "kf" ) sin >> shd.kf;

			else if( field == "user0" ) sin >> shd.user0;
			else if( field == "user1" ) sin >> shd.user1;
			else if( field == "user2" ) sin >> shd.user2;
			else if( field == "user3" ) sin >> shd.user3;
			else if( field == "user4" ) sin >> shd.user4;
			else if( field == "user5" ) sin >> shd.user5;
			else if( field == "user6" ) sin >> shd.user6;
			else if( field == "user7" ) sin >> shd.user7;
			else if( field == "user8" ) sin >> shd.user8;
			else if( field == "user9" ) sin >> shd.user9;

			if( ! sin )
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
