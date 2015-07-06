/* surfsyn.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    integer nsta, nstai;
    char cod[16000];
    real fig[2000], lam[2000];
    char codi[16000];
    real figi[2000], fici[2000], lami[2000];
} stn_;

#define stn_1 stn_

struct {
    real cor[2000000]	/* was [500][2][2000] */;
} trk_;

#define trk_1 trk_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__5 = 5;
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__8 = 8;

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */
/*<       implicit none >*/
/* Main program */ int MAIN__(void)
{
    /* Initialized data */

    static real geo = .993277f;
    static logical key_compr__ = FALSE_;
    static real fix_vel__ = 2.85f;
    static real const1 = 1e-9f;
    static real r0 = 6371.f;
    static real pi = 3.14159f;
    static complex step = {1.f,0.f};
    static integer im = 6;

    /* Format strings */
    static char fmt_1200[] = "(\002 iq=\002,i4,\002 vmax=\002,f10.3\002    t"
	    "ype=\002,a1)";
    static char fmt_1000[] = "(5x,\002aM\002,8x,\002mxx\002,8x,\002myy\002,8"
	    "x,\002mzz\002,8x,\002mxy\002,8x,\002mxz\002,8x,\002myz\002)";

    /* System generated locals */
    address a__1[6], a__2[2], a__3[4], a__4[3];
    integer i__1[6], i__2[2], i__3, i__4, i__5, i__6, i__7, i__8[4], i__9[3];
    real r__1, r__2;
    complex q__1, q__2, q__3, q__4, q__5;
    char ch__1[259];
    cilist ci__1;
    icilist ici__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_rsli(icilist *), e_rsli(void);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    ;
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen), f_open(olist *), s_rsfe(
	    cilist *), e_rsfe(void), f_clos(cllist *), s_rsle(cilist *), 
	    e_rsle(void), s_wsfe(cilist *), e_wsfe(void);
    double tan(doublereal), atan(doublereal);
    integer pow_ii(integer *, integer *);
    double cos(doublereal), sin(doublereal), sqrt(doublereal);

    /* Local variables */
    real azi_back__;
    extern /* Subroutine */ int calspecl_(integer *, integer *, real *, real *
	    , real *, complex *, real *, real *, real *, real *, integer *), 
	    calspecr_(integer *, integer *, real *, real *, real *, complex *,
	     complex *, real *, real *, real *, real *, real *, real *, 
	    integer *), surfread_(char *, char *, char *, char *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, ftnlen, ftnlen, ftnlen, ftnlen);
    char wavetype[2], outseism[80*50];
    integer i__, j, k, l, m;
    real t[2000], v[6000]	/* was [3][2000] */, w, f0, i0[2000];
    integer n1;
    complex ah[2000];
    real df;
    complex al[2000], bl[6];
    real cl[2000], am, fi;
    integer ii, nd;
    complex br[6];
    real cr[2000];
    integer kl;
    real fr[2000], du[3];
    integer iq;
    complex az[2000];
    integer ll;
    real ql[2000];
    integer lq, mm, nf;
    real tm[6], ul[2000], qr[2000], wl[4096];
    integer nt;
    real cs, ur[2000], dt, wr[4096], vu[3], sc;
    integer nk7;
    real del;
    char chn[8];
    real dip, amp[2048];
    integer lin;
    real azi;
    char net[8*2000];
    integer mmm, lso;
    char its[1];
    integer lsy;
    real wvl[2000];
    char str[80];
    real wvn[2];
    extern /* Subroutine */ int syn_(real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, integer *,
	     integer *, logical *, real *, real *, real *, real *, real *, 
	    integer *);
    real wvr[2000], per1, per2, cs_b__, sc_b__, drad;
    char bred[10];
    integer mode;
    real elat;
    integer narg;
    real dels, ampl[2000], dper, elon;
    char sigl[1], neti[8*2000];
    real ampr[2000];
    integer ierr, ncor;
    real slat;
    integer nper[20];
    char sigr[1];
    real dvdz[6000]	/* was [3][2000] */, dist, simr[4096], slon, srer[
	    4096], simv[4096], wvar[2000];
    complex suml;
    real srev[4096], sret[4096], simt[4096];
    char symb[3];
    integer nout;
    complex sumr;
    real slip, tmin, tmax, vmax;
    extern /* Subroutine */ int wsac_(char *, real *, real *, char *, real *, 
	    real *, real *, char *, char *, integer *, real *, ftnlen, ftnlen,
	     ftnlen, ftnlen);
    integer n2pow;
    char fnam11[256], fnam12[256], fnam13[256];
    extern integer iargc_(void);
    real elatc;
    integer nbase;
    extern /* Subroutine */ int force_(char *, char *, real *, real *, real *,
	     complex *, complex *, ftnlen, ftnlen);
    real elonc;
    char model[256];
    extern /* Subroutine */ int azidl_(real *, real *, real *, real *, real *,
	     real *, real *, real *);
    real depth;
    extern /* Subroutine */ int namer_(real *, integer *, integer *, char *, 
	    char *, char *, char *, char *, char *, integer *, ftnlen, ftnlen,
	     ftnlen, ftnlen, ftnlen, ftnlen);
    real ratio[2000];
    integer ndist;
    extern /* Subroutine */ int angles2tensor_(real *, real *, real *, real *)
	    ;
    char symbo[2];
    integer k_spec__;
    char infile__[255];
    extern /* Subroutine */ int getarg_(integer *, char *, ftnlen);
    integer idepth;
    real spread;
    extern integer lnblnk_(char *, ftnlen);
    real seisme[4096], ql_int__[4096], ul_int__[4096], permin, qr_int__[4096],
	     permax, t_curr__[2048], seismn[4096], ur_int__[4096];
    char symbik[2];
    real seismr[4096], strike, seismt[4096];
    char symbol[10];
    extern /* Subroutine */ int source_(char *, char *, real *, real *, real *
	    , real *, real *, complex *, complex *, ftnlen, ftnlen), intpol_(
	    real *, real *, integer *, real *, real *, integer *, real *, 
	    integer *);
    real seismz[4096], tstart;
    extern /* Subroutine */ int system_(char *, ftnlen);
    char outtit[80*50];
    real fi_back__;
    extern /* Subroutine */ int atracer_(char *, real *, real *, integer *, 
	    ftnlen);
    char stafile[255], runfile[255], outfile[255], current[255], outspec[80*
	    50];
    integer npoints;

    /* Fortran I/O blocks */
    static cilist io___17 = { 0, 6, 0, 0, 0 };
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___20 = { 0, 6, 0, 0, 0 };
    static cilist io___21 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___42 = { 0, 6, 0, 0, 0 };
    static cilist io___45 = { 0, 6, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 6, 0, 0, 0 };
    static cilist io___55 = { 0, 6, 0, 0, 0 };
    static cilist io___60 = { 0, 6, 0, 0, 0 };
    static cilist io___61 = { 0, 6, 0, 0, 0 };
    static cilist io___64 = { 0, 20, 1, 0, 0 };
    static cilist io___66 = { 0, 2, 0, 0, 0 };
    static cilist io___74 = { 0, 6, 0, 0, 0 };
    static cilist io___75 = { 0, 6, 0, 0, 0 };
    static cilist io___76 = { 0, 6, 0, 0, 0 };
    static cilist io___78 = { 0, 6, 0, fmt_1200, 0 };
    static cilist io___79 = { 0, 2, 0, 0, 0 };
    static cilist io___82 = { 0, 2, 0, 0, 0 };
    static cilist io___86 = { 0, 2, 0, 0, 0 };
    static cilist io___90 = { 0, 6, 0, 0, 0 };
    static cilist io___93 = { 0, 6, 0, 0, 0 };
    static cilist io___95 = { 0, 6, 0, 0, 0 };
    static cilist io___96 = { 0, 6, 0, fmt_1000, 0 };
    static cilist io___113 = { 0, 6, 0, 0, 0 };
    static cilist io___117 = { 0, 6, 0, 0, 0 };
    static cilist io___123 = { 0, 6, 0, 0, 0 };
    static cilist io___133 = { 0, 6, 0, 0, 0 };
    static cilist io___140 = { 0, 6, 0, 0, 0 };
    static cilist io___180 = { 0, 41, 0, 0, 0 };
    static cilist io___182 = { 0, 41, 0, 0, 0 };
    static cilist io___185 = { 0, 41, 0, 0, 0 };
    static cilist io___191 = { 0, 11, 0, 0, 0 };
    static cilist io___192 = { 0, 12, 0, 0, 0 };
    static cilist io___193 = { 0, 13, 0, 0, 0 };


/*<       character*8 cod(2000),codi(2000),net(2000),neti(2000) >*/
/*<       integer*4   nsta,nstai,nsize >*/
/*<       real*4      fig(2000),lam(2000),figi(2000),fici(2000),lami(2000) >*/
/*<       common /stn/nsta,nstai,cod,fig,lam,codi,figi,fici,lami >*/
/* --- */
/*<       real*4 cor(500,2,2000) >*/
/*<       common /trk/ cor >*/
/* --- */
/*<       character*8 chn >*/
/*<       complex*8 br(6),bl(6),sumr,suml >*/
/*<       complex*8 al(2000),az(2000),ah(2000) >*/
/*<       complex*8 step >*/
/*<       real*4    elat,elatc,elon,elonc,slat,slon >*/
/*<       real*4 v(3,2000),dvdz(3,2000),tm(6) >*/
/*<       real*4 ampr(2000),ampl(2000),ratio(2000),qR(2000),qL(2000) >*/
/*<       real*4 amp(2048),T_curr(2048) >*/
/*<       parameter (nsize=4096) >*/
/*<       real*4 cr(2000),ur(2000),wvr(2000),t(2000),wvar(2000) >*/
/*<       real*4 cl(2000),ul(2000),wvl(2000),fr(2000) >*/
/*<       real*4 wr(nsize),wl(nsize) >*/
/*<       real*4 sreV(nsize),simV(nsize) >*/
/*<       real*4 sreR(nsize),simR(nsize) >*/
/*<       real*4 sreT(nsize),simT(nsize) >*/
/*<       real*4 vu(3),du(3),wvn(2) >*/
/*<       real*4 seismz(nsize),seismn(nsize),seisme(nsize) >*/
/*<       real*4 ul_int(nsize),ur_int(nsize) >*/
/*<       real*4 seismt(nsize),seismr(nsize),qL_int(nsize),qR_int(nsize) >*/
/*<       integer idate(5),nper(20) >*/
/*<       character*2 wavetype >*/
/*<       character*10 symbol,bred >*/
/*<       character*3  symb >*/
/*<       character*2  symbo,symbik >*/
/*<       character*1 sigR,sigL,its >*/
/*<       character*255 runfile,infile,outfile,current,stafile >*/
/*<       character*256 model,fnam11,fnam12,fnam13 >*/
/*<       character*80 str,outseism(50),outtit(50),outspec(50) >*/
/*<       logical key_compr >*/
/* --- */
/*<       integer*4 i,idepth,ierr,ii,im,iq,k,j,k_spec,kl,l,lin,ll,lq >*/
/*<       integer*4 lso,lsy,m,mm,mmm,mode,n1,n2pow,narg,iargc,nbase,ncor >*/
/*<       integer*4 nd,NDIST,nf,nk7,nout,lnblnk,npoints,nt >*/
/*<       real*4    GEO,aM,azi,azi_back,cazz,cince,cinch,cincz,const1 >*/
/*<       real*4    cs,cs_b,DEL,DELS,df,dist,dper,drad,dt,f0,fi,fi_back >*/
/*<       real*4    fix_vel,cincn,depth,strike,dip,slip,per1,per2 >*/
/*<       real*4    PERMAX,PERMIN,pi,R0,sc,sc_b,sdate,spread,sqr_8pi >*/
/*<       real*4    TMIN,TMAX,tstart,vmax,w,I0(2000) >*/
/* --- */
/*<       data GEO/0.993277/ >*/
/*<       data pi/3.14159/,sqr_8pi/5.01325655/,step/(1.0,0.0)/ >*/
/*<       data idate/1999,1,1,1,1/,sdate/0.0/,im/6/ >*/
/*<       data cincz/0.0/,cazz/0.0/,cinch/90.0/,cincn/0.0/,cince/90.0/ >*/
/*<       data key_compr/.false./,fix_vel/2.85/,const1/1.E-09/,R0/6371./ >*/
/* -----------ARGUMENTS----------------------------------------------- */
/* -----runfile: file with parameters for running program------------- */
/* -----infile: output file of SURFLEV-------------------------------- */
/* -----outfiles: output files of SURFLEV----------------------------- */
/* -----wavetype: R/L/RL/LR ------------------------------------------ */
/* -----mode:0,1,2,3,...----------------------------------------------- */
/* -------------INITIATION--------------------------------------------S */
/*<        print*,'to make synthetic seismogram using output of SURF_DEEP' >*/
    s_wsle(&io___17);
    do_lio(&c__9, &c__1, "to make synthetic seismogram using output of SURF_"
	    "DEEP", (ftnlen)54);
    e_wsle();
/*<        print*,'vertical direction is down; displacements as output;' >*/
    s_wsle(&io___18);
    do_lio(&c__9, &c__1, "vertical direction is down; displacements as outpu"
	    "t;", (ftnlen)52);
    e_wsle();
/*<        print*,'step-function source-time function for earthquake;' >*/
    s_wsle(&io___19);
    do_lio(&c__9, &c__1, "step-function source-time function for earthquake;",
	     (ftnlen)50);
    e_wsle();
/*<        print*,'delta-function as a source-function for simple force or' >*/
    s_wsle(&io___20);
    do_lio(&c__9, &c__1, "delta-function as a source-function for simple for"
	    "ce or", (ftnlen)55);
    e_wsle();
/*<        print*,'explosion; (-pi/4) and polarization factors are Included' >*/
    s_wsle(&io___21);
    do_lio(&c__9, &c__1, "explosion; (-pi/4) and polarization factors are In"
	    "cluded", (ftnlen)56);
    e_wsle();
/* ------------------------------------------------------------------------ */
/* ------notations of eigenfunctions follow the Kluwer's book------ */
/* -------------------Initiation-------------------------------------------S */
/*<            narg=iargc() >*/
    narg = iargc_();
/*<            drad=pi/180. >*/
    drad = pi / 180.f;
/*<       if (narg.lt.6.or.narg.gt.8)then >*/
    if (narg < 6 || narg > 8) {
/*<        >*/
	s_wsle(&io___24);
	do_lio(&c__9, &c__1, "usage : surfsyn runfile infile outfiles wavety"
		"pe mode", (ftnlen)53);
	do_lio(&c__9, &c__1, "  depth [-c [fix_vel]]", (ftnlen)22);
	e_wsle();
/*<                     stop >*/
	s_stop("", (ftnlen)0);
/*<             end if >*/
    }
/*<            call GETARG(1,runfile) >*/
    getarg_(&c__1, runfile, (ftnlen)255);
/*<            call GETARG(2,infile) >*/
    getarg_(&c__2, infile__, (ftnlen)255);
/*<            lin=lnblnk(infile) >*/
    lin = lnblnk_(infile__, (ftnlen)255);
/*<            call GETARG(3,outfile) >*/
    getarg_(&c__3, outfile, (ftnlen)255);
/*<            nout=lnblnk(outfile) >*/
    nout = lnblnk_(outfile, (ftnlen)255);
/*<            call GETARG(4,wavetype) >*/
    getarg_(&c__4, wavetype, (ftnlen)2);
/*<            call GETARG(5,symbo) >*/
    getarg_(&c__5, symbo, (ftnlen)2);
/*<            read(symbo,*)mode >*/
    ici__1.icierr = 0;
    ici__1.iciend = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = 2;
    ici__1.iciunit = symbo;
    ici__1.icifmt = 0;
    s_rsli(&ici__1);
    do_lio(&c__3, &c__1, (char *)&mode, (ftnlen)sizeof(integer));
    e_rsli();
/*<            mode=mode+1 >*/
    ++mode;
/*<            symbik='  ' >*/
    s_copy(symbik, "  ", (ftnlen)2, (ftnlen)2);
/*<            if(mode.lt.10)write(symbik,'(i1)')mode >*/
    if (mode < 10) {
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 2;
	ici__1.iciunit = symbik;
	ici__1.icifmt = "(i1)";
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&mode, (ftnlen)sizeof(integer));
	e_wsfi();
    }
/*<            if(mode.ge.10)write(symbik,'(i2)')mode >*/
    if (mode >= 10) {
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 2;
	ici__1.iciunit = symbik;
	ici__1.icifmt = "(i2)";
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&mode, (ftnlen)sizeof(integer));
	e_wsfi();
    }
/*<            lso=lnblnk(symbo) >*/
    lso = lnblnk_(symbo, (ftnlen)2);
/*<            call GETARG(6,symbol) >*/
    getarg_(&c__6, symbol, (ftnlen)10);
/*<            read(symbol,*)depth >*/
    ici__1.icierr = 0;
    ici__1.iciend = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = 10;
    ici__1.iciunit = symbol;
    ici__1.icifmt = 0;
    s_rsli(&ici__1);
    do_lio(&c__4, &c__1, (char *)&depth, (ftnlen)sizeof(real));
    e_rsli();
/*<            idepth=depth >*/
    idepth = depth;
/*<            symb='   ' >*/
    s_copy(symb, "   ", (ftnlen)3, (ftnlen)3);
/*<            if(idepth.lt.10)write(symb,'(I1)')idepth >*/
    if (idepth < 10) {
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 3;
	ici__1.iciunit = symb;
	ici__1.icifmt = "(I1)";
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&idepth, (ftnlen)sizeof(integer));
	e_wsfi();
    }
/*<            if(idepth.ge.10.and.idepth.lt.100)write(symb,'(I2)')idepth >*/
    if (idepth >= 10 && idepth < 100) {
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 3;
	ici__1.iciunit = symb;
	ici__1.icifmt = "(I2)";
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&idepth, (ftnlen)sizeof(integer));
	e_wsfi();
    }
/*<            if(idepth.ge.100)write(symb,'(I3)')idepth >*/
    if (idepth >= 100) {
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 3;
	ici__1.iciunit = symb;
	ici__1.icifmt = "(I3)";
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&idepth, (ftnlen)sizeof(integer));
	e_wsfi();
    }
/*<            lsy=lnblnk(symb) >*/
    lsy = lnblnk_(symb, (ftnlen)3);
/*<            n1=nout+3+lsy+lso >*/
    n1 = nout + 3 + lsy + lso;
/*<       current(1:n1)=outfile(1:nout)//'_'//symb(1:lsy)//'_'//symbo(1:lso)//'_' >*/
/* Writing concatenation */
    i__1[0] = nout, a__1[0] = outfile;
    i__1[1] = 1, a__1[1] = "_";
    i__1[2] = lsy, a__1[2] = symb;
    i__1[3] = 1, a__1[3] = "_";
    i__1[4] = lso, a__1[4] = symbo;
    i__1[5] = 1, a__1[5] = "_";
    s_cat(current, a__1, i__1, &c__6, n1);
/*<        >*/
    if (s_cmp(wavetype, "R ", (ftnlen)2, (ftnlen)2) != 0 && s_cmp(wavetype, 
	    "L ", (ftnlen)2, (ftnlen)2) != 0 && s_cmp(wavetype, "RL", (ftnlen)
	    2, (ftnlen)2) != 0 && s_cmp(wavetype, "LR", (ftnlen)2, (ftnlen)2) 
	    != 0) {
/*<              print*, 'WRONG WAVE TYPE' >*/
	s_wsle(&io___42);
	do_lio(&c__9, &c__1, "WRONG WAVE TYPE", (ftnlen)15);
	e_wsle();
/*<             STOP  >*/
	s_stop("", (ftnlen)0);
/*<             end if >*/
    }
/*<             sigR='-' >*/
    *(unsigned char *)sigr = '-';
/*<             sigL='-' >*/
    *(unsigned char *)sigl = '-';
/*<             if(wavetype(1:1).eq.'R'.or.wavetype(2:2).eq.'R') sigR='+' >*/
    if (*(unsigned char *)wavetype == 'R' || *(unsigned char *)&wavetype[1] ==
	     'R') {
	*(unsigned char *)sigr = '+';
    }
/*<             if(wavetype(1:1).eq.'L'.or.wavetype(2:2).eq.'L') sigL='+' >*/
    if (*(unsigned char *)wavetype == 'L' || *(unsigned char *)&wavetype[1] ==
	     'L') {
	*(unsigned char *)sigl = '+';
    }
/*<             if(narg.gt.6)THEn >*/
    if (narg > 6) {
/*<             call GETARG(7,symbol) >*/
	getarg_(&c__7, symbol, (ftnlen)10);
/*<             if(symbol(1:2).eq.'-c')then >*/
	if (s_cmp(symbol, "-c", (ftnlen)2, (ftnlen)2) == 0) {
/*<             key_compr=.true. >*/
	    key_compr__ = TRUE_;
/*<             print*,'nondispersive seismograms will be created' >*/
	    s_wsle(&io___45);
	    do_lio(&c__9, &c__1, "nondispersive seismograms will be created", 
		    (ftnlen)41);
	    e_wsle();
/*<             if(narg.eq.8)then >*/
	    if (narg == 8) {
/*<             call GETARG(8,symbol) >*/
		getarg_(&c__8, symbol, (ftnlen)10);
/*<             read(symbol,*)fix_vel >*/
		ici__1.icierr = 0;
		ici__1.iciend = 0;
		ici__1.icirnum = 1;
		ici__1.icirlen = 10;
		ici__1.iciunit = symbol;
		ici__1.icifmt = 0;
		s_rsli(&ici__1);
		do_lio(&c__4, &c__1, (char *)&fix_vel__, (ftnlen)sizeof(real))
			;
		e_rsli();
/*<                      endif >*/
	    }
/*<        >*/
	    s_wsle(&io___46);
	    do_lio(&c__9, &c__1, "nondispersive seismograms will be created; "
		    , (ftnlen)43);
	    do_lio(&c__9, &c__1, "fixed velocity=", (ftnlen)15);
	    do_lio(&c__4, &c__1, (char *)&fix_vel__, (ftnlen)sizeof(real));
	    e_wsle();
/*<                                    else >*/
	} else {
/*<             print*,'Wrong argument, try again' >*/
	    s_wsle(&io___47);
	    do_lio(&c__9, &c__1, "Wrong argument, try again", (ftnlen)25);
	    e_wsle();
/*<                                   endif  >*/
	}
/*<                                   ENDif  >*/
    }
/* -------------------Initiation-------------------------------------------E */
/* -----------reading phase velocity info-----S */
/*<             open(1,file=infile(1:lin)//'.phv',status='OLD') >*/
    o__1.oerr = 0;
    o__1.ounit = 1;
    o__1.ofnmlen = lin + 4;
/* Writing concatenation */
    i__2[0] = lin, a__2[0] = infile__;
    i__2[1] = 4, a__2[1] = ".phv";
    s_cat(ch__1, a__2, i__2, &c__2, (ftnlen)259);
    o__1.ofnm = ch__1;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/*<             k=1 >*/
    k = 1;
/*<             m=0    >*/
    m = 0;
/*<             do i=1,1000 >*/
    for (i__ = 1; i__ <= 1000; ++i__) {
/*<             read(1,'(a)',end=9988)bred >*/
	ci__1.cierr = 0;
	ci__1.ciend = 1;
	ci__1.ciunit = 1;
	ci__1.cifmt = "(a)";
	i__3 = s_rsfe(&ci__1);
	if (i__3 != 0) {
	    goto L9988;
	}
	i__3 = do_fio(&c__1, bred, (ftnlen)10);
	if (i__3 != 0) {
	    goto L9988;
	}
	i__3 = e_rsfe();
	if (i__3 != 0) {
	    goto L9988;
	}
/*<             if(i.eq.1)read(bred,*)per1 >*/
	if (i__ == 1) {
	    ici__1.icierr = 0;
	    ici__1.iciend = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = 10;
	    ici__1.iciunit = bred;
	    ici__1.icifmt = 0;
	    s_rsli(&ici__1);
	    do_lio(&c__4, &c__1, (char *)&per1, (ftnlen)sizeof(real));
	    e_rsli();
	}
/*<             if(i.eq.2)read(bred,*)per2 >*/
	if (i__ == 2) {
	    ici__1.icierr = 0;
	    ici__1.iciend = 0;
	    ici__1.icirnum = 1;
	    ici__1.icirlen = 10;
	    ici__1.iciunit = bred;
	    ici__1.icifmt = 0;
	    s_rsli(&ici__1);
	    do_lio(&c__4, &c__1, (char *)&per2, (ftnlen)sizeof(real));
	    e_rsli();
	}
/*<             m=m+1 >*/
	++m;
/*<             if(bred(6:10).eq.'     ')then >*/
	if (s_cmp(bred + 5, "     ", (ftnlen)5, (ftnlen)5) == 0) {
/*<             nper(k)=m-1 >*/
	    nper[k - 1] = m - 1;
/*<             read(1,'(1X)') >*/
	    ci__1.cierr = 0;
	    ci__1.ciend = 0;
	    ci__1.ciunit = 1;
	    ci__1.cifmt = "(1X)";
	    s_rsfe(&ci__1);
	    e_rsfe();
/*<             k=k+1 >*/
	    ++k;
/*<             m=0        >*/
	    m = 0;
/*<                                      endif >*/
	}
/*<             enddo >*/
    }
/*< 9988        print*,k-1,' modes in input' >*/
L9988:
    s_wsle(&io___55);
    i__3 = k - 1;
    do_lio(&c__3, &c__1, (char *)&i__3, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " modes in input", (ftnlen)15);
    e_wsle();
/*<             close(1) >*/
    cl__1.cerr = 0;
    cl__1.cunit = 1;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*<             nt=nper(mode) >*/
    nt = nper[mode - 1];
/*<             dper=per2-per1 >*/
    dper = per2 - per1;
/*<             PERMAX=(nt-1)*dper+per1 >*/
    permax = (nt - 1) * dper + per1;
/*<             PERMIN=per1 >*/
    permin = per1;
/*<        >*/
    s_wsle(&io___60);
    do_lio(&c__3, &c__1, (char *)&nt, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " periods for the mode ", (ftnlen)22);
    do_lio(&c__3, &c__1, (char *)&mode, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " per_incr=", (ftnlen)10);
    do_lio(&c__4, &c__1, (char *)&dper, (ftnlen)sizeof(real));
    do_lio(&c__9, &c__1, " PERMAX=", (ftnlen)8);
    do_lio(&c__4, &c__1, (char *)&permax, (ftnlen)sizeof(real));
    e_wsle();
/*<             if(nt.gt.1998)stop 'too many periods!' >*/
    if (nt > 1998) {
	s_stop("too many periods!", (ftnlen)17);
    }
/* -----------reading phase velocity info-----E */
/* -----------Reading of input parameters--------------S */
/*<            open(2,file=runfile,STATUS='OLD') >*/
    o__1.oerr = 0;
    o__1.ounit = 2;
    o__1.ofnmlen = 255;
    o__1.ofnm = runfile;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/* -----------nt is a number of periods------------ */
/* -----------nd is a number of depth points at SURLEV-------- */
/* -----------npoints is a number of points of the output seismpogram */
/* -----------dt is an output time increment------- */
/* -----------elat,elon  are event coordinates */
/* -----------stlat,stlon  are station coordinates */
/* -----------depth is a source depth (one of given at SURFLEV)----- */
/* -----------TMIN-the lowest untapered period--------------------------- */
/* -----------TMAX-the highest untapered period-------------------------- */
/* -----------iq  is a power of taper          -------------------------- */
/* -----------vmax is a maximal velocity for starting point ------------- */
/* -----------its is type of source:F/single force;Q/quake-tensor;E/explosion */
/* -----------its is type of source:C/quake-double dipole---------------- */
/*<       PRINT*,'Wavetype=',wavetype,' mode=',symbo >*/
    s_wsle(&io___61);
    do_lio(&c__9, &c__1, "Wavetype=", (ftnlen)9);
    do_lio(&c__9, &c__1, wavetype, (ftnlen)2);
    do_lio(&c__9, &c__1, " mode=", (ftnlen)6);
    do_lio(&c__9, &c__1, symbo, (ftnlen)2);
    e_wsle();
/* --- read parameter file  --- */
/*<       read(2,'(a200)')model >*/
    ci__1.cierr = 0;
    ci__1.ciend = 0;
    ci__1.ciunit = 2;
    ci__1.cifmt = "(a200)";
    s_rsfe(&ci__1);
    do_fio(&c__1, model, (ftnlen)256);
    e_rsfe();
/*<       read(2,'(a200)')stafile >*/
    ci__1.cierr = 0;
    ci__1.ciend = 0;
    ci__1.ciunit = 2;
    ci__1.cifmt = "(a200)";
    s_rsfe(&ci__1);
    do_fio(&c__1, stafile, (ftnlen)255);
    e_rsfe();
/* --  read station list  --- */
/*<       open(20,file=stafile,status='old') >*/
    o__1.oerr = 0;
    o__1.ounit = 20;
    o__1.ofnmlen = 255;
    o__1.ofnm = stafile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/*<       nsta = 1 >*/
    stn_1.nsta = 1;
/*<     1 read(20,*,end=98) net(nsta),cod(nsta),fig(nsta),lam(nsta) >*/
L1:
    i__3 = s_rsle(&io___64);
    if (i__3 != 0) {
	goto L98;
    }
    i__3 = do_lio(&c__9, &c__1, net + (stn_1.nsta - 1 << 3), (ftnlen)8);
    if (i__3 != 0) {
	goto L98;
    }
    i__3 = do_lio(&c__9, &c__1, stn_1.cod + (stn_1.nsta - 1 << 3), (ftnlen)8);
    if (i__3 != 0) {
	goto L98;
    }
    i__3 = do_lio(&c__4, &c__1, (char *)&stn_1.fig[stn_1.nsta - 1], (ftnlen)
	    sizeof(real));
    if (i__3 != 0) {
	goto L98;
    }
    i__3 = do_lio(&c__4, &c__1, (char *)&stn_1.lam[stn_1.nsta - 1], (ftnlen)
	    sizeof(real));
    if (i__3 != 0) {
	goto L98;
    }
    i__3 = e_rsle();
    if (i__3 != 0) {
	goto L98;
    }
/*<       nsta = nsta+1 >*/
    ++stn_1.nsta;
/*<       goto 1 >*/
    goto L1;
/*<    98 nsta = nsta-1 >*/
L98:
    --stn_1.nsta;
/*<       close(20) >*/
    cl__1.cerr = 0;
    cl__1.cunit = 20;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*<       read(2,*)nd,npoints,dt,TMIN,TMAX,iq,vmax >*/
    s_rsle(&io___66);
    do_lio(&c__3, &c__1, (char *)&nd, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&npoints, (ftnlen)sizeof(integer));
    do_lio(&c__4, &c__1, (char *)&dt, (ftnlen)sizeof(real));
    do_lio(&c__4, &c__1, (char *)&tmin, (ftnlen)sizeof(real));
    do_lio(&c__4, &c__1, (char *)&tmax, (ftnlen)sizeof(real));
    do_lio(&c__3, &c__1, (char *)&iq, (ftnlen)sizeof(integer));
    do_lio(&c__4, &c__1, (char *)&vmax, (ftnlen)sizeof(real));
    e_rsle();
/*<       if(TMAX.gt.(nt-1)*dper+per1)TMAX=PERMAX-2.*dper >*/
    if (tmax > (nt - 1) * dper + per1) {
	tmax = permax - dper * 2.f;
    }
/*<       if(TMIN.lt.per1)TMIN=PERMIN+2.*dper >*/
    if (tmin < per1) {
	tmin = permin + dper * 2.f;
    }
/*<       print*,' ndepth= ',nd >*/
    s_wsle(&io___74);
    do_lio(&c__9, &c__1, " ndepth= ", (ftnlen)9);
    do_lio(&c__3, &c__1, (char *)&nd, (ftnlen)sizeof(integer));
    e_wsle();
/*<       print*,'npoints= ',npoints, ' dt = ',dt >*/
    s_wsle(&io___75);
    do_lio(&c__9, &c__1, "npoints= ", (ftnlen)9);
    do_lio(&c__3, &c__1, (char *)&npoints, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " dt = ", (ftnlen)6);
    do_lio(&c__4, &c__1, (char *)&dt, (ftnlen)sizeof(real));
    e_wsle();
/*<       print*,'depth= ',depth,' Tmin=',tmin,' Tmax=',Tmax >*/
    s_wsle(&io___76);
    do_lio(&c__9, &c__1, "depth= ", (ftnlen)7);
    do_lio(&c__4, &c__1, (char *)&depth, (ftnlen)sizeof(real));
    do_lio(&c__9, &c__1, " Tmin=", (ftnlen)6);
    do_lio(&c__4, &c__1, (char *)&tmin, (ftnlen)sizeof(real));
    do_lio(&c__9, &c__1, " Tmax=", (ftnlen)6);
    do_lio(&c__4, &c__1, (char *)&tmax, (ftnlen)sizeof(real));
    e_wsle();
/*<            read(2,'(a)') its >*/
    ci__1.cierr = 0;
    ci__1.ciend = 0;
    ci__1.ciunit = 2;
    ci__1.cifmt = "(a)";
    s_rsfe(&ci__1);
    do_fio(&c__1, its, (ftnlen)1);
    e_rsfe();
/*<        >*/
    if (*(unsigned char *)its != 'Q' && *(unsigned char *)its != 'E' && *(
	    unsigned char *)its != 'F' && *(unsigned char *)its != 'C') {
	s_stop(" WRONG TYPE OF SOURCE: Q,E or F are wanted", (ftnlen)42);
    }
/*<       print 1200 ,          iq,vmax,its  >*/
    s_wsfe(&io___78);
    do_fio(&c__1, (char *)&iq, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&vmax, (ftnlen)sizeof(real));
    do_fio(&c__1, its, (ftnlen)1);
    e_wsfe();
/*< 1200   format(' iq=',i4,' vmax=',F10.3'    type=',a1) >*/
/* --if its=Q ---time function is a step function----------- */
/* ---in the following order:xx,yy,zz,xy,xz,yz--- */
/* -from CMT is should be taken as */
/* as del-del, fi-fi, rr, -(fi-del),r-del,-(fi-r) */
/* --if its=E or F time source function is a delta function */
/* --if its=C ---time function is a step function and moment tensor is found from angles */
/* --strike, dip and rake(slip), aM is the best double couple moment; */
/* -----time source function is a step function for earthquake-- */
/*<            if(its.eq.'F'.or.its.eq.'E') im=3  >*/
    if (*(unsigned char *)its == 'F' || *(unsigned char *)its == 'E') {
	im = 3;
    }
/*<            if(its.ne.'C')read(2,*)aM,(tm(m),m=1,im) >*/
    if (*(unsigned char *)its != 'C') {
	s_rsle(&io___79);
	do_lio(&c__4, &c__1, (char *)&am, (ftnlen)sizeof(real));
	i__3 = im;
	for (m = 1; m <= i__3; ++m) {
	    do_lio(&c__4, &c__1, (char *)&tm[m - 1], (ftnlen)sizeof(real));
	}
	e_rsle();
    }
/*<             if(its.eq.'C') then >*/
    if (*(unsigned char *)its == 'C') {
/*<        read(2,*)aM,strike,dip,slip >*/
	s_rsle(&io___82);
	do_lio(&c__4, &c__1, (char *)&am, (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&strike, (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&dip, (ftnlen)sizeof(real));
	do_lio(&c__4, &c__1, (char *)&slip, (ftnlen)sizeof(real));
	e_rsle();
/*<        call angles2tensor(strike,dip,slip,tm) >*/
	angles2tensor_(&strike, &dip, &slip, tm);
/*<                            endif >*/
    }
/* --- read requested station list --- */
/*<        read(2,*) elat,elon >*/
    s_rsle(&io___86);
    do_lio(&c__4, &c__1, (char *)&elat, (ftnlen)sizeof(real));
    do_lio(&c__4, &c__1, (char *)&elon, (ftnlen)sizeof(real));
    e_rsle();
/*<         elatc = atan(GEO*tan(drad*elat))/drad >*/
    elatc = atan(geo * tan(drad * elat)) / drad;
/*<        write(*,*) ' ELAT=',elat,' ELON=',elon >*/
    s_wsle(&io___90);
    do_lio(&c__9, &c__1, " ELAT=", (ftnlen)6);
    do_lio(&c__4, &c__1, (char *)&elat, (ftnlen)sizeof(real));
    do_lio(&c__9, &c__1, " ELON=", (ftnlen)6);
    do_lio(&c__4, &c__1, (char *)&elon, (ftnlen)sizeof(real));
    e_wsle();
/*<        nstai = 1 >*/
    stn_1.nstai = 1;
/*<     3  read(2,'(a80)',end=4) str >*/
L3:
    ci__1.cierr = 0;
    ci__1.ciend = 1;
    ci__1.ciunit = 2;
    ci__1.cifmt = "(a80)";
    i__3 = s_rsfe(&ci__1);
    if (i__3 != 0) {
	goto L4;
    }
    i__3 = do_fio(&c__1, str, (ftnlen)80);
    if (i__3 != 0) {
	goto L4;
    }
    i__3 = e_rsfe();
    if (i__3 != 0) {
	goto L4;
    }
/*<        if(str(1:1).eq.' ') then >*/
    if (*(unsigned char *)str == ' ') {
/*<          read(str,*) figi(nstai),lami(nstai) >*/
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 80;
	ici__1.iciunit = str;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__4, &c__1, (char *)&stn_1.figi[stn_1.nstai - 1], (ftnlen)
		sizeof(real));
	do_lio(&c__4, &c__1, (char *)&stn_1.lami[stn_1.nstai - 1], (ftnlen)
		sizeof(real));
	e_rsli();
/*<          fici(nstai) = atan(GEO*tan(drad*figi(nstai)))/drad >*/
	stn_1.fici[stn_1.nstai - 1] = atan(geo * tan(drad * stn_1.figi[
		stn_1.nstai - 1])) / drad;
/*<          write(codi(nstai),'("Z",i04)') nstai  >*/
	ici__1.icierr = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 8;
	ici__1.iciunit = stn_1.codi + (stn_1.nstai - 1 << 3);
	ici__1.icifmt = "(\"Z\",i04)";
	s_wsfi(&ici__1);
	do_fio(&c__1, (char *)&stn_1.nstai, (ftnlen)sizeof(integer));
	e_wsfi();
/*<          neti(nstai) = 'ZZ' >*/
	s_copy(neti + (stn_1.nstai - 1 << 3), "ZZ", (ftnlen)8, (ftnlen)2);
/*<          do i =1,lnblnk(codi(nstai)) >*/
	i__3 = lnblnk_(stn_1.codi + (stn_1.nstai - 1 << 3), (ftnlen)8);
	for (i__ = 1; i__ <= i__3; ++i__) {
/*<             if(codi(nstai)(i:i).eq.' ') codi(nstai)(i:i)='0' >*/
	    if (*(unsigned char *)&stn_1.codi[(stn_1.nstai - 1 << 3) + (i__ - 
		    1)] == ' ') {
		*(unsigned char *)&stn_1.codi[(stn_1.nstai - 1 << 3) + (i__ - 
			1)] = '0';
	    }
/*<          enddo >*/
	}
/*<        >*/
	s_wsle(&io___93);
	do_lio(&c__3, &c__1, (char *)&stn_1.nstai, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " STA=", (ftnlen)5);
	do_lio(&c__9, &c__1, stn_1.codi + (stn_1.nstai - 1 << 3), (ftnlen)8);
	do_lio(&c__9, &c__1, " NET=", (ftnlen)5);
	do_lio(&c__9, &c__1, neti + (stn_1.nstai - 1 << 3), (ftnlen)8);
	do_lio(&c__9, &c__1, " STLAT=", (ftnlen)7);
	do_lio(&c__4, &c__1, (char *)&stn_1.figi[stn_1.nstai - 1], (ftnlen)
		sizeof(real));
	do_lio(&c__9, &c__1, " STLON=", (ftnlen)7);
	do_lio(&c__4, &c__1, (char *)&stn_1.lami[stn_1.nstai - 1], (ftnlen)
		sizeof(real));
	e_wsle();
/*<        else >*/
    } else {
/*<          read(str,*) codi(nstai) >*/
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = 80;
	ici__1.iciunit = str;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__9, &c__1, stn_1.codi + (stn_1.nstai - 1 << 3), (ftnlen)8);
	e_rsli();
/*<          do ii = 1,nsta >*/
	i__3 = stn_1.nsta;
	for (ii = 1; ii <= i__3; ++ii) {
/*<            if(cod(ii).eq.codi(nstai)) then >*/
	    if (s_cmp(stn_1.cod + (ii - 1 << 3), stn_1.codi + (stn_1.nstai - 
		    1 << 3), (ftnlen)8, (ftnlen)8) == 0) {
/*<              figi(nstai) = fig(ii) >*/
		stn_1.figi[stn_1.nstai - 1] = stn_1.fig[ii - 1];
/*<              fici(nstai) = atan(GEO*tan(drad*fig(ii)))/drad >*/
		stn_1.fici[stn_1.nstai - 1] = atan(geo * tan(drad * stn_1.fig[
			ii - 1])) / drad;
/*<              lami(nstai) = lam(ii) >*/
		stn_1.lami[stn_1.nstai - 1] = stn_1.lam[ii - 1];
/*<              neti(nstai) = net(ii) >*/
		s_copy(neti + (stn_1.nstai - 1 << 3), net + (ii - 1 << 3), (
			ftnlen)8, (ftnlen)8);
/*<        >*/
		s_wsle(&io___95);
		do_lio(&c__3, &c__1, (char *)&stn_1.nstai, (ftnlen)sizeof(
			integer));
		do_lio(&c__9, &c__1, " STA=", (ftnlen)5);
		do_lio(&c__9, &c__1, stn_1.codi + (stn_1.nstai - 1 << 3), (
			ftnlen)8);
		do_lio(&c__9, &c__1, " NET=", (ftnlen)5);
		do_lio(&c__9, &c__1, neti + (stn_1.nstai - 1 << 3), (ftnlen)8)
			;
		do_lio(&c__9, &c__1, " STLAT=", (ftnlen)7);
		do_lio(&c__4, &c__1, (char *)&stn_1.figi[stn_1.nstai - 1], (
			ftnlen)sizeof(real));
		do_lio(&c__9, &c__1, " STLON=", (ftnlen)7);
		do_lio(&c__4, &c__1, (char *)&stn_1.lami[stn_1.nstai - 1], (
			ftnlen)sizeof(real));
		e_wsle();
/*<              goto 5 >*/
		goto L5;
/*<            endif >*/
	    }
/*<          enddo >*/
	}
/*<          goto 3 >*/
	goto L3;
/*<        endif >*/
    }
/*<     5  nstai = nstai+1 >*/
L5:
    ++stn_1.nstai;
/*<        goto 3 >*/
    goto L3;
/*<     4  nstai = nstai-1 >*/
L4:
    --stn_1.nstai;
/*<        if(its.eq.'F')write(*,'(5x,"aM",8x,"fx ",8x,"fy ",8x,"fzz")') >*/
    if (*(unsigned char *)its == 'F') {
	ci__1.cierr = 0;
	ci__1.ciunit = 6;
	ci__1.cifmt = "(5x,\"aM\",8x,\"fx \",8x,\"fy \",8x,\"fzz\")";
	s_wsfe(&ci__1);
	e_wsfe();
    }
/*<        if(its.ne.'F')write(*,1000) >*/
    if (*(unsigned char *)its != 'F') {
	s_wsfe(&io___96);
	e_wsfe();
    }
/*<        write(*,'(7e11.3)') aM,(tm(i),i=1,im) >*/
    ci__1.cierr = 0;
    ci__1.ciunit = 6;
    ci__1.cifmt = "(7e11.3)";
    s_wsfe(&ci__1);
    do_fio(&c__1, (char *)&am, (ftnlen)sizeof(real));
    i__3 = im;
    for (i__ = 1; i__ <= i__3; ++i__) {
	do_fio(&c__1, (char *)&tm[i__ - 1], (ftnlen)sizeof(real));
    }
    e_wsfe();
/*<        call atracer(model,elatc,elon,ncor) >*/
    atracer_(model, &elatc, &elon, &ncor, (ftnlen)256);
/*<        >*/
    surfread_(infile__, sigr, sigl, symbik, &nt, &nd, &depth, t, cr, ur, wvr, 
	    cl, ul, wvl, v, dvdz, ampr, ampl, ratio, qr, ql, i0, (ftnlen)255, 
	    (ftnlen)1, (ftnlen)1, (ftnlen)2);
/*<        PRINT*,nt,' input periods' >*/
    s_wsle(&io___113);
    do_lio(&c__3, &c__1, (char *)&nt, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " input periods", (ftnlen)14);
    e_wsle();
/* -----------Reading of input parameters--------------E */
/* ----------FINDING BASE SIZE ---------------------S */
/*<            do j=1,12 >*/
    for (j = 1; j <= 12; ++j) {
/*<            if(2**j.ge.npoints) go to 2    >*/
	if (pow_ii(&c__2, &j) >= npoints) {
	    goto L2;
	}
/*<            end do >*/
    }
/*< 2          nbase=2**j >*/
L2:
    nbase = pow_ii(&c__2, &j);
/*<            n2pow=j >*/
    n2pow = j;
/*<            PRINT*,'n2pow=',n2pow,'  nbase=',nbase,' npoints=',npoints >*/
    s_wsle(&io___117);
    do_lio(&c__9, &c__1, "n2pow=", (ftnlen)6);
    do_lio(&c__3, &c__1, (char *)&n2pow, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, "  nbase=", (ftnlen)8);
    do_lio(&c__3, &c__1, (char *)&nbase, (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, " npoints=", (ftnlen)9);
    do_lio(&c__3, &c__1, (char *)&npoints, (ftnlen)sizeof(integer));
    e_wsle();
/*<            df=1./dt/nbase >*/
    df = 1.f / dt / nbase;
/*<            f0=0.0  >*/
    f0 = 0.f;
/*<            nf=nbase/2 >*/
    nf = nbase / 2;
/* ----------FINDING BASE SIZE ---------------------E */
/*<            elatc=elatc*drad >*/
    elatc *= drad;
/*<            elonc = elon*drad >*/
    elonc = elon * drad;
/*<       do mm =1,nstai >*/
    i__3 = stn_1.nstai;
    for (mm = 1; mm <= i__3; ++mm) {
/*<            write(*,*) "Station: ",codi(mm) >*/
	s_wsle(&io___123);
	do_lio(&c__9, &c__1, "Station: ", (ftnlen)9);
	do_lio(&c__9, &c__1, stn_1.codi + (mm - 1 << 3), (ftnlen)8);
	e_wsle();
/*<            slat=fici(mm)*drad >*/
	slat = stn_1.fici[mm - 1] * drad;
/*<            slon=lami(mm)*drad >*/
	slon = stn_1.lami[mm - 1] * drad;
/*<            call AZIDL(SLAT,SLON,ELATC,ELONC,DEL,DELS,AZI_BACK,AZI) >*/
	azidl_(&slat, &slon, &elatc, &elonc, &del, &dels, &azi_back__, &azi);
/*<            dist=dels >*/
	dist = dels;
/*<            fi=AZI/drad >*/
	fi = azi / drad;
/*<            fi_back=AZI_BACK/drad >*/
	fi_back__ = azi_back__ / drad;
/*<        >*/
	s_wsle(&io___133);
	do_lio(&c__3, &c__1, (char *)&stn_1.nstai, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " STA= ", (ftnlen)6);
	do_lio(&c__9, &c__1, stn_1.codi + (mm - 1 << 3), (ftnlen)8);
	do_lio(&c__9, &c__1, " NET= ", (ftnlen)6);
	do_lio(&c__9, &c__1, neti + (mm - 1 << 3), (ftnlen)8);
	do_lio(&c__9, &c__1, " STLAT= ", (ftnlen)8);
	do_lio(&c__4, &c__1, (char *)&stn_1.figi[stn_1.nstai - 1], (ftnlen)
		sizeof(real));
	do_lio(&c__9, &c__1, " STLON= ", (ftnlen)8);
	do_lio(&c__4, &c__1, (char *)&stn_1.lami[stn_1.nstai - 1], (ftnlen)
		sizeof(real));
	do_lio(&c__9, &c__1, " DIST= ", (ftnlen)7);
	do_lio(&c__4, &c__1, (char *)&dist, (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, " AZI= ", (ftnlen)6);
	do_lio(&c__4, &c__1, (char *)&fi, (ftnlen)sizeof(real));
	do_lio(&c__9, &c__1, " BACK_AZI= ", (ftnlen)11);
	do_lio(&c__4, &c__1, (char *)&fi_back__, (ftnlen)sizeof(real));
	e_wsle();
/* --- */
/*<            cs=cos(AZI) >*/
	cs = cos(azi);
/*<            sc=sin(AZI) >*/
	sc = sin(azi);
/*<            cs_b=cos(AZI_BACK+pi) >*/
	cs_b__ = cos(azi_back__ + pi);
/*<            sc_b=sin(AZI_BACK+pi) >*/
	sc_b__ = sin(azi_back__ + pi);
/* -----------Reading of spectral data from SURF_DEEP output-----S */
/*<            do j=2,nt+1 >*/
	i__4 = nt + 1;
	for (j = 2; j <= i__4; ++j) {
/*<            fr(j)=1./t(j) >*/
	    fr[j - 1] = 1.f / t[j - 1];
/*<            wvar(j) = pi*2.0/cor(j-1,1,mm)*fr(j) >*/
	    wvar[j - 1] = pi * 2.f / trk_1.cor[j - 1 + ((mm << 1) + 1) * 500 
		    - 1501] * fr[j - 1];
/*          PRint*,t(j),ampr(j),ratio(j),qR(j),wvr(j) */
/*<            enddo >*/
	}
/*<            wvar(1) = wvr(1) >*/
	wvar[0] = wvr[0];
/*<            wvar(nt+2) = wvr(nt+2) >*/
	wvar[nt + 1] = wvr[nt + 1];
/*<            fr(1)=10.0 >*/
	fr[0] = 10.f;
/*<            fr(nt+2)=0.0 >*/
	fr[nt + 1] = 0.f;
/*<            qR(1)=qR(2) >*/
	qr[0] = qr[1];
/*<            qL(1)=qL(2) >*/
	ql[0] = ql[1];
/*<            qR(nt+2)=qR(nt+1) >*/
	qr[nt + 1] = qr[nt];
/*<            qL(nt+2)=qL(nt+1) >*/
	ql[nt + 1] = ql[nt];
/* -----------Reading of spectral data from SURF_DEEP  output-----E */
/* ----------CALCULATION OF SPECTRA--------------------------------------S */
/*<          write(*,*) 'NT ======',nt >*/
	s_wsle(&io___140);
	do_lio(&c__9, &c__1, "NT ======", (ftnlen)9);
	do_lio(&c__3, &c__1, (char *)&nt, (ftnlen)sizeof(integer));
	e_wsle();
/*<            Do j=2,nt+1 >*/
	i__4 = nt + 1;
	for (j = 2; j <= i__4; ++j) {
/* ----------Source term calculations-----------------S */
/*<            vu(1)=v(1,j) >*/
	    vu[0] = v[j * 3 - 3];
/*<            vu(2)=v(2,j) >*/
	    vu[1] = v[j * 3 - 2];
/*<            vu(3)=v(3,j) >*/
	    vu[2] = v[j * 3 - 1];
/*<            du(1)=dvdz(1,j) >*/
	    du[0] = dvdz[j * 3 - 3];
/*<            du(2)=dvdz(2,j) >*/
	    du[1] = dvdz[j * 3 - 2];
/*<            du(3)=dvdz(3,j) >*/
	    du[2] = dvdz[j * 3 - 1];
/*<            w=pi*2.0*fr(j) >*/
	    w = pi * 2.f * fr[j - 1];
/*<            wvn(1)=wvr(j) >*/
	    wvn[0] = wvr[j - 1];
/*<            wvn(2)=wvl(j) >*/
	    wvn[1] = wvl[j - 1];
/*<            if(its.eq.'F') call force(sigR,sigL,cs,sc,vu,br,bl) >*/
	    if (*(unsigned char *)its == 'F') {
		force_(sigr, sigl, &cs, &sc, vu, br, bl, (ftnlen)1, (ftnlen)1)
			;
	    }
/*<          if(its.ne.'F')then >*/
	    if (*(unsigned char *)its != 'F') {
/*<            if(its.eq.'Q'.or.its.eq.'C') step=cmplx(0.,-1./w) >*/
		if (*(unsigned char *)its == 'Q' || *(unsigned char *)its == 
			'C') {
		    r__1 = -1.f / w;
		    q__1.r = 0.f, q__1.i = r__1;
		    step.r = q__1.r, step.i = q__1.i;
		}
/*<            call source(sigR,sigL,cs,sc,wvn,vu,du,br,bl) >*/
		source_(sigr, sigl, &cs, &sc, wvn, vu, du, br, bl, (ftnlen)1, 
			(ftnlen)1);
/*<                       end if >*/
	    }
/*<            suml=(0.0,0.0) >*/
	    suml.r = 0.f, suml.i = 0.f;
/*<            sumr=(0.0,0.0) >*/
	    sumr.r = 0.f, sumr.i = 0.f;
/*<           do m=1,im >*/
	    i__5 = im;
	    for (m = 1; m <= i__5; ++m) {
/*<            sumr= sumr+tm(m)*br(m)*step >*/
		i__6 = m - 1;
		i__7 = m - 1;
		q__3.r = tm[i__6] * br[i__7].r, q__3.i = tm[i__6] * br[i__7]
			.i;
		q__2.r = q__3.r * step.r - q__3.i * step.i, q__2.i = q__3.r * 
			step.i + q__3.i * step.r;
		q__1.r = sumr.r + q__2.r, q__1.i = sumr.i + q__2.i;
		sumr.r = q__1.r, sumr.i = q__1.i;
/*<            suml= suml+tm(m)*bl(m)*step >*/
		i__6 = m - 1;
		i__7 = m - 1;
		q__3.r = tm[i__6] * bl[i__7].r, q__3.i = tm[i__6] * bl[i__7]
			.i;
		q__2.r = q__3.r * step.r - q__3.i * step.i, q__2.i = q__3.r * 
			step.i + q__3.i * step.r;
		q__1.r = suml.r + q__2.r, q__1.i = suml.i + q__2.i;
		suml.r = q__1.r, suml.i = q__1.i;
/*<                   enddo >*/
	    }
/*<            al(j)=suml*ampl(j)*aM*const1 >*/
	    i__5 = j - 1;
	    i__6 = j - 1;
	    q__3.r = ampl[i__6] * suml.r, q__3.i = ampl[i__6] * suml.i;
	    q__2.r = am * q__3.r, q__2.i = am * q__3.i;
	    q__1.r = const1 * q__2.r, q__1.i = const1 * q__2.i;
	    al[i__5].r = q__1.r, al[i__5].i = q__1.i;
/* MB           spread=1./sqrt(R0*sin(DEL)*wvn(1)) */
/*<            spread=1./sqrt(R0*sin(DEL)) >*/
	    spread = 1.f / sqrt(r0 * sin(del));
/*          if(j.eq.2)PRint*,' SPREAD=',spread*sqrt(wvn(1)), ita */
/* MB        az(j)=sumr*ampr(j)*aM*const1*conjg(cor(j-1)) */
/*<            az(j)=sumr*aM*const1*cor(j-1,2,mm)/(2.0*cr(j)*sqrt(ur(j)*I0(j)))/sqrt(2.0*pi) >*/
	    i__5 = j - 1;
	    q__5.r = am * sumr.r, q__5.i = am * sumr.i;
	    q__4.r = const1 * q__5.r, q__4.i = const1 * q__5.i;
	    i__6 = j - 1 + ((mm << 1) + 2) * 500 - 1501;
	    q__3.r = trk_1.cor[i__6] * q__4.r, q__3.i = trk_1.cor[i__6] * 
		    q__4.i;
	    r__1 = cr[j - 1] * 2.f * sqrt(ur[j - 1] * i0[j - 1]);
	    q__2.r = q__3.r / r__1, q__2.i = q__3.i / r__1;
	    r__2 = sqrt(pi * 2.f);
	    q__1.r = q__2.r / r__2, q__1.i = q__2.i / r__2;
	    az[i__5].r = q__1.r, az[i__5].i = q__1.i;
/* MB      write(*,*)t(j),cr(j),ur(j),wvr(j),ampr(j),cor(j-1,1),cor(j-1,2) */
/* MB        az(j)=sumr*ampr(j)*aM*const1 */
/* xx        WRite(33,*)t(j),cabs(az(j))*spread */
/*<            ah(j)=az(j)*ratio(j) >*/
	    i__5 = j - 1;
	    i__6 = j - 1;
	    i__7 = j - 1;
	    q__1.r = ratio[i__7] * az[i__6].r, q__1.i = ratio[i__7] * az[i__6]
		    .i;
	    ah[i__5].r = q__1.r, ah[i__5].i = q__1.i;
/*<                        end do >*/
	}
/* ----------Source term calculations-----------------E */
/* ----------INTERPOLATION OF WAVENUMBERS AND Q-factors--S */
/*<            if (sigR.eq.'+')then >*/
	if (*(unsigned char *)sigr == '+') {
/* MB        call intpol(fr,wvr,nt+2,f0,df,nf,wr,ierr) */
/*<            call intpol(fr,wvar,nt+2,f0,df,nf,wr,ierr) >*/
	    i__4 = nt + 2;
	    intpol_(fr, wvar, &i__4, &f0, &df, &nf, wr, &ierr);
/*<            if(ierr.ne.0)STOP'WRONG INTERPOLATION: RAYLEIGH WAVENUMBER' >*/
	    if (ierr != 0) {
		s_stop("WRONG INTERPOLATION: RAYLEIGH WAVENUMBER", (ftnlen)40)
			;
	    }
/*<            call intpol(fr,qR,nt+2,f0,df,nf,qR_int,ierr) >*/
	    i__4 = nt + 2;
	    intpol_(fr, qr, &i__4, &f0, &df, &nf, qr_int__, &ierr);
/*<            if(ierr.ne.0)STOP'WRONG INTERPOLATION: RAYLEIGH ATTENUATION' >*/
	    if (ierr != 0) {
		s_stop("WRONG INTERPOLATION: RAYLEIGH ATTENUATION", (ftnlen)
			41);
	    }
/*<            call intpol(fr,ur,nt+2,f0,df,nf,ur_int,ierr) >*/
	    i__4 = nt + 2;
	    intpol_(fr, ur, &i__4, &f0, &df, &nf, ur_int__, &ierr);
/*<            if(ierr.ne.0)STOP'WRONG INTERPOLATION: RAYLEIGH GR.VELOCITY' >*/
	    if (ierr != 0) {
		s_stop("WRONG INTERPOLATION: RAYLEIGH GR.VELOCITY", (ftnlen)
			41);
	    }
/*<            DO mmm=1,nf >*/
	    i__4 = nf;
	    for (mmm = 1; mmm <= i__4; ++mmm) {
/*<            if(f0+df*(mmm-1).gt.fr(2))qR_int(mmm)=qR_int(mmm-1) >*/
		if (f0 + df * (mmm - 1) > fr[1]) {
		    qr_int__[mmm - 1] = qr_int__[mmm - 2];
		}
/*<            enddo >*/
	    }
/*<                        endif >*/
	}
/*<            if (sigL.eq.'+')then >*/
	if (*(unsigned char *)sigl == '+') {
/*<            call intpol(fr,wvl,nt+2,f0,df,nf,wl,ierr) >*/
	    i__4 = nt + 2;
	    intpol_(fr, wvl, &i__4, &f0, &df, &nf, wl, &ierr);
/*<            if(ierr.ne.0)STOP'WRONG INTERPOLATION: LOVE WVENUMBER' >*/
	    if (ierr != 0) {
		s_stop("WRONG INTERPOLATION: LOVE WVENUMBER", (ftnlen)35);
	    }
/*<            call intpol(fr,qL,nt+2,f0,df,nf,qL_int,ierr) >*/
	    i__4 = nt + 2;
	    intpol_(fr, ql, &i__4, &f0, &df, &nf, ql_int__, &ierr);
/*<            if(ierr.ne.0)STOP'WRONG INTERPOLATION: LOVE ATTENUATION' >*/
	    if (ierr != 0) {
		s_stop("WRONG INTERPOLATION: LOVE ATTENUATION", (ftnlen)37);
	    }
/*<            call intpol(fr,ul,nt+2,f0,df,nf,ul_int,ierr) >*/
	    i__4 = nt + 2;
	    intpol_(fr, ul, &i__4, &f0, &df, &nf, ul_int__, &ierr);
/*<            if(ierr.ne.0)STOP'WRONG INTERPOLATION: LOVE GR.VELOCITY' >*/
	    if (ierr != 0) {
		s_stop("WRONG INTERPOLATION: LOVE GR.VELOCITY", (ftnlen)37);
	    }
/*<            DO mmm=1,nf >*/
	    i__4 = nf;
	    for (mmm = 1; mmm <= i__4; ++mmm) {
/*<            if(f0+df*(mmm-1).gt.fr(2))then >*/
		if (f0 + df * (mmm - 1) > fr[1]) {
/*<            qL_int(mmm)=qL_int(mmm-1) >*/
		    ql_int__[mmm - 1] = ql_int__[mmm - 2];
/*<            uL_int(mmm)=uL_int(mmm-1) >*/
		    ul_int__[mmm - 1] = ul_int__[mmm - 2];
/*<                                      endif >*/
		}
/*<            enddo >*/
	    }
/*<                          endif >*/
	}
/* ----------INTERPOLATION OF WAVENUMBERS AND Q-factors--E */
/* ----------INTERPOLATION OF SPECTRAL AMPLITUDES-------------S */
/*<        >*/
	if (*(unsigned char *)sigr == '+') {
	    calspecr_(&nt, &nf, &f0, &df, fr, az, ah, srev, simv, srer, simr, 
		    &tmin, &tmax, &iq);
	}
/*<        >*/
	if (*(unsigned char *)sigl == '+') {
	    calspecl_(&nt, &nf, &f0, &df, fr, al, sret, simt, &tmin, &tmax, &
		    iq);
	}
/* ----------INTERPOLATION OF SPECTRAL AMPLITUDES-------------E */
/* ----------CALCULATION OF SPECTRA--------------------------------------E */
/* ----------CALCULATION OF SEISMOGRAMS--------------------S */
/*<            NDIST=1 >*/
	ndist = 1;
/*<       call namer(dels,ndist,n1,current,outseism,outtit,outspec,sigR,sigL,nk7) >*/
	namer_(&dels, &ndist, &n1, current, outseism, outtit, outspec, sigr, 
		sigl, &nk7, (ftnlen)255, (ftnlen)80, (ftnlen)80, (ftnlen)80, (
		ftnlen)1, (ftnlen)1);
/* ----------output filenames are defined--------->>>>>> */
/* ----------DISTANCE LOOP-----------------------S */
/*<            DO l=1,ndist >*/
	i__4 = ndist;
	for (l = 1; l <= i__4; ++l) {
/*<            lq=(l-1)*3 >*/
	    lq = (l - 1) * 3;
/* --------------calculating everything and writing spectra----S */
/*<            if (sigR.eq.'+') THEN >*/
	    if (*(unsigned char *)sigr == '+') {
/*<        >*/
		syn_(wr, &dist, &dt, &f0, &df, srev, simv, ur_int__, qr_int__,
			 seismz, &tstart, &vmax, &n2pow, &npoints, &
			key_compr__, &fix_vel__, amp, t_curr__, &tmin, &tmax, 
			&k_spec__);
/*<       open(41,file=outspec(lq+1)(1:nk7+1)) >*/
		o__1.oerr = 0;
		o__1.ounit = 41;
		o__1.ofnmlen = nk7 + 1;
		o__1.ofnm = outspec + lq * 80;
		o__1.orl = 0;
		o__1.osta = 0;
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
/*<       do kl=2,k_spec >*/
		i__5 = k_spec__;
		for (kl = 2; kl <= i__5; ++kl) {
/*<       write(41,*)T_curr(kl),amp(kl) >*/
		    s_wsle(&io___180);
		    do_lio(&c__4, &c__1, (char *)&t_curr__[kl - 1], (ftnlen)
			    sizeof(real));
		    do_lio(&c__4, &c__1, (char *)&amp[kl - 1], (ftnlen)sizeof(
			    real));
		    e_wsle();
/*<       enddo >*/
		}
/*<       close(41) >*/
		cl__1.cerr = 0;
		cl__1.cunit = 41;
		cl__1.csta = 0;
		f_clos(&cl__1);
/*<        >*/
		syn_(wr, &dist, &dt, &f0, &df, srer, simr, ur_int__, qr_int__,
			 seismr, &tstart, &vmax, &n2pow, &npoints, &
			key_compr__, &fix_vel__, amp, t_curr__, &tmin, &tmax, 
			&k_spec__);
/*<       open(41,file=outspec(lq+2)(1:nk7+1)) >*/
		o__1.oerr = 0;
		o__1.ounit = 41;
		o__1.ofnmlen = nk7 + 1;
		o__1.ofnm = outspec + (lq + 1) * 80;
		o__1.orl = 0;
		o__1.osta = 0;
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
/*<       do kl=2,k_spec >*/
		i__5 = k_spec__;
		for (kl = 2; kl <= i__5; ++kl) {
/*<       write(41,*)T_curr(kl),amp(kl) >*/
		    s_wsle(&io___182);
		    do_lio(&c__4, &c__1, (char *)&t_curr__[kl - 1], (ftnlen)
			    sizeof(real));
		    do_lio(&c__4, &c__1, (char *)&amp[kl - 1], (ftnlen)sizeof(
			    real));
		    e_wsle();
/*<       enddo >*/
		}
/*<       close(41) >*/
		cl__1.cerr = 0;
		cl__1.cunit = 41;
		cl__1.csta = 0;
		f_clos(&cl__1);
/*<             fnam11='w/'//codi(mm)(1:lnblnk(codi(mm)))//'.'//outseism(lq+1)(1:nk7) >*/
/* Writing concatenation */
		i__8[0] = 2, a__3[0] = "w/";
		i__8[1] = lnblnk_(stn_1.codi + (mm - 1 << 3), (ftnlen)8), 
			a__3[1] = stn_1.codi + (mm - 1 << 3);
		i__8[2] = 1, a__3[2] = ".";
		i__8[3] = nk7, a__3[3] = outseism + lq * 80;
		s_cat(fnam11, a__3, i__8, &c__4, (ftnlen)256);
/*<             open(11,file=fnam11) >*/
		o__1.oerr = 0;
		o__1.ounit = 11;
		o__1.ofnmlen = 256;
		o__1.ofnm = fnam11;
		o__1.orl = 0;
		o__1.osta = 0;
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
/*<            ELSEIF (sigL.eq.'+') THEN >*/
	    } else if (*(unsigned char *)sigl == '+') {
/*<        >*/
		syn_(wl, &dist, &dt, &f0, &df, sret, simt, ul_int__, ql_int__,
			 seismt, &tstart, &vmax, &n2pow, &npoints, &
			key_compr__, &fix_vel__, amp, t_curr__, &tmin, &tmax, 
			&k_spec__);
/*<       open(41,file=outspec(lq+3)(1:nk7+1)) >*/
		o__1.oerr = 0;
		o__1.ounit = 41;
		o__1.ofnmlen = nk7 + 1;
		o__1.ofnm = outspec + (lq + 2) * 80;
		o__1.orl = 0;
		o__1.osta = 0;
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
/*<       do kl=2,k_spec >*/
		i__5 = k_spec__;
		for (kl = 2; kl <= i__5; ++kl) {
/*<       write(41,*)T_curr(kl),amp(kl) >*/
		    s_wsle(&io___185);
		    do_lio(&c__4, &c__1, (char *)&t_curr__[kl - 1], (ftnlen)
			    sizeof(real));
		    do_lio(&c__4, &c__1, (char *)&amp[kl - 1], (ftnlen)sizeof(
			    real));
		    e_wsle();
/*<       enddo >*/
		}
/*<       close(41) >*/
		cl__1.cerr = 0;
		cl__1.cunit = 41;
		cl__1.csta = 0;
		f_clos(&cl__1);
/*<            ENDIF >*/
	    }
/* --------------calculating everything and writing spectra----E */
/* --------------writing seismograms----S */
/*<           fnam12='w/'//codi(mm)(1:lnblnk(codi(mm)))//'.'//outseism(lq+2)(1:nk7) >*/
/* Writing concatenation */
	    i__8[0] = 2, a__3[0] = "w/";
	    i__8[1] = lnblnk_(stn_1.codi + (mm - 1 << 3), (ftnlen)8), a__3[1] 
		    = stn_1.codi + (mm - 1 << 3);
	    i__8[2] = 1, a__3[2] = ".";
	    i__8[3] = nk7, a__3[3] = outseism + (lq + 1) * 80;
	    s_cat(fnam12, a__3, i__8, &c__4, (ftnlen)256);
/*<             open(12,file=fnam12) >*/
	    o__1.oerr = 0;
	    o__1.ounit = 12;
	    o__1.ofnmlen = 256;
	    o__1.ofnm = fnam12;
	    o__1.orl = 0;
	    o__1.osta = 0;
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
/*<           fnam13='w/'//codi(mm)(1:lnblnk(codi(mm)))//'.'//outseism(lq+3)(1:nk7) >*/
/* Writing concatenation */
	    i__8[0] = 2, a__3[0] = "w/";
	    i__8[1] = lnblnk_(stn_1.codi + (mm - 1 << 3), (ftnlen)8), a__3[1] 
		    = stn_1.codi + (mm - 1 << 3);
	    i__8[2] = 1, a__3[2] = ".";
	    i__8[3] = nk7, a__3[3] = outseism + (lq + 2) * 80;
	    s_cat(fnam13, a__3, i__8, &c__4, (ftnlen)256);
/*<             open(13,file=fnam13) >*/
	    o__1.oerr = 0;
	    o__1.ounit = 13;
	    o__1.ofnmlen = 256;
	    o__1.ofnm = fnam13;
	    o__1.orl = 0;
	    o__1.osta = 0;
	    o__1.oacc = 0;
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    f_open(&o__1);
/*<            do ll=1,npoints >*/
	    i__5 = npoints;
	    for (ll = 1; ll <= i__5; ++ll) {
/*<            seismn(ll)=0.0 >*/
		seismn[ll - 1] = 0.f;
/*<            seisme(ll)=0.0 >*/
		seisme[ll - 1] = 0.f;
/*<            if (sigR.eq.'+') then >*/
		if (*(unsigned char *)sigr == '+') {
/*          seismz(ll)=-seismz(ll) */
/* -----------Positive 'Z' direction is 'Down'-----!!!!!!!!!!!!!!!!!!!!!! */
/* -----------horizontal components----S */
/*<            seismn(ll)=seismr(ll)*cs_b >*/
		    seismn[ll - 1] = seismr[ll - 1] * cs_b__;
/*<            seisme(ll)=seismr(ll)*sc_b >*/
		    seisme[ll - 1] = seismr[ll - 1] * sc_b__;
/*<            elseif (sigL.eq.'+') then >*/
		} else if (*(unsigned char *)sigl == '+') {
/*<            seismn(ll)=seismn(ll) -seismt(ll)*sc_b >*/
		    seismn[ll - 1] -= seismt[ll - 1] * sc_b__;
/*<            seisme(ll)=seisme(ll) +seismt(ll)*cs_b >*/
		    seisme[ll - 1] += seismt[ll - 1] * cs_b__;
/*<                             endif >*/
		}
/* -----------horizontal components----E */
/* -----------vertical   component-----S */
/*<            if (sigR.eq.'+')write(11,*)tstart+dt*(ll-1),seismz(ll) >*/
		if (*(unsigned char *)sigr == '+') {
		    s_wsle(&io___191);
		    r__1 = tstart + dt * (ll - 1);
		    do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(real));
		    do_lio(&c__4, &c__1, (char *)&seismz[ll - 1], (ftnlen)
			    sizeof(real));
		    e_wsle();
		}
/*<            write(12,*)tstart+dt*(ll-1),seismn(ll) >*/
		s_wsle(&io___192);
		r__1 = tstart + dt * (ll - 1);
		do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&seismn[ll - 1], (ftnlen)sizeof(
			real));
		e_wsle();
/*<            write(13,*)tstart+dt*(ll-1),seisme(ll) >*/
		s_wsle(&io___193);
		r__1 = tstart + dt * (ll - 1);
		do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(real));
		do_lio(&c__4, &c__1, (char *)&seisme[ll - 1], (ftnlen)sizeof(
			real));
		e_wsle();
/* -----------vertical   component-----E */
/*<            enddo >*/
	    }
/* --------------writing seismograms----E */
/*<            close(11) >*/
	    cl__1.cerr = 0;
	    cl__1.cunit = 11;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
/*<            close(12) >*/
	    cl__1.cerr = 0;
	    cl__1.cunit = 12;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
/*<           close(13) >*/
	    cl__1.cerr = 0;
	    cl__1.cunit = 13;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
/* ------------ write sac files --------- */
/*<         fnam11 = 'sac/'//codi(mm)(1:lnblnk(codi(mm)))//'.BHZ.SAC' >*/
/* Writing concatenation */
	    i__9[0] = 4, a__4[0] = "sac/";
	    i__9[1] = lnblnk_(stn_1.codi + (mm - 1 << 3), (ftnlen)8), a__4[1] 
		    = stn_1.codi + (mm - 1 << 3);
	    i__9[2] = 8, a__4[2] = ".BHZ.SAC";
	    s_cat(fnam11, a__4, i__9, &c__3, (ftnlen)256);
/*<         fnam12 = 'sac/'//codi(mm)(1:lnblnk(codi(mm)))//'.BHN.SAC' >*/
/* Writing concatenation */
	    i__9[0] = 4, a__4[0] = "sac/";
	    i__9[1] = lnblnk_(stn_1.codi + (mm - 1 << 3), (ftnlen)8), a__4[1] 
		    = stn_1.codi + (mm - 1 << 3);
	    i__9[2] = 8, a__4[2] = ".BHN.SAC";
	    s_cat(fnam12, a__4, i__9, &c__3, (ftnlen)256);
/*<         fnam13 = 'sac/'//codi(mm)(1:lnblnk(codi(mm)))//'.BHE.SAC' >*/
/* Writing concatenation */
	    i__9[0] = 4, a__4[0] = "sac/";
	    i__9[1] = lnblnk_(stn_1.codi + (mm - 1 << 3), (ftnlen)8), a__4[1] 
		    = stn_1.codi + (mm - 1 << 3);
	    i__9[2] = 8, a__4[2] = ".BHE.SAC";
	    s_cat(fnam13, a__4, i__9, &c__3, (ftnlen)256);
/*<         chn = 'BHZ' >*/
	    s_copy(chn, "BHZ", (ftnlen)8, (ftnlen)3);
/*<        >*/
	    wsac_(fnam11, &elat, &elon, stn_1.codi + (mm - 1 << 3), &
		    stn_1.figi[mm - 1], &stn_1.lami[mm - 1], &dt, chn, neti + 
		    (mm - 1 << 3), &npoints, seismz, (ftnlen)256, (ftnlen)8, (
		    ftnlen)8, (ftnlen)8);
/*<         chn = 'BHN' >*/
	    s_copy(chn, "BHN", (ftnlen)8, (ftnlen)3);
/*<        >*/
	    wsac_(fnam12, &elat, &elon, stn_1.codi + (mm - 1 << 3), &
		    stn_1.figi[mm - 1], &stn_1.lami[mm - 1], &dt, chn, neti + 
		    (mm - 1 << 3), &npoints, seismn, (ftnlen)256, (ftnlen)8, (
		    ftnlen)8, (ftnlen)8);
/*<         chn = 'BHE' >*/
	    s_copy(chn, "BHE", (ftnlen)8, (ftnlen)3);
/*<        >*/
	    wsac_(fnam13, &elat, &elon, stn_1.codi + (mm - 1 << 3), &
		    stn_1.figi[mm - 1], &stn_1.lami[mm - 1], &dt, chn, neti + 
		    (mm - 1 << 3), &npoints, seisme, (ftnlen)256, (ftnlen)8, (
		    ftnlen)8, (ftnlen)8);
/*<         fnam11 = 'mv *.R.spec '//'SPEC/'//codi(mm)(1:lnblnk(codi(mm)))//'.R.sp' >*/
/* Writing concatenation */
	    i__9[0] = 17, a__4[0] = "mv *.R.spec SPEC/";
	    i__9[1] = lnblnk_(stn_1.codi + (mm - 1 << 3), (ftnlen)8), a__4[1] 
		    = stn_1.codi + (mm - 1 << 3);
	    i__9[2] = 5, a__4[2] = ".R.sp";
	    s_cat(fnam11, a__4, i__9, &c__3, (ftnlen)256);
/*<         call system(fnam11) >*/
	    system_(fnam11, (ftnlen)256);
/*<         fnam11 = 'mv *.Z.spec '//'SPEC/'//codi(mm)(1:lnblnk(codi(mm)))//'.Z.sp' >*/
/* Writing concatenation */
	    i__9[0] = 17, a__4[0] = "mv *.Z.spec SPEC/";
	    i__9[1] = lnblnk_(stn_1.codi + (mm - 1 << 3), (ftnlen)8), a__4[1] 
		    = stn_1.codi + (mm - 1 << 3);
	    i__9[2] = 5, a__4[2] = ".Z.sp";
	    s_cat(fnam11, a__4, i__9, &c__3, (ftnlen)256);
/*<         call system(fnam11) >*/
	    system_(fnam11, (ftnlen)256);
/*<            end DO >*/
	}
/* ----------DISTANCE LOOP------ ----------------E */
/* ----------CALCULATION OF SEISMOGRAMS--------------------S */
/*<         enddo >*/
    }
/*<            stop >*/
    s_stop("", (ftnlen)0);
/*< 1000  format(5x,"aM",8x,"mxx",8x,"myy",8x,"mzz",8x,"mxy",8x,"mxz",8x,"myz") >*/
/*<            END >*/
    return 0;
} /* MAIN__ */

