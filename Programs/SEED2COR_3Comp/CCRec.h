#ifndef CCREC_H
#define CCREC_H
#include <iostream>
#include <vector>
#include <string>
//#include <fftw3.h>

struct CC_todo
{
    std::vector < std::string >  infname1, infname2;
    std::string monthdir, sta1, sta2, chan1, chan2;
    int day;
    /* ------------------------------ con/destructors and operators ------------------------------ */
    /* constructors */
    CC_todo( );
    CC_todo(const std::vector <std::string> & IN1, const std::vector <std::string> &IN2, const std::string & MDIR, const std::string & STA1,
            const std::string & STA2, const int & DD);
    /* operators */

    /* destructor */
    ~CC_todo();
};

class CC_output
{
public:
    std::string  sta1,chan1,sta2,chan2, monthdir;
    std::vector< std::string > outfname;
    std::vector< int > day;
public:
    //bool checkfile(const std::string  & STAN1, const std::string  & CHAN1,const std::string  & STAN2, const std::string  & CHAN2, const int & YY, const int & MM);
    //bool checkfile_all(const std::string  & STAN1, const std::string  & CHAN1,const std::string  & STAN2, const std::string  & CHAN2);
    /* ------------------------------ con/destructors and operators ------------------------------ */
    /* constructors */
    CC_output( );
    CC_output(const std::string & OUTname, const std::string & MDIR, const std::string  & STAN1, const std::string  & CHAN1,
              const std::string  & STAN2, const std::string  & CHAN2, const int & DAY);
    /* operators */

    /* destructor */
    ~CC_output();
};
struct CC_pairs
{
    std::string sta1,sta2,chan1,chan2;
    CC_pairs();
    ~CC_pairs();
    CC_pairs(const std::string  & STAN1, const std::string  & CHAN1,
             const std::string  & STAN2, const std::string  & CHAN2);
};
struct sta_pairs
{
    std::string sta1,sta2;
    sta_pairs();
    ~sta_pairs();
    sta_pairs(const std::string  & STAN1,  const std::string  & STAN2);
};

struct Station
{
    std::string sta;
    int CCflag;
    bool checkdoCC(const Station & STA );
    Station();
    ~Station();
    Station(const std::string  & STAN1,  const int & CCFLAG);
};

struct outCCRec
{
    std::vector < SacRec > CCRec;
    int sacRecSize, day, cclistNo;
    std::string monthdir, sta1,sta2,COR;
    std::vector <std::string> outfname;

    outCCRec();
    ~outCCRec();
    outCCRec(std::vector < SacRec > & INRec, int CCSIZE, int DAY, int LISTNo, std::string MDIR,
             std::string STA1, std::string STA2, std::vector <std::string> & OUTFNAME, std::string CORDIR);
};

struct hole_rec
{
    int rec_b, rec_e;
    hole_rec();
    ~hole_rec();
    hole_rec(int BEG, int END);
};


#endif
