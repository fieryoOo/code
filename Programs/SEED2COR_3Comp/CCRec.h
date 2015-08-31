#ifndef CCREC_H
#define CCREC_H
#include <iostream>
#include <vector>
#include <string>


struct outCCRec
{
    std::vector < SacRec > CCRec;
    int sacRecSize, day, cclistNo;
    std::string monthdir, sta1,sta2,COR;
    std::vector <std::string> outfname;
     std::vector < std::string >  chanL1, chanL2;
    outCCRec();
    ~outCCRec();
    outCCRec(std::vector < SacRec > & INRec, int CCSIZE, int DAY, int LISTNo, std::string MDIR,
             std::string STA1,  std::vector < std::string >  ChanL1, std::string STA2, std::vector < std::string >  ChanL2, std::vector <std::string> & OUTFNAME, std::string CORDIR);
};

struct hole_rec
{
    int rec_b, rec_e;
    hole_rec();
    ~hole_rec();
    hole_rec(int BEG, int END);
};


#endif
