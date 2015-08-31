#ifndef STAPAIR_H
#define STAPAIR_H
#include <iostream>
#include <vector>
#include <string>
#include "InfoLists.h"
#include "SacRec.h"
#include "CCDatabase.h"
#include "CCRec.h"

class StaPair
{
public:
    StaInfo sta1, sta2;
    std::string month;
    std::vector < std::string >  chanL1, chanL2;
    int day;
    double dist, az, baz;
public:
    /* constructors */
    StaPair();
    StaPair(const StaInfo Sta1, const StaInfo Sta2);
//    StaPair(const StaInfo Sta1, const StaInfo Sta2, const std::vector < std::string >  ChanL1,const std::vector < std::string > ChanL2);
    /* operators */
    bool GetDistAzBaz();
    void SetDate( const std::string &Month, const int & Day );
    void SetChanLst(const std::vector < std::string >  ChanL1,const std::vector < std::string > ChanL2);

    void staPair2CC(std::vector < SacRec > & D_sacV, std::vector < std::string> & FoutnameD,
                             std::vector <std::string> & M_foutname, CCDatabase & cdb, int LNo );
    /* destructor */
    ~StaPair();
};

#endif
