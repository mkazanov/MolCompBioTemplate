//
//  replicationtime.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 03/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef replicationtime_hpp
#define replicationtime_hpp

#include <set>
#include <string>
#include <map>
#include "ghuman.hpp"
#include "mutation.hpp"

#define RT_NULL 999999
#define RT_NULLBIN_NOVALUE -1
#define RT_NULLBIN_NOBIN -2
#define RT_NULLBIN_CNT  2
#define CHR_NULL -1
#define STRAND_LEADING 1
#define STRAND_LAGGING 0
#define STRAND_NULL -1


using namespace std;

class CRTFileFormat{
public:
    char delimiter;
    int chrNo;
    int startPosNo;
    int endPosNo;
    int RTNo;
    int isHeader;
    CRTFileFormat(char delimiter_,
                   int chrNo_,
                   int startPosNo_,
                   int endPosNo_,
                   int RTNo_,
                   int isHeader);
};


class CReplicationTime {
public:
    //int chrNum;
    unsigned long startpos;
    unsigned long endpos;
    double RTvalue;
    mutable short isForward;
    CReplicationTime(string startpos_, string endpos_, string RTvalue_);
    CReplicationTime(unsigned long pos_);
    CReplicationTime(){};
    bool isRTnull();
    bool operator< (const CReplicationTime &right) const
    {
        if (startpos < right.startpos)
            return true;
        else if (startpos == right.startpos)
            return endpos < right.endpos;
        else
            return false;
    }
};

class CRTBin {
public:
    int binNum;
    double RTleft;
    double RTright;
    CRTBin(int binNum_, double RTleft_, double RTright_);
};

class CReplicationTiming {
    set<CReplicationTime>* RTs;
    vector<CRTFileFormat> fileFormat;
public:
    CReplicationTiming();
    vector<CRTBin> bins;
    void LoadReplicationTiming(string path, int isHeader);
    void LoadRTBins(string path);
    int GetRT(int chrNum, unsigned long pos, double& RTvalue);
    int GetRT(int chrNum, unsigned long pos, CReplicationTime& rtobj);
    int GetRTBin(double RTvalue, vector<CRTBin> bins);
    int GetRTBin(int chrNum, unsigned long pos, vector<CRTBin> bins);
    void GetRTBinStrand(map<string,CReplicationTiming*> rtmap,const CMutation& mut, int& bin, int& strand);
    void ReplicationStrand();
    void SaveToFile(string path);
    int CalculateMotifinRTBins(set<string> motifs, string OUT_PATH, CHumanGenome* phuman_ = NULL);
    static int oppositeStrand(int strand)
    {
        if(strand == STRAND_LAGGING)
            return(STRAND_LEADING);
        else if(strand == STRAND_LEADING)
            return(STRAND_LAGGING);
        else
            return(STRAND_NULL);
    }
    void CalculateRTbinsWidthByMotif(int binsNum, CHumanGenome* phuman, string outputPath);
    void CalculateRTbinsWidth(int binsNum, CHumanGenome* phuman, string outputPath);
};

#endif /* replicationtime_hpp */
