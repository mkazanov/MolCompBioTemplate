//
//  allmotifs.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 14/04/2019.
//  Copyright © 2019 Marat Kazanov. All rights reserved.
//

#ifndef allmotifs_hpp
#define allmotifs_hpp

#include <stdio.h>
#include <string>
#include "mutation.hpp"
#include "replicationtime.hpp"

using namespace std;

class CRTMapKey {
public:
    string motif;
    int RTbin;
    int RTstrand;
    char mutbase;
    CRTMapKey(string motif_, int RTbin_, int RTstrand_, char mutbase_);
    bool operator< (const CRTMapKey &right) const
    {
        if(motif < right.motif)
            return true;
        else if (motif > right.motif)
            return false;
        else
        {
            if(RTbin < right.RTbin)
                return true;
            else if (RTbin > right.RTbin)
                return false;
            else
            {
                if(RTstrand < right.RTstrand)
                    return true;
                else if (RTstrand > right.RTstrand)
                    return false;
                else
                {
                    if(mutbase < right.mutbase)
                        return true;
                    else if(mutbase > right.mutbase)
                        return false;
                    else
                        return false;
                }
            }
        }
    }
};

class CRTMapKeyNoMut {
public:
    string motif;
    int RTbin;
    int RTstrand;
    CRTMapKeyNoMut(string motif_, int RTbin_, int RTstrand_);
    bool operator< (const CRTMapKeyNoMut &right) const
    {
        if(motif < right.motif)
            return true;
        else if (motif > right.motif)
            return false;
        else
        {
            if(RTbin < right.RTbin)
                return true;
            else if (RTbin > right.RTbin)
                return false;
            else
            {
                if(RTstrand < right.RTstrand)
                    return true;
                else if (RTstrand > right.RTstrand)
                    return false;
                else
                    return false;
            }
        }
    }
};


class CAllMotifs{
public:
    void AnalysisRTmutations(CMutations& muts, string dirpath, CReplicationTiming& rti, CHumanGenome* phuman, string cancer, string sample);
    void AnalysisRTtargets(string dirpath, CReplicationTiming& rti, CHumanGenome* phuman);
};

#endif /* allmotifs_hpp */
