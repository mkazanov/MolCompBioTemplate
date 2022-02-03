//
//  g_human.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 03/04/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef ghuman_hpp
#define ghuman_hpp

#include <stdio.h>
#include "biogenome.hpp"
#include "bioread.hpp"
#include "service.h"
#include <assert.h>

class CHumanGenome : public CGenome
{
public:
    void InitializeHuman(string version, string pathPrefix, string pathSuffix, string format, int isLoadFiles = 1);
    static int GetChrNum(string chr)
    {
        int ret;
        size_t chrpos;
        string chr_;
        
        chrpos = chr.find("chr");
        if(chrpos != string::npos && chrpos == 0)
            chr_ = chr.substr(3);
        else
            chr_ = chr;
        
        if (chr_ == "X")
            ret = 22;
        else if (chr_ == "Y")
            ret = 23;
        else
            ret = str2i(chr_) - 1;
        
        if(ret >= 0 && ret < 24)
            return(ret);
        else
            return(-1);
    }
    static int GetChrFromNCBIName(string chr)
    {
        int ret;
        
        if (chr == "NC_000001.10")
            ret = 0;
        else if (chr == "NC_000002.11")
            ret = 1;
        else if (chr == "NC_000003.11")
            ret = 2;
        else if (chr == "NC_000004.11")
            ret = 3;
        else if (chr == "NC_000005.9")
            ret = 4;
        else if (chr == "NC_000006.11")
            ret = 5;
        else if (chr == "NC_000007.13")
            ret = 6;
        else if (chr == "NC_000008.10")
            ret = 7;
        else if (chr == "NC_000009.11")
            ret = 8;
        else if (chr == "NC_000010.10")
            ret = 9;
        else if (chr == "NC_000011.9")
            ret = 10;
        else if (chr == "NC_000012.11")
            ret = 11;
        else if (chr == "NC_000013.10")
            ret = 12;
        else if (chr == "NC_000014.8")
            ret = 13;
        else if (chr == "NC_000015.9")
            ret = 14;
        else if (chr == "NC_000016.9")
            ret = 15;
        else if (chr == "NC_000017.10")
            ret = 16;
        else if (chr == "NC_000018.9")
            ret = 17;
        else if (chr == "NC_000019.9")
            ret = 18;
        else if (chr == "NC_000020.10")
            ret = 19;
        else if (chr == "NC_000021.8")
            ret = 20;
        else if (chr == "NC_000022.10")
            ret = 21;
        else if (chr == "NC_000023.10")
            ret = 22;
        else if (chr == "NC_000024.9")
            ret = 23;
        else
            ret = -1;
        
        return(ret);
    }
};



#endif /* ghuman_hpp */
