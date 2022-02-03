//
//  bioconst.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 02/04/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "biogenome.hpp"
#include "bioread.hpp"

string CGenome::dnaSubstr(int chrNum, unsigned long startpos, unsigned long endpos)
{
    string ret;
    unsigned long i;
    
    for(i=startpos;i<=endpos;i++)
        ret += dna[chrNum][i-1];
    
    return(ret);
}

void CGenome::Read2Memory(string pathPrefix, string pathSuffix, string format)
{
    string path;
    CReader r;
    
    for(int i=0;i<chrCnt;i++)
    {
        path = pathPrefix + chrName[i] + pathSuffix;
        r.ReadFile2Memory(path, format, &dna[i], chrLen[i]);
        string message = "File " + path + " has been loaded.";
        printf("%s\n", message.c_str());
    }
}
