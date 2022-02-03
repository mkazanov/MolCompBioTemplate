//
//  dna.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 01/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef dna_hpp
#define dna_hpp

#include <string>
#include <algorithm>

using namespace std;

class CDNA{
public:
    static string cDNA(string dna)
    {
        string ret;
        
        ret = dna;
        transform(ret.begin(),ret.end(),ret.begin(),::toupper);
        reverse(ret.begin(),ret.end());
        replace(ret.begin(),ret.end(),'A','Y');
        replace(ret.begin(),ret.end(),'T','Z');
        replace(ret.begin(),ret.end(),'Z','A');
        replace(ret.begin(),ret.end(),'Y','T');
        replace(ret.begin(),ret.end(),'C','Y');
        replace(ret.begin(),ret.end(),'G','Z');
        replace(ret.begin(),ret.end(),'Z','C');
        replace(ret.begin(),ret.end(),'Y','G');

        return(ret);
    }
    static bool inACGT(char base)
    {
        if(base == 'A' || base == 'C' || base == 'G' || base == 'T')
            return(true);
        else
            return(false);
    }
};

class CDNAPos{
public:
    int chrNum;
    unsigned long pos;
    int strand;
    CDNAPos(){};
    CDNAPos(int chrNum_, unsigned long pos_);
    int isNull();
};

#endif /* dna_hpp */
