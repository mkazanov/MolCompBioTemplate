//
//  dna.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 01/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "dna.hpp"

int CDNA::END_CHROMOSOME = -2;
int CDNA::END_GENOME = -1;

CDNAPos::CDNAPos(int chrNum_, unsigned long pos_)
{
    chrNum = chrNum_;
    pos = pos_;
    strand = "+";
}

CDNAPos::CDNAPos(int chrNum_, unsigned long pos_, string strand_)
{
    chrNum = chrNum_;
    pos = pos_;
    strand = strand_;
}


int CDNAPos::isNull()
{
    if(chrNum == -1)
        return(1);
    else
        return(0);
}

CDNAInterval::CDNAInterval(int chrNum_, unsigned long startpos_, unsigned long endpos_)
{
    chrNum = chrNum_;
    startpos = startpos_;
    endpos = endpos_;
    strand = "+";
}

CDNAInterval::CDNAInterval(int chrNum_, unsigned long startpos_, unsigned long endpos_, string strand_)
{
    chrNum = chrNum_;
    startpos = startpos_;
    endpos = endpos_;
    strand = strand_;
}

int CDNAInterval::isNull()
{
    if(chrNum == -1)
        return(1);
    else
        return(0);
}

int CDNAInterval::Check()
{
    if(!isNull() && startpos > endpos)
        return(0);
    else
        return(1);
}
