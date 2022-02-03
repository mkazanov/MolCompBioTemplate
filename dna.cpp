//
//  dna.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 01/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "dna.hpp"

CDNAPos::CDNAPos(int chrNum_, unsigned long pos_)
{
    chrNum = chrNum_;
    pos = pos_;
}

int CDNAPos::isNull()
{
    if(chrNum == -1)
        return(1);
    else
        return(0);
}


