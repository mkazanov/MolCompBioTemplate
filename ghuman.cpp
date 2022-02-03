//
//  g_human.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 03/04/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include <string>
#include "ghuman.hpp"

using namespace std;

// Human genome

void CHumanGenome::InitializeHuman(string version_, string pathPrefix, string pathSuffix, string format, int isLoadFiles)
{
    version = version_;
    chrCnt = 24;
    chrName =  new string[chrCnt];
    chrName[0] = "1";
    chrName[1] = "2";
    chrName[2] ="3";
    chrName[3] = "4";
    chrName[4] = "5";
    chrName[5] = "6";
    chrName[6] = "7";
    chrName[7] = "8";
    chrName[8] = "9";
    chrName[9] = "10";
    chrName[10] = "11";
    chrName[11] = "12";
    chrName[12] = "13";
    chrName[13] = "14";
    chrName[14] = "15";
    chrName[15] = "16";
    chrName[16] = "17";
    chrName[17] = "18";
    chrName[18] = "19";
    chrName[19] = "20";
    chrName[20] = "21";
    chrName[21] = "22";
    chrName[22] = "X";
    chrName[23] = "Y";
    
    dna = new char*[chrCnt];
    
    if (version == "37")
    {
        chrLen = new int[chrCnt];
        chrLen[0] = 249250621;
        chrLen[1] = 243199373;
        chrLen[2] = 198022430;
        chrLen[3] = 191154276;
        chrLen[4] = 180915260;
        chrLen[5] = 171115067;
        chrLen[6] = 159138663;
        chrLen[7] = 146364022;
        chrLen[8] = 141213431;
        chrLen[9] = 135534747;
        chrLen[10] = 135006516;
        chrLen[11] = 133851895;
        chrLen[12] = 115169878;
        chrLen[13] = 107349540;
        chrLen[14] = 102531392;
        chrLen[15] = 90354753;
        chrLen[16] = 81195210;
        chrLen[17] = 78077248;
        chrLen[18] = 59128983;
        chrLen[19] = 63025520;
        chrLen[20] = 48129895;
        chrLen[21] = 51304566;
        chrLen[22] = 155270560;
        chrLen[23] = 59373566;
    }
    else if (version == "38")
    {
        chrLen = new int[chrCnt];
        chrLen[0] = 248956422;
        chrLen[1] = 242193529;
        chrLen[2] = 198295559;
        chrLen[3] = 190214555;
        chrLen[4] = 181538259;
        chrLen[5] = 170805979;
        chrLen[6] = 159345973;
        chrLen[7] = 145138636;
        chrLen[8] = 138394717;
        chrLen[9] = 133797422;
        chrLen[10] = 135086622;
        chrLen[11] = 133275309;
        chrLen[12] = 114364328;
        chrLen[13] = 107043718;
        chrLen[14] = 101991189;
        chrLen[15] = 90338345;
        chrLen[16] = 83257441;
        chrLen[17] = 80373285;
        chrLen[18] = 58617616;
        chrLen[19] = 64444167;
        chrLen[20] = 46709983;
        chrLen[21] = 50818468;
        chrLen[22] = 156040895;
        chrLen[23] = 57227415;
    }
    
    if (isLoadFiles)
    {
        Read2Memory(pathPrefix, pathSuffix, format);
    }
}


