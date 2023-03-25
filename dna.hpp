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
#include <map>

using namespace std;

class CDNA{
public:
    static int END_CHROMOSOME;
    static int END_GENOME;
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
    static bool compareMotifDNA(string motif_, string seq_)
    {
        map<char,string> symbols;
        char chs,chm;
        string motif,seq;
        
        motif = motif_;
        seq = seq_;
        
        if(motif.length() != seq.length())
            return(false);
        
        transform(motif.begin(),motif.end(),motif.begin(),::toupper);
        transform(seq.begin(),seq.end(),seq.begin(),::toupper);
        
        symbols['A'] = "A";
        symbols['C'] = "C";
        symbols['G'] = "G";
        symbols['T'] = "T";
        symbols['W'] = "ATW";
        symbols['S'] = "CGS";
        symbols['M'] = "ACM";
        symbols['K'] = "GTK";
        symbols['R'] = "AGR";
        symbols['Y'] = "CTY";
        symbols['B'] = "CGTSKYB";
        symbols['D'] = "AGTRKWD";
        symbols['H'] = "ACTMYWH";
        symbols['V'] = "ACGMSRV";
        symbols['N'] = "ACGTWSMKRYBDHVN";
        
        for(int i=0;i<motif.length();i++)
        {
            chm = motif[i];
            chs = seq[i];
            if (symbols[chm].find(chs) == std::string::npos)
                return(false);
        }
        return(true);
    }
};

class CDNAPos{
public:
    int chrNum;
    unsigned long pos;
    string strand;
    CDNAPos(){chrNum = -1;};
    CDNAPos(int chrNum_, unsigned long pos_);
    CDNAPos(int chrNum_, unsigned long pos_, string strand_);
    int isNull();
};

class CDNAInterval{
public:
    int chrNum;
    unsigned long startpos;
    unsigned long endpos;
    string strand;
    CDNAInterval(){};
    CDNAInterval(int chrNum_);
    CDNAInterval(int chrNum_, unsigned long startpos_, unsigned long endpos_);
    CDNAInterval(int chrNum_, unsigned long startpos_, unsigned long endpos_, string strand_);
    int Check();
    int isNull();
};

#endif /* dna_hpp */
