//
//  genehuman.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 09/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef genehuman_hpp
#define genehuman_hpp

#include <string>
#include <set>
#include <vector>
#include <bitset>

using namespace std;

class CHumanGene {
public:
    //int chrNum;
    unsigned long startpos;
    unsigned long endpos;
    mutable unsigned long maxendpos;
    char strand;
    unsigned long geneID;
    string geneName;
    CHumanGene(string starpos_, string endpos_, string strand_, string info_);
    CHumanGene(unsigned long pos_);
    CHumanGene(unsigned long startpos_, unsigned long endpos_);
    bool operator< (const CHumanGene &right) const
    {
        if (startpos < right.startpos)
            return true;
        else if (startpos == right.startpos)
            return endpos > right.endpos;
        else
            return false;
    }
    //int isNull()
    //{
    //    return((chrNum == -1) ? 1 : 0);
    //}
};

class CHumanGenes {
public:
    CHumanGenes();
    set<CHumanGene>* genes;
    bitset<250000000>* genebits;
    void LoadGenes(string path);
    void PrepareForSearch();
    int GetGenesByPos(int chrNum, unsigned long pos, vector<CHumanGene>& geneList);
    void SaveToFile(string path);
    //void MergeIntervals(vector<CHumanGene>& ret);
};

#endif /* genehuman_hpp */
