//
//  genehuman.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 09/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "genehuman.hpp"
#include <fstream>
#include "service.h"
#include "ghuman.hpp"
#include <iostream>

CHumanGene::CHumanGene(string startpos_, string endpos_, string strand_, string info_)
{
    size_t pos1,pos2;
    
    startpos = str2ul(startpos_);
    endpos = str2ul(endpos_);
    strand = strand_[0];
    pos1 = info_.find("GeneID:");
    if(pos1 != string::npos)
    {
        pos2 = info_.find(";",pos1);
        geneID = str2ul(info_.substr(pos1+7, pos2-pos1-7));
    }
    pos1 = info_.find("gene=");
    if(pos1 != string::npos)
    {
        pos2 = info_.find(";",pos1);
        geneName = info_.substr(pos1+5, pos2-pos1-5);
    }
}

CHumanGene::CHumanGene(unsigned long pos_)
{
    startpos = pos_;
}

CHumanGene::CHumanGene(unsigned long startpos_, unsigned long endpos_)
{
    startpos = startpos_;
    endpos = endpos_;
}

CHumanGenes::CHumanGenes()
{
    genes = new set<CHumanGene>[24];
    genebits = new bitset<250000000>[24];
}

void CHumanGenes::LoadGenes(string path)
{
    cout << "Loading genes" << '\n';
    
    string line;
    int chrNum;
    
    ifstream f(path.c_str());
    if (!f.is_open())
    {
        printf("File not exists\n");
        return;
    }

    vector<string> flds;
    while(getline(f, line))
    {
        if (line.length() != 0)
        {
            if(line[0] == '#')
                continue;
            flds = split(line);
            if(flds[2] != "gene")
                continue;
            chrNum = CHumanGenome::GetChrFromNCBIName(flds[0]);
            if(chrNum < 0 || chrNum > 23)
                continue;
            genes[chrNum].insert(CHumanGene(flds[3],flds[4],flds[6],flds[8]));
        }
    }
    
    f.close();
}

void CHumanGenes::PrepareForSearch()
{
    set<CHumanGene>::iterator it;
    unsigned long maxend;
    int chrNum;

    // set maxpos
    for(chrNum=0;chrNum<24;chrNum++)
    {
        maxend = 0;
        for(it=genes[chrNum].begin();it!=genes[chrNum].end();++it)
        {
            (*it).maxendpos = (maxend < (*it).startpos) ? 0 : maxend;
            maxend = (maxend > (*it).endpos) ? maxend : (*it).endpos;
        }
    }
    
    // set bits
    unsigned long i;
    for(chrNum=0;chrNum<24;chrNum++)
    {
        for(it=genes[chrNum].begin();it!=genes[chrNum].end();++it)
        {
            for(i=it->startpos;i<=it->endpos;i++)
                genebits[chrNum][i-1] = 1;
        }
    }
}

int CHumanGenes::GetGenesByPos(int chrNum, unsigned long pos, vector<CHumanGene>& geneList)
{
    set<CHumanGene>::iterator it;
    CHumanGene g(pos);

    if(genebits[chrNum][pos-1] != 1)
        return(0);
    
    it = genes[chrNum].upper_bound(g);
    if(it == genes[chrNum].begin())
        return(0);
    it--;
    while(1)
    {
        if((*it).startpos <= pos && (*it).endpos >= pos)
            geneList.push_back((*it));
        if(pos > (*it).maxendpos)
            break;
        it--;
    }
    
    if(geneList.empty())
        return(0);
    else
        return(1);
}

/*void CHumanGenes::MergeIntervals(vector<CHumanGene>& ret)
{
    set<CHumanGene>::iterator si;
    CHumanGene prevGene(-1,0);
    unsigned long startpos=-1,endpos=-1;

    for(si=genes.begin();si!=genes.end();si++)
    {
        if((*si).maxendpos == 0 || prevGene.chrNum != (*si).chrNum || prevGene.isNull())
        {
            if(!prevGene.isNull())
                ret.push_back(CHumanGene(prevGene.chrNum,startpos,endpos));
            startpos = (*si).startpos;
            endpos = (*si).endpos;
        }
        else
        {
            endpos = ((*si).endpos > endpos ? (*si).endpos : endpos);
        }
        prevGene = (*si);
    }
    ret.push_back(CHumanGene(prevGene.chrNum,startpos,endpos));
}*/

void CHumanGenes::SaveToFile(string path)
{
    ofstream f;
    f.open(path.c_str());
    set<CHumanGene>::iterator it;
    int chrNum;
    
    for(chrNum=0;chrNum<24;chrNum++)
        for(it=genes[chrNum].begin();it!=genes[chrNum].end();++it)
            f << chrNum << '\t' << (*it).startpos << '\t' << (*it).endpos << '\t' << (*it).maxendpos << '\t' << (*it).strand << '\t' << (*it).geneID << '\t' << (*it).geneName << '\n';
    
    f.close();
}
