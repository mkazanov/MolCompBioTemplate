//
//  replicationtime.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 03/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "replicationtime.hpp"
#include "service.h"
#include <fstream>
#include "ghuman.hpp"
#include "dna.hpp"
#include <map>
#include "options.h"
#include <iostream>
#include <algorithm>
#include "mutsignature.hpp"
#include <cstring>
#include <cmath>

CReplicationTime::CReplicationTime(string startpos_, string endpos_, string RTvalue_)
{
    startpos = str2ul(startpos_);
    endpos = str2ul(endpos_);
    RTvalue = str2d(RTvalue_);
}

CReplicationTime::CReplicationTime(unsigned long pos_)
{
    startpos = pos_;
}

bool CReplicationTime::isRTnull()
{
    if(startpos == -1)
        return true;
    else
        return false;
}

CRTBin::CRTBin(int binNum_, double RTleft_, double RTright_)
{
    binNum = binNum_;
    RTleft = RTleft_;
    RTright = RTright_;
}

CReplicationTiming::CReplicationTiming()
{
    RTs = new set<CReplicationTime>[24];
}

void CReplicationTiming::LoadReplicationTiming(string path, int isHeader)
{
    string line;
    clock_t c1,c2;
    int chrNum;
    
    ifstream f(path.c_str());
    if (!f.is_open())
    {
        printf("File not exists\n");
        exit(1);
    }
    
    if(isHeader)
        getline(f, line);
    
    printf("Replication timing loading ...\n");
    int i=0;
    vector<string> flds;
    c1 = clock();
    while(getline(f, line))
    {
        if (line.length() != 0)
        {
            flds = split(line);
            chrNum = CHumanGenome::GetChrNum(flds[0].substr(3));
            RTs[chrNum].insert(CReplicationTime(flds[1],flds[2],flds[3]));
            i++;
        }
    }
    c2 = clock();
    printf("Replication timing intervals have been loaded\n");
    printf("Executing time: %lu \n", c2 - c1);
}

int CReplicationTiming::GetRT(int chrNum, unsigned long pos, double& RTvalue)
{
    set<CReplicationTime>::iterator it;

    CReplicationTime rt(pos);
    it = RTs[chrNum].upper_bound(rt);
    if(it == RTs[chrNum].begin())
    {
        RTvalue = RT_NULL;
        return(0);
    }
    it--;
    if(it != RTs[chrNum].end() && pos >= it->startpos && pos <= it->endpos)
    {
        RTvalue = it->RTvalue;
        return (1);
    }
    else
    {
        RTvalue = RT_NULL;
        return(0);
    }
}

int CReplicationTiming::GetRT(int chrNum, unsigned long pos, CReplicationTime& rt)
{
    set<CReplicationTime>::iterator it;
    
    CReplicationTime rtquery(pos);
    it = RTs[chrNum].upper_bound(rtquery);
    if(it == RTs[chrNum].begin())
        return(0);
    it--;
    if(it != RTs[chrNum].end() && pos >= it->startpos && pos <= it->endpos)
    {
        rt = (*it);
        return (1);
    }
    else
        return(0);
}

int CReplicationTiming::GetRTBin(double RTvalue, vector<CRTBin> bins)
{
    for(int i=0;i<bins.size();i++)
        if(RTvalue >= bins[i].RTleft && RTvalue < bins[i].RTright)
            return(bins[i].binNum);
    return(RT_NULLBIN_NOBIN);
}

int CReplicationTiming::GetRTBin(int chrNum, unsigned long pos, vector<CRTBin> bins)
{
    double RTvalue;
    int res;
    int bin;
    
    res = GetRT(chrNum,pos,RTvalue);
    if(res)
    {
        bin = GetRTBin(RTvalue,bins);
        return(bin);
    }
    else
        return(RT_NULLBIN_NOVALUE);
}

void CReplicationTiming::ReplicationStrand()
{
    set<CReplicationTime>::iterator it,previt,nextit;
    int chrNum;
    
    for(chrNum=0;chrNum<24;chrNum++)
        for(it=RTs[chrNum].begin();it!=RTs[chrNum].end();++it)
        {
            //cout << it->RTvalue << '\n';
            if(it != RTs[chrNum].begin())
            {
                previt = it;
                previt--;
            }
            nextit = it;
            nextit++;
            if(it == RTs[chrNum].begin())
                it->isForward = -1;
            else
            {
                if(nextit->RTvalue - previt->RTvalue > 0)
                    it->isForward = 0;
                else if(nextit->RTvalue - previt->RTvalue < 0)
                    it->isForward = 1;
                else
                    it->isForward = -1;
            }
        }
}

void CReplicationTiming::SaveToFile(string path)
{
    ofstream f;
    f.open(path.c_str());
    set<CReplicationTime>::iterator it;

    for(int chrNum=0;chrNum<24;chrNum++)
        for(it=RTs[chrNum].begin();it!=RTs[chrNum].end();++it)
            f << chrNum << '\t' << (*it).startpos << '\t' << (*it).endpos << '\t' << (*it).RTvalue << '\t' << (*it).isForward << '\n';
        
    f.close();
}

void CReplicationTiming::GetRTBinStrand(map<string,CReplicationTiming*> rtmap,const CMutation& mut, int& bin, int& strand)
{
    int res;
    CReplicationTime rt;
    
    res = rtmap[string(mut.cancer)]->GetRT(CHumanGenome::GetChrNum(string(mut.chr)), mut.pos, rt);
    if(res)
        bin = rtmap[string(mut.cancer)]->GetRTBin(rt.RTvalue, (rtmap[string(mut.cancer)]->bins));
    else
        bin = RT_NULLBIN_NOVALUE;
    
    if((mut.isForwardStrand == 1 && rt.isForward == 1) || (mut.isForwardStrand == 0 && rt.isForward == 0))
        strand = STRAND_LAGGING;
    else if((mut.isForwardStrand == 1 && rt.isForward == 0) || (mut.isForwardStrand == 0 && rt.isForward == 1))
        strand = STRAND_LEADING;
    else
        strand = STRAND_NULL;

}

int CReplicationTiming::CalculateMotifinRTBins(set<string> motifs, string OUT_PATH, CHumanGenome* phuman)
{
    CMutationSignature msobj;
    CReplicationTime rt;
    int bin;
    map<int,unsigned long> results;
    map<int,unsigned long> leading;
    map<int,unsigned long> lagging;
    map<int,unsigned long>::iterator it;
    set<string> motifsall;
    set<string>::iterator si;
    long motiflen,motifsnum;
    char** motifsarr;
    int* strandarr;
    int i;
    
    msobj.CheckMotifsNotEmpty(motifs);
    msobj.CheckMotifsSameLength(motifs);
    motifsall = msobj.AddcMotifs(motifs);
    motifsnum = (int) motifsall.size();
    motiflen = (motifsall.begin())->length();
    
    // Prepare char array for copying motifs
    motifsarr = new char*[motifsall.size()];
    strandarr = new int[motifsall.size()];
    for(i=0;i<motifsall.size();i++)
    {
        motifsarr[i] = new char[motiflen];
        if(i<motifs.size())
            strandarr[i] = 1;
        else
            strandarr[i] = 0;
    }
        
    // Copy motifs to char array
    i = 0;
    for(si=motifsall.begin();si!=motifsall.end();si++)
    {
        strncpy(motifsarr[i],(*si).c_str(),(*si).length());
        i++;
    }

    results[RT_NULLBIN_NOVALUE] = 0;
    results[RT_NULLBIN_NOBIN] = 0;
    leading[RT_NULLBIN_NOVALUE] = 0;
    leading[RT_NULLBIN_NOBIN] = 0;
    lagging[RT_NULLBIN_NOVALUE] = 0;
    lagging[RT_NULLBIN_NOBIN] = 0;
    for(int i=0;i<bins.size();i++)
    {
        results[bins[i].binNum] = 0;
        leading[bins[i].binNum] = 0;
        lagging[bins[i].binNum] = 0;
    }
    
    CDNAPos pos = CDNAPos(0,0);
    int includeCurPos=1;
    int chrNum = 0;
    int res;
    cout << "Chr:" << pos.chrNum << ", Pos:" << pos.pos << '\n';
    for(pos=msobj.NextMotif(CDNAPos(0,0),motifsarr,strandarr,(int)motifsnum,(int)motiflen,phuman,END_GENOME,includeCurPos);
        !pos.isNull();
        pos=msobj.NextMotif(pos,motifsarr,strandarr,(int)motifsnum,(int)motiflen,phuman,END_GENOME))
    {
        res = GetRT(pos.chrNum, pos.pos+(motiflen/2)+1, rt);
        if(res)
            bin = GetRTBin(rt.RTvalue, bins);
        else
            bin = RT_NULLBIN_NOVALUE;
        results[bin]++;
        
        if(bin != RT_NULLBIN_NOVALUE)
        {
            if((pos.strand == 1 && rt.isForward == 1) || (pos.strand == 0 && rt.isForward == 0))
                lagging[bin]++;
            else if((pos.strand == 1 && rt.isForward == 0) || (pos.strand == 0 && rt.isForward == 1))
                leading[bin]++;
        }
        
        if(chrNum != pos.chrNum)
        {
            cout << "Chr:" << pos.chrNum << ", Pos:" << pos.pos << '\n';
            chrNum = pos.chrNum;
        }
    }
    
    ofstream f;
    f.open(OUT_PATH.c_str());
    
    f << "ReplicationBin" << '\t' << "TargetCnt" << '\t' << "Leading" << '\t' << "Lagging" << '\n';
    for(it=results.begin(); it!=results.end(); ++it)
        f << it->first << '\t' << it->second << '\t' << leading[it->first] << '\t' << lagging[it->first] << '\n';
    
    f.close();
    
    return(1);
}

void CReplicationTiming::CalculateRTbinsWidth(int binsNum, CHumanGenome* phuman, string outputPath)
{
    int i,j,k;
    double rtval;
    int res;
    string motif,motif2;
    map<string, vector<double> > results;
    map<string, vector<double> >::iterator it;
    char nuc[4] = {'A','G','C','T'};
    vector<double> emptyVector;
    string path;

    for(i=0;i<4;i++)
        for(j=0;j<4;j++)
            for(k=0;k<4;k++)
            {
                motif = string(1,nuc[i])+string(1,nuc[j])+string(1,nuc[k]);
                motif2 = CDNA::cDNA(motif);
                motif = (motif < motif2) ? motif : motif2;
                results[motif] = emptyVector;
            }
    
    for(i=0;i<phuman->chrCnt;i++)
    {
        cout << "Chromosome: " << i << '\n';

        for(j=1;j<(phuman->chrLen[i]-1);j++)
        {
            if(!(CDNA::inACGT(phuman->dna[i][j]) &&
                 CDNA::inACGT(phuman->dna[i][j-1]) &&
                 CDNA::inACGT(phuman->dna[i][j+1])))
                continue;
            
            res = GetRT(i, j+1, rtval);
            if(!res)
                continue;

            motif = string(1,phuman->dna[i][j-1]) + string(1,phuman->dna[i][j]) + string(1,phuman->dna[i][j+1]);
            motif2 = CDNA::cDNA(motif);
            motif = (motif < motif2) ? motif : motif2;
            
            results[motif].push_back(rtval);
        }
    }
    
    ofstream f;
    double step;
    f.open(outputPath.c_str());
    for(it=results.begin();it!=results.end();it++)
    {
        sort(it->second.begin(),it->second.end());
        f << it->first << '\t';
        step = it->second.size() / binsNum;
        for(i=1;i<binsNum;i++)
        {
            f << it->second[round(i*step)] << '\t';
        }
        f << '\n';
    }
    f.close();
    
}
