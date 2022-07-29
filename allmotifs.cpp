//
//  allmotifs.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 14/04/2019.
//  Copyright © 2019 Marat Kazanov. All rights reserved.
//

#include <iostream>
#include "allmotifs.hpp"

CRTMapKey::CRTMapKey(string motif_, int RTbin_, int RTstrand_, char mutbase_)
{
    motif = motif_;
    RTbin = RTbin_;
    RTstrand = RTstrand_;
    mutbase = mutbase_;
}

CRTMapKeyNoMut::CRTMapKeyNoMut(string motif_, int RTbin_, int RTstrand_)
{
    motif = motif_;
    RTbin = RTbin_;
    RTstrand = RTstrand_;
}

void CAllMotifs::AnalysisRTmutations(CMutations& muts, string dirpath, CReplicationTiming& rti, CHumanGenome* phuman, string cancer, string sample)
{
    int i;
    int res;
    int RTbin;
    int RTstrand;
    int chrNum;
    string motif;
    
    CReplicationTime rt;
    CMutation* pmut;

    map<CRTMapKey, unsigned long> results, results2strands;
    map<CRTMapKey, unsigned long>::iterator it, it2;
    
    for(i=0;i<muts.mutations.size();i++)
    {
        pmut = &muts.mutations[i];
        
        // RT
        res = rti.GetRT(CHumanGenome::GetChrNum(string(pmut->chr)), pmut->pos, rt);
        if(res)
            RTbin = rti.GetRTBin(rt.RTvalue, (rti.bins));
        else
            RTbin = RT_NULLBIN_NOVALUE;
        
        if(RTbin == RT_NULLBIN_NOVALUE)
            RTstrand = STRAND_NULL;
        else if((pmut->isForwardStrand == 1 && rt.isForward == 0) || (pmut->isForwardStrand == 0 && rt.isForward == 1))
            RTstrand = STRAND_LEADING;
        else if((pmut->isForwardStrand == 1 && rt.isForward == 1) || (pmut->isForwardStrand == 0 && rt.isForward == 0))
            RTstrand = STRAND_LAGGING;
        else
            RTstrand = STRAND_NULL;
        
        // Motif
        chrNum = CHumanGenome::GetChrNum(string(pmut->chr));
        motif = string(1,phuman->dna[chrNum][pmut->pos-2]) + string(1,phuman->dna[chrNum][pmut->pos-1]) + string(1,phuman->dna[chrNum][pmut->pos]);
        
        //cout << mut.chr << '\t' << mut.pos << '\t' << motif << '\t' << RTbin << '\t' << RTstrand << '\t' << expbin << '\t' << expStrand << '\n';
        
        // write results to map
        it = results.find(CRTMapKey(motif,RTbin,RTstrand,pmut->varallele[0]));
        if(it != results.end())
            it->second++;
        else
            results.insert(pair<CRTMapKey, unsigned long>(CRTMapKey(motif,RTbin,RTstrand,pmut->varallele[0]), 1));

    }
    
    results2strands = results;
    
    string cmotif;
    char cvarallele;
    int opRTstrand;
    for(it=results.begin();it!=results.end();it++)
    {
        cmotif = CDNA::cDNA(it->first.motif);
        cvarallele = CDNA::cDNA(string(1,it->first.mutbase))[0];
        opRTstrand = CReplicationTiming::oppositeStrand(it->first.RTstrand);
        it2 = results2strands.find(CRTMapKey(cmotif,it->first.RTbin,opRTstrand,cvarallele));
        if(it2 != results2strands.end())
            it2->second += it->second;
        else
            results2strands.insert(pair<CRTMapKey, unsigned long>(CRTMapKey(cmotif,it->first.RTbin,opRTstrand,cvarallele), it->second));
    }
    
    ofstream f;
    string path;
    path = dirpath + "/" + cancer + "/RT_MUT_" + sample + ".txt";
    f.open(path.c_str());
    if (!f.is_open())
    {
        printf("Directory is not exists\n");
        exit(1);
    }
    for(it=results2strands.begin();it!=results2strands.end();it++)
        f << it->first.motif << '\t' << it->first.RTbin << '\t' << it->first.RTstrand << '\t' << it->first.mutbase << '\t' << it->second << '\n';
    f.close();

}

void CAllMotifs::AnalysisRTtargets(string dirpath, CReplicationTiming& rti, CHumanGenome* phuman)
{
    // Targets
    
    int i,j,k;
    string motif;
    int RTbin;
    int RTstrand;
    int res;
    CReplicationTime rt;
    
    set<string> motifs;
    char nuc[4] = {'A','G','C','T'};
    
    for(i=0;i<4;i++)
        for(j=0;j<4;j++)
            for(k=0;k<4;k++)
            {
                motif = string(1,nuc[i])+string(1,nuc[j])+string(1,nuc[k]);
                motifs.insert(motif);
            }

    map<CRTMapKeyNoMut, unsigned long> resultsT,results2strandsT;
    map<CRTMapKeyNoMut, unsigned long>::iterator itT,it2T;
    
    set<string>::iterator mit;
    for(mit=motifs.begin();mit!=motifs.end();mit++)
        for(RTbin=(-RT_NULLBIN_CNT);RTbin<(int)rti.bins.size();RTbin++)
            for(i=-1;i<2;i++)
                resultsT.insert(pair<CRTMapKeyNoMut, unsigned long>(CRTMapKeyNoMut((*mit),RTbin,i), 0));
    
    for(i=0;i<phuman->chrCnt;i++)
    {
        cout << "Chromosome: " << i << '\n';
        for(j=1;j<(phuman->chrLen[i]-1);j++)
        {
            if(!(CDNA::inACGT(phuman->dna[i][j]) &&
                 CDNA::inACGT(phuman->dna[i][j-1]) &&
                 CDNA::inACGT(phuman->dna[i][j+1])))
                continue;
            
            // RT
            res = rti.GetRT(i, j+1, rt);
            if(res)
                RTbin = rti.GetRTBin(rt.RTvalue, (rti.bins));
            else
                RTbin = RT_NULLBIN_NOVALUE;
            
            if(RTbin == RT_NULLBIN_NOVALUE)
                RTstrand = STRAND_NULL;
            else if(rt.isForward == 0)
                RTstrand = STRAND_LEADING;
            else if(rt.isForward == 1)
                RTstrand = STRAND_LAGGING;
            else
                RTstrand = STRAND_NULL;
            
            // Motif
            motif = string(1,phuman->dna[i][j-1]) + string(1,phuman->dna[i][j]) + string(1,phuman->dna[i][j+1]);
            
            // write results to map
            itT = resultsT.find(CRTMapKeyNoMut(motif,RTbin,RTstrand));
            if(itT != resultsT.end())
                itT->second++;
            else
            {
                cerr << "Error: map key not found";
                exit(1);
            }
        }
    }
    
    results2strandsT = resultsT;
    
    string cmotif;
    int opRTstrand;
    for(itT=resultsT.begin();itT!=resultsT.end();itT++)
    {
        cmotif = CDNA::cDNA(itT->first.motif);
        opRTstrand = CReplicationTiming::oppositeStrand(itT->first.RTstrand);
        it2T = results2strandsT.find(CRTMapKeyNoMut(cmotif,itT->first.RTbin,opRTstrand));
        if(it2T != results2strandsT.end())
            it2T->second += itT->second;
        else
        {
            cerr << "Error: map key not found";
            exit(1);
        }
    }
    
    ofstream f;
    string path;
    path = dirpath + "/RT_TRG.txt";
    f.open(path.c_str());
    if (!f.is_open())
    {
        printf("File not exists\n");
        exit(1);
    }
    for(itT=results2strandsT.begin();itT!=results2strandsT.end();itT++)
        f << itT->first.motif << '\t' << itT->first.RTbin << '\t' << itT->first.RTstrand << '\t' << itT->second << '\n';
    f.close();

}
