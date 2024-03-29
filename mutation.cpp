//
//  mutation.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 28/11/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include "mutation.hpp"
#include "dna.hpp"
#include "ghuman.hpp"
#include <fstream>
#include <iostream>
#include "service.h"

CMutFileFormat::CMutFileFormat(char delimiter_,
               int cancerNo_,
               int sampleNo_,
               int chrNo_,
               int posNo_,
               int refalleleNo_,
               int varalleleNo_,
               int strandNo_,
               int mutationTypeNo_,
               int info1No_,
               int info2No_,
               int isHeader_)
{
    delimiter = delimiter_;
    cancerNo = cancerNo_;
    sampleNo = sampleNo_;
    chrNo = chrNo_;
    posNo = posNo_;
    refalleleNo = refalleleNo_;
    varalleleNo = varalleleNo_;
    strandNo = strandNo_;
    mutationTypeNo = mutationTypeNo_;
    info1No = info1No_;
    info2No = info2No_;
    isHeader = isHeader_;
}

CMutation::CMutation(string cancer_,
                     string sample_,
                     string chr_,
                     string pos_,
                     string refallele_,
                     string varallele_,
                     string isForwardStrand_,
                     string info1_,
                     string info2_)
{
    size_t chrpos;
    
    cancer_.copy(cancer,cancer_.length());
    cancer[cancer_.length()] = '\0';
    sample_.copy(sample,sample_.length());
    sample[sample_.length()] = '\0';
    for(int i=0;i<(STRLEN_CHR+1);i++)
        chr[i] = '\0';
    chrpos = chr_.find("chr");
    if(chrpos != string::npos && chrpos == 0)
        chr_ = chr_.substr(3);
    chr_.copy(chr,chr_.length());
    pos = str2ul(pos_.c_str());
    refallele = refallele_;
    varallele = varallele_;
    if(isForwardStrand_ == "+" || isForwardStrand_ == "1")
        isForwardStrand = 1;
    else if (isForwardStrand_ == "-" || isForwardStrand_ == "0")
        isForwardStrand = 0;
    else
    {
        cerr << "Unknown mutation strand: " << isForwardStrand_ << '\n';
        exit(1);
    }
    info1 = info1_;
    info2 = info2_;
}

CMutations::CMutations()
{
    // Fridriksson file format
    fileFormat.push_back(CMutFileFormat('\t', //separator
                                        1, // cancer field num
                                        0, // sample field num
                                        2, // chromosome field num
                                        3, // position filed num
                                        4, // ref allele num
                                        5, // var allele num
                                        -1, // strand field num
                                        -1, // mutation type num
                                        -1, // info1 field num
                                        -1, // info2 field num
                                        1 // is header
                                        ));
    // PCAWG file format
    fileFormat.push_back(CMutFileFormat('\t',
                                        -1, // cancer field num
                                        8, // sample field num
                                        0, // chromosome field num
                                        1, // position filed num
                                        3, // ref allele num
                                        4, // var allele num
                                        -1, // strand field num
                                        -1, // mutation type num
                                        -1, // info1 field num
                                        -1, // info2 field num
                                        1 // is header
                                        ));
    // Gordenin file format
    fileFormat.push_back(CMutFileFormat('\t',
                                        -1, // cancer field num
                                        20, // sample field num
                                        1, // chromosome field num
                                        2, // position filed num
                                        5, // ref allele num
                                        6, // var allele num
                                        -1, // strand field num
                                        4, // mutation type num
                                        -1, // info1 field num
                                        -1, // info2 field num
                                        1 // is header
                                        ));
    // Gordenin anz4 format
    fileFormat.push_back(CMutFileFormat('\t',
                                        -1, // cancer field num
                                        0, // sample field num
                                        1, // chromosome field num
                                        2, // position filed num
                                        5, // ref allele num
                                        6, // var allele num
                                        -1, // strand field num
                                        4, // mutation type num
                                        -1, // info1 field num
                                        -1, // info2 field num
                                        1 // is header
                                        ));
    
    // Gordenin-Zach anz4 format
    fileFormat.push_back(CMutFileFormat('\t',
                                        -1, // cancer field num
                                        1, // sample field num
                                        2, // chromosome field num
                                        3, // position filed num
                                        4, // ref allele num
                                        5, // var allele num
                                        55, // strand field num
                                        -1, // mutation type num
                                        63, // info1 field num
                                        61, // info2 field num
                                        1 // is header
                                        ));

}

void CMutations::LoadMutations(int fileFormatType /* 0 - Fridriksson, 1 - PCAWG *, 2- Gordenin's*/, string path, vector<string> onlyCancers, vector<string> onlySamples, int onlySubs, CHumanGenome* phuman)
{
    vector<string> flds;
    string line;
    clock_t c1,c2;
    string cancerProject;
    int isForwardStrand;
    
    ifstream f(path.c_str());
    if (!f.is_open())
    {
        printf("File not exists\n");
        exit(1);
    }
    if (fileFormatType == FILE_FORMAT_PCAWG && onlyCancers.empty())
    {
        printf("Cancer project should be specified\n");
        exit(1);
    }
    if (fileFormatType == FILE_FORMAT_PCAWG && onlyCancers.size() != 1)
    {
        printf("Single cancer project should be specified\n");
        exit(1);
    }

    if(fileFormat[fileFormatType].isHeader)
        getline(f, line);
    printf("Calculating number of mutations ...\n");
    mutationsCnt = 0;
    while(1)
    {
        getline(f, line);
        if(f.eof())
            break;
        flds = splitd(line,fileFormat[fileFormatType].delimiter);
        if(fileFormatType == FILE_FORMAT_GORDENIN && flds[fileFormat[fileFormatType].mutationTypeNo] != "SNP")
            continue;
        if (fileFormatType == FILE_FORMAT_FRIDRIKSSON && !onlyCancers.empty() && find(onlyCancers.begin(),onlyCancers.end(),flds[fileFormat[fileFormatType].cancerNo]) == onlyCancers.end())
            continue;
        auto a = flds[fileFormat[fileFormatType].chrNo];
        auto b = phuman->GetChrNum(flds[fileFormat[fileFormatType].chrNo]);
        if(phuman!=NULL && phuman->GetChrNum(flds[fileFormat[fileFormatType].chrNo]) == -1)
            continue;
        if(!onlySamples.empty() && find(onlySamples.begin(),onlySamples.end(),flds[fileFormat[fileFormatType].sampleNo]) == onlySamples.end())
            continue;
        if(onlySubs == 1 && flds[fileFormat[fileFormatType].refalleleNo].size() != flds[fileFormat[fileFormatType].varalleleNo].size())
            continue;
        mutationsCnt++;
    }
    f.clear();
    f.seekg(0, ios::beg);
    
    printf("Allocating space for mutaiton ...\n");
    c1 = clock();
    mutations.reserve(mutationsCnt);
    c2 = clock();
    printf("Executing time: %lu \n", c2 - c1);

    if(fileFormat[fileFormatType].isHeader)
        getline(f, line);
    
    printf("Mutations loading ...\n");
    int i=0;
    c1 = clock();
    flds.clear();
    while(getline(f, line))
    {
        if (line.length() != 0)
        {
            flds = splitd(line,fileFormat[fileFormatType].delimiter);
            if(fileFormatType == FILE_FORMAT_GORDENIN && flds[fileFormat[fileFormatType].mutationTypeNo] != "SNP")
                continue;
            if (fileFormatType == FILE_FORMAT_FRIDRIKSSON && !onlyCancers.empty() && find(onlyCancers.begin(),onlyCancers.end(),flds[fileFormat[fileFormatType].cancerNo]) == onlyCancers.end())
                continue;
            if(phuman!=NULL && phuman->GetChrNum(flds[fileFormat[fileFormatType].chrNo]) == -1)
                continue;
            if(!onlySamples.empty() && find(onlySamples.begin(),onlySamples.end(),flds[fileFormat[fileFormatType].sampleNo]) == onlySamples.end())
                continue;
            if(onlySubs == 1 && flds[fileFormat[fileFormatType].refalleleNo].size() != flds[fileFormat[fileFormatType].varalleleNo].size())
                continue;
            if(fileFormat[fileFormatType].cancerNo == -1)
                cancerProject = onlyCancers[0];
            else
                cancerProject = flds[fileFormat[fileFormatType].cancerNo];
            mutations.push_back(CMutation(cancerProject,
                                              flds[fileFormat[fileFormatType].sampleNo],
                                              flds[fileFormat[fileFormatType].chrNo],
                                              flds[fileFormat[fileFormatType].posNo],
                                              flds[fileFormat[fileFormatType].refalleleNo],
                                              flds[fileFormat[fileFormatType].varalleleNo],
                                              flds[fileFormat[fileFormatType].strandNo],
                                              flds[fileFormat[fileFormatType].info1No],
                                              flds[fileFormat[fileFormatType].info2No]));
            i++;
        }
    }
    c2 = clock();
    printf("Loaded %lu mutations\n", mutationsCnt);
    printf("Executing time: %lu \n", c2 - c1);
    
}

void CMutations::FilterBySample(CMutations &filteredMutations, vector<string> cancers, vector<string> samples)
{
    int i;
    CMutation* pmut;
    
    cout << "Mutation filtering ..." << endl;
    
    for(i=0;i<mutationsCnt;i++)
    {
        pmut = &mutations[i];
        if (string(pmut->chr) == "M")
            continue;
        if (!cancers.empty() && find(cancers.begin(),cancers.end(),pmut->cancer) == cancers.end())
            continue;
        if (!samples.empty() && find(samples.begin(),samples.end(),pmut->sample) == samples.end())
            continue;
        
        filteredMutations.mutations.push_back(*pmut);
    }
    
    cout << filteredMutations.mutations.size() << " mutations selected." << '\n';
    filteredMutations.mutationsCnt = filteredMutations.mutations.size();
    cout << "Done" << endl;
}

void CMutations::FilterBySignature(CMutations& filteredMutations, vector<CMutationSignature>& signatures, CHumanGenome& human,
                                 vector<string> cancers, vector<string> samples, CMutations* pOtherMutations)
{
    int i,j;
    CMutation m;
    CMutationSignature s;
    string mutBase;
    string ss;
    int startPosShift, endPosShift;
    string seq,seqrev;
    int foundMutation;
    
    printf("Mutation filtering ...");
    
    for(i=0;i<mutationsCnt;i++)
    {
        m = mutations[i];
        if (string(m.chr) == "M")
            continue;
        if (!cancers.empty() && find(cancers.begin(),cancers.end(),m.cancer) == cancers.end())
            continue;
        if (!samples.empty() && find(samples.begin(),samples.end(),m.sample) == samples.end())
            continue;
        
        foundMutation = 0;
        for(j=0;j<signatures.size();j++)
        {
            s = signatures[j];
            startPosShift = s.mutationPos - 1;
            endPosShift = s.motif.length() - s.mutationPos;
            mutBase = s.motif[s.mutationPos-1];
            ss = string(m.chr);

            seq = human.dnaSubstr(human.GetChrNum(string(m.chr)),m.pos-startPosShift,m.pos+endPosShift);
            seqrev = CDNA::cDNA(human.dnaSubstr(human.GetChrNum(string(m.chr)),m.pos-endPosShift,m.pos+startPosShift));
                        
            if ((m.refallele == mutBase &&
                 (s.AnyNewBase() || (m.varallele == s.newbase) ) &&
                 CDNA::compareMotifDNA(s.motif,seq)) ||
                (m.refallele == CDNA::cDNA(mutBase) &&
                 (s.AnyNewBase() || m.varallele == CDNA::cDNA(s.newbase) ) &&
                 CDNA::compareMotifDNA(s.motif,seqrev))
                )
                {
                    if (m.refallele == mutBase)
                        m.isForwardStrand = 1;
                    else if(m.refallele == CDNA::cDNA(mutBase))
                        m.isForwardStrand = 0;
                    else
                        m.isForwardStrand = -1;
                    filteredMutations.mutations.push_back(m);
                    foundMutation = 1;
                    break;
                }
        }
        if(foundMutation == 0)
            if(pOtherMutations!=NULL)
            {
                m.isForwardStrand = -1;
                pOtherMutations->mutations.push_back(m);
            }
    }
    cout << "Signature mutations: " << filteredMutations.mutations.size() << '\n';
    filteredMutations.mutationsCnt = filteredMutations.mutations.size();
    if(pOtherMutations!=NULL)
    {
        cout << "Other mutations: " << pOtherMutations->mutations.size() << '\n';
        pOtherMutations->mutationsCnt = pOtherMutations->mutations.size();
    }
    printf("Done\n");
}

void CMutations::CheckRefAllels(CHumanGenome* phuman)
{
    int i;
    int chrNum;
    string dnaAllele;
    CMutation* pmut;
    
    for(i=0;i<mutations.size();i++)
    {
        pmut = &mutations[i];
        chrNum = CHumanGenome::GetChrNum(string(pmut->chr));
        dnaAllele = string(1,phuman->dna[chrNum][pmut->pos - 1]);
        if(dnaAllele != pmut->refallele)
        {
            cerr << "Mutations not correct";
            throw exception();
        }
    }
}

void CMutations::SaveToFile(string path)
{
    ofstream f;
    f.open(path.c_str());
    int i;
    
    for(i=0;i<mutations.size();i++)
        f << string(mutations[i].cancer) << '\t' << string(mutations[i].sample) << '\t' <<
        string(mutations[i].chr) << '\t' << ul2str(mutations[i].pos) << '\t' << mutations[i].refallele <<
        '\t' << mutations[i].varallele << '\t' << (int) mutations[i].isForwardStrand << '\n';
    
    f.close();
}

void CMutations::SaveToFileRTExp(string path)
{
    ofstream f;
    f.open(path.c_str());
    int i;
    
    for(i=0;i<mutations.size();i++)
        f << string(mutations[i].cancer) << '\t' << string(mutations[i].sample) << '\t' <<
        string(mutations[i].chr) << '\t' << ul2str(mutations[i].pos) << '\t' << mutations[i].refallele <<
        '\t' << mutations[i].varallele << '\t' << (int) mutations[i].isForwardStrand << '\t' << i2str(mutations[i].RTbin) << '\t' << i2str(mutations[i].RTstrand) << '\t' << i2str(mutations[i].EXPbin) << '\t' << i2str(mutations[i].EXPstrand) << '\n';
    
    f.close();
}


void CMutations::GetUniqueCancersSamples()
{
    cancerSample.clear();
    for(int i=0;i<mutations.size();i++)
        cancerSample.insert(CCancerSample(mutations[i].cancer,mutations[i].sample));
}

void CMutations::RenameSamples(string renameTablePath, int columnNumOld, int columnNumNew, int newSampleNameLen)
{
    string line;
    vector<string> flds;
    
    ifstream f(renameTablePath.c_str());
    if (!f.is_open())
    {
        printf("File not exists\n");
        exit(1);
    }
    
    getline(f, line);
    while(1)
    {
        getline(f, line);
        if(f.eof())
            break;
        
        flds = splitd(line,',');
        renamemap[flds[columnNumOld]] = flds[columnNumNew];
    }
    
    for(int i=0;i<mutations.size();i++)
    {
        renamemap[mutations[i].sample].copy(mutations[i].sample,STRLEN_SAMPLE);
        mutations[i].sample[newSampleNameLen] = '\0';
    }
}

void CMutations::ClearMutations()
{
    mutations.clear();
    mutationsCnt = 0;
    cancerSample.clear();
    renamemap.clear();
}
