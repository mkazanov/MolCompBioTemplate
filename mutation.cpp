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
               int isHeader_)
{
    delimiter = delimiter_;
    cancerNo = cancerNo_;
    sampleNo = sampleNo_;
    chrNo = chrNo_;
    posNo = posNo_;
    refalleleNo = refalleleNo_;
    varalleleNo = varalleleNo_;
    isHeader = isHeader_;
}

CMutation::CMutation(string cancer_,
                     string sample_,
                     string chr_,
                     string pos_,
                     string refallele_,
                     string varallele_,
                     char isForwardStrand_)
{
    size_t chrpos;
    
    cancer_.copy(cancer,STRLEN_CANCER);
    cancer[STRLEN_CANCER] = '\0';
    sample_.copy(sample,STRLEN_SAMPLE);
    sample[STRLEN_SAMPLE] = '\0';
    for(int i=0;i<(STRLEN_CHR+1);i++)
        chr[i] = '\0';
    chrpos = chr_.find("chr");
    if(chrpos != string::npos && chrpos == 0)
        chr_ = chr_.substr(3);
    chr_.copy(chr,chr_.length());
    pos = str2ul(pos_.c_str());
    refallele = refallele_;
    varallele = varallele_;
    isForwardStrand = isForwardStrand_;
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
                                        1 // is header
                                        ));
}

void CMutations::LoadMutations(int fileFormatType /* 0 - Fridriksson, 1 - PCAWG */, string path, vector<string> onlyCancers, vector<string> onlySamples, int onlySubs, CHumanGenome* phuman)
{
    vector<string> flds;
    string line;
    clock_t c1,c2;
    string cancerProject;
    
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
            if (fileFormatType == FILE_FORMAT_FRIDRIKSSON && !onlyCancers.empty() && find(onlyCancers.begin(),onlyCancers.end(),flds[fileFormat[fileFormatType].cancerNo]) == onlyCancers.end())
                continue;
            if(phuman!=NULL && phuman->GetChrNum(flds[fileFormat[fileFormatType].chrNo]) == -1)
                continue;
            if(!onlySamples.empty() && find(onlySamples.begin(),onlySamples.end(),flds[fileFormat[fileFormatType].sampleNo]) == onlySamples.end())
                continue;
            if(onlySubs == 1 && flds[fileFormat[fileFormatType].refalleleNo].size() != flds[fileFormat[fileFormatType].varalleleNo].size())
                continue;
            if(fileFormatType == FILE_FORMAT_PCAWG)
                cancerProject = onlyCancers[0];
            else
                cancerProject = flds[fileFormat[fileFormatType].cancerNo];
            mutations.push_back(CMutation(cancerProject,
                                              flds[fileFormat[fileFormatType].sampleNo],
                                              flds[fileFormat[fileFormatType].chrNo],
                                              flds[fileFormat[fileFormatType].posNo],
                                              flds[fileFormat[fileFormatType].refalleleNo],
                                              flds[fileFormat[fileFormatType].varalleleNo],1));
            i++;
        }
    }
    c2 = clock();
    printf("Loaded %lu mutations\n", mutationsCnt);
    printf("Executing time: %lu \n", c2 - c1);
    
}

void CMutations::FilterMutations(CMutations& filteredMutations, vector<CMutationSignature>& signatures, CHumanGenome& human,
                                 vector<string> cancers, vector<string> samples, CMutations* pOtherMutations=NULL)
{
    int i,j;
    CMutation m;
    CMutationSignature s;
    string mutBase;
    string ss;
    unsigned long startPosShift, endPosShift;
    string gmotif;
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
            startPosShift = 1 - s.mutationPos;
            endPosShift = s.motif.length() - s.mutationPos;
            mutBase = s.motif[s.mutationPos-1];
            ss = string(m.chr);

            gmotif = human.dnaSubstr(human.GetChrNum(string(m.chr)),m.pos+startPosShift,m.pos+endPosShift);
                        
            if ((m.refallele == mutBase &&
                 (s.AnyNewBase() || (m.varallele == s.newbase) ) &&
                 gmotif == s.motif) ||
                (m.refallele == CDNA::cDNA(mutBase) &&
                 (s.AnyNewBase() || m.varallele == CDNA::cDNA(s.newbase) ) &&
                 gmotif == CDNA::cDNA(s.motif) )
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
    if(pOtherMutations!=NULL)
        cout << "Other mutations: " << pOtherMutations->mutations.size() << '\n';
    printf("Done\n");
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
    int i;
    
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
