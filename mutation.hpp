//
//  mutation.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 28/11/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef mutation_hpp
#define mutation_hpp

#define STRLEN_CANCER 4
#define STRLEN_SAMPLE 36
#define STRLEN_CHR 2
#define STRLEN_REFALLELE 1
#define STRLEN_VARALLELE 1
#define FILE_FORMAT_FRIDRIKSSON 0
#define FILE_FORMAT_PCAWG 1

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "mutsignature.hpp"
#include "ghuman.hpp"
#include "mutation.hpp"
#include <fstream>
#include <map>

using namespace std;

class CMutFileFormat{
public:
    char delimiter;
    int cancerNo;
    int sampleNo;
    int chrNo;
    int posNo;
    int refalleleNo;
    int varalleleNo;
    int isHeader;
    CMutFileFormat(char delimiter_,
              int cancerNo_,
              int sampleNo_,
              int chrNo_,
              int posNo_,
              int refalleleNo_,
              int varalleleNo_,
              int isHeader);
};

class CCancerSample{
public:
    string cancer;
    string sample;
    bool operator< (const CCancerSample &right) const
    {
        if (cancer < right.cancer)
            return true;
        else if (cancer == right.cancer)
            return sample < right.sample;
        else
            return false;
    }
    CCancerSample(string cancer_, string sample_)
    {
        cancer = cancer_;
        sample = sample_;
    }
    static set<CCancerSample> LoadCancerSample(string path)
    {
        string line;
        vector<string> flds;
        set<CCancerSample> ret;
        
        ifstream f(path.c_str());
        if (!f.is_open())
        {
            printf("File not exists\n");
            exit(1);
        }
        
        while(1)
        {
            getline(f, line);
            flds = split(line);
            ret.insert(CCancerSample(flds[0],flds[1]));
            if(f.eof())
                break;
        }
        
        return(ret);
    }
};

class CMutation{
public:
    char cancer[STRLEN_CANCER+1];
    char sample[STRLEN_SAMPLE+1];
    char chr[STRLEN_CHR+1];
    unsigned long pos;
    string refallele;
    string varallele;
    char isForwardStrand;
    int RTbin;
    int RTstrand;
    int EXPbin;
    int EXPstrand;
    CMutation(string cancer_,
              string sample_,
              string chr_,
              string pos_,
              string refallele_,
              string varallele_,
              char isForwardStrand_);
    CMutation(){};
};

class CMutations{
public:
    CMutations();
    vector<CMutFileFormat> fileFormat;
    vector<CMutation> mutations;
    map<string, string> renamemap;
    unsigned long mutationsCnt;
    void LoadMutations(int fileFormatType /* 0 - Fridriksson, 1 - PCAWG */, string path, vector<string> onlyCancers=vector<string>(), vector<string> onlySamples=vector<string>(), int onlySubs=1, CHumanGenome* phuman = NULL);
    void FilterMutations(CMutations& filteredMutations,
                         vector<CMutationSignature>& signatures,
                         CHumanGenome& human,
                         vector<string> cancers,
                         vector<string> samples,
                         CMutations* pOtherMutations);
    set<CCancerSample> cancerSample;
    void GetUniqueCancersSamples();
    void SaveToFile(string path);
    void SaveToFileRTExp(string path);
    void RenameSamples(string renameTablePath, int columnNumOld, int columnNumNew, int newSampleNameLen);
};

#endif /* mutation_hpp */
