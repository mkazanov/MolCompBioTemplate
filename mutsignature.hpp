//
//  mutsignature.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 01/12/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef mutsignature_hpp
#define mutsignature_hpp

#include <string>
#include <stdio.h>
#include "ghuman.hpp"
#include <set>
#include "dna.hpp"

using namespace std;

class CMutationSignature {
public:
    string motif;
    int mutationPos;
    string newbase; //'X' or 'N' - any base
    CMutationSignature(string motif_, int mutationPos_, string newbase_);
    CMutationSignature(){};
    bool AnyNewBase();
    int CheckMotifsNotEmpty(set<string> motifs);
    int CheckMotifsSameLength(set<string> motifs);
    set<string> AddcMotifs(set<string> motifs);
    int NextMotif(CDNAInterval dnainterval, vector<CMutationSignature> ms, string strand /* forward - '+', reverse - '-', both - '+-' or '-+'*/, CHumanGenome* phuman, CDNAPos& res, int direction=1);
    //CDNAPos NextMotif(CDNAPos pos, char** motifsarr, int* strandarr, int motifsnum, int motiflen, CHumanGenome* phuman, int end, int includeCurrentPos=0);
    //unsigned long CountMotifGenome(set<string> motifs, CHumanGenome* phuman);
};

#endif /* mutsignature_hpp */
