//
//  bioconst.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 02/04/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef biogenome_hpp
#define biogenome_hpp

#include <string>

using namespace std;

class CGenome{
public:
    string version;
    int chrCnt;
    string* chrName;
    int* chrLen;
    char** dna;
    CGenome(){};
    string dnaSubstr(int chrNum, unsigned long startpos, unsigned long endpos);
    void Read2Memory(string pathPrefix, string pathSuffix, string format);
};

#endif /* biogenome_hpp */
