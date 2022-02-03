//
//  bioread.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 02/04/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include <fstream>
#include "bioread.hpp"

using namespace std;

void CReader::ReadFile2Memory(string path, string format, char** buffer, int len)
{
    string line;
    int i=0;
    
    ifstream f(path.c_str());
    if (!f.is_open())
    {
        printf("File not exists\n");
        return;
    }
    
    (*buffer) = new char[len+1];
    
    if (format == "FASTA" || format == "fasta")
    {
        getline(f, line);
        if (line[0] != '>')
        {
            printf("Not a FASTA file\n");
            return;
        }
        while(getline(f, line))
        {
            if (line.length() != 0)
            {
                line.copy(&(*buffer)[i], line.length());
                i += line.length();
            }
        }
        printf("Loaded %i symbols\n", i);
        (*buffer)[len] = '\0';
    }
}


