//
//  bioread.hpp
//  APOBECXP
//
//  Created by Marat Kazanov on 02/04/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#ifndef bioread_hpp
#define bioread_hpp

#include <stdio.h>
#include <string>
#include "biogenome.hpp"

class CReader{
public:
    CReader(){};
    void ReadFile2Memory(string path, string format, char** buffer, int len);
};

#endif /* bioread_hpp */
