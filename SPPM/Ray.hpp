//
//  Ray.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/10/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#ifndef Ray_hpp
#define Ray_hpp

#include <stdio.h>

#include "Vector.hpp"

class Ray {
public:
    Vector s;
    Vector d;
    Ray(Vector s, Vector d) {
        this->s = s;
        this->d = d / sqrt(d.norm2());
    }
};

#endif /* Ray_hpp */
