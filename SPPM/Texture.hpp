//
//  Texture.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#ifndef Texture_hpp
#define Texture_hpp

#include <stdio.h>

#include "Vector.hpp"

class Texture {
    Vector** image;
    int width, height;
public:
    Texture(char *filename) {
        import(filename);
    }
    void import(char *filename);
    Vector query(double x, double y);
};

class TextureMapper {
    Texture *texture;
    Vector color;
    double xx, xy, xz, xb, yx, yy, yz, yb;
public:
    TextureMapper(Texture *texture, double xx, double xy, double xz, double xb, double yx, double yy, double yz, double yb) {
        this->texture = texture;
        this->xx = xx;
        this->xy = xy;
        this->xz = xz;
        this->xb = xb;
        this->yx = yx;
        this->yy = yy;
        this->yz = yz;
        this->yb = yb;
    }
    TextureMapper(Vector color) {
        this->texture = nullptr;
        this->color = color;
    }
    Vector query(Vector p);
};


#endif /* Texture_hpp */