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

#include "Vec3d.hpp"

class Texture {
    Vec3d** image;
    int width, height;
public:
    Texture(char *filename) {
        import(filename);
    }
    void import(char *filename);
    Vec3d query(double x, double y);
};

class TextureMapper {
    Texture *texture;
    Vec3d color;
    Vec3d fx, fy;
    double bx, by;
public:
    TextureMapper(Texture *texture, Vec3d fx, Vec3d fy, double bx, double by) {
        this->texture = texture;
        this->bx = bx;
        this->by = by;
        this->fx.assign(fx);
        this->fy.assign(fy);
    }
    TextureMapper(Vec3d color) {
        this->texture = nullptr;
        this->color = color;
    }
    Vec3d query(Vec3d p);
};


#endif /* Texture_hpp */
