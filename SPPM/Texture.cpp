//
//  Texture.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 12/10/18.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#include "Texture.hpp"

void Texture::import(char *filename) {
    FILE *file = fopen(filename, "r");
    fscanf(file, "%*s%d%d%*d", &width, &height);
    image = new Vec3d*[height];
    for (int i = 0; i < height; ++i) {
        image[i] = new Vec3d[width];
        for (int j = 0; j < width; ++j) {
            int r, g, b;
            fscanf(file, "%d%d%d", &r, &g, &b);
            image[i][j] = Vec3d(r / 255., g / 255., b / 255.);
        }
    }
    fclose(file);
}

Vec3d Texture::query(double x, double y) {
    int _x = x * height;
    int _y = y * width;
    if (_x >= 0 && _x < height && _y >= 0 && _y < width)
        return image[_x][_y];
    else
        return Vec3d(0, 0, 0);
}

Vec3d TextureMapper::query(Vec3d p) {
    if (!texture)
        return color;
    double x = dot(this->fx, p) + bx;
    double y = dot(this->fy, p) + by;
    return texture->query(x, y);
}

