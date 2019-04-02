//
//  AABB.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 11/2/18.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#ifndef AABB_hpp
#define AABB_hpp

#include <stdio.h>
#include "Vec3d.hpp"

class AABB;

struct Ray {
    Vec3d s; //source
    Vec3d d; //direction
};


// AABB box
class AABB{
public:
    Vec3d _min;
    Vec3d _max;
    Vec3d _center;
    // triangle mesh
    AABB(){}
    AABB(Vec3d a, Vec3d b, Vec3d c){
        _min = Vec3d(std::min(a.x, std::min(b.x, c.x)), std::min(a.y, std::min(b.y, c.y)), std::min(a.z, std::min(b.z, c.z)));
        _max = Vec3d(std::max(a.x, std::max(b.x, c.x)), std::max(a.y, std::max(b.y, c.y)), std::max(a.z, std::max(b.z, c.z)));
        _center = (_min + _max)/2;
    }
    
    // 3d ball and 2d ball
    AABB(Vec3d c, double r, bool ball){
        if(ball){
            _max = Vec3d(c.x + r, c.y + r, c.z + r);
            _min = Vec3d(c.x - r, c.y - r, c.z - r);
        }
        else{
            _max = Vec3d(c.x + r, c.y + r, c.z);
            _min = Vec3d(c.x - r, c.y - r, c.z);
        }
        _center = (_min + _max)/2;
    }
    
    // 3d cube and 2d square
    AABB(Vec3d min, Vec3d max){
        _max = max;
        _min = min;
        _center = (_min + _max)/2;
    }
    
    bool inside(Vec3d kdmin, Vec3d kdmax) {
        for(int i=0; i<3; i++){
            if(!(_min._p[i] < kdmax._p[i] || (_min._p[i] == kdmax._p[i] && _min._p[i] == _max._p[i])))
                return false;
            if(!(_max._p[i] > kdmin._p[i] || (_max._p[i] == kdmin._p[i] && _min._p[i] == _max._p[i])))
                return false;
        }
        return true;
    }
    
    void print(){
        printf("AABB\n");
        _min.print();
        _max.print();
        _center.print();
    }
    double intersect(const Ray &ray);
};



#endif /* AABB_hpp */
