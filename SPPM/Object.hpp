//
//  Object.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#ifndef Object_hpp
#define Object_hpp

#include <stdio.h>
//#include "Config.h"
#include "Texture.hpp"
#include "BRDF.hpp"

using namespace std;

class Object;
class Mesh;
class TriMesh;
class SphereMesh;
class CircleMesh;
class Sphere;
struct Ray;

// AABB box
class AABB{
public:
    Vec3d _min;
    Vec3d _max;
    Vec3d _center;
    // triangle mesh
    AABB(Vec3d a, Vec3d b, Vec3d c){
        _min = Vec3d(std::min(a.x, std::min(b.x, c.x)), std::min(a.y, std::min(b.y, c.y)), std::min(a.z, std::min(b.z, c.z)));
        _max = Vec3d(std::max(a.x, std::max(b.x, c.x)), std::max(a.y, std::max(b.y, c.y)), std::max(a.z, std::max(b.z, c.z)));
        _center = (_min+_max)/2;
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
         _center = c;
    }
    
    // 3d cube and 2d square
    AABB(Vec3d min, Vec3d max){
        _max = max;
        _min = min;
        _center = (_min+_max)/2;
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
};

class Mesh {
public:
    Object *object;
    TextureMapper *texture;
    int brdf;
    AABB *aabb;
    virtual void scale(
                       double fxx, double fxy, double fxz, double fxb,
                       double fyx, double fyy, double fyz, double fyb,
                       double fzx, double fzy, double fzz, double fzb
                       ) = 0;
    virtual pair<double, Vec3d> intersect(Ray ray) = 0;
    virtual void updateAABB() = 0;
};

class CylinderMesh : public Mesh{
public:
    Vec3d *c;
    double r, h;
    CylinderMesh(Vec3d *c, double r, double h, TextureMapper *image, int brdf=0){
        this->c = c;
        this->r = r;
        this->h = h;
        this->texture = image;
        this->brdf = brdf;
        this->aabb = new AABB(Vec3d(c->x-r, c->y-r, c->z), Vec3d(c->x+r, c->y+r, c->z+h));
    }
    pair<double, Vec3d> intersect(Ray ray);
    double intersectPlane(Ray ray);
    void updateAABB(){
        this->aabb = new AABB(Vec3d(this->c->x - this->r, this->c->y - this->r, this->c->z),
                              Vec3d(this->c->x + this->r, this->c->y + this->r, this->c->z + this->h));
    }
};

class TriMesh : public Mesh {
public:
    Vec3d *a, *b, *c;
    TriMesh(Vec3d *a, Vec3d *b, Vec3d *c, TextureMapper *texture = nullptr, int brdf = 0) {
        this->a = a;
        this->b = b;
        this->c = c;
        this->texture = texture;
        this->brdf = brdf;
        this->aabb = new AABB(*a, *b, *c);
    }
    void scale(
               double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               ){}
    pair<double, Vec3d> intersect(Ray ray);
    double intersectPlane(Ray ray);
    void updateAABB(){
        this->aabb = new AABB(*a, *b, *c);
    }
};

class SphereMesh : public Mesh {
public:
    Vec3d c;
    double r;
    SphereMesh(Vec3d c, double r, TextureMapper *texture, int brdf) {
        this->c = c;
        this->r = r;
        this->texture = texture;
        this->brdf = brdf;
        this->aabb = new AABB(c, r, true);
    }
    void scale(
               double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               ){}
    pair<double, Vec3d> intersect(Ray ray);
    void updateAABB(){
        this->aabb = new AABB(c, r, true);
    }
};

class CircleMesh : public Mesh {
public:
    Vec3d c;
    double r;
    CircleMesh(Vec3d c, double r, TextureMapper *texture, int brdf) {
        this->c = c;
        this->r = r;
        this->texture = texture;
        this->brdf = brdf;
        this->aabb = new AABB(c, r, false);
    }
    void scale(
               double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               ){}
    pair<double, Vec3d> intersect(Ray ray);
    void updateAABB(){
        this->aabb = new AABB(c, r, false);
    }
};

class Object {
public:
    Vec3d** vertexes;
    Vec3d* center;
    Mesh** meshes;
    int numVertexes, numFaces;
    Object() {
        center = nullptr;
    }
    void importPly(char *filename, TextureMapper *texture, int brdf);
    void calcCenter();
    void printBox();
    void scale(double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               );
    void rotXZ(double theta);
};

class Sphere : public Object {
public:
    Sphere(Vec3d c, double r, TextureMapper *texture, int brdf);
};

struct Ray {
    Vec3d s; //source
    Vec3d d; //direction
};


#endif /* Object_hpp */
