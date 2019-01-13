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
    Vec3d _min;
    Vec3d _max;
public:
    // triangle mesh
    AABB(Vec3d a, Vec3d b, Vec3d c){
        _min.x = std::min(a.x, std::min(b.x, c.x));
        _min.y = std::min(a.y, std::min(b.y, c.y));
        _min.z = std::min(a.z, std::min(b.z, c.z));
        _max.x = std::max(a.x, std::max(b.x, c.x));
        _max.y = std::max(a.y, std::max(b.y, c.y));
        _max.z = std::max(a.z, std::max(b.z, c.z));
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
    }
    
    // 3d cube and 2d square
    AABB(Vec3d min, Vec3d max){
        _max = max;
        _min = min;
    }
    
//    bool ObjectKDTreeNode::inside(Mesh *mesh) {
//        Vec3d faceMin = mesh->min();
//        Vec3d faceMax = mesh->max();
//        return (faceMin.x < max.x || (faceMin.x == max.x && faceMin.x == faceMax.x))
//        && (faceMax.x > min.x || (faceMax.x == min.x && faceMin.x == faceMax.x))
//        && (faceMin.y < max.y || (faceMin.y == max.y && faceMin.y == faceMax.y))
//        && (faceMax.y > min.y || (faceMax.y == min.y && faceMin.y == faceMax.y))
//        && (faceMin.z < max.z || (faceMin.z == max.z && faceMin.z == faceMax.z))
//        && (faceMax.z > min.z || (faceMax.z == min.z && faceMin.z == faceMax.z));
//    }

};

class Mesh {
public:
    Object *object;
    TextureMapper *texture;
    int brdf;
    virtual Vec3d min() = 0;
    virtual Vec3d max() = 0;
    virtual Vec3d center() = 0;
    virtual void scale(
                       double fxx, double fxy, double fxz, double fxb,
                       double fyx, double fyy, double fyz, double fyb,
                       double fzx, double fzy, double fzz, double fzb
                       ) = 0;
    virtual pair<double, Vec3d> intersect(Ray ray) = 0;
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
    }
    Vec3d min();
    Vec3d max();
    Vec3d center();
    void scale(
               double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               ){}
    pair<double, Vec3d> intersect(Ray ray);
    double intersectPlane(Ray ray);
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
    }
    Vec3d min();
    Vec3d max();
    Vec3d center();
    void scale(
               double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               ){}
    pair<double, Vec3d> intersect(Ray ray);
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
    }
    Vec3d min();
    Vec3d max();
    Vec3d center();
    void scale(
               double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               ){}
    pair<double, Vec3d> intersect(Ray ray);
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
