//
//  Object.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 11/3/18.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#ifndef Object_hpp
#define Object_hpp

#define INF 1e100

#include <stdio.h>
#include "Texture.hpp"
#include "BRDF.hpp"
#include "AABB.hpp"

using namespace std;

class Object;
class Mesh;
class TriMesh;
class SphereMesh;
class CircleMesh;
class Sphere;

class Mesh {
public:
    Object *object;
    TextureMapper *texture;
    int brdf;
    AABB *aabb;
    virtual pair<double, Vec3d> intersect(Ray ray) = 0;
    virtual void updateAABB() = 0;
    Mesh(){
        aabb = new AABB();
    }
};

class WaterDropMesh : public Mesh{
public:
    Vec3d position;
    double a, b;
    WaterDropMesh(Vec3d pos, double xs, double ys, TextureMapper *texture = nullptr, int brdf = 0){
        this->position = pos;
        this->a = xs;
        this->b = ys;
        this->texture = texture;
        this->brdf = brdf;
    }
    pair<double, Vec3d> intersect(Ray ray);
    void updateAABB(){
        aabb->_min = position + Vec3d(-0.5*a, -0.5*a, 0);
        aabb->_min = position + Vec3d(0.5*a, 0.5*a, b);
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
    pair<double, Vec3d> intersect(Ray ray);
    void updateAABB(){
        this->aabb = new AABB(c, r, true);
    }
};

class Object {
public:
    Vec3d** vertexes;
    Mesh** meshes;
    AABB* aabb;
    int nv, nm;
    Object() {
        aabb = new AABB(Vec3d(-1e100, -1e100, -1e100), Vec3d(1e100, 1e100, 1e100));
    }
    void importPly(char *filename, TextureMapper *texture, int brdf);
    void importObj(char *filename, TextureMapper *texture, int brdf);
    void updateAABB();
    void locate(Vec3d locate_min, Vec3d locate_max);
};

class Sphere : public Object {
public:
    Vec3d c;
    double r;
    Sphere(Vec3d c, double r, TextureMapper *texture, int brdf);
};



#endif /* Object_hpp */
