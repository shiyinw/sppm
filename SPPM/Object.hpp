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
#include "Ray.hpp"

using namespace std;

class Object;

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

enum { TRIANGULAR_FACE, SPHERE_FACE, BEZIER_FACE } FACE_TYPES;

class TriangularFace : public Mesh {
public:
    Vec3d *a, *b, *c;
    TriangularFace(Vec3d *a, Vec3d *b, Vec3d *c, TextureMapper *texture = nullptr, int brdf = 0) {
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

class SphereFace : public Mesh {
public:
    Vec3d c;
    double r;
    SphereFace(Vec3d c, double r, TextureMapper *texture, int brdf) {
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

class DiscFace : public Mesh {
public:
    Vec3d c;
    double r;
    DiscFace(Vec3d c, double r, TextureMapper *texture, int brdf) {
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

class BezierFace : public Mesh {
    int **binom;
    Vec3d P(double u, double v);
    double B(int n, int k, double u);
    double dB(int n, int k, double u);
    Vec3d F(Vec3d x, Ray ray);
    Vec3d d(Vec3d x, Ray ray);
public:
    int n, m;
    Vec3d **p, m_min, m_max, m_center;
    BezierFace(int n, int m, Vec3d **p, TextureMapper *texture, int brdf);
    Vec3d min();
    Vec3d max();
    Vec3d center();
    void scale(double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               );
    pair<double, Vec3d> intersect(Ray ray);
    //Vector intersectNorm(Ray ray);
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


#endif /* Object_hpp */
