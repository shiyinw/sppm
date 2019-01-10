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

class Face {
public:
    Object *object;
    TextureMapper *texture;
    int brdf;
    virtual Vector min() = 0;
    virtual Vector max() = 0;
    virtual Vector center() = 0;
    virtual void scale(
                       double fxx, double fxy, double fxz, double fxb,
                       double fyx, double fyy, double fyz, double fyb,
                       double fzx, double fzy, double fzz, double fzb
                       ) = 0;
    virtual pair<double, Vector> intersect(Ray ray) = 0;
};

enum { TRIANGULAR_FACE, SPHERE_FACE, BEZIER_FACE } FACE_TYPES;

class TriangularFace : public Face {
public:
    Vector *a, *b, *c;
    TriangularFace(Vector *a, Vector *b, Vector *c, TextureMapper *texture = nullptr, int brdf = 0) {
        this->a = a;
        this->b = b;
        this->c = c;
        this->texture = texture;
        this->brdf = brdf;
    }
    Vector min();
    Vector max();
    Vector center();
    void scale(
               double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               ){}
    pair<double, Vector> intersect(Ray ray);
    double intersectPlane(Ray ray);
};

class SphereFace : public Face {
public:
    Vector c;
    double r;
    SphereFace(Vector c, double r, TextureMapper *texture, int brdf) {
        this->c = c;
        this->r = r;
        this->texture = texture;
        this->brdf = brdf;
    }
    Vector min();
    Vector max();
    Vector center();
    void scale(
               double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               ){}
    pair<double, Vector> intersect(Ray ray);
};

class DiscFace : public Face {
public:
    Vector c;
    double r;
    DiscFace(Vector c, double r, TextureMapper *texture, int brdf) {
        this->c = c;
        this->r = r;
        this->texture = texture;
        this->brdf = brdf;
    }
    Vector min();
    Vector max();
    Vector center();
    void scale(
               double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               ){}
    pair<double, Vector> intersect(Ray ray);
};

class BezierFace : public Face {
    int **binom;
    Vector P(double u, double v);
    double B(int n, int k, double u);
    double dB(int n, int k, double u);
    Vector F(Vector x, Ray ray);
    Vector d(Vector x, Ray ray);
public:
    int n, m;
    Vector **p, m_min, m_max, m_center;
    BezierFace(int n, int m, Vector **p, TextureMapper *texture, int brdf);
    Vector min();
    Vector max();
    Vector center();
    void scale(
               double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               );
    pair<double, Vector> intersect(Ray ray);
    //Vector intersectNorm(Ray ray);
};

class Object {
public:
    Vector** vertexes;
    Vector* center;
    Face** faces;
    int numVertexes, numFaces;
    Object() {
        center = nullptr;
    }
    void importPly(char *filename, TextureMapper *texture, int brdf);
    void calcCenter();
    void printBox();
    void scale(
               double fxx, double fxy, double fxz, double fxb,
               double fyx, double fyy, double fyz, double fyb,
               double fzx, double fzy, double fzz, double fzb
               );
    void rotXZ(double theta);
};

class Sphere : public Object {
public:
    Sphere(Vector c, double r, TextureMapper *texture, int brdf);
};


#endif /* Object_hpp */
