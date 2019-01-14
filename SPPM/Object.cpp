//
//  Object.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#define BUFFER_SIZE 1024
#define EPSILON 1e-6

#include "Object.hpp"
#include "Random.hpp"

#include <iostream>
using namespace std;
void Object::importPly(char *filename,  TextureMapper *texture, int brdf) {
    FILE *file = fopen(filename, "r");
    char buffer[BUFFER_SIZE];
    numVertexes = numFaces = 0;
    while (fgets(buffer, BUFFER_SIZE, file)) {
        if (string(buffer) == "end_header\n") break;
        if (string(buffer).substr(0, 14) == "element vertex"){
            string line = string(buffer);
            line = line.substr(15, line.size()-1);
            numVertexes = stoi(line);
        }
        else if (string(buffer).substr(0, 12) == "element face"){
            string line = string(buffer);
            line = line.substr(13, line.size()-1);
            numFaces = stoi(line);
        }
    }
    vertexes = new Vec3d*[numVertexes];
    for (int i = 0; i < numVertexes; ++i) {
        double x, y, z;
        fscanf(file, "%lf%lf%lf", &x, &y, &z);
        fgets(buffer, BUFFER_SIZE, file);
        vertexes[i] = new Vec3d(x, y, z);
    }
    meshes = new Mesh*[numFaces];
    for (int i = 0; i < numFaces; ++i) {
        int n, a, b, c;
        fscanf(file, "%d", &n);
        if (n != 3)
            throw runtime_error("Only trianglar faces are supported!");
        fscanf(file, "%d%d%d", &a, &b, &c);
        meshes[i] = new TriMesh(vertexes[a], vertexes[b], vertexes[c], texture, brdf);
    }
    fclose(file);
    fprintf(stderr, "Imported object %s: %d vertexes and %d faces\n", filename, numVertexes, numFaces);
}


void Object::calcCenter() {
    center = new Vec3d(0, 0, 0);
    for (int i = 0; i < numVertexes; ++i)
        *center = *center + *vertexes[i];
    *center = *center / numVertexes;
}

void Object::printBox() {
    Vec3d min(1e100, 1e100, 1e100);
    Vec3d max = min * -1;
    for (int i = 0; i < numFaces; ++i) {
        min = ::min(min, meshes[i]->aabb->_min);
        max = ::max(max, meshes[i]->aabb->_max);
    }
    min.print();
    max.print();
}

void Object::scale(
                   double fxx, double fxy, double fxz, double fxb,
                   double fyx, double fyy, double fyz, double fyb,
                   double fzx, double fzy, double fzz, double fzb) {
    for (int i = 0; i < numVertexes; ++i) {
        Vec3d ver = *vertexes[i];
        vertexes[i]->x = fxx * ver.x + fxy * ver.y + fxz * ver.z + fxb;
        vertexes[i]->y = fyx * ver.x + fyy * ver.y + fyz * ver.z + fyb;
        vertexes[i]->z = fzx * ver.x + fzy * ver.y + fzz * ver.z + fzb;
    }
    for (int i = 0; i < numFaces; ++i){
        meshes[i]->scale(fxx, fxy, fxz, fxb, fyx, fyy, fyz, fyb, fzx, fzy, fzz, fzb);
        meshes[i]->updateAABB();
    }
    calcCenter();
}

void Object::rotXZ(double theta) {
    calcCenter();
    for (int i = 0; i < numVertexes; ++i) {
        Vec3d _d = *vertexes[i] - *center;
        *vertexes[i] = *center + Vec3d(cos(theta) * _d.x - sin(theta) * _d.z,  _d.y, sin(theta) * _d.x + cos(theta) * _d.z);
    }
    for (int i = 0; i < numFaces; ++i){
        meshes[i]->updateAABB();
    }
}


pair<double, Vec3d> TriMesh::intersect(Ray ray) {
    Vec3d E1 = *a - *b, E2 = *a - *c, S = *a - ray.s;
    double t = det(S, E1, E2);
    double beta = det(ray.d, S, E2);
    double gamma = det(ray.d, E1, S);
    double n = det(ray.d, E1, E2);
    t /= n;
    beta /= n;
    gamma /= n;
    if (!(-EPSILON <= beta && beta <= 1 + EPSILON && -EPSILON <= gamma && gamma <= 1 + EPSILON && beta + gamma <= 1 + EPSILON))
        t = -1;
    
    Vec3d norm = cross(*b - *a, *c - *a);
    if (dot(norm, ray.d) > 0)
        norm = norm * -1;
    norm.normalize();
    return make_pair(t, norm);
}

double TriMesh::intersectPlane(Ray ray) {
    Vec3d E1 = *a - *b, E2 = *a - *c, S = *a - ray.s;
    double t = det(S, E1, E2);
    double n = det(ray.d, E1, E2);
    t /= n;
    return t;
}

pair<double, Vec3d> SphereMesh::intersect(Ray ray) {
    Vec3d l = c - ray.s;
    double dis2 = l.norm2();
    int pos = 0;
    if (dis2 > r * r) pos = 1;
    else if (dis2 < r * r) pos = -1;
    double tp = dot(l, ray.d);
    if (pos > 0 && tp < 0) return make_pair(-1, Vec3d(0, 0, 0));
    double d2 = dis2 - tp * tp;
    if (d2 > r * r || pos == 0) return make_pair(-1, Vec3d(0, 0, 0));
    double t;
    if (pos > 0) t = tp - sqrt(r * r - d2);
    else if (pos < 0) t = tp + sqrt(r * r - d2);
    Vec3d norm = ray.s + ray.d * t - c;
    if (dot(norm, ray.d) > 0) norm = norm * -1;
    norm.normalize();
    return make_pair(t, norm);
}

pair<double, Vec3d> CircleMesh::intersect(Ray ray) {
    double t = 0;
    if (fabs(ray.d.y) < EPSILON && fabs(ray.s.y - c.y) > EPSILON)
        return make_pair(-1, Vec3d(0, 0, 0));
    else t = (c.y - ray.s.y) / ray.d.y;
    Vec3d p = ray.s + ray.d * t;
    if ((p.x - c.x)*(p.x - c.x) + (p.z - c.z)*(p.z - c.z) <= r * r)
        return make_pair(t, Vec3d(0, -1, 0));
    else return make_pair(-1, Vec3d(0, 0, 0));
}


Sphere::Sphere(Vec3d c, double r, TextureMapper *texture, int brdf) {
    numVertexes = 0;
    numFaces = 1;
    meshes = new Mesh*[1] {
        new SphereMesh(c, r, texture, brdf)
    };
    center = new Vec3d(c);
}


pair<double, Vec3d> intersect(Ray ray){
    return make_pair(1, ray.s);
}

double intersectPlane(Ray ray){
    return 1;
}
