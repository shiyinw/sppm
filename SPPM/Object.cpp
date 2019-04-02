//
//  Object.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 11/3/18.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#define BUFFER_SIZE 1024
#define EPSILON 1e-6
#include <algorithm>
#include <fstream>
#include <iostream>

#include "Object.hpp"

#include <iostream>
using namespace std;
void Object::importPly(char *filename,  TextureMapper *texture, int brdf) {
    FILE *file = fopen(filename, "r");
    char buffer[BUFFER_SIZE];
    nv = nm = 0;
    while (fgets(buffer, BUFFER_SIZE, file)) {
        if (string(buffer) == "end_header\n") break;
        if (string(buffer).substr(0, 14) == "element vertex"){
            string line = string(buffer);
            line = line.substr(15, line.size()-1);
            nv = stoi(line);
        }
        else if (string(buffer).substr(0, 12) == "element face"){
            string line = string(buffer);
            line = line.substr(13, line.size()-1);
            nm = stoi(line);
        }
    }
    vertexes = new Vec3d*[nv];
    for (int i = 0; i < nv; ++i) {
        double x, y, z;
        fscanf(file, "%lf%lf%lf", &x, &y, &z);
        fgets(buffer, BUFFER_SIZE, file);
        vertexes[i] = new Vec3d(x, y, z);
    }
    meshes = new Mesh*[nm];
    for (int i = 0; i < nm; ++i) {
        int n, a, b, c;
        fscanf(file, "%d", &n);
        fscanf(file, "%d%d%d", &a, &b, &c);
        meshes[i] = new TriMesh(vertexes[a], vertexes[b], vertexes[c], texture, brdf);
    }
    fclose(file);
    fprintf(stderr, "%s: %d vertexes and %d meshes\n", filename, nv, nm);
    updateAABB();
}

void Object::importObj(char *filename,  TextureMapper *texture, int brdf) {
    FILE *file = fopen(filename, "r");
    nv = nm = 0;
    fscanf(file, "%d %d", &nv, &nm);
    vertexes = new Vec3d*[nv];
    meshes = new Mesh*[nm];
    
    char type;
    double a, b, c;
    int d, e, f;
    
    for (int i = 0; i < nv; i++) {
        fscanf(file, "\n%c %lf %lf %lf", &type, &a, &b, &c);
        vertexes[i] = new Vec3d(a, b, c);
    }
    
    for (int i = 0; i < nm; i++) {
        fscanf(file, "\n%c %d %d %d", &type, &d, &e, &f);
        // obj file start with indice 1
        meshes[i] = new TriMesh(vertexes[d-1], vertexes[e-1], vertexes[f-1], texture, brdf);
    }
    fclose(file);
    fprintf(stderr, "%s: %d vertexes and %d meshes\n", filename, nv, nm);
    updateAABB();
}


void Object::updateAABB() {
    Vec3d min(1e100, 1e100, 1e100);
    Vec3d max = min * -1;
    Vec3d center = Vec3d(0, 0, 0);
    for (int i = 0; i < nv; ++i) {
        min.x = std::min(vertexes[i]->x, min.x);
        min.y = std::min(vertexes[i]->y, min.y);
        min.z = std::min(vertexes[i]->z, min.z);
        max.x = std::max(vertexes[i]->x, max.x);
        max.y = std::max(vertexes[i]->y, max.y);
        max.z = std::max(vertexes[i]->z, max.z);
        center = center + *vertexes[i];
    }
    aabb->_min = min;
    aabb->_max = max;
    aabb->_center = center/nv;
}

void Object::locate(Vec3d locate_min, Vec3d locate_max){
    updateAABB();
    Vec3d ideal_center = (locate_min + locate_max)/2;
    Vec3d current_center = (this->aabb->_max + this->aabb->_min)/2;
    Vec3d current_size = this->aabb->_max - this->aabb->_min;
    Vec3d ideal_size = locate_max - locate_min;
    double scale_x = (current_size.x>0) ? ideal_size.x/current_size.x : 1;
    double scale_y = (current_size.y>0) ? ideal_size.y/current_size.y : 1;
    double scale_z = (current_size.z>0) ? ideal_size.z/current_size.z : 1;
    
    Vec3d b = ideal_center - current_center * Vec3d(scale_x, scale_y, scale_z);
    for (int i = 0; i < nv; ++i) {
        vertexes[i]->x = vertexes[i]->x * scale_x + b.x;
        vertexes[i]->y = vertexes[i]->y * scale_y + b.y;
        vertexes[i]->z = vertexes[i]->z * scale_z + b.z;
    }
    for (int i = 0; i< nm; ++i){
        meshes[i]->updateAABB();
    }
    this->updateAABB();
}

pair<double, Vec3d> TriMesh::intersect(Ray ray) {
    Vec3d E1 = *a - *b, E2 = *a - *c, S = *a - ray.s;
    double n = det(ray.d, E1, E2);
    double t = det(S, E1, E2)/n;
    double beta = det(ray.d, S, E2)/n;
    double gamma = det(ray.d, E1, S)/n;

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
    if (dis2 > r * r)
        pos = 1;
    else if (dis2 < r * r)
        pos = -1;
    double tp = dot(l, ray.d);
    if (pos > 0 && tp < 0)
        return make_pair(-1, Vec3d(0, 0, 0));
    double d2 = dis2 - tp * tp;
    if (d2 > r * r || pos == 0)
        return make_pair(-1, Vec3d(0, 0, 0)); // no intersection
    double t;
    if (pos > 0)
        t = tp - sqrt(r * r - d2);
    else if (pos < 0)
        t = tp + sqrt(r * r - d2);
    Vec3d norm = ray.s + ray.d * t - c;
    if (dot(norm, ray.d) > 0)
        norm = norm * -1;
    norm.normalize();
    return make_pair(t, norm);
}

// x^2=y(y-1)^2
// x^2 + z^2 = y (y-1)^2
// yrange[0, 1]
// xrange approx (-0.5, 0.5)
pair<double, Vec3d> WaterDropMesh::intersect(Ray ray) {
    double aabb_t = aabb->intersect(ray);
    if(aabb_t<-1000 || aabb_t>1000)
        return make_pair(-1, Vec3d(0, 0, 0));
    
    Vec3d relative_s = ray.s - position;
    double dis = b * b * b /a/a;
    double a3 = ray.d.y * ray.d.y * ray.d.y;
    double a2 = 3 * relative_s.y * ray.d.y * ray.d.y - 2 * ray.d.y * ray.d.y - ray.d.x * ray.d.x - ray.d.z * ray.d.z;
    double a1 = 3 * relative_s.y * relative_s.y * ray.d.y - 2 * relative_s.y * ray.d.y + b * b * ray.d.y - 2 * b * relative_s.y * ray.d.y - 2 * dis * relative_s.x * ray.d.x - 2 * dis * ray.d.z * relative_s.z;
    double a0 = relative_s.y * relative_s.y * relative_s.y - 2 * b * relative_s.y * relative_s.y + b * b * relative_s.y - dis * ray.d.x * ray.d.x - ray.d.z * ray.d.z * dis;
    double p = (3 * a3 * a1 - a2 * a2)/3/a3/a3;
    double q = (27 * a3 * a0 - 9 * a3 * a2 * a1 + 2 * a2 * a2)/27/a3/a3/a3;
    double delta = q*q/4 + p*p*p/27;
    if(delta>=0) // no intersection
        return make_pair(-1, Vec3d(0, 0, 0));
    
    double t = aabb_t;
    for(int i=0; i<100; i++){
        t = a3 * t * t * t + a2 * t * t + (a1+1) * t + a0;
    }
    Vec3d res = ray.s + ray.d * t;
    if(!((res.x * res.x + res.z * res.z - res.y * (res.y-1) * (res.y-1) < EPSILON)
         &&(res.x * res.x + res.z * res.z - res.y * (res.y-1) * (res.y-1) > -EPSILON)))
        return make_pair(-1, Vec3d(0, 0, 0));
    
    double slop = -(3 * res.y - 1) * (res.y -1)/2/sqrt(res.x*res.x + res.z*res.z);
    Vec3d inner = Vec3d(position.x, position.y, ray.s.y  - sqrt(relative_s.x * relative_s.x + relative_s.y * relative_s.y)/slop);
    Vec3d norm = res - inner;
    norm.normalize();
    return make_pair(t, norm);
}


Sphere::Sphere(Vec3d c, double r, TextureMapper *texture, int brdf) {
    nv = 0;
    nm = 1;
    meshes = new Mesh*[1] {
        new SphereMesh(c, r, texture, brdf)
    };
    this->aabb->_center = Vec3d(c);
    this->aabb->_max = c + Vec3d(r, r, r);
    this->aabb->_min = c + Vec3d(-r, -r, -r);
    this->c = c;
    this->r = r;
}
