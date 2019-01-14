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
    updateAABB();
    calcCenter();
}

void Object::importObj(char *filename,  TextureMapper *texture, int brdf) {
    FILE *file = fopen(filename, "r");
    numVertexes = numFaces = 0;
    fscanf(file, "%d %d", &numVertexes, &numFaces);
    vertexes = new Vec3d*[numVertexes];
    meshes = new Mesh*[numFaces];
    
    char type;
    double a, b, c;
    int d, e, f;
    
    for (int i = 0; i < numVertexes; i++) {
        fscanf(file, "\n%c %lf %lf %lf", &type, &a, &b, &c);
        vertexes[i] = new Vec3d(a, b, c);
    }
    
    for (int i = 0; i < numFaces; i++) {
        fscanf(file, "\n%c %d %d %d", &type, &d, &e, &f);
        // obj file start with indice 1
        meshes[i] = new TriMesh(vertexes[d-1], vertexes[e-1], vertexes[f-1], texture, brdf);
    }
    fclose(file);
    fprintf(stderr, "Imported object %s: %d vertexes and %d faces\n", filename, numVertexes, numFaces);
    updateAABB();
    calcCenter();
}


void Object::updateAABB() {
    Vec3d min(1e100, 1e100, 1e100);
    Vec3d max = min * -1;
    for (int i = 0; i < numVertexes; ++i) {
        min.x = std::min(vertexes[i]->x, min.x);
        min.y = std::min(vertexes[i]->y, min.y);
        min.z = std::min(vertexes[i]->z, min.z);
        max.x = std::max(vertexes[i]->x, max.x);
        max.y = std::max(vertexes[i]->y, max.y);
        max.z = std::max(vertexes[i]->z, max.z);
    }
    aabb->_min = min;
    aabb->_max = max;
    aabb->_center = (aabb->_min + aabb->_max)/2;
}


void Object::calcCenter() {
    center = new Vec3d(0, 0, 0);
    for (int i = 0; i < numVertexes; ++i)
        *center = *center + *vertexes[i];
    *center = *center / numVertexes;
}

void Object::locate(Vec3d locate_min, Vec3d locate_max){
    updateAABB();
    Vec3d ideal_center = (locate_min + locate_max)/2;
    Vec3d current_center = this->aabb->_center;
    Vec3d current_size = this->aabb->_max - this->aabb->_min;
    Vec3d ideal_size = locate_max - locate_min;
    double scale_x = (current_size.x>0) ? ideal_size.x/current_size.x : 1;
    double scale_y = (current_size.y>0) ? ideal_size.y/current_size.y : 1;
    double scale_z = (current_size.z>0) ? ideal_size.z/current_size.z : 1;
    
    Vec3d b = ideal_center - current_center * Vec3d(scale_x, scale_y, scale_z);
    for (int i = 0; i < numVertexes; ++i) {
        vertexes[i]->x = vertexes[i]->x * scale_x + b.x;
        vertexes[i]->y = vertexes[i]->y * scale_y + b.y;
        vertexes[i]->z = vertexes[i]->z * scale_z + b.z;
    }
    for (int i = 0; i< numFaces; ++i){
        meshes[i]->updateAABB();
    }
    this->updateAABB();
}



void Object::rotXZ(double theta) {
    calcCenter();
    updateAABB();
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
    this->aabb->_center = c;
    this->aabb->_max = c + Vec3d(r, r, r);
    this->aabb->_min = c + Vec3d(-r, -r, -r);
    this->c = c;
    this->r = r;
}
