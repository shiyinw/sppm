//
//  HitPoint.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 11/29/18.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#define ALPHA 0.99
#define EPSILON 1e-6
#define INF 1e100
#define MIN_RADIUS 1e-8

#include "KDTree.hpp"
#include <iostream>
using namespace std;

HitPointKDTreeNode* HitPointKDTree::build(int l, int r, int d) {
    HitPointKDTreeNode *p = new HitPointKDTreeNode;
    p->maxr2 = 0;
    for (int i = l; i <= r; ++i) {
        p->aabb._min = min(p->aabb._min, hitpoints[i]->p);
        p->aabb._max = max(p->aabb._max, hitpoints[i]->p);
        p->maxr2 = max(p->maxr2, hitpoints[i]->r2);
    }
    int m = (l + r) >> 1;
    switch(d){
        case 0:
            nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, [](HitPoint *a, HitPoint *b){return a->p.x < b->p.x;});
            break;
        case 1:
            nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, [](HitPoint *a, HitPoint *b){return a->p.y < b->p.y;});
            break;
        case 2:
            nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, [](HitPoint *a, HitPoint *b){return a->p.z < b->p.z;});
            break;
        default:
            printf("KD axis error\n");
            exit(1);
    }
    p->hitpoint = hitpoints[m];
    p->ls = (l <= m - 1)? build(l, m - 1, (d + 1) % 3):nullptr;
    p->rs = (m + 1 <= r)? build(m + 1, r, (d + 1) % 3):nullptr;
    return p;
}

HitPointKDTree::HitPointKDTree(vector<HitPoint*>* hitpoints) {
    n = (int)hitpoints->size();
    this->hitpoints = new HitPoint*[n];
    for (int i = 0; i < n; ++i)
        this->hitpoints[i] = (*hitpoints)[i];
    root = build(0, n - 1, 0);
}

void HitPointKDTree::del(HitPointKDTreeNode *p) {
    if (p->ls) del(p->ls);
    if (p->rs) del(p->rs);
    delete p;
}

HitPointKDTree::~HitPointKDTree() {
    if (!root) return;
    del(root);
    delete[] hitpoints;
}

void HitPointKDTree::update(HitPointKDTreeNode * p, Vec3d photon, Vec3d color, Vec3d d) {
    if (!p) return;
    double mind = 0;
    for(int i=0; i< 3; i++){
        if (photon._p[i] > p->aabb._max._p[i]) mind += (photon._p[i] - p->aabb._max._p[i])*(photon._p[i] - p->aabb._max._p[i]);
        if (photon._p[i] < p->aabb._min._p[i]) mind += (p->aabb._min._p[i] - photon._p[i])*(p->aabb._min._p[i] - photon._p[i]);
    }
    if (mind > p->maxr2) return;
    if (p->hitpoint->valid && (photon - p->hitpoint->p).norm2() <= p->hitpoint->r2) {
        HitPoint *hp = p->hitpoint;
        double factor = (hp->n + ALPHA) / (hp->n + 1.);
        if(hp->r2<MIN_RADIUS)
            factor = 1;
        Vec3d dr = d - hp->norm * (2 * dot(d, hp->norm));
        double rho = hp->brdf.rho_d + hp->brdf.rho_s * pow(dot(dr, hp->dir), hp->brdf.phong_s);
        if (rho < 0) rho = 0;
        else if (rho > 1) rho = 1;
        hp->n++;
        hp->r2 *= sqrt(factor);
        hp->flux = (hp->flux + hp->color * color * rho) * factor;
    }
    if (p->ls) update(p->ls, photon, color, d);
    if (p->rs) update(p->rs, photon, color, d);
    p->maxr2 = p->hitpoint->r2;
    if (p->ls && p->ls->hitpoint->r2 > p->maxr2)
        p->maxr2 = p->ls->hitpoint->r2;
    if (p->rs && p->rs->hitpoint->r2 > p->maxr2)
        p->maxr2 = p->rs->hitpoint->r2;
}


bool ObjectKDTreeNode::inside(Mesh *mesh) {
    Vec3d faceMin = mesh->aabb->_min;
    Vec3d faceMax = mesh->aabb->_max;
    
    for(int i=0; i<3; i++){
        if(!(faceMin._p[i] < aabb._max._p[i] || (faceMin._p[i] == aabb._max._p[i] && faceMin._p[i] == faceMax._p[i])))
            return false;
        if(!(faceMax._p[i] > aabb._min._p[i] || (faceMax._p[i] == aabb._min._p[i] && faceMin._p[i] == faceMax._p[i])))
            return false;
    }
    return true;
}

ObjectKDTreeNode* ObjectKDTree::build(int depth, int d, vector<Mesh*>* meshes, Vec3d min, Vec3d max) {
    ObjectKDTreeNode *p = new ObjectKDTreeNode;
    p->aabb._min = min;
    p->aabb._max = max;
    Vec3d maxL, minR;
    switch(d){
        case 0:
            maxL = Vec3d((p->aabb._min.x + p->aabb._max.x) / 2, p->aabb._max.y, p->aabb._max.z);
            minR = Vec3d((p->aabb._min.x + p->aabb._max.x) / 2, p->aabb._min.y, p->aabb._min.z);
            break;
        case 1:
            maxL = Vec3d(p->aabb._max.x, (p->aabb._min.y + p->aabb._max.y) / 2, p->aabb._max.z);
            minR = Vec3d(p->aabb._min.x, (p->aabb._min.y + p->aabb._max.y) / 2, p->aabb._min.z);
            break;
        case 2:
            maxL = Vec3d(p->aabb._max.x, p->aabb._max.y, (p->aabb._min.z + p->aabb._max.z) / 2);
            minR = Vec3d(p->aabb._min.x, p->aabb._min.y, (p->aabb._min.z + p->aabb._max.z) / 2);
            break;
    }
    p->meshes = new vector<Mesh*>;
    for (auto face : *meshes)
        if (p->inside(face))
            p->meshes->push_back(face);
    
    const int max_faces = 8;
    const int max_depth = 24;
    
    if (p->meshes->size() > max_faces && depth < max_depth) {
        p->ls = build(depth + 1, (d + 1) % 3, p->meshes, min, maxL);
        p->rs = build(depth + 1, (d + 1) % 3, p->meshes, minR, max);
        
        vector<Mesh*> *meshL = p->ls->meshes, *meshR = p->rs->meshes;
        map<Mesh*, int> cnt;
        for (auto mesh : *meshL) cnt[mesh]++;
        for (auto mesh : *meshR) cnt[mesh]++;
        p->ls->meshes = new vector<Mesh*>;
        p->rs->meshes = new vector<Mesh*>;
        p->meshes->clear();
        for (auto mesh : *meshL)
            if (cnt[mesh] == 1)
                p->ls->meshes->push_back(mesh);
            else
                p->meshes->push_back(mesh);
        for (auto mesh : *meshR)
            if (cnt[mesh] == 1)
                p->rs->meshes->push_back(mesh);
    }
    else
        p->ls = p->rs = nullptr;
    return p;
}

void ObjectKDTree::getFaces(ObjectKDTreeNode *p, vector<Mesh*>* meshes) {
    p->l = (int)meshes->size();
    for (auto mesh : *(p->meshes))
        meshes->push_back(mesh);
    p->r = (int)meshes->size();
    if (p->ls) getFaces(p->ls, meshes);
    if (p->rs) getFaces(p->rs, meshes);
}

ObjectKDTree::ObjectKDTree(vector<Mesh*>* meshes) {
    Vec3d min = Vec3d(1e100, 1e100, 1e100);
    Vec3d max = min * -1;
    for (auto mesh : *meshes) {
        min = ::min(min, mesh->aabb->_min);
        max = ::max(max, mesh->aabb->_max);
    }
    root = build(1, 0, meshes, min, max);
    this->meshes = new vector<Mesh*>;
    getFaces(root, this->meshes);
}

void ObjectKDTree::getIntersection(ObjectKDTreeNode *p, Ray ray, Mesh* &nextMesh, double &distance, Vec3d &norm) {
    for (int i = 0; i < p->meshes->size(); ++i) {
        pair<double, Vec3d> r = (*p->meshes)[i]->intersect(ray);
        double t = r.first;
        if (t > 0 && t < distance) {
            distance = t;
            nextMesh = (*p->meshes)[i];
            norm = r.second;
        }
    }
    double tl = p->ls ? p->ls->aabb.intersect(ray) : 1e100;
    double tr = p->rs ? p->rs->aabb.intersect(ray) : 1e100;
    
    if(distance<=tr && distance<=tl)
        return;
    
    if (tl < tr) {
        getIntersection(p->ls, ray, nextMesh, distance, norm);
        if (tr<distance) getIntersection(p->rs, ray, nextMesh, distance, norm);
    }
    else {
        getIntersection(p->rs, ray, nextMesh, distance, norm);
        if (tl<distance) getIntersection(p->ls, ray, nextMesh, distance, norm);
    }
}
