//
//  HitPoint.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#define ALPHA 0.7
#define EPSILON 1e-6

#include "KDTree.hpp"

HitPointKDTreeNode* HitPointKDTree::build(int l, int r, int d) {
    HitPointKDTreeNode *p = new HitPointKDTreeNode;
    p->min = Vec3d(1e100, 1e100, 1e100);
    p->max = p->min * (-1);
    p->maxr2 = 0;
    for (int i = l; i <= r; ++i) {
        p->min = min(p->min, hitpoints[i]->p);
        p->max = max(p->max, hitpoints[i]->p);
        p->maxr2 = max(p->maxr2, hitpoints[i]->r2);
    }
    int m = (l + r) >> 1;
    if (d == 0)
        nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, cmpHitPointX);
    else if (d == 1)
        nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, cmpHitPointY);
    else
        nth_element(hitpoints + l, hitpoints + m, hitpoints + r + 1, cmpHitPointZ);
    p->hitpoint = hitpoints[m];
    if (l <= m - 1) p->ls = build(l, m - 1, (d + 1) % 3); else p->ls = nullptr;
    if (m + 1 <= r) p->rs = build(m + 1, r, (d + 1) % 3); else p->rs = nullptr;
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

void HitPointKDTree::update(HitPointKDTreeNode * p, Vec3d photon, Vec3d weight, Vec3d d) {
    if (!p) return;
    double mind = 0, maxd = 0;
    if (photon.x > p->max.x) mind += (photon.x - p->max.x)*(photon.x - p->max.x);
    if (photon.x < p->min.x) mind += (p->min.x - photon.x)*(p->min.x - photon.x);
    if (photon.y > p->max.y) mind += (photon.y - p->max.y)*(photon.y - p->max.y);
    if (photon.y < p->min.y) mind += (p->min.y - photon.y)*(p->min.y - photon.y);
    if (photon.z > p->max.z) mind += (photon.z - p->max.z)*(photon.z - p->max.z);
    if (photon.z < p->min.z) mind += (p->min.z - photon.z)*(p->min.z - photon.z);
    if (mind > p->maxr2) return;
    if (p->hitpoint->valid && (photon - p->hitpoint->p).norm2() <= p->hitpoint->r2) {
        HitPoint *hp = p->hitpoint;
        double factor = (hp->n * ALPHA + ALPHA) / (hp->n * ALPHA + 1.);
        Vec3d dr = d - hp->norm * (2 * dot(d, hp->norm));
        double rho = hp->brdf.rho_d + hp->brdf.rho_s * pow(dot(dr, hp->d), hp->brdf.phong_s);
        if (rho < 0) rho = 0;
        else if (rho > 1) rho = 1;
        hp->n++;
        hp->r2 *= factor;
        hp->flux = (hp->flux + hp->weight * weight * rho) * factor;
    }
    if (p->ls) update(p->ls, photon, weight, d);
    if (p->rs) update(p->rs, photon, weight, d);
    p->maxr2 = p->hitpoint->r2;
    if (p->ls && p->ls->hitpoint->r2 > p->maxr2)
        p->maxr2 = p->ls->hitpoint->r2;
    if (p->rs && p->rs->hitpoint->r2 > p->maxr2)
        p->maxr2 = p->rs->hitpoint->r2;
}

bool cmpHitPointX(HitPoint *a, HitPoint *b) {
    return a->p.x < b->p.x;
}

bool cmpHitPointY(HitPoint *a, HitPoint *b) {
    return a->p.y < b->p.y;
}

bool cmpHitPointZ(HitPoint *a, HitPoint *b) {
    return a->p.z < b->p.z;
}


bool ObjectKDTreeNode::inside(Mesh *mesh) {
    Vec3d faceMin = mesh->min();
    Vec3d faceMax = mesh->max();
    return (faceMin.x < max.x || (faceMin.x == max.x && faceMin.x == faceMax.x))
    && (faceMax.x > min.x || (faceMax.x == min.x && faceMin.x == faceMax.x))
    && (faceMin.y < max.y || (faceMin.y == max.y && faceMin.y == faceMax.y))
    && (faceMax.y > min.y || (faceMax.y == min.y && faceMin.y == faceMax.y))
    && (faceMin.z < max.z || (faceMin.z == max.z && faceMin.z == faceMax.z))
    && (faceMax.z > min.z || (faceMax.z == min.z && faceMin.z == faceMax.z));
}

ObjectKDTreeNode* ObjectKDTree::build(int depth, int d, vector<Mesh*>* meshes, Vec3d min, Vec3d max) {
    ObjectKDTreeNode *p = new ObjectKDTreeNode;
    p->min = min;
    p->max = max;
    Vec3d maxL, minR;
    if (d == 0) {
        maxL = Vec3d((p->min.x + p->max.x) / 2, p->max.y, p->max.z);
        minR = Vec3d((p->min.x + p->max.x) / 2, p->min.y, p->min.z);
    }
    else if (d == 1) {
        maxL = Vec3d(p->max.x, (p->min.y + p->max.y) / 2, p->max.z);
        minR = Vec3d(p->min.x, (p->min.y + p->max.y) / 2, p->min.z);
    }
    else {
        maxL = Vec3d(p->max.x, p->max.y, (p->min.z + p->max.z) / 2);
        minR = Vec3d(p->min.x, p->min.y, (p->min.z + p->max.z) / 2);
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
        min = ::min(min, mesh->min());
        max = ::max(max, mesh->max());
    }
    root = build(1, 0, meshes, min, max);
    this->meshes = new vector<Mesh*>;
    getFaces(root, this->meshes);
}

double ObjectKDTree::getCuboidIntersection(ObjectKDTreeNode *p, Ray ray) {
    if (!(ray.s >= p->min && ray.s <= p->max)) { // outside
        double t = -1e100;
        if (fabs(ray.d.x) > 0)
            t = max(t, min((p->min.x - ray.s.x) / ray.d.x, (p->max.x - ray.s.x) / ray.d.x));
        if (fabs(ray.d.y) > 0)
            t = max(t, min((p->min.y - ray.s.y) / ray.d.y, (p->max.y - ray.s.y) / ray.d.y));
        if (fabs(ray.d.z) > 0)
            t = max(t, min((p->min.z - ray.s.z) / ray.d.z, (p->max.z - ray.s.z) / ray.d.z));
        if (t < -EPSILON) return 1e100;
        Vec3d pp = ray.s + ray.d * t;
        if (!(pp >= p->min && pp <= p->max)) return 1e100;
        return t;
    }
    else return -1e100;
}

void ObjectKDTree::getIntersection(ObjectKDTreeNode *p, Ray ray, Mesh* &nextMesh, double &tMin, Vec3d &norm) {
    for (int i = 0; i < p->meshes->size(); ++i) {
        pair<double, Vec3d> r = (*p->meshes)[i]->intersect(ray);
        double t = r.first;
        if (t > 0 && t < tMin) {
            tMin = t;
            nextMesh = (*p->meshes)[i];
            norm = r.second;
        }
    }
    
    double tl = p->ls ? getCuboidIntersection(p->ls, ray) : 1e100;
    double tr = p->rs ? getCuboidIntersection(p->rs, ray) : 1e100;
    if (tl < tr) {
        if (tMin <= tl) return;
        if (p->ls) getIntersection(p->ls, ray, nextMesh, tMin, norm);
        if (tMin <= tr) return;
        if (p->rs) getIntersection(p->rs, ray, nextMesh, tMin, norm);
    }
    else {
        if (tMin <= tr) return;
        if (p->rs) getIntersection(p->rs, ray, nextMesh, tMin, norm);
        if (tMin <= tl) return;
        if (p->ls) getIntersection(p->ls, ray, nextMesh, tMin, norm);
    }
}
