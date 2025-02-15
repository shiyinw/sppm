//
//  HitPoint.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 11/29/18.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#define RADIUS 3e-6
#define INF 1e100

#ifndef KDTree_hpp
#define KDTree_hpp

#include "Vec3d.hpp"
#include "BRDF.hpp"
#include "Object.hpp"

using namespace std;
class Mesh;
struct Ray;

class HitPoint {
public:
    Vec3d p;
    Vec3d color, flux, fluxLight;
    Vec3d dir, norm;
    double n;
    BRDF brdf;
    double r2;
    bool valid;
    HitPoint() {
        flux = fluxLight = Vec3d(0, 0, 0);
        r2 = RADIUS;
        n = 0;
        valid = false;
    }
};

class HitPointKDTreeNode {
public:
    HitPoint *hitpoint;
    AABB aabb;
    double maxr2;
    HitPointKDTreeNode *ls, *rs;
    HitPointKDTreeNode(){
        aabb._min = Vec3d(INF, INF, INF);
        aabb._max = Vec3d(-INF, -INF, -INF);
    }
};

class HitPointKDTree {
    double n;
    HitPoint** hitpoints;
    HitPointKDTreeNode* build(int l, int r, int d);
    void del(HitPointKDTreeNode *p);
public:
    HitPointKDTreeNode *root;
    HitPointKDTree(vector<HitPoint*>* hitpoints);
    ~HitPointKDTree();
    void update(HitPointKDTreeNode *p, Vec3d photon, Vec3d weight, Vec3d d);
};

class ObjectKDTreeNode {
public:
    AABB aabb;
    vector<Mesh*>* meshes;
    ObjectKDTreeNode *ls, *rs;
    int l, r;
    bool inside(Mesh *mesh);
    ObjectKDTreeNode(){
        aabb._min = Vec3d(INF, INF, INF);
        aabb._max = Vec3d(-INF, -INF, -INF);
    }
};

class ObjectKDTree {
    int n;
    Vec3d** vertexes;
    ObjectKDTreeNode* build(int depth, int d, vector<Mesh*>* meshes, Vec3d min, Vec3d max);
    void getFaces(ObjectKDTreeNode *p, vector<Mesh*>* meshes);
public:
    ObjectKDTreeNode* root;
    vector<Mesh*> *meshes;
    ObjectKDTree(vector<Mesh*>* meshes);
    double getCuboidIntersection(ObjectKDTreeNode *p, Ray ray);
    void getIntersection(ObjectKDTreeNode *p, Ray ray, Mesh* &nextMesh, double &tMin, Vec3d &norm);
};



#endif /* HitPoint_hpp */
