//
//  Scene.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#ifndef Scene_hpp
#define Scene_hpp

#include <stdio.h>
#include "Object.hpp"
#include "HitPoint.hpp"
#include "Ray.hpp"

class ObjectKDTreeNode {
public:
    Vec3d min, max;
    vector<Mesh*>* meshes;
    ObjectKDTreeNode *ls, *rs;
    int l, r;
    bool inside(Mesh *mesh);
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

class Scene {
    vector<Object*> objects;
    vector<HitPoint*> hitpoints;
    HitPointKDTree *hitpointsKDTree;
    ObjectKDTree *objectKDTree;
    Vec3d sourceP, sourceN;
    double sourceR;
    Vec3d sampleReflectedRay(Vec3d norm, int depth, long long i, double s = 1);
public:
    void addObject(Object* object);
    Scene(Vec3d _sourceP, double _sourceR, Vec3d _sourceN) :
    sourceP(_sourceP), sourceR(_sourceR), sourceN(_sourceN) { hitpointsKDTree = nullptr; }
    Ray generateRay(long long i);
    void trace(const Ray &ray, const Vec3d &weight, int depth, long long i, HitPoint *hp = nullptr);
    void initializeHitpointKDTree(vector<HitPoint*>* hitpoints);
    void initializeObjectKDTree();
};



#endif /* Scene_hpp */
