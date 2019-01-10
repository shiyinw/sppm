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
    Vector min, max;
    vector<Face*>* faces;
    ObjectKDTreeNode *ls, *rs;
    int l, r;
    bool inside(Face *face);
};

class ObjectKDTree {
    int n;
    Vector** vertexes;
    ObjectKDTreeNode* build(int depth, int d, vector<Face*>* faces, Vector min, Vector max);
    void getFaces(ObjectKDTreeNode *p, vector<Face*>* faces);
public:
    ObjectKDTreeNode* root;
    vector<Face*> *faces;
    ObjectKDTree(vector<Face*>* faces);
    double getCuboidIntersection(ObjectKDTreeNode *p, Ray ray);
    void getIntersection(ObjectKDTreeNode *p, Ray ray, Face* &nextFace, double &tMin, Vector &norm);
};

class Scene {
    vector<Object*> objects;
    vector<HitPoint*> hitpoints;
    HitPointKDTree *hitpointsKDTree;
    ObjectKDTree *objectKDTree;
    Vector sourceP, sourceN;
    double sourceR;
    Vector sampleReflectedRay(Vector norm, int depth, long long i, double s = 1);
public:
    void addObject(Object* object);
    Scene(Vector _sourceP, double _sourceR, Vector _sourceN) :
    sourceP(_sourceP), sourceR(_sourceR), sourceN(_sourceN) { hitpointsKDTree = nullptr; }
    Ray generateRay(long long i);
    void trace(const Ray &ray, const Vector &weight, int depth, long long i, HitPoint *hp = nullptr);
    void initializeHitpointKDTree(vector<HitPoint*>* hitpoints);
    void initializeObjectKDTree();
};



#endif /* Scene_hpp */
