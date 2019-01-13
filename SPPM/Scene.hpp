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
#include "KDTree.hpp"

class Object;
class HitPointKDTree;
class ObjectKDTree;
class HitPoint;
struct Ray;

class Scene {
    std::vector<Object*> objects;
    std::vector<HitPoint*> hitpoints;
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
    void initializeHitpointKDTree(std::vector<HitPoint*>* hitpoints);
    void initializeObjectKDTree();
};



#endif /* Scene_hpp */
