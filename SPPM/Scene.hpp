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
    std::vector<HitPoint*> hitpoints;
    Vec3d sourceP, sourceN;
    double sourceR;
    Vec3d sampleReflectedRay(Vec3d norm, int depth, long long i, double s = 1);
public:
    std::vector<Object*> objects;
    HitPointKDTree *hitpointsKDTree = nullptr;
    ObjectKDTree *objectKDTree = nullptr;
    void addObject(Object* object);
    Scene(Vec3d _sourceP, double _sourceR, Vec3d _sourceN) : sourceP(_sourceP), sourceR(_sourceR), sourceN(_sourceN) {}
    Ray generateRay(long long i);
    void photonTrace(const Ray &ray, const Vec3d &weight, int depth, long long power);
    void rayTrace(const Ray &ray, const Vec3d &weight, int depth, HitPoint *hp);
};


#endif /* Scene_hpp */
