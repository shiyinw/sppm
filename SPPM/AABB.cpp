//
//  AABB.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 11/2/18.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#include "AABB.hpp"
#define INF 1e100
#define EPSILON 1e-6
using namespace std;

double AABB::intersect(const Ray &ray) {
    if (!(ray.s >= _min && ray.s <= _max)) { // outside
        double t = -1e100;
        if (fabs(ray.d.x) > 0)
            t = max(t, min((_min.x - ray.s.x) / ray.d.x, (_max.x - ray.s.x) / ray.d.x));
        if (fabs(ray.d.y) > 0)
            t = max(t, min((_min.y - ray.s.y) / ray.d.y, (_max.y - ray.s.y) / ray.d.y));
        if (fabs(ray.d.z) > 0)
            t = max(t, min((_min.z - ray.s.z) / ray.d.z, (_max.z - ray.s.z) / ray.d.z));
        if (t < -EPSILON) return 1e100;
        Vec3d pp = ray.s + ray.d * t; // intersection point
        if (!(pp >= _min && pp <= _max)) return INF; //inside
        return t;
    }
    else return -INF;
}
