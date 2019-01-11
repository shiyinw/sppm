//
//  HitPoint.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#define RADIUS 1e-5

#ifndef HitPoint_hpp
#define HitPoint_hpp

#include "Vector.hpp"
#include "BRDF.hpp"

using namespace std;

class HitPoint {
public:
    Vec3d p;
    Vec3d weight, flux, fluxLight;
    Vec3d d, norm;
    int n;
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
    Vec3d min, max;
    double maxr2;
    HitPointKDTreeNode *ls, *rs;
};

class HitPointKDTree {
    int n;
    HitPoint** hitpoints;
    HitPointKDTreeNode* build(int l, int r, int d);
    void del(HitPointKDTreeNode *p);
public:
    HitPointKDTreeNode *root;
    HitPointKDTree(vector<HitPoint*>* hitpoints);
    ~HitPointKDTree();
    void update(HitPointKDTreeNode *p, Vec3d photon, Vec3d weight, Vec3d d);
};

bool cmpHitPointX(HitPoint *a, HitPoint *b);
bool cmpHitPointY(HitPoint *a, HitPoint *b);
bool cmpHitPointZ(HitPoint *a, HitPoint *b);


#endif /* HitPoint_hpp */
