//
//  Camera.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#ifndef SPPM_hpp
#define SPPM_hpp

#include <stdio.h>

#include "Scene.hpp"

class SPPM {
    int w, h;
    double fx, fy, aperture, focus;
    Vec3d **canvas, light;
    Scene* scene;
    Vec3d s;
    vector<HitPoint*> *hitpoints;
    void evaluateRadiance(int numRounds, int numPhotons);
public:
    SPPM(int w, int h, Scene *scene, Vec3d s);
    void setLens(double fx, double fy, double aperture, double focus) {
        this->fx = fx;
        this->fy = fy;
        this->aperture = aperture;
        this->focus = focus;
    }
    void render(int numRounds, int numPhotons = 20480);
    void save(char *filename);
    void load(char *filename);
};


#endif /* SPPM_hpp */
