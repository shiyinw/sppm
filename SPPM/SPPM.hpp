//
//  Camera.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 12/13/18.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#ifndef SPPM_hpp
#define SPPM_hpp

#include <stdio.h>

#include "Scene.hpp"
#include "Vec3d.hpp"

class Scene;
class HitPoint;


class SPPM {
    int w, h; // 画布大小
    int round, photon; // 已经迭代的轮数和每一轮的光子数
    Vec3d **canvas, light;
    Scene* scene;
    Vec3d s;
    std::vector<HitPoint*> *hitpoints;
public:
    SPPM(int w, int h, Scene *scene, Vec3d s, int nPhoton=200000){
        this->w = w;
        this->h = h;
        this->scene = scene;
        this->s = s;
        this->round = 0;
        this->photon = nPhoton;
        
        canvas = new Vec3d*[w];
        for (int i = 0; i < w; ++i)
            canvas[i] = new Vec3d[h];
        
        hitpoints = new vector<HitPoint*>;
        for (int u = 0; u < w; ++u)
            for (int v = 0; v < h; ++v)
                hitpoints->push_back(new HitPoint);
    }
    void render(int numRounds);
    void save(char *filename1, char *filename2);
    void load(char *filename, int round);
};

#endif /* SPPM_hpp */
