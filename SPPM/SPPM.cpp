//
//  SPPM.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#include "SPPM.hpp"
#include "Scene.hpp"
#include <iostream>
#include <fstream>
using namespace std;

#define GAMMA 0.5
#define BUFFER_SIZE 1024

void SPPM::evaluateRadiance(int numRounds, int numPhotons) {
    for (int u = 0; u < w; ++u)
        for (int v = 0; v < h; ++v) {
            HitPoint *hp = (*hitpoints)[u * h + v];
            canvas[u][v] = hp->flux / (M_PI * hp->r2 * numPhotons * numRounds) + light * hp->fluxLight / numRounds;
            canvas[u][v].x = pow(canvas[u][v].x, GAMMA);
            canvas[u][v].y = pow(canvas[u][v].y, GAMMA);
            canvas[u][v].z = pow(canvas[u][v].z, GAMMA);
            canvas[u][v] = min(canvas[u][v], Vec3d(1, 1, 1));
            canvas[u][v] = max(canvas[u][v], Vec3d(0, 0, 0));
        }
}

void SPPM::render(int numRounds, int numPhotons) {
    scene->initializeObjectKDTree();
    
    int cx = w / 2, cy = h / 2;
    
    TriangularFace focusPlane(new Vec3d(0, 0, s.z + focus), new Vec3d(0, 1, s.z + focus), new Vec3d(1, 0, s.z + focus));
    
    // initialize hitpoints
    hitpoints = new vector<HitPoint*>;
    for (int u = 0; u < w; ++u)
        for (int v = 0; v < h; ++v)
            hitpoints->push_back(new HitPoint);
    
    light = Vec3d(1, 1, 1);
    const Vec3d weight_init = Vec3d(2.5, 2.5, 2.5);
    
    // SPPM
    for (int round = 0; round < numRounds; ++round) {
        fprintf(stderr, "Round %d/%d:\n", round + 1, numRounds);
        
        // ray tracing pass
        cout<<"Start Ray Tracing......"<<endl;
        for (int u = 0; u < w; ++u) {
            if(u%100==0){
                cout<<"Ray tracing pass "<<u<<"/"<<w<<endl;
            }
            for (int v = 0; v < h; ++v) {
                Vec3d p(double(u - cx) / w / fx, double(v - cy) / h / fy, 0);
                Ray ray(s, p - s);
                double t = focusPlane.intersectPlane(ray);
                Vec3d focusP = ray.s + ray.d * t;
                double theta = Utils::random(0, 2 * M_PI);
                ray.s = ray.s + Vec3d(cos(theta), sin(theta), 0) * aperture;
                ray.d = focusP - ray.s;
                ray.d.normalize();
                (*hitpoints)[u * h + v]->valid = false;
                (*hitpoints)[u * h + v]->d = ray.d * -1;
                scene->trace(ray, Vec3d(1, 1, 1), 1, (long long)round * (numPhotons + w * h) + u * h + v, (*hitpoints)[u * h + v]);
            }
        }
        
        scene->initializeHitpointKDTree(hitpoints);
        
        // photon tracing pass
        for (int i = 0; i < numPhotons; ++i) {
            Ray ray = scene->generateRay((long long)round * numPhotons + (round + 1) * w * h + i);
            scene->trace(ray, weight_init * light, 1, (long long)round * numPhotons + i);
        }
        fprintf(stderr, "\rPhoton tracing pass done\n");
        
        // save checkpoints
        if ((round+1)%1==0) {
            evaluateRadiance(round + 1, numPhotons);
            char filename[100];
            sprintf(filename, "checkpoints/checkpoint-%d.ppm", round + 1);
            save(filename);
        }
    }
    evaluateRadiance(numRounds, numPhotons);
}

void SPPM::save(char *filename) {
    FILE *file = fopen(filename, "w");
    fprintf(file, "P3\n%d %d\n255\n", w, h);
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j)
            fprintf(file, "%d %d %d ", int(canvas[j][h - i - 1].x * 255 + 0.5), int(canvas[j][h - i - 1].y * 255 + 0.5), int(canvas[j][h - i - 1].z * 255 + 0.5));
        fprintf(file, "\n");
    }
    fclose(file);
    fprintf(stderr, "Image saved to %s\n", filename);
}

void SPPM::load(char *filename) {
//    FILE *file = fopen(filename, "r");
//    char buffer[BUFFER_SIZE];
//    int weight, height;
//    fgets(buffer, BUFFER_SIZE, file);
//    fgets(buffer, BUFFER_SIZE, file);
//    assert(string(buffer)==std::to_string(w)+" "+std::to_string(h)+"\n");
//    fgets(buffer, BUFFER_SIZE, file);
//    assert(string(buffer)=="255\n");
//    fgets(buffer, BUFFER_SIZE, file);
//    cout<<buffer<<endl;
    
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }
    string r, g, b;
    inFile>>r;
    assert(r=="P3");
    inFile>>r;
    assert(r==to_string(w));
    inFile>>r;
    assert(r==to_string(h));
    inFile>>r;
    assert(r=="255");
    
    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j){
            inFile>>r>>g>>b;
            canvas[j][h - i - 1].x = stoi(r)/255;
            canvas[j][h - i - 1].y = stoi(g)/255;
            canvas[j][h - i - 1].z = stoi(b)/255;
        }
    }
    
    inFile.close();
}

SPPM::SPPM(int w, int h, Scene *scene, Vec3d s) {
    this->w = w;
    this->h = h;
    this->scene = scene;
    this->s = s;
    
    canvas = new Vec3d*[w];
    for (int i = 0; i < w; ++i)
        canvas[i] = new Vec3d[h];
}
