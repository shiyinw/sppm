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

#define BUFFER_SIZE 1024

struct Ray;

void SPPM::render(int numRounds) {
    
    scene->initializeObjectKDTree();
    int cx = w / 2, cy = h / 2;
    
    TriMesh focusPlane(new Vec3d(0, 0, s.z + focus), new Vec3d(0, 1, s.z + focus), new Vec3d(1, 0, s.z + focus));
    
    light = Vec3d(1, 1, 1);
    const Vec3d weight_init = Vec3d(2.5, 2.5, 2.5);
    
    // SPPM
    for (; round < numRounds; ++round) {
        
        // save checkpoints
        if (round%1==0 && round!=0) {
            char filename1[100], filename2[100];
            sprintf(filename1, "checkpoints/%d_image.ppm", round);
            sprintf(filename2, "checkpoints/%d_hitpoints.txt", round);
            save(filename1, filename2);
        }
        
        printf("Round %d/%d:\n", round + 1, numRounds);
        // PASS1: ray tracing pass
        printf("Start ray tracing......\n");
        for (int u = 0; u < w; ++u) {
            for (int v = 0; v < h; ++v) {
                Vec3d p(double(u - cx) / w / fx, double(v - cy) / h / fy, 0);
                Ray ray;
                ray.s = s;
                ray.d = p - s;
                double t = focusPlane.intersectPlane(ray);
                Vec3d focusP = ray.s + ray.d * t;
                double theta = randomdist(0, 2 * M_PI);
                ray.s = ray.s + Vec3d(cos(theta), sin(theta), 0) * aperture;
                ray.d = focusP - ray.s;
                ray.d.normalize();
                (*hitpoints)[u * h + v]->valid = false;
                (*hitpoints)[u * h + v]->d = ray.d * -1;
                scene->trace(ray, Vec3d(1, 1, 1), 1, (long long)round * (photon + w * h) + u * h + v, (*hitpoints)[u * h + v]);
            }
        }
        
        scene->initializeHitpointKDTree(hitpoints);
        
        // photon tracing pass
        printf("Start photon mapping......\n");
        for (int i = 0; i < photon; ++i) {
            Ray ray = scene->generateRay((long long)round * photon + (round + 1) * w * h + i);
            scene->trace(ray, weight_init * light, 1, (long long)round * photon + i);
        }
    }
}

void SPPM::save(char *filename1, char *filename2) {
    for (int u = 0; u < w; ++u)
        for (int v = 0; v < h; ++v) {
            HitPoint *hp = (*hitpoints)[u * h + v];
            canvas[u][v] = hp->flux / (M_PI * hp->r2 * this->round * this->photon) + light * hp->fluxLight / this->round;
            canvas[u][v].x = sqrt(canvas[u][v].x);
            canvas[u][v].y = sqrt(canvas[u][v].y);
            canvas[u][v].z = sqrt(canvas[u][v].z);
            canvas[u][v] = min(canvas[u][v], Vec3d(1, 1, 1));
            canvas[u][v] = max(canvas[u][v], Vec3d(0, 0, 0));
        }
    
    FILE *file1 = fopen(filename1, "w");
    fprintf(file1, "P3\n%d %d\n255\n", w, h);
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++)
            fprintf(file1, "%d %d %d ", int(canvas[j][h - i - 1].x * 255 + 0.5), int(canvas[j][h - i - 1].y * 255 + 0.5), int(canvas[j][h - i - 1].z * 255 + 0.5));
        fprintf(file1, "\n");
    }
    fclose(file1);
    
    //p, weight, flux, fluxLight, d, norm, n, brdf, r2
    FILE *file2 = fopen(filename2, "w");
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++){
            fprintf(file2, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf\n",
                    (*hitpoints)[i*w+j]->p.x, (*hitpoints)[i*w+j]->p.y, (*hitpoints)[i*w+j]->p.z,
                    (*hitpoints)[i*w+j]->weight.x, (*hitpoints)[i*w+j]->weight.y, (*hitpoints)[i*w+j]->weight.z,
                    (*hitpoints)[i*w+j]->flux.x, (*hitpoints)[i*w+j]->flux.y, (*hitpoints)[i*w+j]->flux.z,
                    (*hitpoints)[i*w+j]->fluxLight.x, (*hitpoints)[i*w+j]->fluxLight.y, (*hitpoints)[i*w+j]->fluxLight.z,
                    (*hitpoints)[i*w+j]->norm.x, (*hitpoints)[i*w+j]->norm.y, (*hitpoints)[i*w+j]->norm.z,
                    (*hitpoints)[i*w+j]->n, (*hitpoints)[i*w+j]->r2);
        }
    }
    fclose(file2);
    fprintf(stderr, "checkpoint %s and %s saved\n", filename1, filename2);
}

void SPPM::load(char *filename, int round) {
    this->round = round;

    FILE *file = fopen(filename, "w");
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++){
            fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf\n",
                   &(*hitpoints)[i*w+j]->p.x, &(*hitpoints)[i*w+j]->p.y, &(*hitpoints)[i*w+j]->p.z,
                   &(*hitpoints)[i*w+j]->weight.x, &(*hitpoints)[i*w+j]->weight.y, &(*hitpoints)[i*w+j]->weight.z,
                   &(*hitpoints)[i*w+j]->flux.x, &(*hitpoints)[i*w+j]->flux.y, &(*hitpoints)[i*w+j]->flux.z,
                   &(*hitpoints)[i*w+j]->fluxLight.x, &(*hitpoints)[i*w+j]->fluxLight.y, &(*hitpoints)[i*w+j]->fluxLight.z,
                   &(*hitpoints)[i*w+j]->norm.x, &(*hitpoints)[i*w+j]->norm.y, &(*hitpoints)[i*w+j]->norm.z,
                   &(*hitpoints)[i*w+j]->n, &(*hitpoints)[i*w+j]->r2);
        }
    }
    fclose(file);
}
