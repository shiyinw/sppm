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
#include <ctime>
using namespace std;

#define BUFFER_SIZE 1024
#define APERTURE 1.5e-3 //1.5mm

struct Ray;

void SPPM::render(int numRounds) {
    vector<Mesh*> *meshes = new vector<Mesh*>;
    for (auto object : scene->objects) {
        for (int i = 0; i < object->numFaces; ++i)
            meshes->push_back(object->meshes[i]);
    }
    scene->objectKDTree = new ObjectKDTree(meshes);
    
    TriMesh focusPlane(new Vec3d(0, 0, s.z + focus), new Vec3d(0, 1, s.z + focus), new Vec3d(1, 0, s.z + focus));
    
    light = Vec3d(1, 1, 1);
    const Vec3d weight_init = Vec3d(2.5, 2.5, 2.5);
    
    // PASS1: ray tracing pass
    printf("PASS1: Ray Tracing, Hitpoint KD Tree......\n");
    for (int u = 0; u < w; ++u) {
        for (int v = 0; v < h; ++v) {
            Vec3d p((double(u)/double(w) - 0.5)/fx, (double(v)/double(h) - 0.5)/fy, 0);
            Ray ray;
            ray.s = s;
            ray.d = p - s;
            
            // depth of focus
            double t = focusPlane.intersectPlane(ray);
            Vec3d focusP = ray.s + ray.d * t;
            double theta = randomdist(0, 2 * M_PI);
            ray.s = ray.s + Vec3d(cos(theta), sin(theta), 0) * APERTURE;
            ray.d = focusP - ray.s;
            ray.d.normalize();
            
            // initialize photons
            (*hitpoints)[u * h + v]->valid = false;
            (*hitpoints)[u * h + v]->dir = ray.d * -1;
            scene->rayTrace(ray, Vec3d(1, 1, 1), 1, (*hitpoints)[u * h + v]);
        }
    }
    scene->hitpointsKDTree = new HitPointKDTree(hitpoints);
    
    // SPPM
    printf("PASS2: Photon tracing......\n");
    clock_t begin = clock();
    for (; round < numRounds; ++round) {
        // save checkpoints
        if (round%1==0 && round!=0) {
            printf("Start saving......\n");
            char filename1[100], filename2[100];
            sprintf(filename1, "checkpoints/%d_image.ppm", round);
            sprintf(filename2, "checkpoints/%d_hitpoints.txt", round);
            if(round%10==0)
                save(filename1, filename2);
            else
                save(filename1, NULL);
        }
        for (int i = 0; i < photon; ++i) {
            Ray ray = scene->generateRay((long long)round * photon + (round + 1) * w * h + i);
            scene->photonTrace(ray, weight_init * light, 1, (long long)round * photon + i);
        }
        printf("Round %d/%d, time %f\n", round + 1, numRounds, float(clock()-begin)/CLOCKS_PER_SEC);
    }
}

void SPPM::save(char *filename1, char *filename2=NULL) {
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
    if(filename2==NULL){
        fprintf(stderr, "image %s saved\n", filename1);
        return;
    }
    //p, weight, flux, fluxLight, d, norm, n, brdf, r2
    FILE *file2 = fopen(filename2, "w");
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++){
            fprintf(file2, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf\n",
                    (*hitpoints)[i*w+j]->p.x, (*hitpoints)[i*w+j]->p.y, (*hitpoints)[i*w+j]->p.z,
                    (*hitpoints)[i*w+j]->color.x, (*hitpoints)[i*w+j]->color.y, (*hitpoints)[i*w+j]->color.z,
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
                   &(*hitpoints)[i*w+j]->color.x, &(*hitpoints)[i*w+j]->color.y, &(*hitpoints)[i*w+j]->color.z,
                   &(*hitpoints)[i*w+j]->flux.x, &(*hitpoints)[i*w+j]->flux.y, &(*hitpoints)[i*w+j]->flux.z,
                   &(*hitpoints)[i*w+j]->fluxLight.x, &(*hitpoints)[i*w+j]->fluxLight.y, &(*hitpoints)[i*w+j]->fluxLight.z,
                   &(*hitpoints)[i*w+j]->norm.x, &(*hitpoints)[i*w+j]->norm.y, &(*hitpoints)[i*w+j]->norm.z,
                   &(*hitpoints)[i*w+j]->n, &(*hitpoints)[i*w+j]->r2);
        }
    }
    fclose(file);
}
