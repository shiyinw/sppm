//
//  Scene.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#define EPSILON 1e-6
#define MAX_DEPTH 10

#include "Scene.hpp"
#include "Random.hpp"
#include <iostream>

struct Params{
    
};

void Scene::addObject(Object* object) {
    for (int i = 0; i < object->numFaces; ++i)
        object->meshes[i]->object = object;
    objects.push_back(object);
}

Ray Scene::generateRay(long long i) {
    double alpha = randomdist(0, 2 * M_PI, 0, i);
    Vec3d s = sourceP +  Vec3d(cos(alpha), 0, sin(alpha)) * sourceR;
    Vec3d d = sampleReflectedRay(sourceN, 0, i);
    Ray ray;
    ray.d = d;
    ray.s = s + ray.d * EPSILON;
    return ray;
}

Vec3d Scene::sampleReflectedRay(Vec3d norm, int depth, long long i, double s) {
    Vec3d u = cross(Vec3d(1, 0, 0), norm);
    if (u.norm2() < EPSILON) u = cross(Vec3d(0, 1, 0), norm);
    u.normalize();
    Vec3d v = cross(norm, u);
    v.normalize();
    double theta = randomdist(0, 2 * M_PI, 2 * depth + 1, i);
    double phi = asin(pow(randomdist(0, 1, 2 * depth + 2, i), 1. / (s + 1)));
    return (norm * cos(phi) + (u * cos(theta) + v * sin(theta)) * sin(phi)).normalize();
}

// scene->trace(ray, Vec3d(1, 1, 1), 1, (long long)round * (numPhotons + w * h) + u * h + v, (*hitpoints)[u * h + v]);
void Scene::trace(const Ray &ray, const Vec3d &weight, int depth, long long i, HitPoint *hp) {
    if (depth > MAX_DEPTH)
        return;
    double tMin = 1e100;
    Mesh* nextMesh = nullptr;
    Vec3d norm;
    
    //search next mesh through KD Tree
    objectKDTree->getIntersection(objectKDTree->root, ray, nextMesh, tMin, norm);
    
    if (!nextMesh || (!hp && tMin < 1e-3)) return;
    Vec3d p = ray.s + ray.d * tMin;
    
    // russian roulette
    double s = BRDFs[nextMesh->brdf].specular + BRDFs[nextMesh->brdf].diffuse + BRDFs[nextMesh->brdf].refraction;
    double action = randomdist(0, 1) * s;
    Vec3d dr = ray.d - norm * (2 * dot(ray.d, norm));
    
    // specular
    if (BRDFs[nextMesh->brdf].specular > 0 && action <= BRDFs[nextMesh->brdf].specular) {
        Ray ray;
        ray.d = dr;
        ray.s = p + ray.d * EPSILON;
        trace(ray, weight * nextMesh->texture->query(p) * s, depth + 1, i, hp);
        return;
    }
    action -= BRDFs[nextMesh->brdf].specular;
    
    // diffuse
    if (BRDFs[nextMesh->brdf].diffuse > 0 && action <= BRDFs[nextMesh->brdf].diffuse) {
        if (hp) {
            hp->p = p;
            hp->weight = weight * nextMesh->texture->query(p) * s;
            hp->fluxLight = hp->fluxLight + hp->weight * (nextMesh->brdf == LIGHT);
            hp->brdf = BRDFs[nextMesh->brdf];
            hp->norm = norm;
            if (nextMesh->brdf == LIGHT) {
                hp->fluxLight = hp->fluxLight + hp->weight;
                hp->valid = false;
            }
            else
                hp->valid = true;
        }
        else {
            double a = randomdist();
            // phong specular
            if (((double) rand() / (RAND_MAX)) <= BRDFs[nextMesh->brdf].rho_s) {
                Ray rayout;
                rayout.d = sampleReflectedRay(dr, depth, i, BRDFs[nextMesh->brdf].phong_s);
                rayout.s = p + rayout.d * EPSILON;
                
                trace(rayout, weight * nextMesh->texture->query(p) * s, depth + 1, i, hp);
            }
            else {
                a -= BRDFs[nextMesh->brdf].rho_s;
                hitpointsKDTree->update(hitpointsKDTree->root, p, weight, ray.d);
                Ray rayout;
                rayout.d = sampleReflectedRay(norm, depth, i);
                if (dot(ray.d, norm) < 0)
                    rayout.d = rayout.d * -1;
                if (a <= BRDFs[nextMesh->brdf].rho_d) {
                    rayout.s = p + rayout.d * EPSILON;
                    trace(ray, weight * nextMesh->texture->query(p) * s, depth + 1, i, hp);
                }
            }
        }
        return;
    }
    action -= BRDFs[nextMesh->brdf].diffuse;
    
    // refraction
    if (BRDFs[nextMesh->brdf].refraction > 0 && action <= BRDFs[nextMesh->brdf].refraction) {
        if (!nextMesh->object->center)
            nextMesh->object->calcCenter();
        bool incoming = dot(*(nextMesh->object)->center - p, norm) < 0;
        double refractiveIndex = BRDFs[nextMesh->brdf].refractiveIndex;
        if (!incoming) refractiveIndex = 1. / refractiveIndex;
        double cosThetaIn = -dot(ray.d, norm);
        double cosThetaOut2 = 1 - (1 - cosThetaIn*cosThetaIn) / refractiveIndex/refractiveIndex;
        
        if (cosThetaOut2 >= -EPSILON) {
            double cosThetaOut = sqrt(cosThetaOut2);
            
            // schlick's approximation
            double R0 = ((1 - refractiveIndex) / (1 + refractiveIndex))*((1 - refractiveIndex) / (1 + refractiveIndex));
            double cosTheta = incoming ? cosThetaIn : cosThetaOut;
            double R = R0 + (1 - R0) * pow(1 - cosTheta, 5);
            
            
            if (((double) rand() / (RAND_MAX)) <= R){
                Ray rayout;
                rayout.d = dr;
                rayout.s = p + rayout.d * EPSILON;
                trace(ray, weight * nextMesh->texture->query(p) * s, depth + 1, i, hp);
            }
            else {
                Ray rayout;
                rayout.d = ray.d / refractiveIndex + norm * (cosThetaIn / refractiveIndex - cosThetaOut);
                rayout.s = p + rayout.d * EPSILON;
                trace(rayout, weight * nextMesh->texture->query(p) * s, depth + 1, i, hp);
            }
        }
        else { // total internal reflection
            Ray rayout;
            rayout.d = dr;
            rayout.s = p + rayout.d * EPSILON;
            trace(rayout, weight * nextMesh->texture->query(p) * s, depth + 1, i, hp);
        }
    }
}

void Scene::initializeHitpointKDTree(vector<HitPoint*>* hitpoints) {
    if (hitpointsKDTree)
        delete hitpointsKDTree;
    hitpointsKDTree = new HitPointKDTree(hitpoints);
    fprintf(stderr, "Hitpoint KD tree built\n");
}

void Scene::initializeObjectKDTree() {
    vector<Mesh*> *meshes = new vector<Mesh*>;
    for (auto object : objects) {
        for (int i = 0; i < object->numFaces; ++i)
            meshes->push_back(object->meshes[i]);
    }
    objectKDTree = new ObjectKDTree(meshes);
    
    fprintf(stderr, "Object KD tree built\n");
}
