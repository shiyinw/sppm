//
//  Scene.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 11/1/18.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#define EPSILON 1e-6
#define MAX_DEPTH 10
#define SEARCH_RADIUS 1e-4

#include "Scene.hpp"
#include <iostream>

struct Params{
    
};

void Scene::addObject(Object* object) {
    for (int i = 0; i < object->nm; ++i)
        object->meshes[i]->object = object;
    objects.push_back(object);
}

Ray Scene::generateRay(long long i) {
    double alpha = double(rand())*2*M_PI/RAND_MAX;
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
    double theta = double(rand())*2.*M_PI/RAND_MAX;
    double phi = asin(double(rand())/RAND_MAX);
    return (norm * cos(phi) + (u * cos(theta) + v * sin(theta)) * sin(phi)).normalize();
}

void Scene::rayTrace(const Ray &ray, const Vec3d &color, int depth, HitPoint *hp) {
    if (depth > MAX_DEPTH)
        return;
    double tMin = 1e100;
    Mesh* nextMesh = nullptr;
    Vec3d norm;
    
    //search next mesh through KD Tree
    objectKDTree->getIntersection(objectKDTree->root, ray, nextMesh, tMin, norm);
    
    if (!nextMesh) return;
    Vec3d p = ray.s + ray.d * tMin;
    
    // russian roulette
    double s = BRDFs[nextMesh->brdf].specular + BRDFs[nextMesh->brdf].diffuse + BRDFs[nextMesh->brdf].refraction;
    double action = double(rand())/RAND_MAX * s;
    Vec3d dr = ray.d - norm * (2 * dot(ray.d, norm));
    
    // specular 反射
    if (BRDFs[nextMesh->brdf].specular > 0 && action <= BRDFs[nextMesh->brdf].specular) {
        Ray ray;
        ray.d = dr;
        ray.s = p + ray.d * EPSILON;
        rayTrace(ray, color * nextMesh->texture->query(p) * s, depth + 1, hp);
        return;
    }
    action -= BRDFs[nextMesh->brdf].specular;
    
    // diffuse
    if (BRDFs[nextMesh->brdf].diffuse > 0 && action <= BRDFs[nextMesh->brdf].diffuse) {
        // stop ray tracing
        hp->p = p;
        hp->color = color * nextMesh->texture->query(p) * s;
        hp->fluxLight = hp->fluxLight + hp->color * (nextMesh->brdf == LIGHT);
        hp->brdf = BRDFs[nextMesh->brdf];
        hp->norm = norm;
        if (nextMesh->brdf == LIGHT) {
            hp->fluxLight = hp->fluxLight + hp->color;
            hp->valid = false;
        }
        else
            hp->valid = true;
        return;
    }
    action -= BRDFs[nextMesh->brdf].diffuse;
    
    // refraction
    if (BRDFs[nextMesh->brdf].refraction > 0 && action <= BRDFs[nextMesh->brdf].refraction) {
        double refractiveIndex = dot(nextMesh->object->aabb->_center - p, norm) < 0? BRDFs[nextMesh->brdf].refractiveIndex: 1. / BRDFs[nextMesh->brdf].refractiveIndex;
        double cosThetaIn = -dot(ray.d, norm);
        double cosThetaOut2 = 1 - (1 - cosThetaIn*cosThetaIn) / refractiveIndex/refractiveIndex;
        
        if (cosThetaOut2 >= -EPSILON) {
            double cosThetaOut = sqrt(cosThetaOut2);
            double R0 = ((1 - refractiveIndex) / (1 + refractiveIndex))*((1 - refractiveIndex) / (1 + refractiveIndex));
            double cosTheta = dot(nextMesh->object->aabb->_center - p, norm) < 0 ? cosThetaIn : cosThetaOut;
            double R = R0 + (1 - R0) * pow(1 - cosTheta, 5);
            
            if (((double) rand() / (RAND_MAX)) <= R){
                Ray rayout;
                rayout.d = dr;
                rayout.s = p + rayout.d * EPSILON;
                rayTrace(ray, color * nextMesh->texture->query(p) * s, depth + 1, hp);
            }
            else {
                Ray rayout;
                rayout.d = ray.d / refractiveIndex + norm * (cosThetaIn / refractiveIndex - cosThetaOut);
                rayout.s = p + rayout.d * EPSILON;
                rayTrace(rayout, color * nextMesh->texture->query(p) * s, depth + 1, hp);
            }
        }
        else {
            Ray rayout;
            rayout.d = dr;
            rayout.s = p + rayout.d * EPSILON;
            rayTrace(rayout, color * nextMesh->texture->query(p) * s, depth + 1, hp);
        }
    }
}

void Scene::photonTrace(const Ray &ray, const Vec3d &color, int depth, long long power) {
    if (depth > MAX_DEPTH)
        return;
    double distance = 1e100;
    Mesh* nextMesh = nullptr;
    Vec3d norm;
    //search next mesh through KD Tree
    objectKDTree->getIntersection(objectKDTree->root, ray, nextMesh, distance, norm);
    
    if (!nextMesh || distance < SEARCH_RADIUS)
        return;
    Vec3d p = ray.s + ray.d * distance;
    
    // russian roulette
    double s = BRDFs[nextMesh->brdf].specular + BRDFs[nextMesh->brdf].diffuse + BRDFs[nextMesh->brdf].refraction;
    double action = double(rand())/RAND_MAX * s;
    Vec3d dr = ray.d - norm * (2 * dot(ray.d, norm));
    
    // specular reflection
    if (BRDFs[nextMesh->brdf].specular > 0 && action <= BRDFs[nextMesh->brdf].specular) {
        Ray ray;
        ray.d = dr;
        ray.s = p + ray.d * EPSILON;
        photonTrace(ray, color * nextMesh->texture->query(p) * s, depth + 1, power);
        return;
    }
    action -= BRDFs[nextMesh->brdf].specular;
    
    // diffuse reflection
    if (BRDFs[nextMesh->brdf].diffuse > 0 && action <= BRDFs[nextMesh->brdf].diffuse) {
        double a = double(rand())/RAND_MAX;
        // phong specular
        if (((double) rand() / (RAND_MAX)) <= BRDFs[nextMesh->brdf].rho_s) {
            Ray rayout;
            rayout.d = sampleReflectedRay(dr, depth, power, BRDFs[nextMesh->brdf].phong_s);
            rayout.s = p + rayout.d * EPSILON;
            photonTrace(rayout, color * nextMesh->texture->query(p) * s, depth + 1, power);
        }
        else {
            a -= BRDFs[nextMesh->brdf].rho_s;
            hitpointsKDTree->update(hitpointsKDTree->root, p, color, ray.d);
            Ray rayout;
            rayout.d = sampleReflectedRay(norm, depth, power);
            if (dot(ray.d, norm) < 0)
                rayout.d = rayout.d * -1;
            if (a <= BRDFs[nextMesh->brdf].rho_d) {
                rayout.s = p + rayout.d * EPSILON;
                photonTrace(ray, color * nextMesh->texture->query(p) * s, depth + 1, power);
            }
        }
        return;
    }
    action -= BRDFs[nextMesh->brdf].diffuse;
    
    // refraction
    if (BRDFs[nextMesh->brdf].refraction > 0 && action <= BRDFs[nextMesh->brdf].refraction) {
        bool incoming = dot(nextMesh->object->aabb->_center - p, norm) < 0;
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
                photonTrace(ray, color * nextMesh->texture->query(p) * s, depth + 1, power);
            }
            else {
                Ray rayout;
                rayout.d = ray.d / refractiveIndex + norm * (cosThetaIn / refractiveIndex - cosThetaOut);
                rayout.s = p + rayout.d * EPSILON;
                photonTrace(rayout, color * nextMesh->texture->query(p) * s, depth + 1, power);
            }
        }
        else { // total internal reflection
            Ray rayout;
            rayout.d = dr;
            rayout.s = p + rayout.d * EPSILON;
            photonTrace(rayout, color * nextMesh->texture->query(p) * s, depth + 1, power);
        }
    }
}
