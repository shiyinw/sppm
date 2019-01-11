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

#include "Utils.hpp"

bool ObjectKDTreeNode::inside(Mesh *mesh) {
    Vec3d faceMin = mesh->min();
    Vec3d faceMax = mesh->max();
    return (faceMin.x < max.x || (faceMin.x == max.x && faceMin.x == faceMax.x))
    && (faceMax.x > min.x || (faceMax.x == min.x && faceMin.x == faceMax.x))
    && (faceMin.y < max.y || (faceMin.y == max.y && faceMin.y == faceMax.y))
    && (faceMax.y > min.y || (faceMax.y == min.y && faceMin.y == faceMax.y))
    && (faceMin.z < max.z || (faceMin.z == max.z && faceMin.z == faceMax.z))
    && (faceMax.z > min.z || (faceMax.z == min.z && faceMin.z == faceMax.z));
}

ObjectKDTreeNode* ObjectKDTree::build(int depth, int d, vector<Mesh*>* meshes, Vec3d min, Vec3d max) {
    ObjectKDTreeNode *p = new ObjectKDTreeNode;
    p->min = min;
    p->max = max;
    Vec3d maxL, minR;
    if (d == 0) {
        maxL = Vec3d((p->min.x + p->max.x) / 2, p->max.y, p->max.z);
        minR = Vec3d((p->min.x + p->max.x) / 2, p->min.y, p->min.z);
    }
    else if (d == 1) {
        maxL = Vec3d(p->max.x, (p->min.y + p->max.y) / 2, p->max.z);
        minR = Vec3d(p->min.x, (p->min.y + p->max.y) / 2, p->min.z);
    }
    else {
        maxL = Vec3d(p->max.x, p->max.y, (p->min.z + p->max.z) / 2);
        minR = Vec3d(p->min.x, p->min.y, (p->min.z + p->max.z) / 2);
    }
    p->meshes = new vector<Mesh*>;
    for (auto face : *meshes)
        if (p->inside(face))
            p->meshes->push_back(face);
    
    const int max_faces = 8;
    const int max_depth = 24;
    
    if (p->meshes->size() > max_faces && depth < max_depth) {
        p->ls = build(depth + 1, (d + 1) % 3, p->meshes, min, maxL);
        p->rs = build(depth + 1, (d + 1) % 3, p->meshes, minR, max);
        
        vector<Mesh*> *meshL = p->ls->meshes, *meshR = p->rs->meshes;
        map<Mesh*, int> cnt;
        for (auto mesh : *meshL) cnt[mesh]++;
        for (auto mesh : *meshR) cnt[mesh]++;
        p->ls->meshes = new vector<Mesh*>;
        p->rs->meshes = new vector<Mesh*>;
        p->meshes->clear();
        for (auto mesh : *meshL)
            if (cnt[mesh] == 1)
                p->ls->meshes->push_back(mesh);
            else
                p->meshes->push_back(mesh);
        for (auto mesh : *meshR)
            if (cnt[mesh] == 1)
                p->rs->meshes->push_back(mesh);
    }
    else
        p->ls = p->rs = nullptr;
    return p;
}

void ObjectKDTree::getFaces(ObjectKDTreeNode *p, vector<Mesh*>* meshes) {
    p->l = (int)meshes->size();
    for (auto mesh : *(p->meshes))
        meshes->push_back(mesh);
    p->r = (int)meshes->size();
    if (p->ls) getFaces(p->ls, meshes);
    if (p->rs) getFaces(p->rs, meshes);
}

ObjectKDTree::ObjectKDTree(vector<Mesh*>* meshes) {
    Vec3d min = Vec3d(1e100, 1e100, 1e100);
    Vec3d max = min * -1;
    for (auto mesh : *meshes) {
        min = ::min(min, mesh->min());
        max = ::max(max, mesh->max());
    }
    root = build(1, 0, meshes, min, max);
    this->meshes = new vector<Mesh*>;
    getFaces(root, this->meshes);
}

double ObjectKDTree::getCuboidIntersection(ObjectKDTreeNode *p, Ray ray) {
    if (!(ray.s >= p->min && ray.s <= p->max)) { // outside
        double t = -1e100;
        if (fabs(ray.d.x) > 0)
            t = max(t, min((p->min.x - ray.s.x) / ray.d.x, (p->max.x - ray.s.x) / ray.d.x));
        if (fabs(ray.d.y) > 0)
            t = max(t, min((p->min.y - ray.s.y) / ray.d.y, (p->max.y - ray.s.y) / ray.d.y));
        if (fabs(ray.d.z) > 0)
            t = max(t, min((p->min.z - ray.s.z) / ray.d.z, (p->max.z - ray.s.z) / ray.d.z));
        if (t < -EPSILON) return 1e100;
        Vec3d pp = ray.s + ray.d * t;
        if (!(pp >= p->min && pp <= p->max)) return 1e100;
        return t;
    }
    else return -1e100;
}

void ObjectKDTree::getIntersection(ObjectKDTreeNode *p, Ray ray, Mesh* &nextMesh, double &tMin, Vec3d &norm) {
    for (int i = 0; i < p->meshes->size(); ++i) {
        pair<double, Vec3d> r = (*p->meshes)[i]->intersect(ray);
        double t = r.first;
        if (t > 0 && t < tMin) {
            tMin = t;
            nextMesh = (*p->meshes)[i];
            norm = r.second;
        }
    }
    
    double tl = p->ls ? getCuboidIntersection(p->ls, ray) : 1e100;
    double tr = p->rs ? getCuboidIntersection(p->rs, ray) : 1e100;
    if (tl < tr) {
        if (tMin <= tl) return;
        if (p->ls) getIntersection(p->ls, ray, nextMesh, tMin, norm);
        if (tMin <= tr) return;
        if (p->rs) getIntersection(p->rs, ray, nextMesh, tMin, norm);
    }
    else {
        if (tMin <= tr) return;
        if (p->rs) getIntersection(p->rs, ray, nextMesh, tMin, norm);
        if (tMin <= tl) return;
        if (p->ls) getIntersection(p->ls, ray, nextMesh, tMin, norm);
    }
}

void Scene::addObject(Object* object) {
    for (int i = 0; i < object->numFaces; ++i)
        object->meshes[i]->object = object;
    objects.push_back(object);
}

Ray Scene::generateRay(long long i) {
    double alpha = Utils::random(0, 2 * M_PI, 0, i);
    Vec3d s = sourceP +  Vec3d(cos(alpha), 0, sin(alpha)) * sourceR;
    Vec3d d = sampleReflectedRay(sourceN, 0, i);
    return Ray(s + d * EPSILON, d);
}

Vec3d Scene::sampleReflectedRay(Vec3d norm, int depth, long long i, double s) {
    Vec3d u = cross(Vec3d(1, 0, 0), norm);
    if (u.norm2() < EPSILON) u = cross(Vec3d(0, 1, 0), norm);
    u.normalize();
    Vec3d v = cross(norm, u);
    v.normalize();
    double theta = Utils::random(0, 2 * M_PI, 2 * depth + 1, i);
    double phi = asin(pow(Utils::random(0, 1, 2 * depth + 2, i), 1. / (s + 1)));
    return (norm * cos(phi) + (u * cos(theta) + v * sin(theta)) * sin(phi)).normalize();
}

void Scene::trace(const Ray &ray, const Vec3d &weight, int depth, long long i, HitPoint *hp) {
    if (depth > MAX_DEPTH)
        return;
    double tMin = 1e100;
    Mesh* nextMesh = nullptr;
    Vec3d norm;
    
    objectKDTree->getIntersection(objectKDTree->root, ray, nextMesh, tMin, norm);
    
    if (!nextMesh || (!hp && tMin < 1e-3)) return;
    Vec3d p = ray.s + ray.d * tMin;
    
    // russian roulette
    double s = BRDFs[nextMesh->brdf].specular + BRDFs[nextMesh->brdf].diffuse + BRDFs[nextMesh->brdf].refraction;
    double action = Utils::random(0, 1) * s;
    
    Vec3d dr = ray.d - norm * (2 * dot(ray.d, norm));
    
    // specular
    if (BRDFs[nextMesh->brdf].specular > 0 && action <= BRDFs[nextMesh->brdf].specular) {
        trace(
              Ray(p + dr * EPSILON, dr),
              weight * nextMesh->texture->query(p) * s, depth + 1, i, hp
              );
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
            double a = Utils::random();
            // phong specular
            if (a <= BRDFs[nextMesh->brdf].rho_s) {
                Vec3d d = sampleReflectedRay(dr, depth, i, BRDFs[nextMesh->brdf].phong_s);
                trace(
                      Ray(p + d * EPSILON, d),
                      weight * nextMesh->texture->query(p) * s,
                      depth + 1, i, hp
                      );
            }
            else {
                a -= BRDFs[nextMesh->brdf].rho_s;
                hitpointsKDTree->update(hitpointsKDTree->root, p, weight, ray.d);
                Vec3d d = sampleReflectedRay(norm, depth, i);
                if (dot(d, norm) < 0) d = d * -1;
                if (a <= BRDFs[nextMesh->brdf].rho_d) {
                    trace(
                          Ray(p + d * EPSILON, d),
                          weight * nextMesh->texture->query(p) * s,
                          depth + 1, i, hp
                          );
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
            
            if (Utils::random() <= R)
                trace(Ray(p + dr * EPSILON, dr), weight * nextMesh->texture->query(p) * s, depth + 1, i, hp);
            else {
                Vec3d d = ray.d / refractiveIndex + norm * (cosThetaIn / refractiveIndex - cosThetaOut);
                trace(Ray(p + d * EPSILON, d), weight * nextMesh->texture->query(p) * s, depth + 1, i, hp);
            }
        }
        else { // total internal reflection
            trace(Ray(p + dr * EPSILON, dr), weight * nextMesh->texture->query(p) * s, depth + 1, i, hp);
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
