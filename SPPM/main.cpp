//
//  main.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 11/25/18.
//  Copyright © 2018 Sherilyn Wankins. All rights reserved.
//

#define EPSILON 1e-6

#include <cstdio>
#include <iostream>
#include "SPPM.hpp"
#include "Scene.hpp"
using namespace std;

#define BOX_X 1.2 // half the size of the box
#define BOX_Y 1  // half the size of the box
#define BOX_Z_F -2.5
#define BOX_Z_B 1.5

Object* genBox(){
    Texture *desert = new Texture((char*)"objects/desert.ppm");
    Object *wall = new Sphere(Vec3d(0, 0, 0), 2, new TextureMapper(desert, Vec3d(0, -1, 0), Vec3d(1/1.2, 0, 0), 0.5, 0.5), GLASS);
    return wall;
}

Object* genWalls() {
    Texture *textureBottom = new Texture((char*)"objects/bottom.ppm");
    Texture *textureBack = new Texture((char*)"objects/desert.ppm");
    Texture *textureTop = new Texture((char*)"objects/top.ppm");
    
    TextureMapper *color_front = new TextureMapper(Vec3d(0.8, 0.2, 0.2));
    TextureMapper *color_left = new TextureMapper(Vec3d(0.2, 0.5, 0.8));
    TextureMapper *color_right = new TextureMapper(Vec3d(0.2, 0.8, 0.5));
    
    Object *walls = new Object;
    walls->numVertexes = 8;
    walls->numFaces = 12;
    walls->vertexes = new Vec3d*[walls->numVertexes]{
        new Vec3d(-BOX_X, BOX_Y, BOX_Z_F),
        new Vec3d(-BOX_X, BOX_Y, BOX_Z_B),
        new Vec3d(BOX_X, BOX_Y, BOX_Z_B),
        new Vec3d(BOX_X, BOX_Y, BOX_Z_F),
        new Vec3d(-BOX_X, -BOX_Y, BOX_Z_F),
        new Vec3d(-BOX_X, -BOX_Y, BOX_Z_B),
        new Vec3d(BOX_X, -BOX_Y, BOX_Z_B),
        new Vec3d(BOX_X, -BOX_Y, BOX_Z_F)
    };
    for (int i = 0; i < walls->numVertexes; ++i)
        *walls->vertexes[i] = *walls->vertexes[i] * 0.5;
    Vec3d** vertexes = walls->vertexes;
    walls->meshes = new Mesh*[walls->numFaces]{
        // back
        new TriMesh(vertexes[1], vertexes[2], vertexes[5], new TextureMapper(textureBack, Vec3d(0, -1, 0), Vec3d(1/1.2, 0, 0), 0.5, 0.5), MARBLE),
        new TriMesh(vertexes[2], vertexes[5], vertexes[6], new TextureMapper(textureBack, Vec3d(0, -1, 0), Vec3d(1/1.2, 0, 0), 0.5, 0.5), MARBLE),
        // top
        new TriMesh(vertexes[0], vertexes[1], vertexes[2], new TextureMapper(textureTop, Vec3d(0, 0, 1./2), Vec3d(1./1.2, 0, 0), 2.5/4, 0.5), DIFFUSE),
        new TriMesh(vertexes[0], vertexes[2], vertexes[3], new TextureMapper(textureTop, Vec3d(0, 0, 1./2), Vec3d(1./1.2, 0, 0), 2.5/4, 0.5), DIFFUSE),
        // bottom
        new TriMesh(vertexes[4], vertexes[5], vertexes[6], new TextureMapper(textureBottom, Vec3d(0, 0, 1./2), Vec3d(1./1.2, 0, 0), 2.5/4, 0.5), FLOOR),
        new TriMesh(vertexes[4], vertexes[6], vertexes[7], new TextureMapper(textureBottom, Vec3d(0, 0, 1./2), Vec3d(1./1.2, 0, 0), 2.5/4, 0.5), FLOOR),
        // left
        new TriMesh(vertexes[0], vertexes[1], vertexes[4], color_left, WALL),
        new TriMesh(vertexes[1], vertexes[4], vertexes[5], color_left, WALL),
        // right
        new TriMesh(vertexes[2], vertexes[3], vertexes[6], color_right, WALL),
        new TriMesh(vertexes[3], vertexes[6], vertexes[7], color_right, WALL),
        // front
        new TriMesh(vertexes[0], vertexes[4], vertexes[7], color_front, WALL),
        new TriMesh(vertexes[0], vertexes[3], vertexes[7], color_front, WALL),
    };
    return walls;
}

Object* genDesk() {
    
    TextureMapper *color = new TextureMapper(Vec3d(1, 1, 1));
    
    Object *desk = new Object;
    desk->numVertexes = 4;
    desk->numFaces = 2;
    const double theta = 45 / 180. * M_PI;
    desk->vertexes = new Vec3d*[desk->numVertexes]{
        new Vec3d(-10, 0 - 10 * sin(theta), -10 * cos(theta)),
        new Vec3d(-10, 0 + 10 * sin(theta), 10 * cos(theta)),
        new Vec3d(10, 0 + 10 * sin(theta), 10 *  cos(theta)),
        new Vec3d(10, 0 - 10 * sin(theta), -10 * cos(theta))
    };
    Vec3d** vertexes = desk->vertexes;
    desk->meshes = new Mesh*[desk->numFaces]{
        // bottom
        new TriMesh(vertexes[0], vertexes[1], vertexes[2], color, DESK),
        new TriMesh(vertexes[0], vertexes[2], vertexes[3], color, DESK),
    };
    return desk;
}

Object* genLight(Vec3d p, double r) {
    Object *light = new Object;
    light->numFaces = 1;
    light->meshes = new Mesh*[1]{
        new SphereMesh(p, r, new TextureMapper(Vec3d(1, 1, 1)), LIGHT)
    };
    return light;
}

Scene *sceneBox() {
    Object *bunny = new Object;
    bunny->importPly((char*)"objects/bunny.ply",  new TextureMapper(Vec3d(0.4, 0.8, 0.8)), STANFORD_MODEL);

    bunny->scale(2.8, 0, 0, 0.24,
                 0, 2.8, 0, -0.09 - 0.5,
                 0, 0, -2.8, 0.12
                 );
    bunny->rotXZ(15 * M_PI / 180);
    bunny->center->print();
    
    Object *water = new Object;
    water->importPly((char*)"objects/water.ply",  new TextureMapper(Vec3d(1, 1, 1)), WATER);
    water->scale(
                 1. / 5.52799 * 1.2, 0, 0, -0.6,
                 0, 0.12 / (1.85354 - 1.34492), 0, -0.31731 + 0.08,
                 0, 0, 1.2 / (5.59200 + 0.00456), -1.19902 + 0.75
                 );
    water->center = new Vec3d(0, -1, 0);
    water->printBox();
    
    Scene *scene = new Scene(Vec3d(0.0, 0.5 - 1e-5, 0.1), 0.2, Vec3d(0, -1, 0));
    scene->addObject(bunny);
    scene->addObject(water);
    scene->addObject(genWalls());
    scene->addObject(genBox());
    scene->addObject(genLight(Vec3d(0, 0.5 - EPSILON, 0.1), 0.05));
    scene->addObject(new Sphere(Vec3d(-0.32, -0.30, 0.3), 0.18, new TextureMapper(Vec3d(1, 1, 1)), GLASS));
    scene->addObject(new Sphere(Vec3d(0.42, 0.20, 0), 0.15, new TextureMapper(Vec3d(1, 1, 1)), MIRROR));
    return scene;
}


int main(int argc, char *argv[]) {
    Scene *scene = sceneBox();
    //Scene *scene = sceneTeapot();
    
    SPPM *camera = new SPPM(1024, 768, scene, Vec3d(0, 0.15, -1), 200000);
    camera->setLens(0.684, 0.811, 1e-3, 1 + 0.09);
    
    
    //camera->load((char*)"checkpoints/11_image.ppm", (char*)"checkpoints/11_hitpoints.txt", 11);
    
    camera -> render(2500);
    
    // 2500, 200000
    
    camera->save((char*)"result_test.ppm", (char*)"hitpoints.txt");
    return 0;
}

