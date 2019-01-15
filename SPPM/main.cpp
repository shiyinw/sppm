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

#define BOX_X 2 // half the size of the box 1.2
#define BOX_Y 1  // half the size of the box
#define BOX_Z_F -2.5
#define BOX_Z_B 1.5

Object* genBox(){
    Texture *desert = new Texture((char*)"objects/desert.ppm");
    Object *wall = new Sphere(Vec3d(0, 0, 0), 1, new TextureMapper(desert, Vec3d(0, -1, 0), Vec3d(1/1.2, 0, 0), 0.5, 0.5), GLASS);
    return wall;
}

Object* genWalls() {
    Texture *textureBottom = new Texture((char*)"objects/bottom.ppm");
    Texture *textureBack = new Texture((char*)"objects/desert.ppm");
    Texture *textureFront = new Texture((char*)"objects/desert_2.ppm");
    Texture *textureTop = new Texture((char*)"objects/top.ppm");
    Texture *textureleft = new Texture((char*)"objects/desert_3.ppm");
    Texture *textureright = new Texture((char*)"objects/desert_1.ppm");
    
    TextureMapper *gold = new TextureMapper(Vec3d(1.0, 0.8, 0.0));
    TextureMapper *white = new TextureMapper(Vec3d(1, 1, 1));
    
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
        new TriMesh(vertexes[1], vertexes[2], vertexes[5], new TextureMapper(textureBack, Vec3d(0, -1/BOX_Y, 0), Vec3d(1/BOX_X, 0, 0), 0.5, 0.5), MARBLE),
        new TriMesh(vertexes[2], vertexes[5], vertexes[6], new TextureMapper(textureBack, Vec3d(0, -1/BOX_Y, 0), Vec3d(1/BOX_X, 0, 0), 0.5, 0.5), MARBLE),
        // top
        new TriMesh(vertexes[0], vertexes[1], vertexes[2], new TextureMapper(textureTop, Vec3d(0, 0, 1./2), Vec3d(1./BOX_X, 0, 0), 2.5/4, 0.5), DIFFUSE),
        new TriMesh(vertexes[0], vertexes[2], vertexes[3], new TextureMapper(textureTop, Vec3d(0, 0, 1./2), Vec3d(1./BOX_X, 0, 0), 2.5/4, 0.5), DIFFUSE),
        // bottom
        new TriMesh(vertexes[4], vertexes[5], vertexes[6], new TextureMapper(textureBottom, Vec3d(0, 0, 1./2), Vec3d(1./BOX_X, 0, 0), -BOX_Z_F/4, 0.5), FLOOR),
        new TriMesh(vertexes[4], vertexes[6], vertexes[7], new TextureMapper(textureBottom, Vec3d(0, 0, 1./2), Vec3d(1./BOX_X, 0, 0), -BOX_Z_F/4, 0.5), FLOOR),
        // left
        new TriMesh(vertexes[0], vertexes[1], vertexes[4], gold, WALL),
        new TriMesh(vertexes[1], vertexes[4], vertexes[5], gold, WALL),
        // right
        new TriMesh(vertexes[2], vertexes[3], vertexes[6], gold, WALL),
        new TriMesh(vertexes[3], vertexes[6], vertexes[7], gold, WALL),
        // front
        new TriMesh(vertexes[0], vertexes[4], vertexes[7], new TextureMapper(textureFront, Vec3d(0, -1/BOX_Y, 0), Vec3d(1/BOX_X, 0, 0), 0.5, 0.5), WALL),
        new TriMesh(vertexes[0], vertexes[3], vertexes[7], new TextureMapper(textureFront, Vec3d(0, -1/BOX_Y, 0), Vec3d(1/BOX_X, 0, 0), 0.5, 0.5), WALL),
    };
    return walls;
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
    bunny->importPly((char*)"objects/bunny.ply",  new TextureMapper(Vec3d(1, 1, 1)), STANFORD_MODEL);
    bunny->locate(Vec3d(-0.3, -0.7, 0.3), Vec3d(0.3, 0.3, 0.6));
//    bunny->rotXZ(15 * M_PI / 180);
    bunny->aabb->print();
    
    Object *relic = new Object;
    relic->importObj((char*)"objects/relic.obj",  new TextureMapper(Vec3d(0.8, 0.8, 0.8)), STANFORD_MODEL);
    relic->locate(Vec3d(0.3, -0.7, 0.3), Vec3d(0.9, 0.3, 0.6));
//    relic->rotXZ(15 * M_PI / 180);
    relic->aabb->print();
    
    Scene *scene = new Scene(Vec3d(0.0, 0.5 - 1e-5, 0.1), 0.2, Vec3d(0, -1, 0));
    scene->addObject(relic);
    scene->addObject(bunny);
    scene->addObject(genWalls());
    scene->addObject(genBox());
    scene->addObject(genLight(Vec3d(0, 0.5 - EPSILON, 0.1), 0.05));
    scene->addObject(new Sphere(Vec3d(-0.32, -0.30, 0.3), 0.08, new TextureMapper(Vec3d(1, 1, 1)), GLASS));
    scene->addObject(new Sphere(Vec3d(0.42, 0.20, 0), 0.05, new TextureMapper(Vec3d(1, 1, 1)), MIRROR));
    return scene;
}


int main(int argc, char *argv[]) {
    Scene *scene = sceneBox();
    
    SPPM *camera = new SPPM(1024, 768, scene, Vec3d(0, 0.15, -1), 200000);
    camera->setLens(0.684, 0.811, 1e-3, 1 + 0.09);
    
    //camera->load((char*)"checkpoints/11_hitpoints.txt", 11);
    
    camera -> render(2500);
    
    // 2500, 200000
    
    camera->save((char*)"result_test.ppm", (char*)"hitpoints.txt");
    return 0;
}

