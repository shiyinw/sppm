//
//  main.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 11/1/18.
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

Object* genLens(){
    Texture *desert = new Texture((char*)"objects/desert.ppm");
    Object *wall = new Sphere(Vec3d(0, 0, 0), 1, new TextureMapper(desert, Vec3d(0, -1/BOX_Y, 0), Vec3d(1/BOX_X, 0, 0), 0.5, 0.5), GLASS);
    return wall;
}

Object* genBox() {
    Texture *textureBottom = new Texture((char*)"objects/floor.ppm");
    Texture *textureBack = new Texture((char*)"objects/desert.ppm");
    Texture *textureFront = new Texture((char*)"objects/desert.ppm");
    Texture *textureTop = new Texture((char*)"objects/top.ppm");
    
    TextureMapper *gold = new TextureMapper(Vec3d(1.0, 0.8, 0.0));
    
    Object *walls = new Object;
    walls->nv = 8;
    walls->nm = 12;
    walls->vertexes = new Vec3d*[walls->nv]{
        new Vec3d(-BOX_X, BOX_Y, BOX_Z_F),
        new Vec3d(-BOX_X, BOX_Y, BOX_Z_B),
        new Vec3d(BOX_X, BOX_Y, BOX_Z_B),
        new Vec3d(BOX_X, BOX_Y, BOX_Z_F),
        new Vec3d(-BOX_X, -BOX_Y, BOX_Z_F),
        new Vec3d(-BOX_X, -BOX_Y, BOX_Z_B),
        new Vec3d(BOX_X, -BOX_Y, BOX_Z_B),
        new Vec3d(BOX_X, -BOX_Y, BOX_Z_F)
    };
    for (int i = 0; i < walls->nv; ++i)
        *walls->vertexes[i] = *walls->vertexes[i] * 0.5;
    Vec3d** vertexes = walls->vertexes;
    walls->meshes = new Mesh*[walls->nm]{
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

Scene *sceneBox() {
    Object *bunny = new Object;
    bunny->importObj((char*)"objects/simplify_bunny.obj",  new TextureMapper(Vec3d(1, 1, 1)), STANFORD_MODEL);
    bunny->locate(Vec3d(0, -0.7, 0.3), Vec3d(0.6, 0.3, 0.6));
    bunny->aabb->print();
    
    Object *relic = new Object;
    relic->importObj((char*)"objects/simplify_relic.obj",  new TextureMapper(Vec3d(0.8, 0.8, 0.8)), GLASS);
    relic->locate(Vec3d(-0.2, -0.7, -0.2), Vec3d(0.2, 0.0, 0.2));
    relic->aabb->print();
    
    
    Object *teardrop = new Object;
    teardrop->nv=0;
    teardrop->nm=1;
    teardrop->meshes = new Mesh*[2]{new WaterDropMesh(Vec3d(0, 0, 0), 1, 1, new TextureMapper(Vec3d(1,0.1 , 0.1)), DIFFUSE)};
    
    
    Scene *scene = new Scene(Vec3d(0.0, 0.5 - 1e-5, 0.1), 0.2, Vec3d(0, -1, 0));
    scene->addObject(relic);
    scene->addObject(bunny);
    scene->addObject(teardrop);
    scene->addObject(genBox());
    //scene->addObject(genLens()); // 镜头前放透镜的效果
    scene->addObject(new Sphere(Vec3d(0, 0.5, 0.1), 0.05, new TextureMapper(Vec3d(1, 1, 1)), LIGHT));
    scene->addObject(new Sphere(Vec3d(-0.2, -0.2, 0.4), 0.10, new TextureMapper(Vec3d(1, 0.8, 0)), DESK));
    scene->addObject(new Sphere(Vec3d(-0.32, -0.30, 0.3), 0.08, new TextureMapper(Vec3d(1, 1, 1)), GLASS));
    scene->addObject(new Sphere(Vec3d(0.32, 0.20, 0), 0.05, new TextureMapper(Vec3d(1, 1, 1)), MIRROR));
    scene->addObject(new Sphere(Vec3d(-0.32, 0.20, 0.3), 0.1, new TextureMapper(Vec3d(1,0.1 , 0.1)), DIFFUSE));
    return scene;
}

int main(int argc, char *argv[]) {
    Scene *scene = sceneBox();

    SPPM *camera = new SPPM(1024, 768, scene, Vec3d(0, 0.15, -1), 200000);

    //camera->load((char*)"checkpoints/50_hitpoints.txt", 50);

    camera -> render(5000);

    camera->save((char*)"result_test.ppm", (char*)"hitpoints.txt");
    return 0;
}
