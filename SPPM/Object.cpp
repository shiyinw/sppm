//
//  Object.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#define BUFFER_SIZE 1024
#define EPSILON 1e-6

#include "Object.hpp"
#include "Utils.hpp"

#include <iostream>
using namespace std;
void Object::importPly(char *filename,  TextureMapper *texture, int brdf) {
    FILE *file = fopen(filename, "r");
    char buffer[BUFFER_SIZE];
    numVertexes = numFaces = 0;
    while (fgets(buffer, BUFFER_SIZE, file)) {
        if (string(buffer) == "end_header\n") break;
        if (string(buffer).substr(0, 14) == "element vertex"){
            string line = string(buffer);
            line = line.substr(15, line.size()-1);
            numVertexes = stoi(line);
        }
        else if (string(buffer).substr(0, 12) == "element face"){
            string line = string(buffer);
            line = line.substr(13, line.size()-1);
            numFaces = stoi(line);
        }
    }
    vertexes = new Vec3d*[numVertexes];
    for (int i = 0; i < numVertexes; ++i) {
        double x, y, z;
        fscanf(file, "%lf%lf%lf", &x, &y, &z);
        fgets(buffer, BUFFER_SIZE, file);
        vertexes[i] = new Vec3d(x, y, z);
    }
    meshes = new Mesh*[numFaces];
    for (int i = 0; i < numFaces; ++i) {
        int n, a, b, c;
        fscanf(file, "%d", &n);
        if (n != 3)
            throw runtime_error("Only trianglar faces are supported!");
        fscanf(file, "%d%d%d", &a, &b, &c);
        meshes[i] = new TriangularFace(vertexes[a], vertexes[b], vertexes[c], texture, brdf);
    }
    fclose(file);
    fprintf(stderr, "Imported object %s: %d vertexes and %d faces\n", filename, numVertexes, numFaces);
}


void Object::calcCenter() {
    center = new Vec3d(0, 0, 0);
    for (int i = 0; i < numVertexes; ++i)
        *center = *center + *vertexes[i];
    *center = *center / numVertexes;
}

void Object::printBox() {
    Vec3d min(1e100, 1e100, 1e100);
    Vec3d max = min * -1;
    for (int i = 0; i < numFaces; ++i) {
        min = ::min(min, meshes[i]->min());
        max = ::max(max, meshes[i]->max());
    }
    min.print();
    max.print();
}

void Object::scale(
                   double fxx, double fxy, double fxz, double fxb,
                   double fyx, double fyy, double fyz, double fyb,
                   double fzx, double fzy, double fzz, double fzb) {
    for (int i = 0; i < numVertexes; ++i) {
        Vec3d ver = *vertexes[i];
        vertexes[i]->x = fxx * ver.x + fxy * ver.y + fxz * ver.z + fxb;
        vertexes[i]->y = fyx * ver.x + fyy * ver.y + fyz * ver.z + fyb;
        vertexes[i]->z = fzx * ver.x + fzy * ver.y + fzz * ver.z + fzb;
    }
    for (int i = 0; i < numFaces; ++i)
        meshes[i]->scale(fxx, fxy, fxz, fxb, fyx, fyy, fyz, fyb, fzx, fzy, fzz, fzb);
    calcCenter();
}

void Object::rotXZ(double theta) {
    calcCenter();
    for (int i = 0; i < numVertexes; ++i) {
        Vec3d _d = *vertexes[i] - *center;
        *vertexes[i] = *center + Vec3d(cos(theta) * _d.x - sin(theta) * _d.z,  _d.y, sin(theta) * _d.x + cos(theta) * _d.z);
    }
}

Vec3d TriangularFace::min() {
    return ::min(*a, ::min(*b, *c));
}

Vec3d TriangularFace::max() {
    return ::max(*a, ::max(*b, *c));
}

Vec3d TriangularFace::center() {
    return (*a + *b + *c) / 3;
}

pair<double, Vec3d> TriangularFace::intersect(Ray ray) {
    Vec3d E1 = *a - *b, E2 = *a - *c, S = *a - ray.s;
    double t = det(S, E1, E2);
    double beta = det(ray.d, S, E2);
    double gamma = det(ray.d, E1, S);
    double n = det(ray.d, E1, E2);
    t /= n;
    beta /= n;
    gamma /= n;
    if (!(-EPSILON <= beta && beta <= 1 + EPSILON &&
          -EPSILON <= gamma && gamma <= 1 + EPSILON && beta + gamma <= 1 + EPSILON))
        t = -1;
    
    Vec3d norm = cross(*b - *a, *c - *a);
    if (dot(norm, ray.d) > 0) norm = norm * -1;
    norm.normalize();
    
    return make_pair(t, norm);
}

double TriangularFace::intersectPlane(Ray ray) {
    Vec3d E1 = *a - *b, E2 = *a - *c, S = *a - ray.s;
    double t = det(S, E1, E2);
    double n = det(ray.d, E1, E2);
    t /= n;
    return t;
}

Vec3d SphereFace::min() {
    return c - Vec3d(r, r, r);
}

Vec3d SphereFace::max() {
    return c + Vec3d(r, r, r);
}

Vec3d SphereFace::center() {
    return c;
}

pair<double, Vec3d> SphereFace::intersect(Ray ray) {
    Vec3d l = c - ray.s;
    double dis2 = l.norm2();
    int pos = 0;
    if (dis2 > r * r) pos = 1;
    else if (dis2 < r * r) pos = -1;
    double tp = dot(l, ray.d);
    if (pos > 0 && tp < 0) return make_pair(-1, Vec3d(0, 0, 0));
    double d2 = dis2 - tp * tp;
    if (d2 > r * r || pos == 0) return make_pair(-1, Vec3d(0, 0, 0));
    double t;
    if (pos > 0) t = tp - sqrt(r * r - d2);
    else if (pos < 0) t = tp + sqrt(r * r - d2);
    Vec3d norm = ray.s + ray.d * t - c;
    if (dot(norm, ray.d) > 0) norm = norm * -1;
    norm.normalize();
    return make_pair(t, norm);
}

Vec3d DiscFace::min() {
    return c - Vec3d(r, 0, r);
}

Vec3d DiscFace::max() {
    return c + Vec3d(r, 0, r);
}

Vec3d DiscFace::center() {
    return c;
}

pair<double, Vec3d> DiscFace::intersect(Ray ray) {
    double t = 0;
    if (fabs(ray.d.y) < EPSILON && fabs(ray.s.y - c.y) > EPSILON)
        return make_pair(-1, Vec3d(0, 0, 0));
    else t = (c.y - ray.s.y) / ray.d.y;
    Vec3d p = ray.s + ray.d * t;
    if ((p.x - c.x)*(p.x - c.x) + (p.z - c.z)*(p.z - c.z) <= r * r)
        return make_pair(t, Vec3d(0, -1, 0));
    else return make_pair(-1, Vec3d(0, 0, 0));
}

Vec3d BezierFace::min() {
    return m_min;
}

Vec3d BezierFace::max() {
    return m_max;
}

Vec3d BezierFace::center() {
    return m_center;
}

BezierFace::BezierFace(int n, int m, Vec3d **p, TextureMapper *texture, int brdf) {
    this->n = n;
    this->m = m;
    this->p = new Vec3d*[n + 1];
    this->texture = texture;
    this->brdf = brdf;
    m_min = Vec3d(1e100, 1e100, 1e100);
    m_max = m_min * -1;
    m_center = Vec3d(0, 0, 0);
    for (int i = 0; i <= n; ++i) {
        this->p[i] = new Vec3d[m + 1];
        for (int j = 0; j <= m; ++j) {
            this->p[i][j] = p[i][j];
            m_min = ::min(m_min, p[i][j]);
            m_max = ::max(m_max, p[i][j]);
            m_center = m_center + p[i][j];
        }
    }
    m_center = m_center / ((n + 1) * (m + 1));
    int nm = std::max(n, m);
    binom = new int*[nm + 1];
    for (int i = 0; i <= nm; ++i)
        binom[i] = new int[nm + 1];
    binom[0][0] = 1;
    for (int i = 1; i <= n || i <= m; ++i) {
        binom[i][0] = 1;
        for (int j = 1; j <= i; ++j)
            binom[i][j] = binom[i - 1][j] + binom[i - 1][j - 1];
    }
}

double BezierFace::B(int n, int k, double u) {
    return binom[n][k] * pow(u, k) * pow(1 - u, n - k);
}

double BezierFace::dB(int n, int k, double u) {
    return binom[n][k] * (k * pow(u, k - 1) * pow(1 - u, n - k) - (n - k) * pow(u, k) * pow(1 - u, n - k - 1));
}

Vec3d BezierFace::P(double u, double v) {
    Vec3d res;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= m; ++j)
            res = res + p[i][j] * B(n, i, u) * B(m, j, v);
    return res;
}

Vec3d BezierFace::F(Vec3d x, Ray ray) {
    return ray.s + ray.d * x.x - P(x.y, x.z);
}

Vec3d BezierFace::d(Vec3d x, Ray ray) {
    double jac[3][3], jac_inv[3][3];
    jac[0][0] = ray.d.x;
    jac[1][0] = ray.d.y;
    jac[2][0] = ray.d.z;
    Vec3d du, dv;
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= m; ++j) {
            du = du - p[i][j] * dB(n, i, x.y) * B(m, j, x.z);
            dv = dv - p[i][j] * B(n, i, x.y) * dB(m, j, x.z);
        }
    jac[0][1] = du.x;
    jac[1][1] = du.y;
    jac[2][1] = du.z;
    jac[0][2] = dv.x;
    jac[1][2] = dv.y;
    jac[2][2] = dv.z;
    
    jac_inv[0][0] = jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1];
    jac_inv[0][1] = jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2];
    jac_inv[0][2] = jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1];
    jac_inv[1][0] = jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2];
    jac_inv[1][1] = jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0];
    jac_inv[1][2] = jac[1][0] * jac[0][2] - jac[0][0] * jac[1][2];
    jac_inv[2][0] = jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0];
    jac_inv[2][1] = jac[0][1] * jac[2][0] - jac[0][0] * jac[2][1];
    jac_inv[2][2] = jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1];
    double d = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
    - jac[1][0] * (jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1])
    + jac[2][0] * (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]);
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j)
            jac_inv[i][j] /= d;
    
    Vec3d f = F(x, ray);
    return Vec3d(jac_inv[0][0] * f.x + jac_inv[0][1] * f.y + jac_inv[0][2] * f.z, jac_inv[1][0] * f.x + jac_inv[1][1] * f.y + jac_inv[1][2] * f.z, jac_inv[2][0] * f.x + jac_inv[2][1] * f.y + jac_inv[2][2] * f.z);
}

void BezierFace::scale(double fxx, double fxy, double fxz, double fxb,
                       double fyx, double fyy, double fyz, double fyb,
                       double fzx, double fzy, double fzz, double fzb ) {
    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= m; ++j) {
            Vec3d _p = p[i][j];
            p[i][j].x = fxx * _p.x + fxy * _p.y + fxz * _p.z + fxb;
            p[i][j].y = fyx * _p.x + fyy * _p.y + fyz * _p.z + fyb;
            p[i][j].z = fzx * _p.x + fzy * _p.y + fzz * _p.z + fzb;
        }
    m_min = Vec3d(1e100, 1e100, 1e100);
    m_max = m_min * -1;
    m_center = Vec3d(0, 0, 0);
    for (int i = 0; i <= n; ++i) {
        for (int j = 0; j <= m; ++j) {
            m_min = ::min(m_min, p[i][j]);
            m_max = ::max(m_max, p[i][j]);
            m_center = m_center + p[i][j];
        }
    }
    m_center = m_center / ((n + 1) * (m + 1));
}

pair<double, Vec3d> BezierFace::intersect(Ray ray) {
    double mint = 1e100;
    double resu, resv;
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= m; ++j) {
            Vec3d min = Vec3d(1e100, 1e100, 1e100);
            Vec3d max = min * -1;
            
            for (int _i = i - 1; _i <= i; ++_i)
                for (int _j = j - 1; _j <= j; ++_j) {
                    Vec3d p = this->p[_i][_j];
                    min = ::min(min, p);
                    max = ::max(max, p);
                }
            
            double t = -1e100;
            if (!(ray.s >= min && ray.s <= max)) { // outside
                if (fabs(ray.d.x) > 0)
                    t = std::max(t, std::min((min.x - ray.s.x) / ray.d.x, (max.x - ray.s.x) / ray.d.x));
                if (fabs(ray.d.y) > 0)
                    t = std::max(t, std::min((min.y - ray.s.y) / ray.d.y, (max.y - ray.s.y) / ray.d.y));
                if (fabs(ray.d.z) > 0)
                    t = std::max(t, std::min((min.z - ray.s.z) / ray.d.z, (max.z - ray.s.z) / ray.d.z));
                if (t < 0) continue;
                Vec3d pp = ray.s + ray.d * t;
                if (!(pp >= min && pp <= max)) continue;
            }
            else t = 0;
            
            for (int __ = 0; __ < 20; ++__) {
                double u = Utils::random(0, 1), v = Utils::random(0, 1);
                Vec3d x = Vec3d(t, u, v);
                double last = F(x, ray).norm2();
                for (int _ = 0; _ < 50; ++_) {
                    x = x - d(x, ray);
                    double cost = F(x, ray).norm2();
                    if (!(x.norm2() <= 1e3)) {
                        x.y = x.z = 1e100;
                        break;
                    }
                    if (last - cost < 1e-8) break;
                    last = cost;
                }
                t = x.x; u = x.y; v = x.z;
                if (0 <= u && u <= 1 && 0 <= v && v <= 1 && F(x, ray).norm2() < 1e-5) {
                    mint = std::min(mint, t);
                    resu = u;
                    resv = v;
                    break;
                }
            }
        }
    
    if (mint < 1e100) {
        Vec3d du, dv;
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= m; ++j) {
                du = du - p[i][j] * dB(n, i, resu) * B(m, j, resv);
                dv = dv - p[i][j] * B(n, i, resu) * dB(m, j, resv);
            }
        
        return make_pair(mint, cross(du, dv));
    }
    else return make_pair(-1, Vec3d());
}

Sphere::Sphere(Vec3d c, double r, TextureMapper *texture, int brdf) {
    numVertexes = 0;
    numFaces = 1;
    meshes = new Mesh*[1] {
        new SphereFace(c, r, texture, brdf)
    };
    center = new Vec3d(c);
}
