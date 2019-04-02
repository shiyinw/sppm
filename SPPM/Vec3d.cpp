//
//  Vector.cpp
//  SPPM
//
//  Created by Sherilyn Wankins on 11/1/18.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//


#define EPSILON 1e-6


#include "Vec3d.hpp"


Vec3d operator+(const Vec3d &a, const Vec3d &b) {
    return Vec3d(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vec3d operator-(const Vec3d &a, const Vec3d &b) {
    return Vec3d(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vec3d operator*(const Vec3d &a, const Vec3d &b) {
    return Vec3d(a.x * b.x, a.y * b.y, a.z * b.z);
}

Vec3d operator/(const Vec3d &a, const Vec3d &b) {
    return Vec3d(a.x / b.x, a.y / b.y, a.z / b.z);
}

Vec3d operator*(const Vec3d &a, const double &b) {
    return Vec3d(a.x * b, a.y * b, a.z * b);
}

Vec3d operator/(const Vec3d &a, const double &b) {
    return Vec3d(a.x / b, a.y / b, a.z / b);
}

Vec3d min(const Vec3d &a, const Vec3d &b) {
    return Vec3d(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
}

Vec3d max(const Vec3d &a, const Vec3d &b) {
    return Vec3d(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
}

Vec3d cross(const Vec3d &a, const Vec3d &b) {
    return Vec3d(a.y * b.z - b.y * a.z, b.x * a.z - a.x * b.z, a.x * b.y - b.x * a.y);
}

bool operator<=(const Vec3d &a, const Vec3d &b) {
    return a.x <= b.x + EPSILON && a.y <= b.y + EPSILON && a.z <= b.z + EPSILON;
}

bool operator>=(const Vec3d &a, const Vec3d &b) {
    return a.x + EPSILON >= b.x && a.y + EPSILON >= b.y && a.z + EPSILON >= b.z;
}

double dot(const Vec3d &a, const Vec3d &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double det(const Vec3d &a, const Vec3d &b, const Vec3d &c) {
    return a.x * (b.y * c.z - b.z * c.y)
    - b.x * (a.y * c.z - a.z * c.y)
    + c.x * (a.y * b.z - a.z * b.y);
}

