//
//  Vector.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#ifndef Vec3d_hpp
#define Vec3d_hpp
#include "Random.hpp"

class Vec3d {
public:
    union
    {
        struct
        { double _p[3]; };
        struct
        { double x,y,z; };
        struct
        { double r,g,b; };
    };
    Vec3d() {
        x = y = z = 0;
    }
    Vec3d(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    double norm2() const {
        return x * x + y * y + z * z;
    }
    Vec3d normalize() {
        double norm = sqrt(norm2());
        x /= norm; y /= norm; z /= norm;
        return *this;
    }
    double max() {
        return std::max(x, std::max(y, z));
    }
    double min() {
        return std::min(x, std::min(y, z));
    }
    void print() const {
        printf("%.5lf %.5lf %.5lf\n", x, y, z);
    }
};

Vec3d operator+(const Vec3d &a, const Vec3d &b);
Vec3d operator-(const Vec3d &a, const Vec3d &b);
Vec3d operator*(const Vec3d &a, const Vec3d &b);
Vec3d operator/(const Vec3d &a, const Vec3d &b);
Vec3d operator*(const Vec3d &a, const double &b);
Vec3d operator/(const Vec3d &a, const double &b);
Vec3d min(const Vec3d &a, const Vec3d &b);
Vec3d max(const Vec3d &a, const Vec3d &b);
Vec3d cross(const Vec3d &a, const Vec3d &b);
bool operator<=(const Vec3d &a, const Vec3d &b);
bool operator>=(const Vec3d &a, const Vec3d &b);
double dot(const Vec3d &a, const Vec3d &b);

// compute det([a b c])
double det(const Vec3d &a, const Vec3d &b, const Vec3d &c);


#endif /* Vector_hpp */
