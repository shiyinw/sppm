//
//  Vector.hpp
//  SPPM
//
//  Created by Sherilyn Wankins on 1/9/19.
//  Copyright © 2019 Sherilyn Wankins. All rights reserved.
//

#ifndef Vector_hpp
#define Vector_hpp
#include "Utils.hpp"

class Vector {
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
    Vector() {
        x = y = z = 0;
    }
    Vector(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    double norm2() const {
        return x * x + y * y + z * z;
    }
    Vector normalize() {
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

class Matrix {
public:
    int n;
    double **v;
    Matrix(int n) {
        this->n = n;
        v = new double*[n];
        for (int i = 0; i < n; ++i)
            v[i] = new double[n];
    }
    ~Matrix() {
        for (int i = 0; i < n; ++i)
            delete[] v[i];
        delete[] v;
    }
    Matrix inv();
    double& at(int i, int j) {
        return v[i - 1][j - 1];
    }
};

double sqr(double x);
Vector operator*(const Matrix &a, const Vector &b);
Vector operator+(const Vector &a, const Vector &b);
Vector operator-(const Vector &a, const Vector &b);
Vector operator*(const Vector &a, const Vector &b);
Vector operator/(const Vector &a, const Vector &b);
Vector operator*(const Vector &a, const double &b);
Vector operator/(const Vector &a, const double &b);
Vector min(const Vector &a, const Vector &b);
Vector max(const Vector &a, const Vector &b);
Vector cross(const Vector &a, const Vector &b);
bool operator<=(const Vector &a, const Vector &b);
bool operator>=(const Vector &a, const Vector &b);
double dot(const Vector &a, const Vector &b);

// compute det([a b c])
double det(const Vector &a, const Vector &b, const Vector &c);


#endif /* Vector_hpp */
