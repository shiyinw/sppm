#ifndef VEC_HPP
#define VEC_HPP

#include <cstring>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

class Vec3f{
public:
    
    //Constructors
    Vec3f();
    Vec3f(float x,float y, float z);
    Vec3f(const Vec3f& v);
    //Deconstructor
    virtual ~Vec3f();
    
    static Vec3f multi_out(Vec3f const &a, Vec3f const &b) {
        return Vec3f(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    }
    void normal() {
        float t = sqrt(x * x + y * y + z * z);
        x /= t;
        y /= t;
        z /= t;
    }
    
    friend std::ostream &operator <<(ostream &, Vec3f const &);
    
public:
    //Operators
    
    //Operator []
    inline float& operator [](int index){
        assert(index>=0&&index<3);
        return _p[index];
    }
    inline const float& operator [](int index) const{
        assert(index >= 0 && index<3);
        return _p[index];
    }
    
    //Operator =
    Vec3f& operator = (const Vec3f& v);
    
    //Operators +=,-=, *=, /=
    void operator +=(const Vec3f& v);
    void operator +=(float f);
    void operator -=(const Vec3f& v);
    void operator -=(float f);
    void operator *=(const Vec3f& v);
    void operator *=(float f);
    void operator /=(const Vec3f& v);
    void operator /=(float f);
    
    //Operators +,-.*,/
    Vec3f operator +(const Vec3f&v) const;
    Vec3f operator +(float f) const;
    Vec3f operator -(const Vec3f&v) const;
    Vec3f operator -(float f) const;
    Vec3f operator *(const Vec3f&v) const;
    Vec3f operator *(float f) const;
    Vec3f operator /(const Vec3f&v) const;
    Vec3f operator /(float f) const;
    
    Vec3f operator -() const;
    
public:
    void Normalize();
    float L2Norm_Sqr();
    
public:
    union
    {
        struct
        { float _p[3]; };
        struct
        { float x,y,z; };
        struct
        { float r,g,b; };
    };
    
};

#endif // VEC3F_H
