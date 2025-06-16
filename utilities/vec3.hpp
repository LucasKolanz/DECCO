#pragma once

#ifndef VEC3_HPP
#define VEC3_HPP

#include <iostream>

class vec3
{
public:
    double x;
    double y;
    double z;

    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3();

    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3(const double newx, const double newy, const double newz);

    // vec3(const vec3& v) : x(v.x), y(v.y), z(v.z)
    //{
    //}

    // vec3& operator=(const vec3& v)
    //{
    //	x = v.x;
    //	y = v.y;
    //	z = v.z;
    //	return *this;
    //}

    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3 operator-() const ;
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3 operator+(const vec3& v) const;
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3 operator+=(const vec3& v);
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3 operator-(const vec3& v) const;
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3 operator-=(const vec3& v);
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3 operator*(const double scalar) const;
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3 operator*=(const double scalar);
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3 operator/(const double scalar) const;
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3 operator/=(const double scalar);
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    double& operator[](const int i);
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    double operator[](const int i) const;

    // bool operator==(const vec3& v) const
    //{
    //	return ((x == v.x) && (y == v.y) &&
    //		(z == v.z));
    //}
    // bool operator!=(const vec3& v) const
    //{
    //	return !(*this == v);
    //}

    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    [[nodiscard]] double dot(const vec3& v) const;
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    [[nodiscard]] vec3 cross(const vec3& v) const;
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    double norm() const;
    double normsquared() const;
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3 normalized() const;
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    vec3 normalized_safe() const;
    void print() const;
    #ifdef GPU_ENABLE
        #pragma acc routine seq
    #endif
    [[nodiscard]] vec3 rot(char axis, double angle) const;
    // [[nodiscard]] vec3 rot(char axis, double angle);
    vec3 arbitrary_orthogonal() const;
};

inline vec3
operator*(const double scalar, const vec3& v)
{
    return v * scalar;
}

// Output vec3 to console easily.
inline std::ostream&
operator<<(std::ostream& s, const vec3& v)
{
    return s << v.x << ',' << v.y << ',' << v.z;
}

// class rotation
// {
// public:
//     rotation()
//     {
//         w = 1.0;
//         x = y = z = 0;
//     }
//     rotation(const rotation& q)
//     {
//         w = q.w;
//         x = q.x;
//         y = q.y;
//         z = q.z;
//     }
//     rotation(const double angle, const vec3& axis)
//     {
//         w = cos(angle / 2.0);
//         const vec3 goodaxis = axis.normalized();
//         const double sinangle_over2 = sin(angle / 2.0);
//         x = goodaxis.x * sinangle_over2;
//         y = goodaxis.y * sinangle_over2;
//         z = goodaxis.z * sinangle_over2;
//     }

//     rotation& operator=(const rotation& q)
//     {
//         w = q.w;
//         x = q.x;
//         y = q.y;
//         z = q.z;
//         return *this;
//     }

//     rotation operator*(const rotation& q) const
//     {
//         return rotation(
//                    w * q.w - x * q.x - y * q.y - z * q.z,
//                    w * q.x + x * q.w + y * q.z - z * q.y,
//                    w * q.y - x * q.z + y * q.w + z * q.x,
//                    w * q.z + x * q.y - y * q.x + z * q.w)
//             .normalized();
//     }
//     rotation operator*=(const rotation& q)
//     {
//         *this = q * (*this);
//         return *this;
//     }

//     bool operator==(const rotation& q) const
//     {
//         return ((w == q.w) && (x == q.x) && (y == q.y) && (z == q.z));
//     }
//     bool operator!=(const rotation& q) const { return !(*this == q); }

//     rotation conj() const { return rotation(w, -x, -y, -z); }

//     vec3 rotate_vector(const vec3& v) const
//     {
//         const rotation p(
//             v.x * x + v.y * y + v.z * z,
//             v.x * w - v.y * z + v.z * y,
//             v.x * z + v.y * w - v.z * x,
//             -v.x * y + v.y * x + v.z * w);
//         const rotation product(
//             w * p.w - x * p.x - y * p.y - z * p.z,
//             w * p.x + x * p.w + y * p.z - z * p.y,
//             w * p.y - x * p.z + y * p.w + z * p.x,
//             w * p.z + x * p.y - y * p.x + z * p.w);
//         return vec3(product.x, product.y, product.z);
//     }

//     // void tostr(char str[]) const {
//     //	const double theta = 2 * acos(w);
//     //	const double fac = 1.0 / sin(theta / 2.0);
//     //	sfprintf(stderr,str, "[%6.2f, (%6.2f, %6.2f, %6.2f)]", theta, x*fac, y*fac, z*fac);
//     //}

// private:
//     double w;
//     double x;
//     double y;
//     double z;

//     rotation(const double neww, const double newx, const double newy, const double newz)
//     {
//         w = neww;
//         x = newx;
//         y = newy;
//         z = newz;
//     }

//     rotation operator/(const double scalar) const
//     {
//         return rotation(w / scalar, x / scalar, y / scalar, z / scalar);
//     }
//     rotation operator/=(const double scalar)
//     {
//         w /= scalar;
//         x /= scalar;
//         y /= scalar;
//         z /= scalar;
//         return *this;
//     }

//     rotation normalized() const { return *this / sqrt(w * w + x * x + y * y + z * z); }
// };




#endif