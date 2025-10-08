#pragma once

#ifndef VEC3_HPP
#define VEC3_HPP

#include <iostream>
#include "MPI_utilities.hpp"

class vec3
{
public:
    double x;
    double y;
    double z;

    vec3();
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

    vec3 operator-() const noexcept;
    vec3 operator+(const vec3& v) const noexcept;
    vec3& operator+=(const vec3& v);
    vec3 operator-(const vec3& v) const noexcept;
    vec3& operator-=(const vec3& v);
    vec3 operator*(const double scalar) const noexcept;
    vec3& operator*=(const double scalar);
    vec3 operator/(const double scalar) const noexcept;
    vec3& operator/=(const double scalar);
    double& operator[](const int i);
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


    [[nodiscard]] double dot(const vec3& v) const noexcept;
    [[nodiscard]] vec3 cross(const vec3& v) const noexcept;
    double norm() const noexcept;
    double normsquared() const noexcept;
    vec3 normalized() const noexcept;
    vec3 normalized_safe() const;
    void print() const;
    [[nodiscard]] vec3 rot(char axis, double angle) const;
    // [[nodiscard]] vec3 rot(char axis, double angle);
    vec3 arbitrary_orthogonal() const;
};

inline vec3 operator*(const double s, const vec3& v) { return v * s; }
inline vec3 operator*(const int s,    const vec3& v) { return v * static_cast<double>(s); }

// Output vec3 to console easily.
inline std::ostream&
operator<<(std::ostream& s, const vec3& v)
{
    return s << v.x << ',' << v.y << ',' << v.z;
}

class rotation
{
public:
    rotation();
    rotation(const rotation& q);
    rotation(const double scalar, const vec3& vector);
    // rotation(const double angle, const vec3& axis);
    rotation& operator=(const rotation& q);
    rotation operator*(const rotation& q) const noexcept;
    rotation& operator*=(const rotation& q);
    bool operator==(const rotation& q) const noexcept;
    bool operator!=(const rotation& q) const noexcept;
    rotation conj() const noexcept;
    vec3 rotate_vector(const vec3& v) const;
    vec3 worldToLocal(const vec3& vecWorld) const noexcept;
    vec3 localToWorld(const vec3& vecLocal) const noexcept;
    void exponential_integrate(const vec3 w_body, const double dt);
    vec3 get_vec() const noexcept;
    double get_scalar() const noexcept;


    // void tostr(char str[]) const;
    double w;
    double x;
    double y;
    double z;
private:

    rotation(const double neww, const double newx, const double newy, const double newz);
    rotation operator/(const double scalar) const;
    rotation operator/=(const double scalar);
    rotation normalized() const;
};

//THIS ONLY WORKS FOR UNIT VEC QUATERNIONS
vec3 quatRotate(const rotation& q, const vec3& vec);

// Output vec3 to console easily.
inline std::ostream&
operator<<(std::ostream& s, const rotation& r)
{
    return s << r.w << ',' << r.x << ',' << r.y << ',' << r.z;
}

#endif