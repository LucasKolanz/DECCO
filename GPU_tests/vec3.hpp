#pragma once

#include <iostream>
#include <cassert>
#include <string>
#include <random>
#include <sstream>
#include <iomanip>
#include <omp.h>

// #pragma omp declare target
// #pragma acc declare
class vec3
{
public:

    double x;
    double y;
    double z;

    vec3()
        : x(0)
        , y(0)
        , z(0)
    {
    }


    // #pragma acc routine
    vec3(const double newx, const double newy, const double newz)
        : x(newx)
        , y(newy)
        , z(newz)
    {
    }

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


    vec3 operator-() const { return vec3(-x, -y, -z); }

    vec3 operator+(const vec3& v) const { return vec3(x + v.x, y + v.y, z + v.z); }

    // #pragma acc routine
    vec3 operator+=(const vec3& v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    vec3 operator-(const vec3& v) const { return vec3(x - v.x, y - v.y, z - v.z); }

    // #pragma acc routine
    vec3 operator-=(const vec3& v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    vec3 operator*(const double scalar) const { return vec3(scalar * x, scalar * y, scalar * z); }

    vec3 operator*=(const double scalar)
    {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    vec3 operator/(const double scalar) const { return vec3(x / scalar, y / scalar, z / scalar); }

    vec3 operator/=(const double scalar)
    {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }

    // bool operator==(const vec3& v) const
    //{
    //	return ((x == v.x) && (y == v.y) &&
    //		(z == v.z));
    //}
    // bool operator!=(const vec3& v) const
    //{
    //	return !(*this == v);
    //}

    // #pragma acc routine seq
    double& operator[](const int i)
    {
        switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        }
        assert(0);
    }
    // #pragma acc routine seq
    double operator[](const int i) const
    {
        switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        }
        assert(0);
    }

    [[nodiscard]] double dot(const vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    [[nodiscard]] vec3 cross(const vec3& v) const
    {
        return vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }
    double norm() const { return sqrt(x * x + y * y + z * z); }
    double normsquared() const { return x * x + y * y + z * z; }
    vec3 normalized() const { return *this / this->norm(); }

    vec3 normalized_safe() const
    {
        if (fabs(this->norm()) < 1e-13) {
            std::cerr << "dividing by zero in unit vector calculation!!!!!!!!!!!!!!";
            system("pause");
        }
        return *this / this->norm();
    }

    void print() const { std::cout << x << ',' << y << ',' << z << "\n"; }

    [[nodiscard]] vec3 rot(char axis, double angle) const
    {
        double rotx[3][3] = {{1, 0, 0}, {0, cos(angle), -sin(angle)}, {0, sin(angle), cos(angle)}},
               roty[3][3] = {{cos(angle), 0, sin(angle)}, {0, 1, 0}, {-sin(angle), 0, cos(angle)}},
               rotz[3][3] = {{cos(angle), -sin(angle), 0}, {sin(angle), cos(angle), 0}, {0, 0, 1}};
        vec3 newVec;
        switch (axis) {
        case 'x':
            newVec[0] = rotx[0][0] * x + rotx[0][1] * y + rotx[0][2] * z;
            newVec[1] = rotx[1][0] * x + rotx[1][1] * y + rotx[1][2] * z;
            newVec[2] = rotx[2][0] * x + rotx[2][1] * y + rotx[2][2] * z;
            break;
        case 'y':
            newVec[0] = roty[0][0] * x + roty[0][1] * y + roty[0][2] * z;
            newVec[1] = roty[1][0] * x + roty[1][1] * y + roty[1][2] * z;
            newVec[2] = roty[2][0] * x + roty[2][1] * y + roty[2][2] * z;
            break;
        case 'z':
            newVec[0] = rotz[0][0] * x + rotz[0][1] * y + rotz[0][2] * z;
            newVec[1] = rotz[1][0] * x + rotz[1][1] * y + rotz[1][2] * z;
            newVec[2] = rotz[2][0] * x + rotz[2][1] * y + rotz[2][2] * z;
            break;
        default:
            std::cerr << "Must choose x, y, or z rotation axis.";
            break;
        }
        return newVec;
    }

    vec3 arbitrary_orthogonal() const
    {
        bool b0 = (x < y) && (x < z);
        bool b1 = (y <= x) && (y < z);
        bool b2 = (z <= x) && (z <= y);

        return this->cross(vec3(int(b0), int(b1), int(b2)));
    }
};
// #pragma omp end declare target

#pragma acc routine seq
inline vec3
operator*(const double scalar, const vec3& v)
{
    return v * scalar;
}


class rotation
{
public:
    rotation()
    {
        w = 1.0;
        x = y = z = 0;
    }
    rotation(const rotation& q)
    {
        w = q.w;
        x = q.x;
        y = q.y;
        z = q.z;
    }
    rotation(const double angle, const vec3& axis)
    {
        w = cos(angle / 2.0);
        const vec3 goodaxis = axis.normalized();
        const double sinangle_over2 = sin(angle / 2.0);
        x = goodaxis.x * sinangle_over2;
        y = goodaxis.y * sinangle_over2;
        z = goodaxis.z * sinangle_over2;
    }

    rotation& operator=(const rotation& q)
    {
        w = q.w;
        x = q.x;
        y = q.y;
        z = q.z;
        return *this;
    }

    rotation operator*(const rotation& q) const
    {
        return rotation(
                   w * q.w - x * q.x - y * q.y - z * q.z,
                   w * q.x + x * q.w + y * q.z - z * q.y,
                   w * q.y - x * q.z + y * q.w + z * q.x,
                   w * q.z + x * q.y - y * q.x + z * q.w)
            .normalized();
    }
    rotation operator*=(const rotation& q)
    {
        *this = q * (*this);
        return *this;
    }

    bool operator==(const rotation& q) const
    {
        return ((w == q.w) && (x == q.x) && (y == q.y) && (z == q.z));
    }
    bool operator!=(const rotation& q) const { return !(*this == q); }

    rotation conj() const { return rotation(w, -x, -y, -z); }

    vec3 rotate_vector(const vec3& v) const
    {
        const rotation p(
            v.x * x + v.y * y + v.z * z,
            v.x * w - v.y * z + v.z * y,
            v.x * z + v.y * w - v.z * x,
            -v.x * y + v.y * x + v.z * w);
        const rotation product(
            w * p.w - x * p.x - y * p.y - z * p.z,
            w * p.x + x * p.w + y * p.z - z * p.y,
            w * p.y - x * p.z + y * p.w + z * p.x,
            w * p.z + x * p.y - y * p.x + z * p.w);
        return vec3(product.x, product.y, product.z);
    }

    // void tostr(char str[]) const {
    //	const double theta = 2 * acos(w);
    //	const double fac = 1.0 / sin(theta / 2.0);
    //	sfprintf(stderr,str, "[%6.2f, (%6.2f, %6.2f, %6.2f)]", theta, x*fac, y*fac, z*fac);
    //}

private:
    double w;
    double x;
    double y;
    double z;

    rotation(const double neww, const double newx, const double newy, const double newz)
    {
        w = neww;
        x = newx;
        y = newy;
        z = newz;
    }

    rotation operator/(const double scalar) const
    {
        return rotation(w / scalar, x / scalar, y / scalar, z / scalar);
    }
    rotation operator/=(const double scalar)
    {
        w /= scalar;
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }

    rotation normalized() const { return *this / sqrt(w * w + x * x + y * y + z * z); }
};

// Output vec3 to console easily.
std::ostream&
operator<<(std::ostream& s, const vec3& v)
{
    return s << v.x << ',' << v.y << ',' << v.z;
}
