
#include <cmath>
#include <iostream>
#include <cassert>
#include <sstream>
#include "vec3.hpp"

vec3::vec3()
    : x(0)
    , y(0)
    , z(0)
{
}

vec3::vec3(const double newx, const double newy, const double newz)
    : x(newx)
    , y(newy)
    , z(newz)
{
}

vec3 vec3::operator-() const noexcept { return vec3(-x, -y, -z); }

vec3 vec3::operator+(const vec3& v) const noexcept { return vec3(x + v.x, y + v.y, z + v.z); }

vec3& vec3::operator+=(const vec3& v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

vec3 vec3::operator-(const vec3& v) const noexcept { return vec3(x - v.x, y - v.y, z - v.z); }

vec3& vec3::operator-=(const vec3& v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

vec3 vec3::operator*(const double scalar) const noexcept { return vec3(scalar * x, scalar * y, scalar * z); }

vec3& vec3::operator*=(const double scalar)
{
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}

vec3 vec3::operator/(const double scalar) const noexcept { return vec3(x / scalar, y / scalar, z / scalar); }

vec3& vec3::operator/=(const double scalar)
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

double& vec3::operator[](const int i)
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
double vec3::operator[](const int i) const
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

[[nodiscard]] double vec3::dot(const vec3& v) const noexcept {return x * v.x + y * v.y + z * v.z;}
[[nodiscard]] vec3 vec3::cross(const vec3& v) const noexcept {return vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);}
double vec3::norm() const noexcept { return sqrt(x * x + y * y + z * z); }
double vec3::normsquared() const noexcept { return x * x + y * y + z * z; }
vec3 vec3::normalized() const noexcept { return *this / this->norm(); }

vec3 vec3::normalized_safe() const
{
    const double n = this->norm();
    if (n < 1e-13) {
        std::cerr << "dividing by zero in unit vector calculation!!!!!!!!!!!!!!\n";
        MPIsafe_exit(-1);
    }
    return *this / n;
}

void vec3::print() const { std::cout << x << ',' << y << ',' << z << "\n"; }

[[nodiscard]] vec3 vec3::rot(char axis, double angle) const
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

vec3 vec3::arbitrary_orthogonal() const
{
    bool b0 = (x < y) && (x < z);
    bool b1 = (y <= x) && (y < z);
    bool b2 = (z <= x) && (z <= y);

    return this->cross(vec3(int(b0), int(b1), int(b2)));
}


// inline vec3
// operator*(const double scalar, const vec3& v)
// {
// return v * scalar;
// }


rotation::rotation()
{
    w = 1.0;
    x = y = z = 0;
}
rotation::rotation(const rotation& q)
{
    w = q.w;
    x = q.x;
    y = q.y;
    z = q.z;
}
rotation::rotation(const double scalar, const vec3& vector)
{
    w = scalar;
    x = vector.x;
    y = vector.y;
    z = vector.z;
}
// rotation::rotation(const double angle, const vec3& axis)
// {
//     w = cos(angle / 2.0);
//     const vec3 goodaxis = axis.normalized();
//     const double sinangle_over2 = sin(angle / 2.0);
//     x = goodaxis.x * sinangle_over2;
//     y = goodaxis.y * sinangle_over2;
//     z = goodaxis.z * sinangle_over2;
// }

rotation& rotation::operator=(const rotation& q)
{
    w = q.w;
    x = q.x;
    y = q.y;
    z = q.z;
    return *this;
}

rotation rotation::operator*(const rotation& q) const noexcept
{
    return rotation(
               w * q.w - x * q.x - y * q.y - z * q.z,
               w * q.x + x * q.w + y * q.z - z * q.y,
               w * q.y - x * q.z + y * q.w + z * q.x,
               w * q.z + x * q.y - y * q.x + z * q.w)
        .normalized();
}
rotation& rotation::operator*=(const rotation& q)
{
    *this = q * (*this);
    return *this;
}

bool rotation::operator==(const rotation& q) const noexcept
{
    return ((w == q.w) && (x == q.x) && (y == q.y) && (z == q.z));
}
bool rotation::operator!=(const rotation& q) const noexcept { return !(*this == q); }

vec3 rotation::get_vec() const noexcept {return vec3(x,y,z);}
double rotation::get_scalar() const noexcept {return w;}

rotation rotation::conj() const noexcept { return rotation(w, -x, -y, -z); }


// void tostr(char str[]) const {
//	const double theta = 2 * acos(w);
//	const double fac = 1.0 / sin(theta / 2.0);
//	sfprintf(stderr,str, "[%6.2f, (%6.2f, %6.2f, %6.2f)]", theta, x*fac, y*fac, z*fac);
//}

//We integrate attitude on SO(3) with a Lie-group midpoint (exponential map)
//  update implemented in unit quaternions (Δq right-multiplication with 
//  half-step body-frame angular velocity). Might have energy drift over very long
//  time periods. If so, try Lie-group variational integrators.
// --- build Δq from body-frame ω at half step (use your w_body) ---
void rotation::exponential_integrate(const vec3 w_body, const double dt)
{
    const double wb = w_body.norm();
    // double theta = wn * dt;
    // double half = 0.5 * theta;

    // double c, s_over;
    // if (half*half < 1e-12) {  // small angle
    //     // cos(half)
    //     double h2 = half*half;
    //     c = 1.0 - 0.5*h2 + (1.0/24.0)*h2*h2;
    //     // sin(half)/‖ω‖
    //     s_over = 0.5*dt * (1.0 - (theta*theta)/24.0 + (theta*theta*theta*theta)/1920.0);
    // } else {
    //     c = std::cos(half);
    //     s_over = std::sin(half) / wn;
    // }
    const double half = 0.5 * wb * dt;
    double s_over = (wb > 0.0) ? std::sin(half) / wb : 0.5 * dt;   // small-angle safe
    const vec3  dqv = s_over * w_body;
    const double dq0 = std::cos(half);

    // --- right multiply: q_new = q_old ⊗ Δq ---
    const vec3 qv = {x,y,z};
    const double q0 = w;

    w = q0*dq0 - qv.dot(dqv);
    const vec3 newv = q0*dqv + dq0*qv + qv.cross(dqv);
    x = newv.x;
    y = newv.y;
    z = newv.z;

    this->normalized();
    // normalize
    // double L = std::sqrt(new0*new0 + newv.normsquared());
    // Eu0[Ball] =  new0 / L;
    // Eu[Ball]  =  newv / L;
}

rotation::rotation(const double neww, const double newx, const double newy, const double newz)
{
    w = neww;
    x = newx;
    y = newy;
    z = newz;
}

rotation rotation::operator/(const double scalar) const
{
    return rotation(w / scalar, x / scalar, y / scalar, z / scalar);
}
rotation rotation::operator/=(const double scalar)
{
    w /= scalar;
    x /= scalar;
    y /= scalar;
    z /= scalar;
    return *this;
}

rotation rotation::normalized() const { return *this / sqrt(w * w + x * x + y * y + z * z); }

// Output vec3 to console easily.
// std::ostream&
// operator<<(std::ostream& s, const vec3& v)
// {
//     return s << v.x << ',' << v.y << ',' << v.z;
// // }


vec3 rotation::rotate_vector(const vec3& v) const
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

vec3 rotation::worldToLocal(const vec3& vecWorld) const noexcept {return quatRotate(conj(),vecWorld);}
vec3 rotation::localToWorld(const vec3& vecLocal) const noexcept {return quatRotate(*this,vecLocal);}

//THIS ONLY WORKS FOR UNIT VEC QUATERNIONS
//I think this should be faster than using the 
//defined quaternion multiplication and returning
//q*vec*(q^-1), but this should be tested
vec3 quatRotate(const rotation& q, const vec3& vec)
{
    // t = 2 q_vec × vec
    const vec3 v = q.get_vec();
    const double s = q.get_scalar();
    vec3 t = 2.0 * v.cross(vec);

    // return vec + s t + q_vec × t
    return vec + s * t + v.cross(t);
}

