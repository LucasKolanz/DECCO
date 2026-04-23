#pragma once

#ifndef MATH_HPP
#define MATH_HPP


mat3 inverse3x3(const mat3& A);
vec3 matTimesVec(const mat3& A, const vec3& v);

//THIS ONLY WORKS FOR UNIT VEC QUATERNIONS
vec3 quatRotate(const rotation& q, const vec3& vec);

#endif