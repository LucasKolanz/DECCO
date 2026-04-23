#include "vec3.hpp"
#include "Utils.hpp"

#include <array>
using mat3 = std::array<vec3, 3>;




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


mat3 inverse3x3(const mat3& A)
{
    double A_det = A[0][0]*(A[1][1]*A[2][2] - A[1][2]*A[2][1]) -
                     A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) +
                     A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]);

    if (A_det <= 1e-100)
    {
        MPIsafe_print(std::cerr,"ERROR: matrix determinent is "+dToSci(A_det)+". Exiting to avoid divide by zero. Now exiting. . .\n");
        MPIsafe_exit(-1);
    }

    double inv_A_det = 1.0/A_det;
    mat3 A_inverse = {{
        vec3{(A[1][1]*A[2][2] - A[1][2]*A[2][1]) * inv_A_det, (A[0][2]*A[2][1] - A[0][1]*A[2][2]) * inv_A_det, (A[0][1]*A[1][2] - A[0][2]*A[1][1]) * inv_A_det},
        vec3{(A[1][2]*A[2][0] - A[1][0]*A[2][2]) * inv_A_det, (A[0][0]*A[2][2] - A[0][2]*A[2][0]) * inv_A_det, (A[0][2]*A[1][0] - A[0][0]*A[1][2]) * inv_A_det},
        vec3{(A[1][0]*A[2][1] - A[1][1]*A[2][0]) * inv_A_det, (A[0][1]*A[2][0] - A[0][0]*A[2][1]) * inv_A_det, (A[0][0]*A[1][1] - A[0][1]*A[1][0]) * inv_A_det}
    }};

    return A_inverse;
}

vec3 matTimesVec(const mat3& A, const vec3& v)
{
    return 
    {
        A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2],
        A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2],
        A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2],
    };
}