#include <cmath>
#include "pch.h"
#include "../ball_group.hpp"

// TEST(TestAddProjectile, OutsideClusterRadius) {
//	Ball_group group(1);
//
//	for (size_t i = 0; i < 1000; i++) {
//		group = group.add_projectile();
//	}
//
//	for (size_t i = 1; i < group.num_particles; i++) {
//		ASSERT_GT(group.pos[i].norm(), 3);
//	}
//}

TEST(TestRandomDistribution, GaussianCbrt)
{
    for (size_t i = 0; i < 1000; i++) {
        double radius = cbrt(rand_between(0, 5));
        std::cout << rand_vec3(5) << '\n';
    }
}
