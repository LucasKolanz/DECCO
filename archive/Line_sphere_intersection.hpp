#pragma once

#include <iostream>
#include <cmath>
#include "vector3d.hpp"

// I got this from https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
bool line_sphere_intersect(const vector3d& position, const vector3d& velocity, const vector3d& center, const double radius) {
	vector3d direction = velocity.normalized();

	double grad =
		direction.dot((position - center)) *
		direction.dot((position - center)) -
		(
			(position - center).normsquared() -
			radius * radius
			);

	// If I ever want to distance to the two intersections:
	//double d1 = -u_vec.dot(origin - center) + sqrtf(grad);
	//double d2 = -u_vec.dot(origin - center) - sqrtf(grad);
	//std::cout << '\t' << d1 << '\t' << d2 << '\n';

	if (grad < 0)
	{
		//std::cout << "Position"; position.print();
		//std::cout << "Direction"; direction.print();
		//std::cout << "Target"; center.print();
		//std::cout << "Radius " << radius << '\n';
		//std::cout << "Grad " << grad << '\n';
		return false;
	}
	else
	{
		//std::cout << "Intersection exists.";
		return true;
	}
}