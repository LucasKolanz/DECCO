#include <limits>
#include <vector>
#include <cmath>        // std::abs
#include "vec3.hpp"

struct g
{
	bool has_balls = false;
	std::vector<int> balls;
};


class grid
{
public:
	vec3 *pos;
	int num_particles;
	std::vector<g> grids;
	int num_cells = -1;
	double cell_size;

	grid() = default;
	grid(vec3 *p, int n, double csize): pos(p), num_particles(n), cell_size(csize)
	{
		num_cells = get_num_cells();
	}

	void construct_grid();
	double min_max();
	inline int get_num_cells();
	inline int one_to_three(int index, int xyz);
	inline int three_to_one(int x, int y, int z);

private:

};

double grid::min_max()
{
	double max_val = -1;
	for (int i = 0; i < num_particles; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			if (max_val < std::abs(pos[i][j]))
			{
				max_val = std::abs(pos[i][j]);
			}
		}
	}
	return max_val;
}

inline int grid::get_num_cells(){return std::pow(2*ceil(min_max()/cell_size),3);}

inline int grid::three_to_one(int x, int y, int z)
{
	return x+(y*num_cells)+(z*num_cells*num_cells);
}

inline int grid::one_to_three(int index, int xyz)
{
	if (xyz == 2)
	{
		int z = index/(num_cells*num_cells);
		return z;
	}
	if (xyz == 1)
	{
		int z = index/(num_cells*num_cells);
		int y = (index - z * num_cells * num_cells)/num_cells;
		return y;
	}
	if (xyz == 0)
	{
		int z = index/(num_cells*num_cells);
		int y = (index - z * num_cells * num_cells)/num_cells;
		int x = index - z * num_cells * num_cells - y * num_cells;
		return x;
	}
	else
	{
		std::cerr<<"ERROR in grid: int xyz in one_to_three can be 0, 1, or 2. Not "<<xyz<<std::endl;
		return -1;
	}
}

void grid::construct_grid()
{
	grids = std::vector<g> (num_cells);
	return;
}