#include <iostream>
#include <vector>
// #include <fstream>
#include <unordered_map>
#include <cmath>
#include <climits>
// #include <queue>
// #include <numeric>
#include "../vec3.hpp"
#include "../timing/timing.hpp"

#define MIN_FOR_GRID 0

class grid
{
public:
	// bool initialized = false;
	int numBalls = -1;
	double rad = -1.0;
	double gridSize = -1.0;
	double tolerance = -1.0;
	// int maxGridIndex = -1;
	std::vector<std::vector<int>> gridIDs;
	std::unordered_map<std::string,std::vector<int>> IDToGrid;
	vec3 *pos = nullptr;

	// grid() = default;

	grid(int num_balls,double radius,double grid_size,vec3 *positions,double tol=-1.0)
	{
		numBalls = num_balls;
		rad = radius;
		pos = positions;
		gridSize = grid_size;
		if (tol == -1.0)
		{
			tolerance = grid_size/2;
		}
		else
		{
			tolerance = tol;
		}
		IDToGrid.reserve(num_balls); ///need a way to predict number of grid spaces
		IDToGrid.max_load_factor(0.25);

		gridIDs.resize(numBalls);


		init();
		return;
	}

	void init()
	{
		findGroups();
		mapGroups();
		return;
	}

 	grid(grid &&) = default;


	void findGroups()
	{
		#pragma omp parallel
		{
			#pragma omp parallel for default(none) shared(gridIDs)
			for (int i = 0; i < numBalls; ++i)
			{	
				gridIDs[i] = std::vector<int>(3);
				for (int j = 0; j < 3; j++)
				{
					gridIDs[i][j] = floor(pos[i][j]/gridSize);
				}
			}
		}
  		return;
	}

	void mapGroups()
	{

		for (int i = 0; i < numBalls; ++i)
		{
			std::string key = getKey(gridIDs[i]);
			if (IDToGrid.find(key) == IDToGrid.end()) // key not present
			{
				std::vector<int> indices{i};
				IDToGrid[key] = indices;
				// std::cout<<"key: val  " <<key<<": "<<IDToGrid[key][0]<<std::endl;
			}
			else//key present, add to vector
			{
				IDToGrid[key].push_back(i);
			}
		}
		// printMap();
		return;
	}

	inline std::string getKey(std::vector<int> v)
	{
		return std::to_string(v[0]) + std::to_string(v[1]) + std::to_string(v[2]);
	}

	inline std::string getKey(int x, int y, int z)
	{
		return std::to_string(x) + std::to_string(y) + std::to_string(z);
	}


	std::vector<int> getBalls(int ballIndex)
	{

		int currx, curry, currz;
		currx = gridIDs[ballIndex][0];
		curry = gridIDs[ballIndex][1]; 
		currz = gridIDs[ballIndex][2]; 


		std::vector<int> ball_indicies;
		for (int x = -1; x < 2; ++x)
		{
			for (int y = -1; y < 2; ++y)
			{
				for (int z = -1; z < 2; ++z)
				{
					std::string key = getKey(currx+x,curry+y,currz+z);
					if (IDToGrid.find(key) != IDToGrid.end())
					{
						for (auto ind_it = begin(IDToGrid[key]); ind_it != end(IDToGrid[key]); ++ind_it)
						{
							// if (ball_indicies.find(*ind_it) == ball_indicies.end())
							if(std::find(ball_indicies.begin(), ball_indicies.end(), *ind_it) == ball_indicies.end())
							{
								ball_indicies.push_back(*ind_it);
							}
						}
					}		
				}
			}
		}		


		return ball_indicies;
	}

	
	void printGroups()
	{
		std::cout<<"================Groups================"<<std::endl;
		for (int i = 0; i < numBalls; ++i)
		{
			std::cout<<"(ball; group): "<<i<<"; {"<<gridIDs[i][0]<<","<<gridIDs[i][1]<<","<<gridIDs[i][2]<<"}"<<std::endl;
		}
		std::cout<<"======================================"<<std::endl;
		return;
	}

  	void printVector(std::vector<int> vec)
  	{
		std::cout<<"================Vector================"<<std::endl;

  		for (auto it = begin(vec); it != end(vec); ++it) 
		{
			std::cout<<*it<<", ";
		}	
		std::cout<<std::endl;
		std::cout<<"======================================"<<std::endl;
		return;
  	}


};