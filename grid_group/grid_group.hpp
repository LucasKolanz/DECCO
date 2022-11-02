#include <iostream>
#include <vector>
// #include <fstream>
#include <unordered_map>
#include <cmath>
#include <climits>
// #include <queue>
// #include <numeric>
#include "../vec3.hpp"

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
		// for (int i = 0; i < numBalls; ++i)
		// {
		// 	gridIDs[i] = std::vector<int>(3);	
		// }

		// std::cout<<"grid num_particles: "<<numBalls<<'\n';
  //       std::cout<<"grid gridSize: "<<gridSize<<'\n';
  //       std::cout<<"grid tolerance: "<<tolerance<<'\n';
  //       std::cout<<"grid radius: "<<rad<<'\n';

		init();
		return;
	}

	void init()
	{
		findGroups();
		mapGroups();
		return;
	}

	// ~grid() 
	// {
	// 	// delete[] gridIDs;
 // 	}

 	grid(grid &&) = default;
	// //Copy constructor
	// grid(const grid& copyme)
	// {
	// 	numBalls = copyme.numBalls;
	// 	gridIDs = copyme.gridIDs;
	// 	pos = copyme.pos;
	// 	gridSize = copyme.gridSize;
	// 	rad = copyme.rad;
	// }

	// //Assignment operator
	// grid& operator=(const grid& src)
	// {
	// 	if (this != &src)
	// 	{
	// 		numBalls = src.numBalls;
	// 		gridIDs = src.gridIDs;
	// 		pos = src.pos;
	// 		gridSize = src.gridSize;
	// 		rad = src.rad;
	// 	}
	// 	return *this;
	// }

	void findGroups()
	{
		for (int i = 0; i < numBalls; ++i)
		{	
			gridIDs[i] = std::vector<int>(3);
			for (int j = 0; j < 3; j++)
			{
				gridIDs[i][j] = floor(pos[i][j]/gridSize);
			}
		}
  		return;
	}

	void mapGroups()
	{
		// if (IDToGrid.size() != 0)
		// {
		// 	IDToGrid.clear();
		// }

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

	//@brief finds the necessary size of the overall grid
	// void setMaxGridIndex()
	// {
	// 	// static vec3 extrema[2];
	// 	// extrema[0] = pos[0];
	// 	// extrema[1] = pos[0];
	// 	double max = -1;
	// 	for (int i = 0; i < numBalls; i++)
	// 	{
	// 		for (int j = 0; j < 3; j++)
	// 		{
	// 			if (abs(pos[i][j]) > max)
	// 			{
	// 				max = abs(pos[i][j]);	
	// 			}
	// 		}
	// 	}
	// 	// maxGridIndex = ceil(max/gridSize);
	// }

	std::vector<int> getBalls(int ballIndex)
	{
		// for (int i = 0; i < numBalls; ++i)
		// {
		// 	std::cout<<pos[i]<<'\n';
		// }
		// std::vector<std::vector<int>> neededGroups(27);
		//check if ball is too close to a boarder
		// double tolerancesq = tolerance * tolerance;
		// bool addVertex = false;
		int currx, curry, currz;
		currx = gridIDs[ballIndex][0];
		curry = gridIDs[ballIndex][1]; 
		currz = gridIDs[ballIndex][2]; 

		// vec3 v1(gridIDs[ballIndex][0]*gridSize, gridIDs[ballIndex][1]*gridSize, gridIDs[ballIndex][2]*gridSize); 

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

		// //add all points in the groups in neededGroups to a vector 
		// //of ints and return that
		// for (auto it = begin(neededGroups); it != end(neededGroups); ++it)
		// {
			
		// }

		return ball_indicies;
	}

	// void resetGroups()
	// {
	// 	for (int i = 0; i < numBalls; i++)
	// 	{
	// 		gridIDs[i] = std::vector<int> {INT_MAX,INT_MAX,INT_MAX};
	// 	}
	// 	// numGroups = -1;
	// 	// delete[] groups;
	// 	return;
	// }
	
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

	// template<typename K, typename V>
	// void printMap(std::unordered_map<K, V> const &m)
	// {
	//     for (auto const &pair: m) {
	//         std::cout << "{" << pair.first << ": ";
	//         for (auto it = begin(pair.second); it != end(pair.second); ++it)
	//         {
	//         	std::cout<<pair.second[*it];
	//         	if (it+1 != end(pair.second))
	//         	{
	//         		std::cout<<",";
	//         	}
	//         }
	//         std::cout<<"}" << std::endl;
	//     }
	// }

	// void printMap()
	// {
	//     for (auto mit = IDToGrid.cbegin(); mit != IDToGrid.cend(); mit++) {
	//         std::cout << "{" << (*mit).first << ": ";
	//         for (auto vit = begin((*mit).second); vit != end((*mit).second); ++vit)
	//         {
	//         	std::cout<<(*mit).second[*vit];
	//         	if (vit+1 != end((*mit).second))
	//         	{
	//         		std::cout<<",";
	//         	}
	//         }
	//         std::cout<<"}" << std::endl;
	//     }
	// }


	// template <typename T>
	// void printArray(T *arr, int size)
	// {
	// 	std::cout<<"================Array================"<<std::endl;
	// 	for (int i = 0; i < size; i++)
	// 	{
	// 		std::cout<<arr[i]<<", ";
	// 	}
	// 	std::cout<<std::endl;
	// 	std::cout<<"====================================="<<std::endl;
	// 	return;
	// }

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