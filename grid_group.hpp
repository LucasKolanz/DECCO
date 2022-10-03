#include <iostream>
#include <vector>
// #include <fstream>
#include <cmath>
#include <climits>
// #include <queue>
// #include <numeric>
#include "vec3.hpp"

class grid
{
public:
	// bool initialized = false;
	int numBalls = -1;
	double rad = -1.0;
	double gridSize = -1.0;
	int maxGridIndex = -1;
	std::vector<std::vector<int>> gridIDs;
	vec3 *pos = nullptr;

	grid() = default;

	grid(int num_balls,double radius,double grid_size,vec3 *positions)
	{
		numBalls = num_balls;
		rad = radius;
		pos = positions;
		gridSize = grid_size;

		std::vector<int> index = {INT_MAX,INT_MAX,INT_MAX};
		for (int i = 0; i < numBalls; i++)
		{

			gridIDs.push_back(index); 	
		}

		///TAKE THIS OUT FOR REAL THING
		//This is just meant to make multiple groups out of a single group
		// srand(0);
		// for (int i = 0; i < 4; i++)
		// {
		// 	double x,y,z;
		// 	x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		// 	y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		// 	z = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		// 	vec3 offset = {x,y,z};
		// 	for (int j = 0; j < 25; j++)
		// 	{
		// std::cout<<"pre pos call"<<std::endl;
		// 		pos[i*25+j] += offset;					
		// std::cout<<"post pos call"<<std::endl;
		// 	}
		// }

		// pos[10] = {-100,0,0};
		// pos[20] = {0,-100,0};
		// pos[30] = {0,0,-100};

		// pos[100] = {100,0,0};
		// pos[120] = {0,100,0};
		// pos[130] = {0,0,100};

		// pos[10] = {0,11,0};
		// pos[100] = {-12,11,0};

		///

		//To see what group a ball belongs to
		//index into group_ids with ball index and 
		//the value is what group that ball belongs to.
		//Groups also acts as a visited list.
		//If group_id[index] = -1 then it hasnt been visited
		// gridIDs = new int[num_balls];
		findGroups();
		return;
	}

	~grid() 
	{
		// delete[] gridIDs;
 	}

	//Copy constructor
	grid(const grid& copyme)
	{
		numBalls = copyme.numBalls;
		gridIDs = copyme.gridIDs;
		pos = copyme.pos;
		gridSize = copyme.gridSize;
		rad = copyme.rad;
	}

	//Assignment operator
	grid& operator=(const grid& src)
	{
		if (this != &src)
		{
			numBalls = src.numBalls;
			gridIDs = src.gridIDs;
			pos = src.pos;
			gridSize = src.gridSize;
			rad = src.rad;
		}
		return *this;
	}

	void findGroups()
	{
  		// resetGroups();
		// vec3 *extrema;
		setMaxGridIndex();
		
		for (int i = 0; i < numBalls; i++)
		{	
			std::vector<int> iID;
			iID.push_back(pos[i][0]/gridSize);
			iID.push_back(pos[i][1]/gridSize);
			iID.push_back(pos[i][2]/gridSize);
			gridIDs[i] = iID;
		}
		printGroups();

  		return;
	}

	int findGroup(int ball_pos)
	{
		return -12;
	}

	//@brief finds the necessary size of the overall grid
	void setMaxGridIndex()
	{
		// static vec3 extrema[2];
		// extrema[0] = pos[0];
		// extrema[1] = pos[0];
		double max = -1;
		for (int i = 0; i < numBalls; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (abs(pos[i][j]) > max)
				{
					max = abs(pos[i][j]);	
				}
			}
		}
		maxGridIndex = ceil(max/gridSize);
	}

	void resetGroups()
	{
		for (int i = 0; i < numBalls; i++)
		{
			gridIDs[i] = std::vector<int> {INT_MAX,INT_MAX,INT_MAX};
		}
		// numGroups = -1;
		// delete[] groups;
		return;
	}
	
	void printGroups()
	{
		std::cout<<"================Groups================"<<std::endl;
		for (int i = 0; i < numBalls; i++)
		{
			std::cout<<"(ball; group): "<<i<<"; {"<<gridIDs[i][0]<<","<<gridIDs[i][1]<<","<<gridIDs[i][2]<<"}"<<std::endl;
		}
		std::cout<<"======================================"<<std::endl;
		return;
	}

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

 //  	void printVector(std::vector<int> vec)
 //  	{
	// 	std::cout<<"================Vector================"<<std::endl;

 //  		for (auto it = begin(vec); it != end(vec); ++it) 
	// 	{
	// 		std::cout<<*it<<", ";
	// 	}	
	// 	std::cout<<std::endl;
	// 	std::cout<<"======================================"<<std::endl;
	// 	return;
 //  	}


};