#include <iostream>
#include <vector>
// #include <fstream>
#include <unordered_map>
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
	// int maxGridIndex = -1;
	std::vector<std::vector<int>> gridIDs;
	std::unordered_map<std::string,std::vector<int>> IDToGrid;
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
		mapGroups();
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
		// setMaxGridIndex();
		
		for (int i = 0; i < numBalls; i++)
		{	
			std::vector<int> iID;
			iID.push_back(floor(pos[i][0]/gridSize));
			iID.push_back(floor(pos[i][1]/gridSize));
			iID.push_back(floor(pos[i][2]/gridSize));
			gridIDs[i] = iID;
		}
  		return;
	}

	void mapGroups()
	{
		for (int i = 0; i < numBalls; i++)
		{
			std::string key = getKey(gridIDs[i]);
			if (IDToGrid.find(key) == IDToGrid.end()) // key not present
			{
				std::vector<int> indices{i};
				IDToGrid[key] = indices;
				std::cout<<"key: val  " <<key<<": "<<IDToGrid[key][0]<<std::endl;
			}
			else//key present, add to vector
			{
				IDToGrid[key].push_back(i);
			}
		}
		printMap();
		return;
	}

	std::string getKey(std::vector<int> v)
	{
		return std::to_string(v[0]) + std::to_string(v[1]) + std::to_string(v[2]);
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
		// maxGridIndex = ceil(max/gridSize);
	}

	std::vector<int> getBalls(int ballIndex)
	{
		std::vector<std::vector<int>> neededGroups{gridIDs[ballIndex]};
		//check if ball is too close to a boarder
		double tolerance = radius*4;


		//check 4 verticies first

		//check x value
		double gridXmin = gridIDs[ballIndex][0]*gridSize; 
		double gridXmax = (gridIDs[ballIndex][0]+1)*gridSize; 
		if (abs(pos[ballIndex][0] - gridXmin) <= tolerance) // check if too close to min
		{
			//need to add group to neededGroups
			std::vector<int> tempID = gridIDs[ballIndex]
			tempID[0]--;
			neededGroups.push_back(tempID);
		}



		return IDToGrid[getKey(in)];
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