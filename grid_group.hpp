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
		double tolerancesq = rad*4 * rad*4;


		//check 8 verticies first
		// vec3 v1,v2,v3,v4,v5,v6,v7,v8;
		// v1.push_back(gridIDs[ballIndex][0]*gridSize);
		// v1.push_back(gridIDs[ballIndex][1]*gridSize);
		// v1.push_back(gridIDs[ballIndex][2]*gridSize);
		vec3 v1(gridIDs[ballIndex][0]*gridSize, gridIDs[ballIndex][1]*gridSize, gridIDs[ballIndex][2]*gridSize); 

		vec3 v2 = v1;
		v2[0] += gridSize;

		vec3 v3 = v1;
		v3[1] += gridSize;

		vec3 v4 = v1;
		v4[2] += gridSize;

		vec3 v5 = v1;
		v5[0] += gridSize;
		v5[1] += gridSize;

		vec3 v6 = v1;
		v6[1] += gridSize;
		v6[2] += gridSize;

		vec3 v7 = v1;
		v7[0] += gridSize;
		v7[2] += gridSize;

		vec3 v8 = v1 * gridSize;

		///check vertex1
		if (v1.dot(pos[ballIndex]) <= tolerancesq)
		{
			//need 7 surrounding groups
			std::vector<int> g1 = gridIDs[ballIndex];
			g1[0]--;
			neededGroups.push_back(g1);

			std::vector<int> g2 = gridIDs[ballIndex];
			g2[0]--;
			g2[1]--;
			neededGroups.push_back(g2);

			std::vector<int> g3 = gridIDs[ballIndex];
			g3[1]--;
			neededGroups.push_back(g3);

			std::vector<int> g4 = gridIDs[ballIndex];
			g4[2]--;
			neededGroups.push_back(g4);

			std::vector<int> g5 = gridIDs[ballIndex];
			g5[0]--;
			g5[2]--;
			neededGroups.push_back(g5);

			std::vector<int> g6 = gridIDs[ballIndex];
			g6[0]--;
			g6[1]--;
			g6[2]--;
			neededGroups.push_back(g6);

			std::vector<int> g7 = gridIDs[ballIndex];
			g7[1]--;
			g7[2]--;
			neededGroups.push_back(g7);
		} 

		///check vertex 2
		if (v2.dot(pos[ballIndex]) <= tolerancesq)
		{
			//need 7 surrounding groups
			std::vector<int> g1 = gridIDs[ballIndex];
			g1[1]--;
			neededGroups.push_back(g1);

			std::vector<int> g2 = gridIDs[ballIndex];
			g2[1]--;
			g2[2]--;
			neededGroups.push_back(g2);

			std::vector<int> g3 = gridIDs[ballIndex];
			g3[2]--;
			neededGroups.push_back(g3);

			std::vector<int> g4 = gridIDs[ballIndex];
			g4[0]++;
			neededGroups.push_back(g4);

			std::vector<int> g5 = gridIDs[ballIndex];
			g5[1]--;
			g5[0]++;
			neededGroups.push_back(g5);

			std::vector<int> g6 = gridIDs[ballIndex];
			g6[1]--;
			g6[2]--;
			g6[0]++;
			neededGroups.push_back(g6);

			std::vector<int> g7 = gridIDs[ballIndex];
			g7[0]++;
			g7[2]--;
			neededGroups.push_back(g7);
		}

		///check vertex 3
		if (v3.dot(pos[ballIndex]) <= tolerancesq)
		{
			//need 7 surrounding groups
			std::vector<int> g1 = gridIDs[ballIndex];
			g1[2]--;
			neededGroups.push_back(g1);

			std::vector<int> g2 = gridIDs[ballIndex];
			g2[0]--;
			g2[2]--;
			neededGroups.push_back(g2);

			std::vector<int> g3 = gridIDs[ballIndex];
			g3[0]--;
			neededGroups.push_back(g3);

			std::vector<int> g4 = gridIDs[ballIndex];
			g4[1]++;
			neededGroups.push_back(g4);

			std::vector<int> g5 = gridIDs[ballIndex];
			g5[2]--;
			g5[1]++;
			neededGroups.push_back(g5);

			std::vector<int> g6 = gridIDs[ballIndex];
			g6[0]--;
			g6[2]--;
			g6[1]++;
			neededGroups.push_back(g6);

			std::vector<int> g7 = gridIDs[ballIndex];
			g7[1]++;
			g7[0]--;
			neededGroups.push_back(g7);
		}

		///check vertex 4
		if (v4.dot(pos[ballIndex]) <= tolerancesq)
		{
			//need 7 surrounding groups
			std::vector<int> g1 = gridIDs[ballIndex];
			g1[0]--;
			neededGroups.push_back(g1);

			std::vector<int> g2 = gridIDs[ballIndex];
			g2[0]--;
			g2[1]--;
			neededGroups.push_back(g2);

			std::vector<int> g3 = gridIDs[ballIndex];
			g3[1]--;
			neededGroups.push_back(g3);

			std::vector<int> g4 = gridIDs[ballIndex];
			g4[2]++;
			neededGroups.push_back(g4);

			std::vector<int> g5 = gridIDs[ballIndex];
			g5[0]--;
			g5[2]++;
			neededGroups.push_back(g5);

			std::vector<int> g6 = gridIDs[ballIndex];
			g6[0]--;
			g6[1]--;
			g6[2]++;
			neededGroups.push_back(g6);

			std::vector<int> g7 = gridIDs[ballIndex];
			g7[2]++;
			g7[1]--;
			neededGroups.push_back(g7);
		}

		///check vertex 5
		if (v5.dot(pos[ballIndex]) <= tolerancesq)
		{
			//need 7 surrounding groups
			std::vector<int> g1 = gridIDs[ballIndex];
			g1[0]++;
			neededGroups.push_back(g1);

			std::vector<int> g2 = gridIDs[ballIndex];
			g2[0]++;
			g2[1]++;
			neededGroups.push_back(g2);

			std::vector<int> g3 = gridIDs[ballIndex];
			g3[1]++;
			neededGroups.push_back(g3);

			std::vector<int> g4 = gridIDs[ballIndex];
			g4[2]--;
			neededGroups.push_back(g4);

			std::vector<int> g5 = gridIDs[ballIndex];
			g5[2]--;
			g5[0]++;
			neededGroups.push_back(g5);

			std::vector<int> g6 = gridIDs[ballIndex];
			g6[2]--;
			g6[0]++;
			g6[1]++;
			neededGroups.push_back(g6);

			std::vector<int> g7 = gridIDs[ballIndex];
			g7[1]++;
			g7[2]--;
			neededGroups.push_back(g7);
		}

		///check vertex 6
		if (v6.dot(pos[ballIndex]) <= tolerancesq)
		{
			//need 7 surrounding groups
			std::vector<int> g1 = gridIDs[ballIndex];
			g1[1]++;
			neededGroups.push_back(g1);

			std::vector<int> g2 = gridIDs[ballIndex];
			g2[1]++;
			g2[0]--;
			neededGroups.push_back(g2);

			std::vector<int> g3 = gridIDs[ballIndex];
			g3[0]--;
			neededGroups.push_back(g3);

			std::vector<int> g4 = gridIDs[ballIndex];
			g4[0]++;
			neededGroups.push_back(g4);

			std::vector<int> g5 = gridIDs[ballIndex];
			g5[1]++;
			g5[2]++;
			neededGroups.push_back(g5);

			std::vector<int> g6 = gridIDs[ballIndex];
			g6[1]++;
			g6[0]--;
			g6[2]++;
			neededGroups.push_back(g6);

			std::vector<int> g7 = gridIDs[ballIndex];
			g7[2]++;
			g7[0]--;
			neededGroups.push_back(g7);
		}
		///check vertex 7
		if (v7.dot(pos[ballIndex]) <= tolerancesq)
		{
			//need 7 surrounding groups
			std::vector<int> g1 = gridIDs[ballIndex];
			g1[]--;
			neededGroups.push_back(g1);

			std::vector<int> g2 = gridIDs[ballIndex];
			g2[]--;
			g2[]--;
			neededGroups.push_back(g2);

			std::vector<int> g3 = gridIDs[ballIndex];
			g3[]--;
			neededGroups.push_back(g3);

			std::vector<int> g4 = gridIDs[ballIndex];
			g4[]++;
			neededGroups.push_back(g4);

			std::vector<int> g5 = gridIDs[ballIndex];
			g5[]--;
			g5[]++;
			neededGroups.push_back(g5);

			std::vector<int> g6 = gridIDs[ballIndex];
			g6[]--;
			g6[]--;
			g6[]++;
			neededGroups.push_back(g6);

			std::vector<int> g7 = gridIDs[ballIndex];
			g7[]++;
			g7[]--;
			neededGroups.push_back(g7);
		}
		///check vertex 8
		if (v8.dot(pos[ballIndex]) <= tolerancesq)
		{
			//need 7 surrounding groups
			std::vector<int> g1 = gridIDs[ballIndex];
			g1[]--;
			neededGroups.push_back(g1);

			std::vector<int> g2 = gridIDs[ballIndex];
			g2[]--;
			g2[]--;
			neededGroups.push_back(g2);

			std::vector<int> g3 = gridIDs[ballIndex];
			g3[]--;
			neededGroups.push_back(g3);

			std::vector<int> g4 = gridIDs[ballIndex];
			g4[]++;
			neededGroups.push_back(g4);

			std::vector<int> g5 = gridIDs[ballIndex];
			g5[]--;
			g5[]++;
			neededGroups.push_back(g5);

			std::vector<int> g6 = gridIDs[ballIndex];
			g6[]--;
			g6[]--;
			g6[]++;
			neededGroups.push_back(g6);

			std::vector<int> g7 = gridIDs[ballIndex];
			g7[]++;
			g7[]--;
			neededGroups.push_back(g7);
		}	

		//check x value
		// double gridXmin = gridIDs[ballIndex][0]*gridSize; 
		// double gridXmax = (gridIDs[ballIndex][0]+1)*gridSize; 
		// if (abs(pos[ballIndex][0] - gridXmin) <= tolerance) // check if too close to min
		// {
		// 	//need to add group to neededGroups
		// 	std::vector<int> tempID = gridIDs[ballIndex]
		// 	tempID[0]--;
		// 	neededGroups.push_back(tempID);
		// }



		return IDToGrid[getKey(gridIDs[ballIndex])];
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