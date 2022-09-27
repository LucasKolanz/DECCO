#include <iostream>
#include <vector>
#include <fstream>
#include <queue>
#include <numeric>
#include "vec3.hpp"

//Facilitates a SINGLE group of balls
class group
{
public:
	double init_rad = -0.1; //radius of spheres
	int id = -1;
	vec3 com = {-1,-1,-1};
	int numBalls = -1;
	std::vector<vec3> pos;
	vec3 center = {-1,-1,-1};
	double radius = -1.0;
	std::vector<int> ball_indices;

	group() = default;



	group(double rad, int group_id, vec3 center_of_mass, 
		int num_balls, std::vector<vec3> positions)
	{
		init_rad = rad;
		id = group_id;
		com = center_of_mass;
		numBalls = num_balls;
		// pos = new vec3[numBalls];
		pos = positions;
		enclosingSphere();
	}

	void enclosingSphere();
	
	~group()
	{
		// delete[] pos;
		// delete[] com;
		// delete[] id;
	}
	
};

class graph
{
public:
	bool initialized = false;
	double rad = -1.0;
	int numVerts;
	int numGroups = -1;
	std::vector<int> *adj = nullptr;


	int *group_ids = nullptr;
	std::vector<group> groups;

	vec3 *pos = nullptr;
	// bool *visited = nullptr;
	// std::


	graph() = default;

	graph(int num_balls, double radius,vec3 *positions)
	{
		numVerts = num_balls;
		adj = new std::vector<int>[num_balls];
		rad = radius;
		pos = positions;

		///TAKE THIS OUT FOR REAL THING
		//This is just meant to make multiple groups out of a single group
		srand(0);
		for (int i = 0; i < 4; i++)
		{
			double x,y,z;
			x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			z = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
			vec3 offset = {x,y,z};
			for (int j = 0; j < 25; j++)
			{
				pos[i*25+j] += offset;					
			}
		}
		///

		graph_init();
		//To see what group a ball belongs to
		//index into group_ids with ball index and 
		//the value is what group that ball belongs to.
		//Groups also acts as a visited list.
		//If group_id[index] = -1 then it hasnt been visited
		group_ids = new int[num_balls];
		findGroups();

		initialized = true;
		return;
	}

	void graph_init()
	{
		for (int A = 1; A < numVerts; A++)
		{
			for (int B = 0; B < A; B++)
			{
				const double sumRaRb = 2*rad;
                const double dist = (pos[A] - pos[B]).norm();
                const double overlap = dist - sumRaRb;
                //TODO: check overlap conditions
                if (overlap <= 0)
                {
                	addEdge(A,B);
                }
                // else
                // {
                // 	std::cout<<"EDGE not ADDED: "<<overlap<<std::endl;
                // }
			}
		}
	}

	~graph() 
	{
		delete[] adj;
		delete[] group_ids;
		delete[] groups;
 	}

	//Copy constructor
	graph(const graph& copyme)
	{
		numVerts = copyme.numVerts;
		adj = copyme.adj;
		group_ids = copyme.group_ids;
		initialized = copyme.initialized;
	}

	//Assignment operator
	graph& operator=(const graph& src)
	{
		if (this != &src)
		{
			group_ids = src.group_ids;
			adj = src.adj;
			numVerts = src.numVerts;
			initialized = src.initialized;
		}
		return *this;
	}

	// int sumVisited()
	// {
	// 	int accum = 0;
	// 	for (int i = 0; i < numVerts; i++)
	// 	{
	// 		accum += visited[i];
	// 	}
	// 	return accum;
	// }

	int notVisited()
	{
		for (int i = 0; i < numVerts; i++)
		{
			if (group_ids[i] == -1)
			{
				return i;
			}
		}
		return -1;
	}

	void findGroups(int start=0)
	{
  		resetGroups();
		int curr = start;
		int group_id = 0;
  		do {
  			getWalk(curr, group_id);
  		// std::cout<<group_id<<std::endl;
  			group_id++;
  			curr = notVisited();
  			// printGroups();
  		}while(curr != -1);
  		// std::cout<<group_id<<std::endl;
  		numGroups = group_id;
  		groups = new group[numGroups];

  		// int lengths[numGroups] = {0};
  		// // printArray(lengths,numGroups);
  		// for (int i = 0; i < numVerts; i++)
  		// {
  		// 	lengths[group_ids[i]]++;
  		// }

  		// printArray(lengths,numGroups);
  		// printArray(group_ids,numVerts);


  		int gpid;
  		for (int groupnum = 0; groupnum < numGroups; groupnum++)
  		{
  			// std::vector<vec3> group_pos[lengths[groupnum]];
  			std::vector<vec3> group_pos;
  			vec3 center_of_mass = {0,0,0};
  			gpid = 0;
	  		for (int i = 0; i < numVerts; i++)
	  		{
	  			if (group_ids[i] == groupnum)
	  			{
	  				center_of_mass += pos[i];
	  				group_pos.push_back(pos[i]);
	  				gpid++;
	  			}
	  		}
	  		groups[groupnum] = group(rad,groupnum,center_of_mass, numVerts,group_pos);
	  	}

  		return;
	}


	

	void resetGroups()
	{
		for (int i = 0; i < numVerts; i++)
		{
			group_ids[i] = -1;
		}
		numGroups = -1;
		// delete[] groups;
		return;
	}

	void getWalk(int start,int group_id)
	{
		std::vector<vec3> positions;
		vec3 center_of_mass = {0,0,0};
		int group_size = 0;
  		std::queue<int> q;
		q.push(start); 
  		while (!q.empty())
  		{
  			int front = q.front();
  			q.pop();
  			if (adj[front].size() == 0)
  			{
  				group_ids[front] = group_id;
  			}
  			else
  			{
	  			for (auto it = begin(adj[front]); it != end(adj[front]); ++it) 
				{	

					if (group_ids[*it] == -1) 
					{
						q.push(*it);
						group_size++;
						positions.push_back(pos[*it])
						ball_indices.push_back(*it)
						// group_ids[*it] = group_id;

					}
				}
			}
		}
		group new_group(rad,group_id,center_of_mass,group_size,positions);
		groups.push_back(new_group);
	}


	void addEdge(int n1, int n2)
	{
		if (std::find(adj[n1].begin(), adj[n1].end(), n2) == adj[n1].end()) 
		{
			adj[n1].push_back(n2);
		}
		if (std::find(adj[n2].begin(), adj[n2].end(), n1) == adj[n2].end()) 
		{
			adj[n2].push_back(n1);
		}
		return;
	}

	void deleteEdge(int n1, int n2)
	{
		auto find_me = std::find(adj[n1].begin(), adj[n1].end(), n2);
		if (find_me != adj[n1].end()) 
		{
			adj[n1].erase(find_me);
		}
		find_me = std::find(adj[n2].begin(), adj[n2].end(), n1);
		if (find_me != adj[n2].end()) 
		{
			adj[n2].erase(find_me);
		}
		return;
	}
	
	void printGroups()
	{
		std::cout<<"================Groups================"<<std::endl;
		for (int i = 0; i < numVerts; i++)
		{
			std::cout<<"(ball; group): "<<i<<"; "<<group_ids[i]<<std::endl;
		}
		std::cout<<"======================================"<<std::endl;
		return;
	}

	template <typename T>
	void printArray(T arr[], int size)
	{
		std::cout<<"================Array================"<<std::endl;
		for (int i = 0; i < size; i++)
		{
			std::cout<<arr[i]<<", ";
		}
		std::cout<<std::endl;
		std::cout<<"====================================="<<std::endl;
		return;
	}

	void printGraph() 
	{
		std::cout<<"================Graph================"<<std::endl;
		if (initialized)
		{
	    	for (int d = 0; d < numVerts; ++d) 
	    	{
	    		std::cout << "\n Vertex "<< d << ":";
	    		for (auto x : adj[d])
	      		{
	      			std::cout << "-> " << x;
	  			}
	    		std::cout<<std::endl;
	  		}
	  	}
	  	else
	  	{
	  		std::cout<<"ERROR: graph not initalized"<<std::endl;
	  	}
		std::cout<<"====================================="<<std::endl;
  		return;
  	}

  	void printGraph(std::string file) 
	{
		if (initialized)
		{
			std::string to_file = "";
	    	for (int d = 0; d < numVerts; ++d) 
	    	{
	      		// std::cout<<"d: "<<d<<std::endl;
	      		// std::cout<<"what is this: "<< "Vertex " + std::to_string(d) + ":"<<std::endl;
	    		to_file += "\n Vertex " + std::to_string(d) + ":";
	    		for (auto x : adj[d])
	      		{
	      			to_file += "-> " + std::to_string(x);
	  			}
	    		to_file += '\n';
	  		}
			std::ofstream write;
			write.open(file,std::ofstream::out);
			write << to_file;
			write.close();
		}
	  	else
	  	{
	  		std::cout<<"ERROR: graph not initalized"<<std::endl;
	  	}
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

  	// void mergeGraphs(const graph src)

  	std::vector<int> nNearestNeighbors(int start, int numNeighbor)
  	{

  		//If there are no balls touching the start ball
  		//TODO: return group information for non-contact forces
  		//At the moment: return all balls besides start ball 
  		std::vector<int> NN;
  		if (adj[start].size() == 0)
  		{
  			NN.resize(numVerts);
  			std::iota (std::begin(NN),std::end(NN),0);
  			NN.erase(NN.begin() + start);
  			// std::cout<<start<<": "<<numVerts<<std::endl;
  			// printVector(NN);
  			// exit(0);
  		}
  		else
  		{
  			for (auto it = begin(adj[start]); it != end(adj[start]); ++it)
  			{
  				if (it != begin(adj[start]))
  				{
  					NN.push_back(*it);
  				}
  			}
  		// 	int n = 0;
	  	// 	std::queue<int> q;
	  	// 	q.push(start);
	  	// 	while (!q.empty())
	  	// 	{
	  	// 		int front = q.front();
	  	// 		q.pop();
	  	// 		for (auto it = begin(adj[front]); it != end(adj[front]); ++it) 
				// {	
				// 	if (*it != front and *it != start)
				// 	{
				// 		if (std::find(NN.begin(), NN.end(), *it) == NN.end()) 
				// 		{
				// 			q.push(*it);

				// 			NN.push_back(*it);
				// 			n++;
				// 			if (n >= numNeighbor)
				// 			{
				// 				return NN;
				// 			}
				// 		}
				// 	}
				// }
	  		// }
  		}
  		return NN;
  	}

};


//This function uses Ritter Bounding Sphere
//Ritter, Jack. "An efficient bounding sphere." Graphics gems 1 (1990): 301-303.
void group::enclosingSphere()
{
	double minx, miny, minz, maxx, maxy, maxz;
	int minxi, minyi, minzi, maxxi, maxyi, maxzi;
	minxi = -1;
	minyi = -1;
	minzi = -1;
	maxxi = -1;
	maxyi = -1;
	maxzi = -1;

	minx = 9999;
	miny = 9999;
	minz = 9999;
	maxx = -1;
	maxy = -1;
	maxz = -1;

	//find max and min x,y,z
	for (int i = 0; i < numBalls; i++)
	{
		if (pos[i][0] < minx)
		{
			minx = pos[i][0];
			minxi = i;
		}
		else if (pos[i][0] > maxx)
		{
			maxx = pos[i][0];
			maxxi = i;
		}

		if (pos[i][1] < miny)
		{
			miny = pos[i][1];
			minyi = i;
		}
		else if (pos[i][1] > maxy)
		{
			maxy = pos[i][1];
			maxyi = i;
		}

		if (pos[i][2] < minz)
		{
			minz = pos[i][2];
			minzi = i;
		}
		else if (pos[i][2] > maxz)
		{
			maxz = pos[i][2];
			maxzi = i;
		}
	}

	//which two of these six points are farthest apart?
	int indices[6] = {maxxi, minxi, maxyi, minyi, maxzi, minzi};
	double max_dist = -1;
	double dist;
	int max_indi = -1;
	int max_indj = -1;

	for (int i = 1; i < 6; i++)
	{
		for (int j = 0; j < i; j++)
		{
			dist = (pos[i]-pos[j]).norm();
			if (dist > max_dist)
			{
				max_dist = dist;
				max_indi = i;
				max_indj = j;
			}
		}
	}
	//distance between farthest points is the initial diameter
	radius = dist/2;
	center = (pos[max_indi] + pos[max_indj]) / 2;
	
	//Check if points are inside circle and update circle if not
	double old_to_p_sq,old_to_p,radius_sq,old_to_new;
	vec3 dxyz;
	radius_sq = radius*radius;
	for (int i = 0; i < numBalls; i++)
	{
		dxyz = (pos[i] - center);
		old_to_p_sq = dxyz.normsquared();
		if (old_to_p_sq > radius_sq)
		{
			old_to_p = sqrt(old_to_p_sq);
			radius = (radius + old_to_p)/2.0;
			radius_sq = radius * radius;
			old_to_new = old_to_p - radius;
			center = (radius*center + old_to_new*pos[i])/old_to_p;
		}
	}
	radius += init_rad;
	return;

}