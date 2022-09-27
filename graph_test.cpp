#include <iostream>
#include <cmath>
#include <fstream>
#include "graph_group.hpp"
#include "vec3.hpp"

class wrapper
{
public:
	graph *g = nullptr;
	vec3 *pos = nullptr;
	int num_particles = -1;
	double init_rad = 1e-5;
	
	wrapper() = default;

	// wrapper(int num)
	// {
	// 	num_particles = num;
	// 	g = new graph(num);	
	// }

	

	wrapper(int num,std::string file)
	{
		num_particles = num;
		pos = new vec3[num];
		parseSimData(getLastLine(file));
		std::cout<<"TEST: "<<pos[0]<<std::endl;
		g = new graph(num, init_rad, pos);
		std::cout<<"TEST1: "<<g->pos[0]<<std::endl;
		// g -> graph_init(init_rad);
	}


	void writeFile(std::string file, vec3 center, double radius)
	{
		std::ofstream myfile;
		myfile.open(file, std::ios::out);
		for (int i = 0; i < num_particles; i++)
		{
			myfile<<pos[i]<<"\n";
		}
		myfile<<"center: "<<center<<"\n";
		myfile<<"radius: "<<radius<<"\n";
	}

	//This function uses Ritter Bounding Sphere
	//Ritter, Jack. "An efficient bounding sphere." Graphics gems 1 (1990): 301-303.
	// void enclosingSphere(vec3 &center, double &radius)
	// {
	// 	double minx, miny, minz, maxx, maxy, maxz;
	// 	int minxi, minyi, minzi, maxxi, maxyi, maxzi;
	// 	minxi = -1;
	// 	minyi = -1;
	// 	minzi = -1;
	// 	maxxi = -1;
	// 	maxyi = -1;
	// 	maxzi = -1;

	// 	minx = 9999;
	// 	miny = 9999;
	// 	minz = 9999;
	// 	maxx = -1;
	// 	maxy = -1;
	// 	maxz = -1;

	// 	//find max and min x,y,z
	// 	for (int i = 0; i < num_particles; i++)
	// 	{
	// 		if (pos[i][0] < minx)
	// 		{
	// 			minx = pos[i][0];
	// 			minxi = i;
	// 		}
	// 		else if (pos[i][0] > maxx)
	// 		{
	// 			maxx = pos[i][0];
	// 			maxxi = i;
	// 		}

	// 		if (pos[i][1] < miny)
	// 		{
	// 			miny = pos[i][1];
	// 			minyi = i;
	// 		}
	// 		else if (pos[i][1] > maxy)
	// 		{
	// 			maxy = pos[i][1];
	// 			maxyi = i;
	// 		}

	// 		if (pos[i][2] < minz)
	// 		{
	// 			minz = pos[i][2];
	// 			minzi = i;
	// 		}
	// 		else if (pos[i][2] > maxz)
	// 		{
	// 			maxz = pos[i][2];
	// 			maxzi = i;
	// 		}
	// 	}

	// 	//which two of these six points are farthest apart?
	// 	int indices[6] = {maxxi, minxi, maxyi, minyi, maxzi, minzi};
	// 	double max_dist = -1;
	// 	double dist;
	// 	int max_indi = -1;
	// 	int max_indj = -1;

	// 	for (int i = 1; i < 6; i++)
	// 	{
	// 		for (int j = 0; j < i; j++)
	// 		{
	// 			dist = (pos[i]-pos[j]).norm();
	// 			if (dist > max_dist)
	// 			{
	// 				max_dist = dist;
	// 				max_indi = i;
	// 				max_indj = j;
	// 			}
	// 		}
	// 	}
	// 	//distance between farthest points is the initial diameter
	// 	radius = dist/2;
	// 	center = (pos[max_indi] + pos[max_indj]) / 2;
	// 	// std::cout<<"init rad = "<<radius<<std::endl;
	// 	// std::cout<<"init center = {"<<center<<"}"<<std::endl;

	// 	// std::cout<<"min: "<<minx<<", "<<miny<<", "<<minz<<std::endl;
	// 	// std::cout<<"mini: "<<minxi<<", "<<minyi<<", "<<minzi<<std::endl;
	// 	// std::cout<<"max: "<<maxx<<", "<<maxy<<", "<<maxz<<std::endl;
	// 	// std::cout<<"maxi: "<<maxxi<<", "<<maxyi<<", "<<maxzi<<std::endl;
		
	// 	//Check if points are inside circle and update circle if not
	// 	double old_to_p_sq,old_to_p,radius_sq,old_to_new;
	// 	vec3 dxyz;
	// 	radius_sq = radius*radius;
	// 	for (int i = 0; i < num_particles; i++)
	// 	{
	// 		dxyz = (pos[i] - center);
	// 		old_to_p_sq = dxyz.normsquared();
	// 		if (old_to_p_sq > radius_sq)
	// 		{
	// 			old_to_p = sqrt(old_to_p_sq);
	// 			radius = (radius + old_to_p)/2.0;
	// 			radius_sq = radius * radius;
	// 			old_to_new = old_to_p - radius;
	// 			center = (radius*center + old_to_new*pos[i])/old_to_p;
	// 		}
	// 	}
	// 	radius += init_rad;
	// 	return;

	// }

	void parseSimData(std::string line)
    {
        std::string lineElement;
        // Get number of balls in file
        int count = 269;
        int properties = 11;
        
        // int count = std::count(line.begin(), line.end(), ',') / properties + 1;
        // allocate_group(count);
        std::stringstream chosenLine(line);  // This is the last line of the read file, containing all data
                                             // for all balls at last time step
        // Get position and angular velocity data:
        for (int A = 0; A < num_particles; A++) {
            for (int i = 0; i < 3; i++)  // Position
            {
                std::getline(chosenLine, lineElement, ',');
                pos[A][i] = std::stod(lineElement);
            }
            for (int i = 0; i < 3; i++)  // Angular Velocity
            {
                std::getline(chosenLine, lineElement, ',');
                // w[A][i] = std::stod(lineElement);
            }
            std::getline(chosenLine, lineElement, ',');  // Angular velocity magnitude skipped
            for (int i = 0; i < 3; i++)                  // velocity
            {
                std::getline(chosenLine, lineElement, ',');
                // vel[A][i] = std::stod(lineElement);
            }
            for (int i = 0; i < properties - 10; i++)  // We used 10 elements. This skips the rest.
            {
                std::getline(chosenLine, lineElement, ',');
            }
        }
    }

	/// Get last line of previous simData by filename.
    [[nodiscard]] static std::string getLastLine(const std::string& filename)
    {
        std::string simDataFilepath = filename;
        if (auto simDataStream = std::ifstream(simDataFilepath, std::ifstream::in)) {
            std::cerr << "\nParsing last line of data.\n";

            simDataStream.seekg(-1, std::ios_base::end);  // go to one spot before the EOF

            bool keepLooping = true;
            while (keepLooping) {
                char ch = ' ';
                simDataStream.get(ch);  // Get current byte's data

                if (static_cast<int>(simDataStream.tellg()) <=
                    1) {                     // If the data was at or before the 0th byte
                    simDataStream.seekg(0);  // The first line is the last line
                    keepLooping = false;     // So stop there
                } else if (ch == '\n') {     // If the data was a newline
                    keepLooping = false;     // Stop at the current position.
                } else {                     // If the data was neither a newline nor at the 0 byte
                    simDataStream.seekg(-2, std::ios_base::cur);  // Move to the front of that data, then to
                                                                  // the front of the data before it
                }
            }
            std::string line;
            std::getline(simDataStream, line);  // Read the current line
            return line;
        } else {
            std::cerr << "Could not open simData file: " << simDataFilepath << "... Existing program."
                      << '\n';
            exit(EXIT_FAILURE);
        }
        return "-1";
    }

	void addEdge()
	{
		g -> addEdge(0,1);
		g -> addEdge(1,2);
		g -> addEdge(1,3);
		g -> addEdge(1,4);
		g -> addEdge(1,9);
		g -> addEdge(3,4);
		g -> addEdge(3,5);
		g -> addEdge(3,6);
		g -> addEdge(3,7);
		g -> addEdge(5,8);
		g -> addEdge(5,9);

		g -> addEdge(10,11);
		// g -> addEdge(11,12);
		g -> addEdge(11,13);
		g -> addEdge(11,14);
		g -> addEdge(11,19);
		g -> addEdge(13,14);
		g -> addEdge(13,15);
		g -> addEdge(13,16);
		g -> addEdge(13,17);
		g -> addEdge(15,18);
		g -> addEdge(15,19);

		g -> findGroups();
	}

	~wrapper()
	{
		delete g;
	}
	
};

// void test_NN()
// {
// 	graph g(10);

// 	g.addEdge(0,1);
// 	g.addEdge(1,2);
// 	g.addEdge(1,3);
// 	g.addEdge(1,4);
// 	g.addEdge(1,9);
// 	g.addEdge(3,4);
// 	g.addEdge(3,5);
// 	g.addEdge(3,6);
// 	g.addEdge(3,7);
// 	g.addEdge(5,8);
// 	g.addEdge(5,9);

// 	g.printVector(g.nNearestNeighbors(3,4));

// }

// void test_sparse_NN()
// {
// 	graph g(10);

// 	// g.addEdge(0,1);
// 	// g.addEdge(1,2);
// 	// g.addEdge(1,3);
// 	// g.addEdge(1,4);
// 	// g.addEdge(1,9);
// 	// g.addEdge(3,4);
// 	// g.addEdge(3,5);
// 	// g.addEdge(3,6);
// 	// g.addEdge(3,7);
// 	// g.addEdge(5,8);
// 	// g.addEdge(5,9);

// 	g.printVector(g.nNearestNeighbors(3,4));

// }

// void test_delete()
// {
// 	graph g(10);

// 	g.addEdge(0,1);
// 	g.addEdge(1,2);
// 	g.addEdge(1,3);
// 	g.addEdge(1,4);
// 	g.addEdge(1,9);
// 	g.addEdge(3,4);
// 	g.addEdge(3,5);
// 	g.addEdge(3,6);
// 	g.addEdge(3,7);
// 	g.addEdge(5,8);
// 	g.addEdge(5,9);


// 	// std::vector<int> test = g.nNearestNeighbors(1,5);

// 	// g.printVector(test);
// 	g.printGraph();


// 	g.deleteEdge(1,9);
// 	std::cout<<"=========================================="<<std::endl;

// 	g.printGraph();

// }

// void test_groups()
// {
// 	graph *g = nullptr;

// 	g = new graph(20);

// 	g -> addEdge(0,1);
// 	g -> addEdge(1,2);
// 	g -> addEdge(1,3);
// 	g -> addEdge(1,4);
// 	g -> addEdge(1,9);
// 	g -> addEdge(3,4);
// 	g -> addEdge(3,5);
// 	g -> addEdge(3,6);
// 	g -> addEdge(3,7);
// 	g -> addEdge(5,8);
// 	g -> addEdge(5,9);

// 	g -> addEdge(10,11);
// 	// g -> addEdge(11,12);
// 	g -> addEdge(11,13);
// 	g -> addEdge(11,14);
// 	g -> addEdge(11,19);
// 	g -> addEdge(13,14);
// 	g -> addEdge(13,15);
// 	g -> addEdge(13,16);
// 	g -> addEdge(13,17);
// 	g -> addEdge(15,18);
// 	g -> addEdge(15,19);

// 	g -> findGroups();
// 	g -> printGraph();
// 	g -> printGroups();
// }

// void test_wrapper_methods()
// {
// 	std::cout<<"Start wrapper test"<<std::endl;
// 	wrapper w(20);
// 	w.addEdge();

	
// 	w.g -> printGraph();
// 	w.g -> printGroups();
// 	std::cout<<"END TEST"<<std::endl;
// }

void test_sbs()
{
	std::string file = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/large_aggregate/N_1000/198_2_R2e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_simData.csv";
	wrapper w(200, file);

	// for (int i = 0; i < 20; i++)
	// {
	// 	std::cout<< w.pos[i] << std::endl;
	// }

	vec3 center{0,0,0};
	double radius = 0; 

	// vec3 *center = &c;
	// double *radius = &r;

	// w.g -> enclosingSphere(center, radius);
	w.writeFile("enclosingSphere.txt", center, radius);
}

int main(int argc, char const *argv[])
{
	// test_wrapper_methods();
	// test_sparse_NN();
	// test_NN();
	test_sbs();

	return 0;
}

