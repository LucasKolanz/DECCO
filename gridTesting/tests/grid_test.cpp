#include "../grid.hpp"
#include "../vec3.hpp"
#include "../Utils.hpp"


int main()
{
	int length = 100;
	vec3 *pos = new vec3[length];

	srand(100);
	for (int i = 0; i < length; i++)
	{
		pos[i] = {rand_between(-0.999,0.999),rand_between(-0.999,0.999),rand_between(-0.999,0.999)};
	}
	double max_val = 0.9999;
	pos[length-1] = {max_val,max_val,max_val};

	grid g(pos,length,0.1);

	//Test maxmin
	double mm = g.min_max(); 
	if (mm == max_val)
	{
		std::cout<<"maxmin passes test."<<std::endl;
	}
	else
	{
		std::cout<<"maxmin FAILS test, with value of "<<mm<<std::endl;
	}

	//Test grid_length
	int gl = g.get_num_cells();
	if (gl == std::pow(20,3))
	{
		std::cout<<"grid_length passes test."<<std::endl;
	}
	else
	{
		std::cout<<"grid_length FAILS test, with value of "<<gl<<std::endl;
	}

	//Test one_to_three and three_to_one
	int idx = 9; //should be the last grid space in x and the most negative y and z spaces
	int idxx = g.one_to_three(idx,0); 
	int idxy = g.one_to_three(idx,1); 
	int idxz = g.one_to_three(idx,2); 
	if (idxx == 9 && idxy == 0 && idxz == 0)
	{
		std::cout<<"one_to_three (1) passes test."<<std::endl;
	}
	else
	{
		std::cout<<"one_to_three FAILS test, with values of "<<idxx<<", "<<idxy<<", "<<idxz<<std::endl;
	}

	idx = 21; //should be the first grid space in y and smallest in x and z
	idxx = g.one_to_three(idx,0); 
	idxy = g.one_to_three(idx,1); 
	idxz = g.one_to_three(idx,2); 
	if (idxx == 0 && idxy == 1 && idxz == 0)
	{
		std::cout<<"one_to_three passes test."<<std::endl;
	}
	else
	{
		std::cout<<"one_to_three (2) FAILS test, with values of "<<idxx<<", "<<idxy<<", "<<idxz<<std::endl;
	}
}