#include <iostream>
#include <chrono>
#include <unordered_map>
#include <fstream>


struct time_event
{
	bool running = false;
	std::string name;
	std::chrono::high_resolution_clock::time_point start;
	std::chrono::high_resolution_clock::time_point end;
	double duration = 0.0;
	time_event(std::string n)
	{
		name = n;
	}
	time_event() = default;
	
};

class timey
{
public:

	std::unordered_map<std::string,time_event> events;
	timey() = default;

	std::chrono::high_resolution_clock::time_point now()
	{
		return std::chrono::high_resolution_clock::now();
	}

	void start_event(std::string name)
	{
		if (events.find(name) == events.end()) // key not present
		{
			events[name] = time_event(name);
			events[name].start = now();
			events[name].running = true;
		}
		else
		{
			events[name].start = now();
			events[name].running = true;
		}
	}

	void end_event(std::string name)
	{
		if (events.find(name) == events.end()) // key not present
		{
			std::cout<<"Key not present"<<std::endl;
		}	
		else
		{
			events[name].end = now();
			std::chrono::duration<double, std::milli> ms_double = events[name].end - events[name].start;
			events[name].duration += std::chrono::duration<double>(ms_double).count();
			events[name].running = false;
		}
	}

	void print_events()
	{
		for (auto i : events)
		{
			std::cout<<"Event '"<<i.first<<"' took "<<i.second.duration<<" ms\n";
		}
		std::cout<<std::endl;
	}

	void save_events(std::string save_name)
	{
		std::ofstream myfile;
		myfile.open(save_name, std::ios::app);
		for (auto i : events)
		{
			myfile<<"Event '"<<i.first<<"' took "<<i.second.duration<<" ms\n";
		}
		myfile<<std::endl;

	}

};