#include <iostream>
#include "timing.hpp"
#include <math.h> 

void func1()
{
	int n;
	double c;
	n = 10000000;
	c = 0.0;
	for (int i = 0; i < n; ++i)
	{
		c += i;
		if (i > 100)
		{
			c /= (i*0.5);
			c = sqrt(c);
		}
	}
}


int main()
{
	timey t;
	t.start_event("func1");
	t.start_event("func2");
	t.start_event("func1");
	func1();
	t.end_event("func1");
	t.print_events();
	t.start_event("func1");
	func1();
	t.end_event("func1");
	t.end_event("func2");
	t.print_events();
	t.save_events("timing.txt");

	return 0;
} 