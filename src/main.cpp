#include <iostream>

//#include "drawGraph.h"
#include "Ara.h"
#include "Ad.h"

int main(int argc, char** argv)
{
	//g::argc = argc;
	//g::argv = argv;
	/*    TESTE ARA*
	Ara a("dimacs//USA-road-d.BAY.co","dimacs//USA-road-d.BAY.gr");
	Timer timer;
	timer.start();
	a.araSpeedUpTest(101);
	timer.stop();
	cout << "tempo total:" << timer << endl;
	*/
	//Ad a("../data/dimacs/USA-road-d.NY.co", "../data/dimacs/USA-road-d.NY.gr");
	
	
	if(argc < 4) {
	  std::cout << "wrong options : use ARA co_file gr_file #iterations" << std::endl; 
	  return -2;
	}
	Ad a(argv[1], argv[2]);
	//Ad a("dimacs//grid_c.gr", "dimacs//grid_d.gr");
	
#ifdef DRAW_GRAPH
	drawGraph(&a);
#else
	a.test(atoi(argv[3]));
#endif
	std::cout << "fim"<< std::endl;

	return 0;
}