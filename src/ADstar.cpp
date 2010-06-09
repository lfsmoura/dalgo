
#include <iostream>
#include <fstream>
#include <limits>
#include <queue>
#include <list>
#include <algorithm>


#include <boost/graph/random.hpp>
#include <boost/random.hpp>
//#include <boost/foreach.hpp>

#include "Timer.h"
#include "ADstar.h"

#define test(nome,var) cout << nome << " = " << var << endl

const double INFINITO = std::numeric_limits<double>::max();

ADstar::ADstar(G* graph){
	numNodes_ = num_vertices(*graph);
	//x = new double[numNodes_]();
	//y = new double[numNodes_]();

	f_ = new double[numNodes_]();
	g_ = new double[numNodes_]();
	v_ = new double[numNodes_]();
	backPointer_ = new int[numNodes_]();
	set_ = new Sets[numNodes_]();

	open_ = new PriorityQueue(numNodes_, f_);
}