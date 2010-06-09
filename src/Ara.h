#ifndef ARA
#define ARA

#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <limits>
#include <queue>
#include <list>
#include <algorithm>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/astar_search.hpp>

#include <boost/graph/random.hpp>
#include <boost/random.hpp>
//#include <boost/foreach.hpp>

#include "PriorityQueue.h"
#include "Timer.h"

#define DEBUG
//#define DIST_EUCLIDIANA

using namespace boost;
using namespace std;

class Ara{
  typedef int Edge_Weight;
  typedef adjacency_list < listS, vecS, directedS,
		  no_property, property <edge_weight_t, Edge_Weight> > G;
  typedef graph_traits < G >::vertex_descriptor vertex_descriptor;
  typedef graph_traits < G >::edge_descriptor edge_descriptor;
  typedef graph_traits < G >::vertex_iterator vertex_iterator;
  typedef graph_traits < G >::out_edge_iterator out_edge_iterator;
  typedef graph_traits < G >::adjacency_iterator adjacency_iterator;

  typedef std::pair<int, int> Edge;

  enum Sets { OPEN, CLOSED, INCONS, NONE };


/********************************************************************
*        VARIABLES AND CONSTANTS                                    *
********************************************************************/
  const double PI;
  const double DG_TO_RAD;	
  const Edge_Weight INFINITO;
  int num_nodes, num_edges;
	
  Edge* edge_array;
  Edge_Weight* weights;
  double* x;
  double* y;

  Edge_Weight* f;
  Edge_Weight* g;
  Edge_Weight* v;
  int* backPointer;

  Sets* set;
  G* graph;
  PriorityQueue<Edge_Weight>* open;
  list<int> incons;

#ifdef DEBUG
  int estadosChecados;
  int iter;
  Timer lTimer, gTimer; 
#endif


    void load(string co_filename, string dist_filename);
	void initARA();
	void initialize();
	void initializeEdges();

	double key(int id, const int goal, const double &e);
	void update_open(const int goal, const double &e);
	void empty_closed();
	void move_incons_open();
	void greedyPath(const int &start, const int &goal);
	void computePath(const double &e, int start, int goal);
public:
	double h(int id, const int goal);

	Ara(string co_filename, string dist_filename);
	void publishSolution(int &iter, int et, int start, int goal, int max = 1);
    void araSpeedUpTest(int repeat = 1);

	void main(int start, int goal, double e0, double delta);

	Edge_Weight astar(int start, int goal, bool verbose = false);
	void dijkstra(int start,int goal);
};

#endif
