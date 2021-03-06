#ifndef AD
#define AD

#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <limits>
#include <set>
#include <stack>
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

//#define DEBUG
//#define DEBUG_HEAP
//#define DIST_EUCLIDIANA
//#define SHOW_OPEN
//#define SHOW_PATH
//#define ONLY_CHANGE_PATH
#define ONLY_CHANGE_REGION
//#define MUDANCAS_PASSO 200
//#define DRAW_GRAPH


using namespace boost;
using namespace std;

class Ad{
public:
  typedef long Edge_Weight;
  typedef adjacency_list < listS, vecS, directedS,
		  no_property, property <edge_weight_t, Edge_Weight> > G;
  typedef graph_traits < G >::vertex_descriptor vertex_descriptor;
  typedef graph_traits < G >::edge_descriptor edge_descriptor;
  typedef graph_traits < G >::vertex_iterator vertex_iterator;
  typedef graph_traits < G >::out_edge_iterator out_edge_iterator;
  typedef graph_traits < G >::adjacency_iterator adjacency_iterator;

  typedef std::pair<int, int> Edge;

  enum Sets { OPEN, CLOSED, INCONS, NONE };
  /* GRAY  = not visited
   * BLACK = ???
   * RED   = visited once 
   * YELLOW = visited more than once
   */
  enum Color{ GRAY, BLACK, RED, YELLOW };
  typedef pair<int,Edge> Change;

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

	double* drawx;
	double* drawy;
	Color* color, *colorastar;

	Edge_Weight* f;
	Edge_Weight* g;
	Edge_Weight* v;
	int* backPointer;
	vector<Change> changes;

	Sets* setAd;

	G* graph;

	PriorityQueue<Edge_Weight>* open;

	int estadosChecados;

	Timer lTimer, gTimer; 

    void load(string co_filename, string dist_filename);
	void init_adstar(int start, int goal, double e);
	void initialize();
	void initializeEdges();

	bool inconsistent(int s);
	bool overconsistent(int s);
	bool underconsistent(int s);
	void printOpen();

	Edge_Weight key(int id, const int goal, double e = 1.0);
	void move_incons_open();
	void update_open(const int goal, const double &e);
	void clearClosed();

	void updateSetMembership(int s, int goal, double e);
	int argmin(int u, Edge_Weight* funcs, int goal);

	Edge_Weight c(int u, int v);

	bool changeGraph(Change change);
	vector<Change> getChanges();
	vector<Change> oldGetChanges();
	void fixInitialization(int vtarget, int start, int goal, double e);
	Edge_Weight pathCost(const int &start, const int &goal);

	double h(int id, const int goal);
	Ad(string co_filename, string dist_filename);

	//void publishSolution(int &iter, int et, int start, int goal, int max = 1);
    //void araSpeedUpTest(int repeat = 1);
	void computePath(int start, int goal, double e);
	void main(int start, int goal, double e0, double delta, int maxIter);
	void publishSolution(int &iter, double e, int start, int goal);

	void test(int repeat);
	void graphicTest();

	Edge_Weight astar(int start, int goal, bool verbose = false);
	Ad::Edge_Weight dijkstra(int start,int goal);
};
#ifdef DRAW_GRAPH
#include "drawGraph.h"
#endif

#endif
