#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include "PriorityQueue.h"

using namespace boost;

struct coord_t {
	typedef vertex_property_tag kind;
};
struct Pos{
	int x;
	int y;
	Pos() : x(0), y(0) { }
	Pos(int x1, int y1) : x(x1), y(y1) { }
};

typedef adjacency_list < listS, vecS, undirectedS,
		  property <coord_t, Pos >, property <edge_weight_t, int> > G;

typedef graph_traits < G >::vertex_descriptor vertex_descriptor;
typedef graph_traits < G >::edge_descriptor edge_descriptor;
typedef graph_traits < G >::vertex_iterator vertex_iterator;
typedef graph_traits < G >::out_edge_iterator adjacency_iterator;

property_map<G, coord_t>::type coord;

typedef std::pair<int, int> Edge;

#define PI 3.14159265
const double DG_TO_RAD = PI / 180.0;	

class ADstar {
  enum Sets { OPEN, CLOSED, INCONS, NONE };
  
  int numNodes_,
	  numEdges_;

  int estadosChecados;

  double* x;
  double* y;

  double* f_;
  double* g_;
  double* v_;
  int* backPointer_;
  
  Sets* set_;

  G* graph_;

  PriorityQueue* open_;

public:
	ADstar(G* graph);
};