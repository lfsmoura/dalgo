
#include <iostream>
#include <fstream>
#include <limits>
#include <queue>
#include <list>
#include <algorithm>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/astar_search.hpp>

#include <boost/graph/random.hpp>
#include <boost/random.hpp>
//#include <boost/foreach.hpp>

#include "Timer.h"
#include "ADstar.h"

using namespace std;

//int main(){
//
//	return 0;
//}

G* load(string co_filename, string dist_filename){
	using namespace boost;
	G* graph;
	int num_nodes,
		num_edges;
	Edge* edge_array;
	int* weights;

	double x, y;

	cout << "carregando arquivo de distancias" << endl;
	ifstream dist_file(dist_filename.c_str());

	char c[100];
	int i = 0;
	int a,b;
	while( dist_file.good() ){
		dist_file >> c;
		if( c[0] == 'p'){
			dist_file >> c >> num_nodes >> num_edges;
			
			edge_array = new Edge[num_edges]();;
			weights = new int[num_edges]();

		}
		else if( c[0] == 'a'){
			dist_file >> edge_array[i].first >> edge_array[i].second >> weights[i];			
			i++;

			//grafo nao-direcionado
			dist_file >> c >> c >> c;
		}
		else{
			dist_file.getline(c,100);
		}
	}
	dist_file.close();
	
	graph = new G(edge_array, edge_array + num_edges, weights, num_nodes);

	cout << "carregando arquivo de coordenadas" << endl;
	ifstream co_file(co_filename.c_str());

	while( co_file.good() ){
		co_file >> c;
		if( c[0] == 'v'){
			co_file >> i >> a >> b;
			// transforma em graus, depois em radianos
			x = ((double) a) / 1000000.0; 
			x *= DG_TO_RAD;
			y = ((double) b) / 1000000.0;
			y *= DG_TO_RAD;

			put(coord,i,Pos(x,y));
		}
		else{
			co_file.getline(c,100);
		}
	}
	co_file.close();
}