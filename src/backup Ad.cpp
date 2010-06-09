//#include <boost/config.hpp>
//#include <iostream>
//#include <fstream>
//#include <limits>
//#include <queue>
//#include <list>
//#include <algorithm>
//
//#include <boost/graph/graph_traits.hpp>
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/dijkstra_shortest_paths.hpp>
//#include <boost/graph/astar_search.hpp>
//
//#include <boost/graph/random.hpp>
//#include <boost/random.hpp>
////#include <boost/foreach.hpp>
//
//#include "PriorityQueue.h"
//#include "Timer.h"
//
////#define DEBUG
////#define DEBUG_HEAP
////#define CASO_SIMPLES
////#define SHOW_OPEN
////#define SHOW_PATH
//#define ONLY_CHANGE_PATH
//#define MUDANCAS_PASSO 200
//
//namespace ara_star
//{
//	using namespace boost;
//	using namespace std;
//
///********************************************************************
//*        TYPE DEFINITIONS                                           *
//********************************************************************/
//
//  #define test(nome,var) cout << nome << " = " << var << endl
//  #define PI 3.14159265
//  const double DG_TO_RAD = PI / 180.0;	
//
//  typedef long Edge_Weight;
//
//  const Edge_Weight INFINITO = std::numeric_limits<Edge_Weight>::max() / 2;
//
//  typedef adjacency_list < listS, vecS, undirectedS,
//		  no_property, property <edge_weight_t, Edge_Weight> > G;
//
//  typedef graph_traits < G >::vertex_descriptor vertex_descriptor;
//  typedef graph_traits < G >::edge_descriptor edge_descriptor;
//  typedef graph_traits < G >::vertex_iterator vertex_iterator;
//  typedef graph_traits < G >::out_edge_iterator out_edge_iterator;
//  typedef graph_traits < G >::adjacency_iterator adjacency_iterator;
//
//  typedef std::pair<int, int> Edge;
//
///********************************************************************
//*        FUNCTIONS                                                  *
//********************************************************************/
//
//Edge_Weight astar(int start, int goal, bool verbose = false);
//
/////********************************************************************
////*        VARIABLES                                                  *
////********************************************************************/
//  int num_nodes, num_edges;
//	
//  Edge* edge_array;
//  Edge_Weight* weights;
//  double* x;
//  double* y;
//
//  Edge_Weight* f;
//  Edge_Weight* g;
//  Edge_Weight* v;
//  int* backPointer;
//
//  enum Sets { OPEN, CLOSED, INCONS, NONE };
//  Sets* set;
//
//  G* graph;
//
//  PriorityQueue<Edge_Weight>* open;
// 
//  list<int> incons;
//
//  int estadosChecados;
//
//  Timer lTimer, gTimer; 
//
///********************************************************************
//*        FUNCTIONS                                                  *
//********************************************************************/
//
//void initialize(){
//	edge_array = new Edge[num_edges]();;
//	weights = new Edge_Weight[num_edges]();
//	x = new double[num_nodes+1]();
//	y = new double[num_nodes+1]();
//
//	f = new Edge_Weight[num_nodes+1]();
//	g = new Edge_Weight[num_nodes+1]();
//	v = new Edge_Weight[num_nodes+1]();
//	backPointer = new int[num_nodes+1]();
//	set = new Sets[num_nodes+1]();
//	for(int i = 0; i <= num_nodes; i++)
//	{
//		f[i] = g[i] = v[i] = INFINITO;
//		backPointer[i] = NULL;
//		set[i] = NONE;
//	}
//	open = new PriorityQueue<Edge_Weight>(num_nodes+1, f);
//}
//void initializeEdges(){
//	graph = new G(edge_array, edge_array + num_edges, weights, num_nodes);
//}
//void load(string co_filename, string dist_filename){
//
//	cout << "carregando arquivo de distancias" << endl;
//	ifstream dist_file(dist_filename.c_str());
//
//	char c[100];
//	int i = 0;
//	int a,b;
//	int u,v;
//	vector< vector < int >> edges;
//
//	while( dist_file.good() ){
//		dist_file >> c;
//		if( c[0] == 'p'){
//			dist_file >> c >> num_nodes >> num_edges;
//			initialize();
//			edges.resize(num_edges);
//		}
//		else if( c[0] == 'a'){
//			bool adiciona = true;
//			dist_file >> u >> v >> weights[i];
//			
//			// checa se a aresta ja foi adicionada
//			// se sim, nao a adiciona
//			for(unsigned int k = 0; k < edges[u].size(); k++)
//			{
//				if( edges[u][k] == v )
//					adiciona = false;
//			}
//
//			if( adiciona ){
//				edge_array[i].first = u;
//				edge_array[i].second = v; 
//				i++;
//				edges[u].push_back(v);
//				edges[v].push_back(u);
//			}
//
//			//grafo nao-direcionado
//			dist_file >> c >> c >> c;
//		}
//		else{
//			dist_file.getline(c,100);
//		}
//	}
//	dist_file.close();
//	initializeEdges();
//	cout << "carregando arquivo de coordenadas" << endl;
//	ifstream co_file(co_filename.c_str());
//
//	while( co_file.good() ){
//		co_file >> c;
//		if( c[0] == 'v'){
//			co_file >> i >> a >> b;
//			// transforma em graus, depois em radianos
//#ifndef CASO_SIMPLES
//			x[i] = ((double) a) / 1000000.0; 
//			x[i] *= DG_TO_RAD;
//			y[i] = ((double) b) / 1000000.0;
//			y[i] *= DG_TO_RAD;
//#else
//			x[i] = a;
//			y[i] = b;
//#endif
//		}
//		else{
//			co_file.getline(c,100);
//		}
//	}
//	co_file.close();
//}
//
//double h(int id, const int goal){
//#ifndef CASO_SIMPLES
//	double d = acos((sin(y[id])*sin(y[goal]))+(cos(y[id])*cos(y[goal])*cos(x[id] - x[goal])));
//	return 63787000.0 * d; 
//	//return 6378700.0 * d;
//#else
//	return sqrt( pow(x[id] - x[goal],2) + pow(y[id] - y[goal],2));
//#endif
//}
///********************************************************************
//*        ANYTIME DYNAMIC A* (AD*)                                *
//********************************************************************/
//// no ad* eu não usei uma lista INCONS, ao inves disso eu criei uma array que diz o conjunto
//// que cada nodo pertence
//namespace ad
//{
//	// variaveis
//	vector<int> path;
//	double eg;
//	int inCount = 0;
//
//	int times = 0;
//	int aNodesVisited = 0;
//	float	aTime = 0;
//	int adNodesVisited = 0;
//	float   adTime = 0;
//	//
//	bool inconsistent(int s){
//		return v[s] != g[s];
//	}
//	bool overconsistent(int s){      // confirmar isso
//		return v[s] > g[s];
//	}
//	bool underconsistent(int s){
//		return v[s] < g[s];
//	}
//
//	void printOpen()
//	{
//#ifdef SHOW_OPEN
//		cout << "open: ";
//		for(int i = 0; i <= num_nodes; i++)
//		{
//			if( set[i] == OPEN)
//			{
//				cout << i << " (f=" << f[i] << ", v+h="<<v[i]+h(i,24) <<", g="<< g[i]+h(i,24)<< ") . ";
//			}
//		}
//		cout << endl;
//#endif
//	}
//
//	Edge_Weight key(int id, const int goal,double e = 1.0){
//		Edge_Weight anterior = f[id];
//
//#ifdef DEBUG_HEAP
//		if( v[id] >= g[id] ){
//			cout << id << " g " << f[id] << " -> "; 
//			f[id] = max(g[id] + (Edge_Weight)( e * h(id,goal) ),g[id]);
//			cout << f[id] << endl;
//		}else{
//			cout << id << " v " << f[id] << " -> "; 
//			f[id] = max(v[id] + (Edge_Weight) h(id,goal),v[id]);
//			cout << f[id] << endl;
//		}
//#else
//		if( v[id] >= g[id] ){
//			f[id] = max(g[id] + (Edge_Weight)( e * h(id,goal) ),g[id]);
//		}else{
//			f[id] = max(v[id] + (Edge_Weight) h(id,goal),v[id]);
//		}
//#endif
//
//
//		if(set[id] == OPEN)
//			if( f[id] < anterior)
//				open->decrease(id);
//			else
//				open->increase(id);
//		
//		return f[id];
//	}
//	void move_incons_open(){
//		for(int i = 0; i <= num_nodes; i++)
//		{
//			if( set[i] == INCONS)
//			{
//				open->push(i);
//				set[i] = OPEN;
//				}
//		}
//	}
//	void update_open(const int goal, const double &e){
//		for(int i = 0; i <= num_nodes; i++)
//		{
//			if(	set[i] == OPEN){
//				key(i,goal,e);
//			} 
//		}
//	}
//	void clearClosed(){
//		for(int i = 0; i <= num_nodes; i++)
//		{
//			if(set[i] == CLOSED)
//				set[i] = NONE;
//		} 
//	}
//	void updateSetMembership(int s, int goal, double e){
//		key(s, goal, e);
//		if( inconsistent(s) )
//		{
//			if( set[s] != CLOSED && set[s] != INCONS)
//			{
//				if( set[s] != OPEN )
//					open->push(s);
//				set[s] = OPEN;
//			}
//			else{
//				
//				if(set[s] == CLOSED)
//				{
//					//inCount++;
//					//cout << "%%%%%%%%%%%%%%%INCONS: "<< s << " f[s] = " << f[s] << endl;
//					open->push(s);
//					set[s] = OPEN;
//				}else
//					set[s] = INCONS;
//			}
//		}
//		else
//		{
//			if( set[s] == OPEN)
//			{
//				set[s] = NONE;
//				open->remove(s);
//			}
//			else if( set[s] == INCONS )
//			{
//				set[s] = NONE;
//			}
//		}
//	}
//	int argmin(int u, Edge_Weight* funcs, int goal){
//		
//		int argMin = -15;
//		Edge_Weight minCost = INFINITO;
//		Edge_Weight cost = 0;
//
//		out_edge_iterator ki,kend;
//		for(tie(ki,kend) = out_edges(u, *graph); ki != kend; ki++)
//		{
//			vertex_iterator p1 = target(*ki, *graph);		
//			int neighbour = *p1;
//
//			cost = get(edge_weight, *graph, *ki);
//#ifdef DEBUG
//			cout << "  a " << *ki << ") peso = " << funcs[neighbour] + cost << " " << funcs[neighbour] << " + " << cost << " bp = " << backPointer[neighbour] << " v[neighbour] = " << funcs[neighbour] << " - v[u] = " << funcs[u] << " minCost " << minCost << " minArg " << argMin << endl;
//			//cout << "aresta: " << *ki << " h nodo = " << h(neighbour,goal) << " h u = " << h(u,goal) << " x u =" << x[neighbour] << " y u = " << y[neighbour] << endl;
//#endif			
//			if( ((funcs[neighbour] + cost) < minCost )   && ( backPointer[neighbour] != u ) )
//			{
//				argMin = neighbour;
//				minCost = funcs[neighbour] + cost;
//			}
//		}
//		
//		if( argMin < 1 ){
//			//exit(-12);
//			argMin = backPointer[u];
//			//cout << " usado o negativo " << endl;
//		}
//#ifdef DEBUG
//		cout << "  - o argmin eh " << argMin << endl;
//#endif
//		return argMin;
//	}
//	Edge_Weight c(int u, int v){
//		edge_descriptor e;
//		bool edgeExists;
//
//		tie(e,edgeExists) = edge(v, u, *graph);
//		if(edgeExists)
//			return get(edge_weight, *graph, e);
//		else
//			throw;
//	}
//
//	void computePath(int start, int goal, double e)
//	{
//		int min_id;
//			if(!open->empty())
//				min_id = open->key_top();
//			else
//				return;
//#ifdef DEBUG
//		std::cout << "Compute Path ===> min id=" << min_id << endl;
//		printOpen();
//#endif
//
//		while( ( key(goal, goal, e) >= key(min_id, goal, e) ) ||
//			   ( underconsistent(goal) ) ||
//			   ( underconsistent(min_id) )
//			)
//		{
//			int s = min_id;
//			open->remove(min_id);
//
//			if( ! underconsistent(s) )
//			{
//#ifdef DEBUG			
//			std::cout << "  - - - expandindo estado ( overconsistent ) " << s << endl;
//#endif
//				v[s] = g[s];
//				set[s] = CLOSED;
//				updateSetMembership(s,goal,e);
//
//				estadosChecados++;
//
//				out_edge_iterator vi,vend;
//				for(tie(vi,vend) = out_edges(s, *graph); vi != vend; vi++)
//				{
//					int custo = get(edge_weight, *graph, *vi);
//					vertex_iterator p1 = target(*vi, *graph);		
//					int neighbour = *p1;
//
//					if( g[neighbour] > g[s] + custo)
//					{
//						backPointer[neighbour] = s;
//						g[neighbour] = g[s] + custo;
//						updateSetMembership(neighbour,goal,e);
//					}
//				}
//				
//			}
//			else // propagating underconsistency
//			{
//#ifdef DEBUG			
//			std::cout << "  - - - expandindo estado ( underconsistent ) " << s << endl;
//#endif
//				v[s] = INFINITO;
//				set[s] = NONE;
//				updateSetMembership(s,goal,e);
//
//				adjacency_iterator vi,vend;
//				for(tie(vi,vend) = adjacent_vertices(s, *graph); vi != vend; vi++)
//				{	
//					int neighbour = *vi;
//					if( backPointer[neighbour] == s)
//					{
//						backPointer[neighbour] = argmin(neighbour,v,goal);
//						g[neighbour] = v[backPointer[neighbour]] + c(backPointer[neighbour],neighbour);
//
//						updateSetMembership(neighbour,goal,e);
//					}
//				}
//			}
//			printOpen();
//			if(!open->empty()){
//				min_id = open->key_top();
//#ifdef DEBUG
//				cout << " o proximo a ser expandido eh " << min_id << endl;
//#endif
//			}else
//				return;
//		} // wend
//		v[goal] = 0;
//	}
//	typedef pair<int,Edge> Change;
//	vector<Change> getChanges(){
//		vector<Change> changes;
//		static mt19937 gen((boost::uint32_t)time(0));
//
//#ifdef ONLY_CHANGE_PATH
//
//		int change = path[(gen() % (path.size())) ];
//		//int change = path[gen() % 10 ]; // limita para o final as mudanças
//		if( change == 1 )
//			change++;
//		Edge_Weight custo = c(change, backPointer[change]);
//		changes.push_back(make_pair(( gen() % 150) + custo ,make_pair(change, backPointer[change])));
//
//#else
//		vector<pair<int,int> > jaFeitas;
//
//		int numChanges = (gen() % (MUDANCAS_PASSO)) + 1;
//
//		for(int i = 0; i < numChanges; ++i)
//		{
//
//			int v = (gen() % num_nodes) + 1;
//
//			int d = out_degree(v, *graph);
//			
//			int nu = gen() % d;
//			int u;
//			//int* o = new int[d];
//			int j = 0;
//			out_edge_iterator vi,vend;
//			Edge_Weight custo = 0;
//			for(tie(vi,vend) = out_edges(v, *graph); vi != vend; vi++){
//				vertex_iterator p1 = target(*vi, *graph);	
//				if( j == nu ){
//					u = *p1;
//					custo = get(edge_weight, *graph, *vi);
//					break;
//				}
//				j++;
//			}
//			
//			// nao muda mais de uma vez a mesma aresta
//			jaFeitas.push_back(make_pair<int,int>(u,v));
//			jaFeitas.push_back(make_pair<int,int>(v,u));
//			vector< pair<int,int> >::iterator it = find(jaFeitas.begin(),jaFeitas.end(),make_pair<int,int>(u,v));
//			bool naoAchou = it != jaFeitas.end();
//
//			if(naoAchou)
//				changes.push_back(make_pair(( gen() % 400)+ custo,make_pair(v, u)));
//		}
//
//#endif
//		return changes;
//	}
//	bool changeGraph(Change change){
//		
//		int newCost = change.first,
//			u = change.second.first,
//			v = change.second.second;
//#ifdef DEBUG
//		cout << u << " -> " << v << " (" << newCost << ")" << endl;
//#endif
//		edge_descriptor e;
//		bool edgeExists;
//
//		tie(e,edgeExists) = edge(u, v, *graph);
//
//		if(edgeExists){
//			put(edge_weight, *graph,  e, newCost);
//			return true;
//		}
//		else
//			return false;
//	}
//	void fixInitialization(int vtarget, int start, int goal, double e){
//		if( vtarget != start )
//		{
//			backPointer[vtarget] = argmin(vtarget,v,goal);
//
//			if(backPointer[vtarget] == NULL) // quando o nodo esta fora da area de interesse, descarta
//				return;
//
//			g[vtarget] = max(v[backPointer[vtarget]] + c(vtarget, backPointer[vtarget]), v[backPointer[vtarget]]);
//
//			updateSetMembership(vtarget,goal,e);
//
//		}
//	}
//
//
//	Edge_Weight pathCost(const int &start, const int &goal)
//	{
//		Edge_Weight custo = 0;
//
//		path.clear();
//		path.push_back(goal);
//
//#ifdef SHOW_PATH 
//		cout << endl << "ARA*: " << goal;
//#endif
//		for( int v = goal; v != start; v = backPointer[v])
//		{
//			edge_descriptor e;
//			bool edgeExists;
//#ifdef SHOW_PATH
//			cout << " -> " << backPointer[v];
//#endif
//
//			// tirar a linha abaixo
//			if( v == backPointer[backPointer[backPointer[v]]] || v == backPointer[backPointer[v]] ){
//				cout << "erro!";
//				cout << v << endl;
//				cout << backPointer[v];
//				cout << endl;
//			}
//
//			tie(e,edgeExists) = edge(v, backPointer[v], *graph);
//			if(edgeExists){
//				custo += get(edge_weight, *graph, e);
//
//				if( find(path.begin(), path.end(), backPointer[v]) != path.end() ){
//					cout << v << " <- " << backPointer[v] << endl; 				
//					return 0;
//				}
//				path.push_back(backPointer[v]);
//			}else{
//				//cout << endl << "999 " << e << endl;
//				return -1;
//			}
//			//cout << v << " e " << backPointer[v]  << " ->";
//		}
//#ifdef SHOW_PATH
//		cout << endl;
//#endif
//		return custo;
//	}
//bool erro = false;
//	void publishSolution(int &iter, double e, int start, int goal)
//	{
//
//		lTimer.stop();
//		static double anterior = 0;
//		Edge_Weight ara, ast;
//
//		eg = e;
//		ara = pathCost(start,goal);
//
//		//ast = ( iter % 100 == 0) ? astar(start,goal,false) : 0; 
//		ast = astar(start,goal,false);
//		erro =  ara != ast;
//		
//		adTime += lTimer.show();
//		adNodesVisited += estadosChecados;
//		times += 1;
//
//		if(times == 10){
//			cout  << endl <<  e << ";" << iter << ";" << (adNodesVisited / 10)  << ";" << (adTime)/10 << ";" << ara << ";" ; // << inCount;
//			astar(start,goal,true);
//		
//			times = 0;
//			adNodesVisited = 0;
//			adTime = 0;
//		}
//
//		if(erro && (e <= 1)){
//			cout  << e <</* ";" << iter <<*/ ";" << estadosChecados << ";" << lTimer << ";" << ara << "==" << ast << endl;
//			cout << "erro 5"<<endl;
//
//			//move_incons_open();
//			//update_open(goal,e);
//			//clearClosed();
//			//computePath(start,goal,e);
//			//publishSolution(iter,e,start,goal);
//
//			//cout << "timex: " << timex << endl;
//			//exit(0);
//		}
//		estadosChecados = 0;
//		anterior = ara;
//		lTimer.start();
//
//
//
//		inCount = 0;
//
//	}
//	void init_adstar(int start, int goal, double e)
//	{
//		for(int i = 0; i <= num_nodes; i++)
//		{
//			f[i] = g[i] = v[i] = INFINITO;
//			backPointer[i] = 0;
//			set[i] = NONE;
//		}
//
//		g[goal] = v[goal] = v[start] = INFINITO;
//		backPointer[start] = NULL;
//		g[start] = 0;
//
//		open->clear();
//		for(int i = 0; i <= num_nodes; i++)
//		{
//			set[i] = NONE;
//		} 
//				
//		key(start,goal,e);
//		open->push(start);
//		set[start] = OPEN;
//	}
//vector<Change> changes;
//	void adstar(int start, int goal, double e0, double delta){
//
//		bool neverBefore = true;
//		double e = e0;
//		init_adstar(start,goal,e);
//
//		lTimer.start();	
//
//		for(int iter = 0; iter < 500000; iter++)
//		{
//			computePath(start,goal,e);
//			publishSolution(iter,e,start,goal);
//
//			if ( e <= 1){
//				if(neverBefore) { gTimer.stop(); neverBefore = false; }
//				
//				//
//				//if(erro){
//				//	cout << "mundancas:" << endl;
//				//	for(vector<Change>::iterator c = changes.begin(); c != changes.end(); c++){
//				//		cout << c->second.first << ", " << c->second.second << " -> " << c->first << endl;
//				//	}
//
//				//	cout << "erros fim" << endl;
//				//}
//				//
//
//				//vector<Change> changes = getChanges();
//				changes.clear();
//
//				changes = getChanges();
//
//				//cout << " *** mudancas(" << changes.size() << ") : ";
//
//				for(vector<Change>::iterator c = changes.begin(); c != changes.end(); c++){
//					if( ! changeGraph(*c) )
//						exit(0);
//				}
//
//				for(vector<Change>::iterator c = changes.begin(); c != changes.end(); c++){
//					//cout << c->first << " [ " << c->second.first << "->" << c->second.second << " ] " << endl;
//
//					int u = c->second.first;
//					int p = c->second.second;
//
//					if( v[u] > v[p]){					
//						fixInitialization(u,start,goal,e);
//						//fixInitialization(p,start,goal,e);
//					}
//					else{
//						fixInitialization(p,start,goal,e);
//						//fixInitialization(u,start,goal,e);
//					}
//				}
//
//				//if( iter % 100 == 0)
//				//{
//				//	e = e0;
//				//	init_adstar(start,goal,e);
//				//}
//			}
//			else
//			{
//				e -= delta;
//			}
//			move_incons_open();
//			update_open(goal,e);
//			clearClosed();
//		}
//	}
//}
///********************************************************************
//*        DIJKSTRA                                                   *
//********************************************************************/
//void dijkstra(int start, int goal){
//	std::vector<vertex_descriptor> p(num_vertices(*graph));
//	std::vector<int> d(num_vertices(*graph)+1);
//	vertex_descriptor s = vertex(start,*graph);
//	
//	Timer timer;
//	timer.start();
//	dijkstra_shortest_paths(*graph, s, predecessor_map(&p[0]).distance_map(&d[0]));
//	timer.stop();
//
//	//int custo = 0;
//	//for( int v = goal; v != start;)
//	//{
//	//	out_edge_iterator vi,vend;
//	//	int minG = numeric_limits<int>::max();
//	//	int argMinG;
//	//	int minCost;
//
//	//	for(tie(vi,vend) = out_edges(v, *graph); vi != vend; vi++)
//	//	{		
//
//	//		if((d[(target(*vi,*graph))] + get(edge_weight, *graph, *vi) < minG))
//	//		{
//	//			minCost = get(edge_weight, *graph, *vi);
//	//			minG = d[(target(*vi,*graph))] + minCost;
//	//			argMinG = (target(*vi,*graph));	
//	//		}
//	//	}
//	//	custo += minCost;
//	//	v = argMinG;
//	//}
//	cout << endl <<"dijkstra: d = " << d[goal] << endl;
//	cout << "tempo total " << timer << endl;
//}
////********************************************************************
////*        A*                                                         *
////********************************************************************/
//
//int aNodesVisisted;
//
//template <class Graph, class CostType>
//class distance_heuristic : public astar_heuristic<Graph, CostType>
//{
//public:
//  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
//  distance_heuristic(Vertex goal)
//    : m_goal(goal) {}
//  CostType operator()(Vertex u)
//  {
//    return h(u,m_goal);
//  }
//private:
//  int m_goal;
//};
//
//struct found_goal {}; // exception for termination
//
//// visitor that terminates when we find the goal
//template <class Vertex>
//class astar_goal_visitor : public boost::default_astar_visitor
//{
//public:
//  astar_goal_visitor(Vertex goal) : m_goal(goal) {}
//  template <class Graph>
//  void examine_vertex(Vertex u, Graph& g) {
//	aNodesVisisted++;
//	if(u == m_goal)
//      throw found_goal();
//  }
//private:
//  Vertex m_goal;
//};
//
//
//Edge_Weight astar(int start, int goal, bool verbose){
//	std::vector<vertex_descriptor> p(num_vertices(*graph));
//	std::vector<Edge_Weight> d(num_vertices(*graph)+1);
//	vertex_descriptor s = vertex(start,*graph),
//					  g = vertex(goal,*graph);
//	
//	Timer timer;
//	timer.start();
//
//	aNodesVisisted = 0;
//
//	try{
//	astar_search
//      (*graph, s, distance_heuristic<G, Edge_Weight> (goal),
//       predecessor_map(&p[0]).distance_map(&d[0]).
//       visitor(astar_goal_visitor<vertex_descriptor>(g)));
//	  } catch(found_goal) {
//#ifdef SHOW_PATH		
//		cout << "A*  : " ;
//
//		bool detec = false;
//		for(int t = goal, i = 0; t != start; t = p[t])
//		{
//			if( ! detec &&  (t != ad::path[i++])){
//				cout << endl << "erro => " << t << " != " << ad::path[i-1]<< endl;
//				detec = true;
//			}
//			cout << t << " -> ";
//		}
//#endif  
//
//		//if( ad::eg <= 1 ){
//		//	for(int i = 0; i <= num_nodes; i++)
//		//	{
//		//		if( h(start,i) > d[i] )
//		//			exit(-9);
//		//	}
//		//}
//		timer.stop();		
//		
//		ad::aTime += timer.show();
//		ad::aNodesVisited += aNodesVisisted;
//
//		if(verbose && ad::times == 10 )
//		{
//			cout << ";" << (int) (ad::aNodesVisited / 10) << ";" << (ad::aTime / 10.0); // << ";" << d[goal];
//			ad::aNodesVisited = 0;
//			ad::aTime = 0;
//		}
//	}
//	return d[goal];
//}
//}
//
/////********************************************************************
////*        MAIN                                                       *
////********************************************************************/
//int main()
//{
//#ifndef CASO_SIMPLES
//	//ara_star::load("dimacs//grid_c.gr","dimacs//grid_d.gr");
//
//	ara_star::load("dimacs//USA-road-d.BAY.co","dimacs//USA-road-d.BAY.gr");
//
//	//ara_star::gTimer.start();
//	//ara_star::ara(1,264345,2.0,0.05);
//	//ara_star::gTimer.stop();
//	//std::cout << "ARA* tempo total: " << ara_star::gTimer << std::endl;
//
//	ara_star::gTimer.start();
//	ara_star::ad::adstar(1,264345,4.8,0.2);
//	std::cout << "AD* tempo total: " << ara_star::gTimer << std::endl;
//	//
//	////ara_star::gTimer.start();
//	////ara_star::weightedAStar(1,264345);
//	////ara_star::gTimer.stop();
//	////std::cout << "ARA* tempo total: " << ara_star::gTimer << std::endl;
//	
//	ara_star::dijkstra(1,264345);
//	
//	std::cout << " astar "<< 
//	ara_star::astar(1,264345,true);
//
//#else
//	ara_star::load("dimacs//grid_c.gr","dimacs//grid_d.gr");
//	ara_star::gTimer.start();
//	ara_star::ad::adstar(1,24,4.8,0.2);
//	std::cout << "AD* tempo total: " << ara_star::gTimer << std::endl;
//	//
//	//
//	ara_star::dijkstra(1,24);
//	//
//	std::cout << " astar "<< 
//	ara_star::astar(1,24,true);
//#endif
//	std::cin.get();
//	std::cin.get();
//  return 0;
//}
