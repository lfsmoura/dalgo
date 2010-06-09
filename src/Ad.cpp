#include "Ad.h"
#include <iomanip>


///********************************************************************
//*        ANYTIME DYNAMIC A* (AD*)                                *
//********************************************************************/

// variaveis
namespace Debug{
vector<int> path;
double eg;
int inCount = 0;

int aNodesVisited = 0;
int times = 0;
int adOverNodesVisited = 0,
	adUnderNodesVisited;

int nvIncA, nvDecA, nvIncOverAd, nvDecOverAd, nvIncUnderAd, nvDecUnderAd;

double adTime = 0;
double aTime;

double adTimeInc,adTimeDec,aTimeInc,aTimeDec;
double dijTime;
double adijt;

unsigned long mudancas_passo;
double porcentagem;
double mp[] = { 0.01, 0.025, 0.05, 0.1, 0.2, 0.5 };
bool only_path_vertices;
stack<Ad::Change> changes;
int maxIter, totIter, iter;
int regionSize;
bool firstPubl;
}

Ad* ad;
//

bool Ad::inconsistent(int s){
	return v[s] != g[s];
}
bool Ad::overconsistent(int s){      // confirmar isso
	return v[s] > g[s];
}
bool Ad::underconsistent(int s){
	return v[s] < g[s];
}

void Ad::printOpen()
{
#ifdef SHOW_OPEN
	cout << "open: ";
	for(int i = 0; i <= num_nodes; i++)
	{
		if( setAd[i] == OPEN)
		{
			cout << i << " (f=" << f[i] << ", v+h="<<v[i]+h(i,24) <<", g="<< g[i]+h(i,24)<< ") . ";
		}
	}
	cout << endl;
#endif
}

Ad::Edge_Weight Ad::key(int id, const int goal,double e){
	Ad::Edge_Weight anterior = f[id];

#ifdef DEBUG_HEAP
		if( v[id] >= g[id] ){
			cout << id << " g " << f[id] << " -> "; 
			f[id] = max(g[id] + (Edge_Weight)( e * h(id,goal) ),g[id]);
			cout << f[id] << endl;
		}else{
			cout << id << " v " << f[id] << " -> "; 
			f[id] = max(v[id] + (Edge_Weight) h(id,goal),v[id]);
			cout << f[id] << endl;
		}
#else
		if( v[id] >= g[id] ){
			f[id] = max(g[id] + (Edge_Weight)( e * h(id,goal) ),g[id]);
		}else{
			f[id] = max(v[id] + (Edge_Weight) h(id,goal),v[id]);
		}
#endif
		if(setAd[id] == OPEN)
			if( f[id] < anterior)
				open->decrease(id);
			else
				open->increase(id);
		
		return f[id];
	}
void Ad::move_incons_open(){
		for(int i = 0; i <= num_nodes; i++)
		{
			if( setAd[i] == INCONS)
			{
				open->push(i);
				setAd[i] = OPEN;
				}
		}
	}
void Ad::update_open(const int goal, const double &e){
		for(int i = 0; i <= num_nodes; i++)
		{
			if(	setAd[i] == OPEN){
				key(i,goal,e);
			} 
		}
	}
void Ad::clearClosed(){
		for(int i = 0; i <= num_nodes; i++)
		{
			if(setAd[i] == CLOSED)
				setAd[i] = NONE;
		} 
	}
void Ad::updateSetMembership(int s, int goal, double e){
		key(s, goal, e);
		if( inconsistent(s) )
		{
			if( setAd[s] != CLOSED && setAd[s] != INCONS)
			{
				if( setAd[s] != OPEN )
					open->push(s);
				setAd[s] = OPEN;
			}
			else{
				
				//if(setAd[s] == CLOSED)
				//{
					open->push(s);
					setAd[s] = OPEN;
				//}else
				//	setAd[s] = INCONS;
			}
		}
		else
		{
			if( setAd[s] == OPEN)
			{
				setAd[s] = NONE;
				open->remove(s);
			}
			else if( setAd[s] == INCONS )
			{
				setAd[s] = NONE;
			}
		}
	}
int Ad::argmin(int u, Edge_Weight* funcs, int goal){
		
		int argMin = -15;
		Edge_Weight minCost = INFINITO;
		Edge_Weight cost = 0;

		out_edge_iterator ki,kend;
		for(tie(ki,kend) = out_edges(u, *graph); ki != kend; ki++)
		{
			vertex_iterator p1 = target(*ki, *graph);		
			int neighbour = *p1;

			cost = get(edge_weight, *graph, *ki);
#ifdef DEBUG
			cout << "  a " << *ki << ") peso = " << funcs[neighbour] + cost << " " << funcs[neighbour] << " + " << cost << " bp = " << backPointer[neighbour] << " v[neighbour] = " << funcs[neighbour] << " - v[u] = " << funcs[u] << " minCost " << minCost << " minArg " << argMin << endl;
			//cout << "aresta: " << *ki << " h nodo = " << h(neighbour,goal) << " h u = " << h(u,goal) << " x u =" << x[neighbour] << " y u = " << y[neighbour] << endl;
#endif			
			if( ((funcs[neighbour] + cost) < minCost )   && ( backPointer[neighbour] != u ) )
			{
				argMin = neighbour;
				minCost = funcs[neighbour] + cost;
			}
		}
		
		if( argMin < 1 ){
			//exit(-12);
			argMin = backPointer[u];
			//cout << " usado o negativo " << endl;
		}
#ifdef DEBUG
		cout << "  - o argmin eh " << argMin << endl;
#endif
		return argMin;
	}
Ad::Edge_Weight Ad::c(int u, int v){
		edge_descriptor e;
		bool edgeExists;

		tie(e,edgeExists) = edge(v, u, *graph);
		if(edgeExists)
			return get(edge_weight, *graph, e);
		else{
			cout << u <<"-"<<v<<endl;
			throw;
		}
	}

void Ad::computePath(int start, int goal, double e)
	{
		int min_id;
			if(!open->empty())
				min_id = open->key_top();
			else
				return;
#ifdef DEBUG
		std::cout << "Compute Path ===> min id=" << min_id << endl;
		printOpen();
#endif

		while( ( key(goal, goal, e) > key(min_id, goal, e) ) ||
			  // ( underconsistent(goal) ) ||
			( underconsistent(min_id) ) ) {
		int s = min_id;
			open->remove(min_id);

#ifdef DRAW_GRAPH
	if( Debug::iter > 0)
		color[s] = YELLOW;
	else
		color[s] = RED;
#endif

			if( ! underconsistent(s) )
			{
				Debug::adUnderNodesVisited++;
#ifdef DEBUG			
			std::cout << "  - - - expandindo estado ( overconsistent ) " << s << endl;
#endif
				v[s] = g[s];
				setAd[s] = CLOSED;
				updateSetMembership(s,goal,e);

				out_edge_iterator vi,vend;
				for(tie(vi,vend) = out_edges(s, *graph); vi != vend; vi++)
				{
					int custo = get(edge_weight, *graph, *vi);
					vertex_iterator p1 = target(*vi, *graph);		
					int neighbour = *p1;

					if( g[neighbour] > g[s] + custo)
					{
						backPointer[neighbour] = s;
						g[neighbour] = g[s] + custo;
						updateSetMembership(neighbour,goal,e);
					}
				}
				
			}
			else // propagating underconsistency
			{
#ifdef DEBUG			
			std::cout << "  - - - expandindo estado ( underconsistent ) " << s << endl;
#endif
				Debug::adOverNodesVisited++;
				v[s] = INFINITO;
				setAd[s] = NONE;
				updateSetMembership(s,goal,e);

				adjacency_iterator vi,vend;
				for(tie(vi,vend) = adjacent_vertices(s, *graph); vi != vend; vi++)
				{	
					int neighbour = *vi;
					if( backPointer[neighbour] == s)
					{
						backPointer[neighbour] = argmin(neighbour,v,goal);
						g[neighbour] = v[backPointer[neighbour]] + c(backPointer[neighbour],neighbour);

						updateSetMembership(neighbour,goal,e);
					}
				}
			}
			printOpen();
			if(!open->empty()){
				min_id = open->key_top();
#ifdef DEBUG
				cout << " o proximo a ser expandido eh " << min_id << endl;
#endif
			}else
				return;
		} // wend
		v[goal] = 0;
	}
vector<Ad::Change> Ad::getChanges(){
	vector<Change> dchanges;
	set< pair<int,int> > cjMud;
	static mt19937 gen((boost::uint32_t)time(0));

#ifdef ONLY_CHANGE_REGION
	if(Debug::porcentagem > 0.0)
			Debug::mudancas_passo = (unsigned long)(Debug::regionSize * Debug::porcentagem);

	while(cjMud.size() < Debug::mudancas_passo){
			pair<int, int> uv;
			out_edge_iterator vi,vend;
			Edge_Weight custo;

			uv.first = (gen() % num_nodes) + 1;
			int d = out_degree(uv.first, *graph);
			int nu = gen() % d;
			int j = 0;
			for(tie(vi,vend) = out_edges(uv.first, *graph); vi != vend; vi++){
				vertex_iterator p1 = target(*vi, *graph);	
				if( j == nu ){
					uv.second = *p1;
					custo = c(uv.first,uv.second);
					break;
				}
				j++;
			}
		
			if(g[nu] >= INFINITO && custo > 16)
				cjMud.insert(uv);
	}
#else
	if(Debug::only_path_vertices){
		if(Debug::porcentagem > 0.0)
			Debug::mudancas_passo = (unsigned long)(Debug::path.size() * Debug::porcentagem);

		while(cjMud.size() < Debug::mudancas_passo){

			int change = Debug::path[(gen() % (Debug::path.size()/2))];
			if(backPointer[change] > 0)
				cjMud.insert(make_pair(change, backPointer[change]));

		}

	}else{//mudancas aleatorias
		
		while(cjMud.size() < Debug::mudancas_passo){
			pair<int, int> uv;
			out_edge_iterator vi,vend;
			Edge_Weight custo;

			uv.first = (gen() % num_nodes) + 1;
			int d = out_degree(uv.first, *graph);
			int nu = gen() % d;
			int j = 0;
			for(tie(vi,vend) = out_edges(uv.first, *graph); vi != vend; vi++){
				vertex_iterator p1 = target(*vi, *graph);	
				if( j == nu ){
					uv.second = *p1;
					custo = c(uv.first,uv.second);
					break;
				}
				j++;
			}
		
			if(custo > 16)
				cjMud.insert(uv);
		}
	}
#endif
	Debug::iter++;
	//decrementos
	if(Debug::iter % 2 == 1 ){

		for( set< pair<int,int> >::iterator it = cjMud.begin(); it != cjMud.end(); it++){
			Edge_Weight custo = c(it->first,it->second);
			if( custo > 16){
				dchanges.push_back( make_pair( custo - ((gen() % 16) + 1) , make_pair(it->first,it->second)));
				Debug::changes.push( make_pair( custo + ((gen() % 16) + 1) , make_pair(it->first,it->second)));
			}
		}
	}else{ // incrementos
		
		while( !Debug::changes.empty()){
			dchanges.push_back(Debug::changes.top());
			Debug::changes.pop();
		}
		
	}
	return dchanges;
}

vector<Ad::Change> Ad::oldGetChanges(){
		vector<Change> changes;
		static mt19937 gen((boost::uint32_t)time(0));

#ifdef ONLY_CHANGE_PATH

		int change = Debug::path[(gen() % (Debug::path.size())) ];
		//int change = path[gen() % 10 ]; // limita para o final as mudan�as
		if( change == 1 )
			change++;
		Edge_Weight custo = c(change, backPointer[change]);
		changes.push_back(make_pair(( gen() % 150) + custo ,make_pair(change, backPointer[change])));

#else
		vector<pair<int,int> > jaFeitas;

		int numChanges = (gen() % (Debug::mudancas_passo)) + 1;

		for(int i = 0; i < numChanges; ++i)
		{

			int v = (gen() % num_nodes) + 1;

			int d = out_degree(v, *graph);
			
			int nu = gen() % d;
			int u;
			//int* o = new int[d];
			int j = 0;
			out_edge_iterator vi,vend;
			Edge_Weight custo = 0;
			for(tie(vi,vend) = out_edges(v, *graph); vi != vend; vi++){
				vertex_iterator p1 = target(*vi, *graph);	
				if( j == nu ){
					u = *p1;
					custo = get(edge_weight, *graph, *vi);
					break;
				}
				j++;
			}
			
			// nao muda mais de uma vez a mesma aresta
			jaFeitas.push_back(make_pair<int,int>(u,v));
			jaFeitas.push_back(make_pair<int,int>(v,u));
			vector< pair<int,int> >::iterator it = find(jaFeitas.begin(),jaFeitas.end(),make_pair<int,int>(u,v));
			bool naoAchou = it != jaFeitas.end();

			if(naoAchou)
				changes.push_back(make_pair(( gen() % 400)+ custo,make_pair(v, u)));
		}

#endif
		return changes;
	}

bool Ad::changeGraph(Change change){
	
	int newCost = change.first,
		u = change.second.first,
		v = change.second.second;
#ifdef DRAW_GRAPH
//color[u] = colorastar[u] = BLACK;
//color[v] = colorastar[v] = BLACK;
#endif
#ifdef DEBUG
	cout << u << " -> " << v << " (" << newCost << ")" << endl;
#endif
	edge_descriptor e;
	bool edgeExists;

	tie(e,edgeExists) = edge(u, v, *graph);

	if(edgeExists){
		put(edge_weight, *graph,  e, newCost);
		return true;
	}
	else
		return false;
}

void Ad::fixInitialization(int vtarget, int start, int goal, double e){
	if( vtarget != start )
	{
		backPointer[vtarget] = argmin(vtarget,v,goal);

		if(backPointer[vtarget] == 0) // quando o nodo esta fora da area de interesse, descarta
			return;

		g[vtarget] = max(v[backPointer[vtarget]] + c(vtarget, backPointer[vtarget]), v[backPointer[vtarget]]);

		updateSetMembership(vtarget,goal,e);
	}
}


Ad::Edge_Weight Ad::pathCost(const int &start, const int &goal)
{
	Edge_Weight custo = 0;

	Debug::path.clear();
	Debug::path.push_back(goal);

#ifdef SHOW_PATH 
	cout << endl << "ARA*: " << goal;
#endif
	for( int v = goal; v != start; v = backPointer[v])
	{
		edge_descriptor e;
		bool edgeExists;
#ifdef SHOW_PATH
		cout << " -> " << backPointer[v];
#endif

		// tirar a linha abaixo
		/*if( v == backPointer[backPointer[backPointer[v]]] || v == backPointer[backPointer[v]] ){
			cout << "erro!";
			cout << v << endl;
			cout << backPointer[v];
			cout << endl;
		}*/

		tie(e,edgeExists) = edge(v, backPointer[v], *graph);
		if(edgeExists){
			custo += get(edge_weight, *graph, e);

			if( find(Debug::path.begin(), Debug::path.end(), backPointer[v]) != Debug::path.end() ){
				cout << v << " <- " << backPointer[v] << endl; 				
				return 0;
			}
			Debug::path.push_back(backPointer[v]);
		}else{
			//cout << endl << "999 " << e << endl;
			return -1;
		}
		//cout << v << " e " << backPointer[v]  << " ->";
	}
#ifdef SHOW_PATH
	cout << endl;
#endif
	//cout << "t: " << g[goal] << " c:" << custo << endl;
	return custo;
}
void Ad::publishSolution(int &iter, double e, int start, int goal)
	{
		//adTimeInc,adTimeDec,aTimeInc,aTimeDec;
		//nvIncA, nvDecA, nvIncAd, nvIncA;
		lTimer.stop();
		
		pathCost(start, goal);

		Edge_Weight ad = g[goal], as = astar(start,goal,false);//, dj = dijkstra(start,goal);
		/*if( ad < as)
			cout << ad << " <> " << as << " dijk:" << dj << endl;
		else
			cout << ad << " == " << as << " dijk:" << dj << endl;
		*/
		if( iter == 0){
			Debug::regionSize = Debug::adUnderNodesVisited;
			Debug::adTimeInc = Debug::adTimeDec = 0.0001; 
			Debug::aTimeInc = Debug::aTimeDec = 0.0;
			Debug::nvIncA = Debug::nvDecA = Debug::nvIncOverAd = Debug::nvDecOverAd = Debug::nvIncUnderAd = Debug::nvDecUnderAd = 0;
			cout << "orig\t0\t" << Debug::aTime*1000 << "\t" << (lTimer.show() + 0.0001)*1000 << "\t" << ( Debug::aTime / (lTimer.show() + 0.0001))
				 << "\t" << Debug::aNodesVisited << "\t" << Debug::adUnderNodesVisited << "\t" << Debug::adOverNodesVisited <<  endl; 
		}else{ 
			if(Debug::iter % 2 == 1 ){
				Debug::aTimeDec += Debug::aTime;
				Debug::adTimeDec += lTimer.show();
				Debug::nvDecA += Debug::aNodesVisited;
				Debug::nvDecUnderAd += Debug::adUnderNodesVisited;
				Debug::nvDecOverAd += Debug::adOverNodesVisited;
				if((Debug::iter+1) == Debug::totIter){
					cout << "dec\t";
					cout << Debug::porcentagem << "\t" << Debug::aTimeDec/Debug::totIter*1000 << "\t" << Debug::adTimeDec/Debug::totIter*1000 << "\t" << /*std::setprecision(3) <<*/( Debug::aTimeDec/Debug::adTimeDec)
						 << "\t" << Debug::nvDecA/Debug::totIter   << "\t" << Debug::nvDecUnderAd/Debug::totIter  << "\t" << Debug::nvDecOverAd/Debug::totIter  << endl; 
				}
			}else{
				Debug::aTimeInc += Debug::aTime;
				Debug::adTimeInc += lTimer.show();
				Debug::nvIncA += Debug::aNodesVisited;
				Debug::nvIncUnderAd += Debug::adUnderNodesVisited;
				Debug::nvIncOverAd += Debug::adOverNodesVisited;

				if((Debug::iter) == Debug::totIter){
					cout << "inc\t";
					cout << Debug::porcentagem << "\t" << Debug::aTimeInc/Debug::totIter*1000 << "\t" << Debug::adTimeInc/Debug::totIter*1000 << "\t" << /*std::setprecision(3) <<*/( Debug::aTimeInc /Debug::adTimeInc)
						 << "\t" << Debug::nvIncA/Debug::totIter   << "\t" << Debug::nvIncUnderAd/Debug::totIter  << "\t" << Debug::nvIncOverAd/Debug::totIter  << endl; 
				}
			}
		}
		
		Debug::adUnderNodesVisited = Debug::adOverNodesVisited = 0;
		
		lTimer.start();

	}
void Ad::init_adstar(int start, int goal, double e)
	{
		for(int i = 0; i <= num_nodes; i++)
		{
			f[i] = g[i] = v[i] = INFINITO;
			backPointer[i] = 0;
			setAd[i] = NONE;
			color[i] = GRAY;
		}

		g[goal] = v[goal] = v[start] = INFINITO;
		backPointer[start] = NULL;
		g[start] = 0;

		open->clear();
		//for(int i = 0; i <= num_nodes; i++)
		//{
		//	setAd[i] = NONE;
		//} 
				
		key(start,goal,e);
		open->push(start);
		setAd[start] = OPEN;
	}

void Ad::main(int start, int goal, double e0, double delta, int maxIter){
	double e = e0;
	init_adstar(start,goal,e);

	lTimer.start();
	for(int iter = 0; iter <= maxIter; iter++)
	{
		computePath(start,goal,e);
		publishSolution(iter,e,start,goal);

		if ( e <= 1){
			changes.clear();
			changes = getChanges();

			for(vector<Change>::iterator c = changes.begin(); c != changes.end(); c++){
				if( ! changeGraph(*c) )
					exit(-97);
			}
			lTimer.start();
			for(vector<Change>::iterator c = changes.begin(); c != changes.end(); c++){
				int u = c->second.first;
				int p = c->second.second;

				if( v[u] > v[p]){					
					fixInitialization(u,start,goal,e);
					fixInitialization(p,start,goal,e);
				}
				else{
					fixInitialization(p,start,goal,e);
					fixInitialization(u,start,goal,e);
				}
			}
		}
		else
		{
			e -= delta;
		}
		move_incons_open();
		update_open(goal,e);
		clearClosed();
	}
}
























Ad::Ad(string co_filename, string dist_filename) : PI(3.14159265), DG_TO_RAD(PI / 180.0), INFINITO(std::numeric_limits<Edge_Weight>::max() / 2)
{
	ad = this;
	load(co_filename, dist_filename);
}

void Ad::initialize(){
	edge_array = new Edge[num_edges]();;
	weights = new Edge_Weight[num_edges]();
	x = new double[num_nodes+1]();
	y = new double[num_nodes+1]();
	drawx = new double[num_nodes+1]();
	drawy = new double[num_nodes+1]();
	color = new Color[num_nodes+1]();
	colorastar = new Color[num_nodes+1]();

	f = new Edge_Weight[num_nodes+1]();
	g = new Edge_Weight[num_nodes+1]();
	v = new Edge_Weight[num_nodes+1]();
	backPointer = new int[num_nodes+1]();
	setAd = new Sets[num_nodes+1]();
	open = new PriorityQueue<Edge_Weight>(num_nodes+1, f);
}
void Ad::initializeEdges(){
	graph = new G(edge_array, edge_array + num_edges, weights, num_nodes);
}
void Ad::load(string co_filename, string dist_filename)
{
	cout << "carregando arquivo de distancias " << dist_filename << endl;
	ifstream dist_file(dist_filename.c_str());

	char c[100];
	int i = 0;
	int a,b;
	int u,v;
	vector< vector < int > > edges;

	if( ! dist_file.good() ){
	  cout << "erro de arquivo" << endl;
	  exit(-89);
	}
	while( dist_file.good() ){
		dist_file >> c;
		if( c[0] == 'p'){
			dist_file >> c >> num_nodes >> num_edges;
			initialize();
			edges.resize(num_edges);
		}
		else if( c[0] == 'a'){
			bool adiciona = true;
			dist_file >> u >> v >> weights[i];

			for(unsigned int k = 0; k < edges[u].size(); k++)
			{
				if( edges[u][k] == v )
					adiciona = false;
			}

			if( adiciona ){
				edge_array[i].first = u;
				edge_array[i].second = v; 
				i++;
				edges[u].push_back(v);
			}
		}
		else{
			dist_file.getline(c,100);
		}
	}
	cout << " nodes: " << num_nodes << " edges: " << num_edges << endl;
	dist_file.close();
	initializeEdges();
	cout << "carregando arquivo de coordenadas" << endl;
	ifstream co_file(co_filename.c_str());

	while( co_file.good() ){
		co_file >> c;
		if( c[0] == 'v'){
			co_file >> i >> a >> b;
			
#ifdef DIST_EUCLIDIANA
			drawx[i] = x[i] = a;
			drawy[i] = y[i] = b;
#else       // transforma em graus, depois em radianos
			drawx[i] = a;
			drawy[i] = b;
			x[i] = ((double) a) / 1000000.0; 
			x[i] *= DG_TO_RAD;
			y[i] = ((double) b) / 1000000.0;
			y[i] *= DG_TO_RAD;
#endif
		}
		else{
			co_file.getline(c,100);
		}
	}
	co_file.close();
}


double Ad::h(int id, const int goal){
#ifdef DIST_EUCLIDIANA
	return sqrt( ((x[id] - x[goal])*(x[id] - x[goal])) + ((y[id] - y[goal])*(y[id] - y[goal])));
#else
	double d = acos((sin(y[id])*sin(y[goal]))+(cos(y[id])*cos(y[goal])*cos(x[id] - x[goal])));
	return 63787000.0 * d;
#endif
}





/********************************************************************
*        DIJKSTRA                                                   *
********************************************************************/
Ad::Edge_Weight Ad::dijkstra(int start, int goal){
	std::vector<vertex_descriptor> p(num_vertices(*graph));
	std::vector<int> d(num_vertices(*graph)+1);
	vertex_descriptor s = vertex(start,*graph);
	
	Timer dTimer;
	dTimer.start();
	dijkstra_shortest_paths(*graph, s, predecessor_map(&p[0]).distance_map(&d[0]));
	dTimer.stop();

	Debug::dijTime = dTimer.show();
	return d[goal];
}
/********************************************************************
*        A*                                                         *
********************************************************************/
template <class Graph, class CostType>
class distance_heuristic : public astar_heuristic<Graph, CostType>
{
public:
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  distance_heuristic(Vertex goal)
    : m_goal(goal) {}
  CostType operator()(Vertex u)
  {
    return (CostType) ad->h(u,m_goal);
  }
private:
  int m_goal;
};

struct found_goal {}; // exception for termination

// visitor that terminates when we find the goal
template <class Vertex>
class astar_goal_visitor : public boost::default_astar_visitor
{
public:
  astar_goal_visitor(Vertex goal) : m_goal(goal) {}
  template <class Graph>
  void examine_vertex(Vertex u, Graph& g) {
	Debug::aNodesVisited += 1;

#ifdef DRAW_GRAPH
	if(ad->colorastar[u] == Ad::GRAY)
		ad->colorastar[u]=Ad::RED;
#endif
	if(u == m_goal)
      throw found_goal();
  }
private:
  Vertex m_goal;
};


Ad::Edge_Weight Ad::astar(int start, int goal, bool verbose){
	vector<vertex_descriptor> p(num_vertices(*graph)+100);
	vector<Edge_Weight> d(num_vertices(*graph)+100);

	vertex_descriptor s = vertex(start,*graph),
					  g = vertex(goal,*graph);
	
	Timer timer;
	timer.start();

	Debug::aNodesVisited = 0;
#ifdef DRAW_GRAPH
	for(int i = 1; i <= num_nodes; i++){
		colorastar[i]=GRAY;
	}
#endif
	try{
	astar_search
      (*graph, s, distance_heuristic<G, Edge_Weight> (goal),
       predecessor_map(&p[0]).distance_map(&d[0]).
       visitor(astar_goal_visitor<vertex_descriptor>(g)));
	  } catch(found_goal) {

#ifdef SHOW_PATH		
		cout << "A*  : " ;

		bool detec = false;
		for(int t = goal, i = 0; t != start; t = p[t])
		{
			if( ! detec &&  (t != ad::path[i++])){
				cout << endl << "erro => " << t << " != " << ad::path[i-1]<< endl;
				detec = true;
			}
			cout << t << " -> ";
		}
#endif  

		timer.stop();	
		Debug::aTime = timer.show();
		if(verbose)
		{
			cout << Debug::aNodesVisited << ";" << timer << ";";
		}
	}
	return d[goal];
}


void Ad::test(int repeat){
	int start, goal;
	static mt19937 gen((boost::uint32_t)time(0));
	Debug::totIter = repeat/2;

	cout << "Mudancas aleatorias" << endl << "type\t#ch\ttA*\ttAD*\tspeedup\tvA*\tvAD* exp\tvAD* cor" << endl;
	Debug::only_path_vertices = false;
	Debug::iter = 0;
	Debug::mudancas_passo = 1;

	//for(int iter = 1; iter <= repeat; iter++){
		do{
			start = (gen()%num_nodes) + 1;
			goal = (gen()%num_nodes) + 1;
		}while( start == goal );
		Debug::iter = 0;
		main(start, goal, 1.0, 0.2, repeat);
	//}
//cout << "-----" << endl;
	for(int i = 0; i < 6; i++)
	{
		//Debug::mudancas_passo = (long) (num_nodes * Debug::mp[i]);
		Debug::porcentagem = Debug::mp[i];
		//for(int iter = 1; iter <= repeat; iter++){
			do{
				start = (gen()%num_nodes) + 1;
				goal = (gen()%num_nodes) + 1;
			}while( start == goal );
			Debug::iter = 0;
			main(start, goal, 1.0, 0.2, repeat);
		//	cout << "-----" << endl;
		//}
	}
/*
	cout << "Mudancas somente no caminho mininimo" << endl << "type\t#ch\ttA*\ttAD*\tspeedup\tvA*\tvAD* exp\tvAD* cor" << endl;
	Debug::only_path_vertices = true;
	Debug::mudancas_passo = 1;
	Debug::porcentagem = 0.0;
	//for(int iter = 1; iter <= repeat; iter++){
		do{
			start = (gen()%num_nodes) + 1;
			goal = (gen()%num_nodes) + 1;
		}while( start == goal );
		Debug::iter = 0;
		main(start, goal, 1.0, 0.2, repeat);
	//}
cout << "-----" << endl;
	for(int i = 2; i < 5; i++)
	{
		Debug::porcentagem = Debug::mp[i];
		for(int iter = 1; iter <= repeat; iter++){
			do{
				start = (gen()%num_nodes) + 1;
				goal = (gen()%num_nodes) + 1;
			}while( start == goal );
			Debug::iter = 0;
			main(start, goal, 1.0, 0.2, repeat);
			cout << "-----" << endl;
		}
	}
*/
}
//#ifdef DRAW_GRAPH
void Ad::graphicTest(){
	int start, goal;
	Debug::only_path_vertices = false;
	Debug::mudancas_passo = (long) (num_nodes * 0.01);
	Debug::iter = 0;
	Debug::totIter = 4;

	static mt19937 gen((boost::uint32_t)time(0));
	do{
			start = (gen()%num_nodes) + 1;
			goal = (gen()%num_nodes) + 1;
	}while( start == goal );

	main(start, goal, 1.0, 0.2, 1);

	color[start] = color[goal] = colorastar[start] = colorastar[goal] = BLACK;
}
//#endif