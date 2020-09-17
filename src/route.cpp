#include <iostream> //cout
#include <sstream>  //input string
#include <fstream>  
#include <string>   //string
#include <unordered_map>
#include <vector>
#include <math.h> //max
#include <algorithm> //sort
#include <time.h> //time
#include <png.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#define PNG_NO_SETJMP
using namespace std;
struct BIN;
struct NET;
struct EDGE;


struct BIN {
	int id;
	int X;
	int Y;
	EDGE* edgeLft = nullptr;
	EDGE* edgeRgt = nullptr;
	EDGE* edgeUp = nullptr;
	EDGE* edgeDn = nullptr;
};

struct NET {
	int id;
	string name;
	int pinNum;
	int binId_s;
	int binId_t;
	int absX;
	int absY;
	int overflow;
	vector<int> edgesId;
};

struct EDGE {
	int id;
	double capacity;
	double demand;
	double history;
	BIN* binLft = nullptr;
	BIN* binRgt = nullptr;
	BIN* binUp = nullptr;
	BIN* binDn = nullptr;
	vector<int> netsId;
};

struct GRID {
	int binId;
	double cost;
	BIN* parentBin = nullptr;
	EDGE* parentEdge = nullptr;
};

unordered_map<int, BIN*> bins;
unordered_map<int, EDGE*> edges;
unordered_map<int, NET*> nets;

int cp_netsBoundingBox_ascending(const void *elem1,const void *elem2 )
{
  return ((NET*)elem1)->absX + ((NET*)elem1)->absY < ((NET*)elem2)->absX + ((NET*)elem2)->absY;
}

int cp_netsBoundingBox_descending(const void *elem1,const void *elem2 )
{
  return ((NET*)elem1)->absX + ((NET*)elem1)->absY > ((NET*)elem2)->absX + ((NET*)elem2)->absY;
}


int cp_gridsCost_descending(const void *elem1,const void *elem2 )
{
  return ((GRID*)elem1)->cost > ((GRID*)elem2)->cost;
}




void first_rout(NET* net){
	BIN* bin_n = bins[net->binId_s];
	BIN* bin_t = bins[net->binId_t];
	unordered_map<int, GRID*> grids;
	vector<GRID*> grids_order;
	
	GRID* grid_s = new GRID();
	grid_s->binId = net->binId_s;
	grid_s->cost = 0;
	grids.insert (std::make_pair(net->binId_s, grid_s));
	while(grids.find(net->binId_t) == grids.end()){
		for(int i = 0; i < 4; ++i){
			int binFind_id;
			EDGE* edge;
			BIN* bin_tmp;
			if(i%4 == 0 && bin_n->edgeUp) {binFind_id = bin_n->edgeUp->binUp->id; edge = bin_n->edgeUp;}
			else if(i%4 == 1 && bin_n->edgeDn) {binFind_id = bin_n->edgeDn->binDn->id; edge = bin_n->edgeDn;}
			else if(i%4 == 2 && bin_n->edgeLft) {binFind_id = bin_n->edgeLft->binLft->id; edge = bin_n->edgeLft;}
			else if(i%4 == 3 && bin_n->edgeRgt) {binFind_id = bin_n->edgeRgt->binRgt->id; edge = bin_n->edgeRgt;}
			else continue;
			
			double cost;
			bin_tmp = bins[binFind_id];
			if(grids.find(binFind_id) == grids.end()){
				GRID* grid = new GRID();
				grid->binId = binFind_id;
				//cost = pow(max(abs(bin_t->X - bin_tmp->X),abs(bin_t->Y - bin_tmp->Y)) * (edge->demand+1) / edge->capacity, 5.0);
				cost = grids[bin_n->id]->cost + pow(max(abs(bin_t->X - bin_tmp->X),abs(bin_t->Y - bin_tmp->Y)), 5.0) + pow((edge->demand+1) / edge->capacity, 5.0);
				grid->cost = cost;
				grid->parentBin = bin_n;
				grid->parentEdge = edge;
				grids.insert(std::make_pair(binFind_id, grid));
				grids_order.push_back(grid);
			}else{
				GRID* grid = grids.find(binFind_id)->second;
				//cost = grid->cost + pow(max(abs(bin_t->X - bin_tmp->X),abs(bin_t->Y - bin_tmp->Y)) * (edge->demand+1) / edge->capacity, 5.0);
				cost = grids[bin_n->id]->cost + pow(max(abs(bin_t->X - bin_tmp->X),abs(bin_t->Y - bin_tmp->Y)), 5.0) + pow((edge->demand+1) / edge->capacity, 5.0);
				if(cost < grid->cost){
					grid->cost = cost;
					grid->parentBin = bin_n;
					grid->parentEdge = edge;
					if(find(grids_order.begin(), grids_order.end(), grid) == grids_order.end())
						grids_order.push_back(grid);
					//grids_order.push_back(grid);
				}
			}
		}
		sort(grids_order.begin(), grids_order.end(), cp_gridsCost_descending);
		bin_n = bins[grids_order.back()->binId];
		grids_order.pop_back();
	}
	/*****trace back*****/
	GRID* grid_n = grids[net->binId_t];
	while(grid_n != grid_s){
		EDGE* edge = grid_n->parentEdge;
		net->edgesId.push_back(edge->id);
		++edge->demand;
		edge->netsId.push_back(net->id);
		grid_n = grids[grid_n->parentBin->id];
	}
	reverse(net->edgesId.begin(), net->edgesId.end());
}


void rip_up(NET* netRerout, int round){
	//parameter
	double k1 = 5.0, k2 = 1.5, k3 = 3.0, B = 5.0, r = 0.1;
	BIN* bin_n = bins[netRerout->binId_s];
	BIN* bin_t = bins[netRerout->binId_t];
	unordered_map<int, GRID*> grids;
	vector<GRID*> grids_order;
	
	GRID* grid_s = new GRID();
	grid_s->binId = netRerout->binId_s;
	grid_s->cost = 0;
	grids.insert (std::make_pair(netRerout->binId_s, grid_s));
	while(grids.find(netRerout->binId_t) == grids.end()){
		for(int i = 0; i < 4; ++i){
			int binFind_id;
			EDGE* edge;
			if(i%4 == 0 && bin_n->edgeUp) {binFind_id = bin_n->edgeUp->binUp->id; edge = bin_n->edgeUp;}
			else if(i%4 == 1 && bin_n->edgeDn) {binFind_id = bin_n->edgeDn->binDn->id; edge = bin_n->edgeDn;}
			else if(i%4 == 2 && bin_n->edgeLft) {binFind_id = bin_n->edgeLft->binLft->id; edge = bin_n->edgeLft;}
			else if(i%4 == 3 && bin_n->edgeRgt) {binFind_id = bin_n->edgeRgt->binRgt->id; edge = bin_n->edgeRgt;}
			else continue;
			double adj = k3 * (1 - exp((0 - B)*exp((0-r)*round)));
			double f = round*(k2+adj)/(round*(k2+adj)-(edge->history-1));
			double cost;
			BIN* bin_tmp = bins[binFind_id];
			if(grids.find(binFind_id) == grids.end()){
				GRID* grid = new GRID();
				grid->binId = binFind_id;
				// original cost = grids[bin_n->id]->cost + edge->history * pow((edge->demand+1) / edge->capacity, k1);
				// paper function cost = grids[bin_n->id]->cost + f * pow((edge->demand+1) / edge->capacity, k1);
				// distance cost = grids[bin_n->id]->cost + max(abs(bin_t->X - bin_tmp->X),abs(bin_t->Y - bin_tmp->Y))*(1 - exp((0 - B)*exp((0-r)*round))) + edge->history * pow((edge->demand+1) / edge->capacity, k1);
				// distance cost better = grids[bin_n->id]->cost + max(abs(bin_t->X - bin_tmp->X),abs(bin_t->Y - bin_tmp->Y))*(1 - exp((0 - B)*exp((0-r)*round))) + pow(f *(edge->demand+1) / edge->capacity, k1);
				cost = grids[bin_n->id]->cost + (1 - exp((0 - B)*exp((0-r)*round))) + pow(f *(edge->demand+1) / edge->capacity, k1);
				grid->cost = cost;
				grid->parentBin = bin_n;
				grid->parentEdge = edge;
				grids.insert(std::make_pair(binFind_id, grid));
				grids_order.push_back(grid);
			}else{
				GRID* grid = grids.find(binFind_id)->second;
				// original cost = grids[bin_n->id]->cost + edge->history * pow((edge->demand+1) / edge->capacity, k1);
				// paper function cost = grids[bin_n->id]->cost + f * pow((edge->demand+1) / edge->capacity, k1);
				// distance cost = grids[bin_n->id]->cost + max(abs(bin_t->X - bin_tmp->X),abs(bin_t->Y - bin_tmp->Y))*(1 - exp((0 - B)*exp((0-r)*round))) + edge->history * pow((edge->demand+1) / edge->capacity, k1);
				// distance cost better = grids[bin_n->id]->cost + max(abs(bin_t->X - bin_tmp->X),abs(bin_t->Y - bin_tmp->Y))*(1 - exp((0 - B)*exp((0-r)*round))) + pow(f *(edge->demand+1) / edge->capacity, k1);
				cost = grids[bin_n->id]->cost + (1 - exp((0 - B)*exp((0-r)*round))) + pow(f *(edge->demand+1) / edge->capacity, k1);
				if(cost < grid->cost){
					grid->cost = cost;
					grid->parentBin = bin_n;
					grid->parentEdge = edge;
					if(find(grids_order.begin(), grids_order.end(), grid) == grids_order.end())
						grids_order.push_back(grid);
					//grids_order.push_back(grid);
				}
			}
		}
		sort(grids_order.begin(), grids_order.end(), cp_gridsCost_descending);
		bin_n = bins[grids_order.back()->binId];
		grids_order.pop_back();
	}
	//trace back*****
	GRID* grid_n = grids[netRerout->binId_t];
	while(grid_n != grid_s){
		EDGE* edge = grid_n->parentEdge;
		netRerout->edgesId.push_back(edge->id);
		++edge->demand;
		edge->netsId.push_back(netRerout->id);
		grid_n = grids[grid_n->parentBin->id];
	}
	reverse(netRerout->edgesId.begin(), netRerout->edgesId.end());
}


int main(int argc, char const *argv[])
{		
	unsigned long start = clock();
	unsigned long now = clock();
	/****************/
	/**input region**/
	/****************/
	string in_str;
	int tilesNum_h, tilesNum_v;
	double capaNum_h, capaNum_v;
	int netNum;
	int x1,x2, y1, y2;
	fstream fs(argv[1], fstream::in);
	fs >> in_str >> tilesNum_v >> tilesNum_h;
	fs >> in_str >> in_str >> capaNum_h;
	fs >> in_str >> in_str >> capaNum_v;
	fs >> in_str >> in_str >> netNum;
	
	int edge_v = 0, edge_h = (tilesNum_v - 1)* tilesNum_h;
	for(int i = 0; i < tilesNum_h; ++i){
		for(int j = 0; j < tilesNum_v; ++j){
			BIN* bin = new BIN();
			bin->id = i*tilesNum_v+j;
			bin->X = j;
			bin->Y = i;
			bins.insert (std::make_pair(bin->id, bin));
			if(j < tilesNum_v - 1){
				EDGE* edge = new EDGE();
				edge->id = edge_v;
				edge->demand = 0;
				edge->capacity = capaNum_v;
				edge->history = 1;
				edge->binLft = bin;
				++edge_v;
				bin->edgeRgt = edge;
				edges.insert (std::make_pair(edge->id, edge));
			}if(j > 0){
				bin->edgeLft = edges[bin->id - bin->Y - 1];
				bin->edgeLft->binRgt = bin;
			}
			if(i < tilesNum_h - 1){
				EDGE* edge = new EDGE();
				edge->id = bin->id + edge_h;
				edge->demand = 0;
				edge->capacity = capaNum_h;
				edge->history = 1;
				edge->binDn = bin;
				bin->edgeUp = edge;
				edges.insert (std::make_pair(edge->id, edge));
			}if(i > 0){
				bin->edgeDn = edges[bin->id + edge_h - tilesNum_v];
				bin->edgeDn->binUp = bin;
			}
		}
	}
	
	for(int i = 0; i < netNum; ++i){
		NET* net = new NET();
		fs >> net->name >> net->id >> net->pinNum;
		fs >> x1 >> y1;
		net->binId_s = y1*tilesNum_v + x1;
		fs >> x2 >> y2;
		net->binId_t = y2*tilesNum_v + x2;
		net->absX = abs(x1 - x2);
		net->absY = abs(y1 - y2);
		net->overflow = 0;
 		nets.insert (std::make_pair(net->id, net));
	}
	
	/*****************/
	/**first routing**/
	/*****************/
	vector<NET*> nets_order;
	for(auto &net : nets)
		nets_order.push_back(net.second);
	sort(nets_order.begin(), nets_order.end(), cp_netsBoundingBox_ascending);

	for(auto &net : nets_order){
		first_rout(net);
	}
	
	/*****************/
	/***Debug stage***/
	/*****************/
	int edge_overflow = 0;
	for(auto edge : edges){
		if(edge.second->demand > edge.second->capacity){
			edge_overflow+=edge.second->demand-edge.second->capacity;
			++edge.second->history;
			for(auto netId : edge.second->netsId){
				++nets[netId]->overflow;
			}
		}
	}//cout << "edge_overflow : " << edge_overflow << endl;

	/******************/
	/***Rerout stage***/
	/******************/
	int round = 1;
	int best_edge_overflow = edge_overflow;
	int bad = 0;
	while((now-start)/CLOCKS_PER_SEC< 550){
		vector<NET*> nets_rerout;
		for(auto &net : nets){
			if(net.second->overflow>0) nets_rerout.push_back(net.second);
		}
		for(auto &netRerout : nets_rerout){
			for(auto &edgeId : netRerout->edgesId){
				EDGE* edge = edges[edgeId];
				--edge->demand;
				int Q;
				for(Q = 0;Q<edge->netsId.size();Q++){
					if(edge->netsId[Q]==netRerout->id)
						break;
				}
				edge->netsId.erase(edge->netsId.begin() + Q);
			}
			netRerout->edgesId.clear();	
		}
		
		sort(nets_rerout.begin(), nets_rerout.end(), cp_netsBoundingBox_ascending);
		for(auto &netRerout : nets_rerout){
			rip_up(netRerout, round);

		}
		
		for(auto &net : nets){
			net.second->overflow = 0;
			for(auto &edgeId : net.second->edgesId){
				EDGE* edge= edges[edgeId];
				if(edge->demand > edge->capacity) ++net.second->overflow;
			}
		}
		
		edge_overflow = 0;
		for(auto &edge : edges){
			if(edge.second->demand > edge.second->capacity){
				++edge.second->history;
				//debug
				edge_overflow+=edge.second->demand-edge.second->capacity;
			}
		}//cout << "edge_overflow : " << edge_overflow << endl;
		
		if(edge_overflow < best_edge_overflow) {bad = 0; best_edge_overflow = edge_overflow;}
		else bad++;
		if(bad > 10) break;
		if(edge_overflow < 1) break;
		now = clock();
		++round;
	}
	
	
	//debug
	/*int net_not_overflow = 0;
	for(auto &net : nets){
		net.second->overflow = 0;
		for(auto &edgeId : net.second->edgesId){
			EDGE* edge= edges[edgeId];
			if(edge->demand > edge->capacity) ++net.second->overflow;
		}
		if(net.second->overflow == 0) ++net_not_overflow;
	}cout << "net_not_overflow : " << net_not_overflow << endl;*/

	
	/******************/
	/***Output Block***/
	/******************/
	int j;
	std::ofstream ofs (argv[2], std::ofstream::out);
	for(int i = 0; i < netNum; ++i){
		NET* net = nets[i];
		ofs << net->name << " " <<  i << endl;
		BIN* bin_n = bins[net->binId_s];
		if(bin_n != bins[net->binId_s]) cout << "s error\n";
		for(auto &edgeId : net->edgesId){
			EDGE* edge = edges[edgeId];
			BIN* bin_tmp = bin_n;
			ofs<<"(" << bin_n->X << ", " << bin_n->Y << ", 1)-" ;
			if(bin_n->edgeUp == edge) bin_n = edge->binUp;
			else if(bin_n->edgeDn == edge) bin_n = edge->binDn;
			else if(bin_n->edgeLft == edge) bin_n = edge->binLft;
			else if(bin_n->edgeRgt == edge) bin_n = edge->binRgt;
			else{bin_n = NULL; cout << "error\n"; cin >> j;}
			ofs<<"(" << bin_n->X << ", " << bin_n->Y << ", 1)" << endl;

			if(abs(bin_tmp->X-bin_n->X)+abs(bin_tmp->Y-bin_n->Y)!=1) cout << "break\n";
		}
		if(bin_n != bins[net->binId_t]) cout << "t error\n";
		ofs<<"!" << endl; 
	}
    ofs.close();
	now = clock();
	cout << "time :" << (now-start)/CLOCKS_PER_SEC << endl;
	
}