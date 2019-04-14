#pragma once
#include <string>
#include<vector>

class EDGE
{
public:
	// Data Members 
	int node_0;
	int node_1;
	int t1;
	int e1;
	int t2;             // if not exist -1
	int e2;             // if not exist -1
	double norm_vec[2];
	double length;
	double s;           // wave speed
	string type;        // Left / Right / Bottom / Top / Interior
	double midpt[2];
};

class ADJ_CELL
{
public:
	// Data Members
	int loc_edge;
	int adj_edge;
	int adj_tri;
};

class CELL
{
public:
	// Data Members 
	int vertex[3];
	int edge[3];
	int tri;                   // the index of triangle (cell)
	vector<ADJ_CELL> adj_cell;
	// double state[4];           // state vector
	// double R[4];               // residual vector
	double dt;                 // time step
	double A;                  // area
	double centroid[2];        // cnetroid
};

class BOUNDARY
{
public:
	// Data Members 
	int nBFace;   // number of faces on that boundary
	int nf;
	string Title;
	int *B;       // the address of B[0], indeices of edges are stored in B ? will that be overwrite?
};

class MESH
{
public:
	// Data Members 
	int nNode;
	int nElem;
	double **node_pos;
	int nBGroup;
	BOUNDARY *boundary;
	CELL *Elems;
	int nEdge;
	vector<EDGE> Edges;
};

class TEMP
{
public:
	// Data Members
	int ***B_temp;
	int ***edges_temp;
};

class PROP               // key properties of the state
{
public:
	double F[2][4];      // Euler flux
	double rho;          // density
	double v[2];         // velocity
	double c;            // speed of sound
	double E;
	double p;
	double H;
};