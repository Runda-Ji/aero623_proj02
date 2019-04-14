/* How to use gsl library
https://blog.csdn.net/raodotcong/article/details/8998379
https://zoujiemeng.github.io/post/2018/03/gsl-vs-linux.html
https://www.cnblogs.com/flipped/p/9314461.html

copy .dll files from ../../gsl-cmake-x64/bin/Release/ to working directory

properties -> VC++ Directories -> Include Directories
Add ../gsl-my-files/header/

properties -> VC++ Directories -> Library Directories
Add ../gsl-my-files/lib/

properties -> Linker -> Input -> Additional Dependencies
Add gsl.lib; gslcblas.lib;
*/

#include "gsl/gsl_spmatrix.h"

void create_adjacent_cell(int t1, int e1, int t2, int e2, MESH *mesh)
{
	ADJ_CELL new_info_1, new_info_2;
	new_info_1.loc_edge = e1;
	new_info_1.adj_edge = e2;
	new_info_1.adj_tri = t2;
	mesh->Elems[t1].adj_cell.push_back(new_info_1);
	new_info_2.loc_edge = e2;
	new_info_2.adj_edge = e1;
	new_info_2.adj_tri = t1;
	mesh->Elems[t2].adj_cell.push_back(new_info_2);
}

void check(int *edge, gsl_spmatrix *conn_vertex, MESH *mesh)
{
	int n0, n1, t, e;
	n0 = edge[0];
	n1 = edge[1];
	t = edge[2];
	e = edge[3];
	int gloal_edge_no;
	if (gsl_spmatrix_get(conn_vertex, n0, n1) == 0)
	{
		gloal_edge_no = mesh->nEdge;
		EDGE new_edge;
		new_edge.node_0 = n0;
		new_edge.node_1 = n1;
		new_edge.t1 = t;
		new_edge.e1 = e;
		new_edge.t2 = -1; //will be updated when visit again
		new_edge.e2 = -1;
		mesh->Edges.push_back(new_edge);
		mesh->nEdge = mesh->nEdge + 1;
		gsl_spmatrix_set(conn_vertex, n0, n1, gloal_edge_no + 1);
		gsl_spmatrix_set(conn_vertex, n1, n0, gloal_edge_no + 1);
		mesh->Elems[t].edge[e] = gloal_edge_no;
		//the index in sparse matrix conn_vertex starts from 1
	}
	else
	{
		//if the edge already in the sparse matrix
		int t1, e1, t2, e2;
		t2 = t;
		e2 = e;
		//decide which edge already in the list
		gloal_edge_no = gsl_spmatrix_get(conn_vertex, n0, n1) - 1;
		mesh->Elems[t2].edge[e2] = gloal_edge_no;
		t1 = mesh->Edges[gloal_edge_no].t1;
		e1 = mesh->Edges[gloal_edge_no].e1;
		//store the new edge
		mesh->Edges[gloal_edge_no].t2 = t2;
		mesh->Edges[gloal_edge_no].e2 = e2;
		mesh->Edges[gloal_edge_no].type = "Interior";
		create_adjacent_cell(t1, e1, t2, e2, mesh);
	}
}

void find_adjecent_cell(MESH *mesh, TEMP temp, gsl_spmatrix *conn_vertex)
{
	int nNode = mesh->nNode;
	int nElem = mesh->nElem;
	int ***edges_temp = temp.edges_temp;
	int i, j;
	for (i = 0; i < nElem; i++)
		for (j = 0; j < 3; j++)
			check(edges_temp[i][j], conn_vertex, mesh);
}

/*-------------------------------------------------------------------------------*/

void furnish_boundary(MESH *mesh, TEMP temp, gsl_spmatrix *conn_vertex)
{
	int i, j, n0, n1, global_edge_no;
	string Title;
	int ***B_temp = temp.B_temp;
	for (i = 0;i < mesh->nBGroup;i++)
	{//4 boundary groups
		Title = mesh->boundary[i].Title;
		mesh->boundary[i].B = new int[mesh->boundary[i].nBFace]; // allocate nBFace integers in B
		for (j = 0;j < mesh->boundary[i].nBFace;j++)
		{
			n0 = B_temp[i][j][0];
			n1 = B_temp[i][j][1];
			global_edge_no = gsl_spmatrix_get(conn_vertex, n0, n1) - 1;
			mesh->Edges[global_edge_no].type = Title;
			mesh->boundary[i].B[j] = global_edge_no;
		}
	}
}

/*-------------------------------------------------------------------------------*/

void furnish_edge(MESH *mesh)
{
	int i, A, B;
	double xA, yA, xB, yB, vec_len;
	for (i = 0;i < mesh->nEdge;i++)
	{
		A = mesh->Edges[i].node_0;
		B = mesh->Edges[i].node_1;
		xA = mesh->node_pos[A][0];
		yA = mesh->node_pos[A][1];
		xB = mesh->node_pos[B][0];
		yB = mesh->node_pos[B][1];
		//find length
		vec_len = sqrt((xA - xB)*(xA - xB) + (yA - yB)*(yA - yB));
		mesh->Edges[i].length = vec_len;
		//find normal vector;
		mesh->Edges[i].norm_vec[0] = (yB - yA) / vec_len;
		mesh->Edges[i].norm_vec[1] = (xA - xB) / vec_len;
		//find mid point
		mesh->Edges[i].midpt[0] = 0.5*(xA + xB);
		mesh->Edges[i].midpt[1] = 0.5*(yA + yB);
	}
}

void furnish_cell(MESH *mesh)
{
	int i, j, k, global_edge_no, global_node_no;
	double len[3], S, vertex_pos[3][2];
	for (i=0;i < mesh->nElem;i++)
	{
		for (j = 0;j < 3;j++)
		{
			global_edge_no = mesh->Elems[i].edge[j];
			len[j] = mesh->Edges[global_edge_no].length;
		}
		S = 0.5*(len[0] + len[1] + len[2]);
		mesh->Elems[i].A = sqrt(S*(S - len[0])*(S - len[1])*(S - len[2]));
		for (j = 0;j < 3;j++)
		{
			global_node_no = mesh->Elems[i].vertex[j];
			for (k = 0;k < 2;k++)
				vertex_pos[j][k] = mesh->node_pos[global_node_no][k];
		}
		for (k = 0;k < 2;k++)
			mesh->Elems[i].centroid[k] = (vertex_pos[0][k] + vertex_pos[1][k] + vertex_pos[2][k]) / 3.0;
	}
}

void pre_calculation(MESH *mesh, TEMP temp)
{
	gsl_spmatrix *conn_vertex = gsl_spmatrix_alloc(mesh->nNode, mesh->nNode);
	find_adjecent_cell(mesh, temp, conn_vertex);
	furnish_boundary(mesh, temp, conn_vertex);
	furnish_edge(mesh);
	furnish_cell(mesh);
}