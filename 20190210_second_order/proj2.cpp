#include "stdafx.h"
#include "read_data.cpp"
#include "pre_calculation.cpp"
#include "initialization.cpp"
#include "fluxes.cpp"
#include "time_iteration.cpp"

void free_stream_test(MESH *mesh)
{
	// task 3 (a)
	cout << "task 3 (a)" << endl;
	double R_L_inf[5000];
	string fname = "data/task_3_a/free_stream_test";
	double **state;
	state = new double*[mesh->nElem];
	for (int i = 0;i < mesh->nElem;i++)
		state[i] = new double[4];
	initialization(mesh, state);
	time_iteration(mesh, state, 1, fname, R_L_inf);
	save_history(R_L_inf, 1, fname + "_history" + ".txt");
	// task 3 (b)
	cout << "task 3 (b)" << endl;
	fname = "data/task_3_b/free_stream_preservation_test";
	initialization(mesh, state);
	time_iteration(mesh, state, 5000, fname, R_L_inf);
	save_history(R_L_inf, 5000, fname + "_history" + ".txt");
}

void FVM(string case_name)
{
	MESH mesh;
	TEMP temp;
	long long int nTime_step = 100000;
	double R_L_inf[100000];
	read_data("mesh/" + case_name + ".gri", &mesh, &temp);
	pre_calculation(&mesh, temp);
	// free_stream_test(&mesh);
	double **state;
	state = new double*[mesh.nElem];
	for (int i = 0; i < mesh.nElem; i++)
		state[i] = new double[4];
	load_state(&mesh, "data/first_order_solution/" + case_name + "_state_final.txt", state);
	time_iteration(&mesh, state, nTime_step, "data/task_3_c/" + case_name + "/" + case_name, R_L_inf);
}

/*
void boundary_reconstruct(MESH *mesh, double **state, int boundary_edge_no, double *U_star)
{
	CELL cell;
	EDGE edge;
	int adj_cell_no, i, j, k, t1, t2;
	double grad[4][2], U_hat[4];
	adj_cell_no = mesh->Edges[boundary_edge_no].t1;
	cell = mesh->Elems[adj_cell_no];
	for (j = 0; j < 4; j++)
		for (k = 0; k < 2; k++)
			grad[j][k] = 0.0;
	// find grad
	for (i = 0; i < 3; i++)
	{
		edge = mesh->Edges[cell.edge[i]];
		if (edge.type == "Interior")
		{
			t1 = edge.t1;
			t2 = edge.t2;
			for (j = 0; j<4; j++)
				U_hat[j] = 0.5*(state[t1][j] + state[t2][j]);
			if (adj_cell_no == t1) // adj_cell is the Left cell
			{
				for (j = 0;j<4;j++)
					for (k=0;k<2;k++)
						grad[j][k] = grad[j][k] + U_hat[j] * edge.norm_vec[k] * edge.length;
			}
			else // adj_cell == t2, i.e. adj_cell is the Right cell
			{
				for (j = 0; j<4; j++)
					for (k = 0; k<2; k++)
						grad[j][k] = grad[j][k] - U_hat[j] * edge.norm_vec[k] * edge.length;
			}
		}
		else // edge.type !=  "Interior"
		{
			for (j = 0; j<4; j++)
				U_hat[j] = state[adj_cell_no][j];
			for (j = 0; j<4; j++)
				for (k = 0; k<2; k++)
					grad[j][k] = grad[j][k] + U_hat[j] * edge.norm_vec[k] * edge.length;
		}
	}
	for (j = 0; j < 4; j++)
		for (k = 0; k < 2; k++)
			grad[j][k] = grad[j][k] / cell.A;
	// find U_star
	double centroid[2], midpt[2];
	for (k = 0; k < 2; k++)
	{
		centroid[k] = mesh->Elems[adj_cell_no].centroid[k];
		midpt[k] = mesh->Edges[boundary_edge_no].midpt[k];
	}
	for (j = 0; j < 4; j++)
		U_star[j] = state[adj_cell_no][j] + (midpt[0] - centroid[0])*grad[j][0] + (midpt[1] - centroid[1])*grad[j][1];
}*/

void outputs(string case_name)
{
	MESH mesh;
	TEMP temp;
	PROP prop;
	ofstream output;
	string fname = "data/task_3_c/" + case_name + "/" + case_name + "_cp_distribution.txt";

	int i, j, adj_cell;
	double **state;
	double n[2], U_star[4], v_b[2], p_b, cl = 0, cd = 0;
	double p_inf = 1.0, M_inf = 0.5, h = 0.0625;
	double x_i, cp_i;

	cout << case_name << " ";
	cout << scientific;
	read_data("mesh/" + case_name + ".gri", &mesh, &temp);
	printf("%5d ", mesh.nElem);
	pre_calculation(&mesh, temp);
	
	state = new double*[mesh.nElem];
	for (int i = 0; i < mesh.nElem; i++)
		state[i] = new double[4];
	load_state(&mesh, "data/task_3_c/" + case_name + "/" + case_name + "_state_final.txt",state);
	output.open(fname);
	/*-------------------------------------------------------------------------------*/
	for (i = 0;i < mesh.nEdge;i++)
	{
		if (mesh.Edges[i].type == "Bottom")
		{
			n[0] = mesh.Edges[i].norm_vec[0];
			n[1] = mesh.Edges[i].norm_vec[1];
			// boundary_reconstruct(&mesh, state, i, U_star);
			// Euler_flux(U_star, &prop);
			adj_cell = mesh.Edges[i].t1;
			Euler_flux(state[adj_cell], &prop);
			for (j = 0; j<2; j++)
				v_b[j] = prop.v[j] - (prop.v[0] * n[0] + prop.v[1] * n[1]) * n[j];
			p_b = (gamma - 1)*(prop.rho*prop.E - 0.5*prop.rho*(v_b[0] * v_b[0] + v_b[1] * v_b[1]));
			cl = cl + (p_b - p_inf) * n[1] * mesh.Edges[i].length;
			cd = cd + (p_b - p_inf) * n[0] * mesh.Edges[i].length;
			// find cp
			x_i = mesh.Edges[i].midpt[0];
			cp_i = (p_b - p_inf) / (0.5*gamma*p_inf*M_inf*M_inf);
			output << x_i << " " << cp_i << endl;
		}
	}
	output.close();
	// normalize cl and cd
	cl = cl / (0.5 * gamma * p_inf * M_inf*M_inf*h);
	cd = cd / (0.5 * gamma * p_inf * M_inf*M_inf*h);
	cout << cl << " " << cd << " ";
	/*-------------------------------------------------------------------------------*/
	// find Es
	double T_t = 1 + ((gamma - 1) / 2) * M_inf * M_inf;
	double p_t = pow(T_t, gamma / (gamma - 1));
	double rho_t = p_t / (R*T_t);
	double s_t = p_t / pow(rho_t, gamma);
	double sums = 0, suma = 0, s, Es;
	for (i = 0;i < mesh.nElem;i++)
	{
		Euler_flux(state[i], &prop);
		s = prop.p / pow(prop.rho, gamma);
		sums = sums + (s / s_t - 1)*(s / s_t - 1)*mesh.Elems[i].A;
		suma = suma + mesh.Elems[i].A;
	}
	Es = sqrt(sums / suma);
	cout << Es << endl;
}

int main()
{
	/*FVM("bump0");
	FVM("bump1");
	FVM("bump2");
	FVM("bump3");
	FVM("bump4");
	*/
	cout << "bump | DOF |     cl     |     cd     |     Es     |" << endl;
	outputs("bump0");
	outputs("bump1");
	outputs("bump2");
	outputs("bump3");
	outputs("bump4");
	
	system("pause");
	return 0;
}