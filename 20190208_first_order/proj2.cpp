#include "stdafx.h"
#include "read_data.cpp"
#include "pre_calculation.cpp"
#include "initialization.cpp"
#include "fluxes.cpp"
#include "time_iteration.cpp"

void Roe_flux_test()
{
	// double gamma = 1.4;                           // specific heat ratio for air
	double U_L[4], U_R[4];
	define_U(0.5, 0, gamma, U_L);                 // Mach, alpha, gamma, U
	define_U(1.0, 0, gamma, U_R);
	/*-------------------------------------------------------------------------------*/
	// task 1 (a)
	cout << "task 1 (a)" << endl;
	double n[2] = { 1, 1 };
	double F_Roe[4], F_L[4];
	double s;
	int i;
	Roe_flux(U_L, U_L, n, F_Roe, &s); // F and s are return values
	cout << "Roe flux: ";
	for (int i = 0; i < 4; i++)
		cout << F_Roe[i] << " ";
	cout << endl;
	PROP prop_L;
	Euler_flux(U_L, &prop_L);
	for (i = 0; i < 4; i++)
		F_L[i] = prop_L.F[0][i] * n[0] + prop_L.F[1][i] * n[1];
	cout << "Euler flux: ";
	for (int i = 0; i < 4; i++)
		cout << F_L[i] << " ";
	cout << endl;
	/*-------------------------------------------------------------------------------*/
	// task 1 (b)
	cout << "task 1 (b)" << endl;
	double n_1[2] = { 1, 1 }, F_Roe_1[4];
	Roe_flux(U_L, U_R, n_1, F_Roe_1, &s); // F and s are return values
	cout << "Roe flux: ";
	for (int i = 0; i < 4; i++)
		cout << F_Roe_1[i] << " ";
	cout << endl;
	double n_2[2] = { -1, -1 }, F_Roe_2[4];
	Roe_flux(U_R, U_L, n_2, F_Roe_2, &s);
	cout << "Roe flux flipped: ";
	for (int i = 0; i < 4; i++)
		cout << F_Roe_2[i] << " ";
	cout << endl;
	/*-------------------------------------------------------------------------------*/
	// task 1 (c)
	cout << "task 1 (c)" << endl;
	define_U(1.5, 0, gamma, U_L);                 // Mach, alpha, gamma, U
	define_U(1.4, 0, gamma, U_R);
	Roe_flux(U_L, U_R, n, F_Roe, &s);
	cout << "Roe flux: ";
	for (int i = 0; i < 4; i++)
		cout << F_Roe[i] << " ";
	cout << endl;
	Euler_flux(U_L, &prop_L);
	for (i = 0; i < 4; i++)
		F_L[i] = prop_L.F[0][i] * n[0] + prop_L.F[1][i] * n[1];
	cout << "Euler flux: ";
	for (int i = 0; i < 4; i++)
		cout << F_L[i] << " ";
	cout << endl;
}

void free_stream_test(MESH *mesh)
{
	// task 2 (a)
	cout << "task 2 (a)" << endl;
	double R_L_inf[5000];
	string fname = "data/task_2_a/free_stream_test";
	initialization(mesh);
	time_iteration(mesh, 1, fname, R_L_inf);
	save_history(R_L_inf, 1, fname + "_history" + ".txt");
	// task 2 (b)
	cout << "task 2 (b)" << endl;
	fname = "data/task_2_b/free_stream_preservation_test";
	initialization(mesh);
	time_iteration(mesh, 5000, fname, R_L_inf);
	save_history(R_L_inf, 5000, fname + "_history" + ".txt");
}

void FVM(string case_name)
{
	MESH mesh;
	TEMP temp;
	long long int nTime_step = 1000000;
	string fname;
	double R_L_inf[100000];

	read_data("mesh/" + case_name + ".gri", &mesh, &temp);
	pre_calculation(&mesh, temp);
	// Roe_flux_test();
	// free_stream_test(&mesh);
	initialization(&mesh);
	fname = "data/task_2_c/" + case_name + "/" + case_name;
	time_iteration(&mesh, nTime_step, fname, R_L_inf);
}

void load_state(MESH *mesh, string fname)
{
	ifstream input(fname);
	for (int i = 0;i < mesh->nElem;i++)
		for (int j = 0;j < 4;j++)
			input >> mesh->Elems[i].state[j];
}

void outputs(string case_name)
{
	MESH mesh;
	TEMP temp;
	PROP prop;
	ofstream output;
	string fname = "data/task_2_c/" + case_name + "/" + case_name + "_cp_distribution.txt";

	int i, j, adj_cell, node_0, node_1;
	double p_inf = 1.0, M_inf = 0.5, h = 0.0625;
	double n[2], v_b[2], p_b, cl = 0, cd = 0;
	double x_i, cp_i;

	cout << case_name << " ";
	cout << scientific;
	read_data("mesh/" + case_name + ".gri", &mesh, &temp);
	printf("%5d ", mesh.nElem);
	pre_calculation(&mesh, temp);
	load_state(&mesh, "data/task_2_c/" + case_name + "/" + case_name + "_state_final.txt");
	output.open(fname);
	/*-------------------------------------------------------------------------------*/
	for (i = 0;i < mesh.nEdge;i++)
	{
		if (mesh.Edges[i].type == "Bottom")
		{
			n[0] = mesh.Edges[i].norm_vec[0];
			n[1] = mesh.Edges[i].norm_vec[1];
			adj_cell = mesh.Edges[i].t1;
			Euler_flux(mesh.Elems[adj_cell].state, &prop);
			for (j = 0; j<2; j++)
				v_b[j] = prop.v[j] - (prop.v[0] * n[0] + prop.v[1] * n[1]) * n[j];
			p_b = (gamma - 1)*(prop.rho*prop.E - 0.5*prop.rho*(v_b[0] * v_b[0] + v_b[1] * v_b[1]));
			cl = cl + (p_b - p_inf) * n[1] * mesh.Edges[i].length;
			cd = cd + (p_b - p_inf) * n[0] * mesh.Edges[i].length;
			// find cp
			node_0 = mesh.Edges[i].node_0;
			node_1 = mesh.Edges[i].node_1;
			x_i = 0.5*(mesh.node_pos[node_0][0] + mesh.node_pos[node_1][0]);
			cp_i = (p_b - p_inf) / (0.5*gamma*p_inf*M_inf*M_inf);
			output << x_i << " " << cp_i << endl;
		}
	}
	output.close();
	// normalize cl and cd
	cl = cl / (0.5 * gamma * p_inf * M_inf*M_inf*h);
	cd = cd / (0.5 * gamma * p_inf * M_inf*M_inf*h);
	cout << cl << " ";
	cout << cd << " ";
	/*-------------------------------------------------------------------------------*/
	// find Es
	double T_t = 1 + ((gamma - 1) / 2) * M_inf * M_inf;
	double p_t = pow(T_t, gamma / (gamma - 1));
	double rho_t = p_t / (R*T_t);
	double s_t = p_t / pow(rho_t, gamma);
	double sums = 0, suma = 0, s, Es;
	for (i = 0;i < mesh.nElem;i++)
	{
		Euler_flux(mesh.Elems[i].state, &prop);
		s = prop.p / pow(prop.rho, gamma);
		sums = sums + (s / s_t - 1)*(s / s_t - 1)*mesh.Elems[i].A;
		suma = suma + mesh.Elems[i].A;
	}
	Es = sqrt(sums / suma);
	cout << Es << endl;
}

int main()
{
	/*
	FVM("bump0");
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