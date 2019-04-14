void define_U(double M, double alpha, double gamma, double U[4])
{
	
	U[0] = 1;
	U[1] = M*cos(alpha);
	U[2] = M*sin(alpha);
	U[3] = 1 / (gamma - 1) + 0.5*M*M*gamma;
	// U[3] = 1 / ((gamma - 1)*gamma) + 0.5*M*M; // <--- AE 523 Proj 2 ?
}

void initialization(MESH *mesh, double **state)
{
	double alpha = 0;                                 // attack angle (radians)
	double M_inf = 0.5;                               // Mach 0.5 flow
	double U_inf[4];
	define_U(M_inf, alpha, gamma, U_inf);
	/*-------------------------------------------------------------------------------*/
	int i, j;
	for (i = 0;i < mesh->nElem;i++)
		for (j = 0;j < 4;j++)
			state[i][j] = U_inf[j];
}

void load_state(MESH *mesh, string fname, double **state)
{
	ifstream input(fname);
	for (int i = 0;i < mesh->nElem;i++)
		for (int j = 0;j < 4;j++)
			input >> state[i][j];
}