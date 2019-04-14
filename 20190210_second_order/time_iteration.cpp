void refresh_R(MESH *mesh, double **state, double **R)
{
	int i, j, k;
	// initialize grad to zero for all cells
	double ***grad;
	grad = new double**[mesh->nElem];
	for (i = 0;i < mesh->nElem;i++)
	{
		grad[i] = new double*[4];
		for (j = 0;j < 4;j++)
			grad[i][j] = new double[2];
	}
	for (i = 0;i < mesh->nElem;i++)
		for (j = 0;j < 4;j++)
			for (k = 0;k < 2;k++)
				grad[i][j][k] = 0;
	// find grad
	EDGE edge;
	double n[2], delta_l, U_hat[4];
	int t, t1, t2;
	for (i = 0;i < mesh->nEdge;i++)
	{
		edge = mesh->Edges[i];
		n[0] = edge.norm_vec[0];
		n[1] = edge.norm_vec[1];
		delta_l = edge.length;
		if (edge.type == "Interior")
		{
			// for interior cells
			t1 = edge.t1;
			t2 = edge.t2;
			for (j = 0;j < 4;j++)
			{
				U_hat[j] = 0.5*(state[t1][j] + state[t2][j]);
				for (k = 0;k < 2;k++)
				{
					grad[t1][j][k] = grad[t1][j][k] + U_hat[j] * n[k] * delta_l;
					grad[t2][j][k] = grad[t2][j][k] - U_hat[j] * n[k] * delta_l;
				}
			}
		}
		else
		{
			// for boundary cells
			t = edge.t1;
			for (j = 0;j < 4;j++)
			{
				U_hat[j] = state[t][j];
				for (k = 0;k < 2;k++)
					grad[t][j][k] = grad[t][j][k] + U_hat[j] * n[k] * delta_l;
			}
		}
	}
	for (i = 0; i < mesh->nElem; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 2; k++)
				grad[i][j][k] = grad[i][j][k] / mesh->Elems[i].A;
	// initialize residual and wave speed on each cell to zero
	for (i = 0;i < mesh->nElem;i++)
		for (j = 0;j < 4;j++)
			R[i][j] = 0.0;
	for (i = 0;i < mesh->nEdge;i++)
		mesh->Edges[i].s = 0.0;
	// find flux
	double U_L_star[4], U_R_star[4], F[4], s, midpt[2], centroid[2];
	double M_inf = 0.5;
	double T_t = 1 + ((gamma - 1) / 2) * M_inf * M_inf;
	double p_t = pow(T_t, gamma / (gamma - 1));
	double alpha = 0;
	double p_inf = 1.0;
	for (i = 0; i < mesh->nEdge; i++)
	{
		edge = mesh->Edges[i];
		delta_l = edge.length;
		n[0] = edge.norm_vec[0];
		n[1] = edge.norm_vec[1];
		midpt[0] = edge.midpt[0];
		midpt[1] = edge.midpt[1];
		if (edge.type == "Interior")
		{
			// for interior cells
			t1 = edge.t1;
			centroid[0] = mesh->Elems[t1].centroid[0];
			centroid[1] = mesh->Elems[t1].centroid[1];
			for (j = 0; j < 4; j++)
				U_L_star[j] = state[t1][j] + (midpt[0] - centroid[0])*grad[t1][j][0] + (midpt[1] - centroid[1])*grad[t1][j][1];
			t2 = edge.t2;
			centroid[0] = mesh->Elems[t2].centroid[0];
			centroid[1] = mesh->Elems[t2].centroid[1];
			for (j = 0; j < 4; j++)
				U_R_star[j] = state[t2][j] + (midpt[0] - centroid[0])*grad[t2][j][0] + (midpt[1] - centroid[1])*grad[t2][j][1];
			Roe_flux(U_L_star, U_R_star, n, F, &s); // F and s are return values
			for (j = 0; j < 4; j++)
			{
				R[t1][j] = R[t1][j] + F[j] * delta_l;
				R[t2][j] = R[t2][j] - F[j] * delta_l;
			}
			mesh->Edges[i].s = s;
		}
		if (edge.type == "Bottom" || edge.type == "Top")
		{
			// for wall
			t = edge.t1;
			centroid[0] = mesh->Elems[t].centroid[0];
			centroid[1] = mesh->Elems[t].centroid[1];
			for (j = 0; j < 4; j++)
				U_L_star[j] = state[t][j] + (midpt[0] - centroid[0])*grad[t][j][0] + (midpt[1] - centroid[1])*grad[t][j][1];
			wall_flux(U_L_star, n, F, &s);
			/*// ----- free stream test -----
			double U_inf[4];
			define_U(0.5, 0, gamma, U_inf);
			Roe_flux(U_L_star, U_inf, n, F, &s);
			// ----- free stream test ----- */
			for (j = 0; j < 4; j++)
				R[t][j] = R[t][j] + F[j] * delta_l;
			mesh->Edges[i].s = s;
		}
		if (edge.type == "Left")
		{
			// for inflow
			t = edge.t1;
			centroid[0] = mesh->Elems[t].centroid[0];
			centroid[1] = mesh->Elems[t].centroid[1];
			for (j = 0; j < 4; j++)
				U_L_star[j] = state[t][j] + (midpt[0] - centroid[0])*grad[t][j][0] + (midpt[1] - centroid[1])*grad[t][j][1];
			inlet_flux(U_L_star, n, T_t, p_t, alpha, F, &s);
			/*// ----- free stream test -----
			double U_inf[4];
			define_U(0.5, 0, gamma, U_inf);
			Roe_flux(U_L_star, U_inf, n, F, &s);
			// ----- free stream test ----- */
			for (j = 0; j < 4; j++)
				R[t][j] = R[t][j] + F[j] * delta_l;
			mesh->Edges[i].s = s;
		}
		if (edge.type == "Right")
		{
			// for outflow
			t = edge.t1;
			centroid[0] = mesh->Elems[t].centroid[0];
			centroid[1] = mesh->Elems[t].centroid[1];
			for (j = 0; j < 4; j++)
				U_L_star[j] = state[t][j] + (midpt[0] - centroid[0])*grad[t][j][0] + (midpt[1] - centroid[1])*grad[t][j][1];
			// flux evaluation
			outflow_flux(U_L_star, n, p_inf, F, &s);
			/*// ----- free stream test -----
			double U_inf[4];
			define_U(0.5, 0, gamma, U_inf);
			Roe_flux(U_L_star, U_inf, n, F, &s);
			// ----- free stream test -----*/
			for (j = 0; j < 4; j++)
				R[t][j] = R[t][j] + F[j] * delta_l;
			mesh->Edges[i].s = s;
		}
	}
	// free memory
	for (i = 0; i < mesh->nElem; i++)
	{
		for (j = 0; j < 4; j++)
			delete[] grad[i][j];
		delete[] grad[i];
	}
	delete[] grad;
}

void compute_time_step(MESH *mesh, double CFL)
{
	int i, j, global_edge_no;
	double SUM, s, delta_l;
	for (i = 0; i < mesh->nElem; i++)
	{
		SUM = 0;
		for (j = 0; j < 3; j++)
		{
			global_edge_no = mesh->Elems[i].edge[j];
			s = mesh->Edges[global_edge_no].s;
			delta_l = mesh->Edges[global_edge_no].length;
			SUM = SUM + s*delta_l;
		}
		mesh->Elems[i].dt = 2 * CFL / SUM;
	}
}

double RK2_step_1(MESH *mesh, double **state, double **state_FE)
{
	int i, j;
	double **R;
	R = new double*[mesh->nElem];
	for (i = 0;i < mesh->nElem;i++)
		R[i] = new double[4];
	double CFL = 0.85;
	refresh_R(mesh, state, R);
	compute_time_step(mesh, CFL);
	double max_R = 0.0, abs_R;
	for (i = 0; i < mesh->nElem; i++)
	{
		for (j = 0; j < 4; j++)
		{
			state_FE[i][j] = state[i][j] - mesh->Elems[i].dt * R[i][j];
			abs_R = abs(R[i][j]);
			if (abs_R > max_R)
				max_R = abs_R;
		}
	}
	// free memory
	for (i = 0; i < mesh->nElem; i++)
		delete[] R[i];
	delete[] R;
	return max_R;
}

void RK2_step_2(MESH *mesh, double **state, double **state_FE)
{
	int i, j;
	double **R_FE;
	R_FE = new double*[mesh->nElem];
	for (i = 0;i < mesh->nElem;i++)
		R_FE[i] = new double[4];
	refresh_R(mesh, state_FE, R_FE);
	for (i = 0; i < mesh->nElem; i++)
		for (j = 0; j < 4; j++)
			state[i][j] = 0.5*(state[i][j] + state_FE[i][j] - mesh->Elems[i].dt * R_FE[i][j]);
	// free memory
	for (i = 0; i < mesh->nElem; i++)
		delete[] R_FE[i];
	delete[] R_FE;
}

void auto_save(MESH *mesh, double **state, string fname)
{
	ofstream output;
	output.open(fname);
	for (int i = 0; i < mesh->nElem; i++)
	{
		for (int j = 0; j < 4; j++)
			output << state[i][j] << " ";
		output << endl;
	}
	output.close();
}

void save_history(double *R_L_inf, int nTime_step_curr, string fname)
{
	ofstream output;
	output.open(fname);
	output << nTime_step_curr << endl;
	for (int i = 0; i < nTime_step_curr; i++)
		output << R_L_inf[i] << endl;
	output.close();
}

void time_iteration(MESH *mesh, double **state, int nTime_step, string fname, double *R_L_inf)
{
	int i,j;
	double **state_FE;
	state_FE = new double*[mesh->nElem];
	for (i = 0;i < mesh->nElem;i++)
		state_FE[i] = new double[4];

	for (i = 0;i < nTime_step;i++)
	{
		// cout << "KR2_STEP1" << endl;
		R_L_inf[i] = RK2_step_1(mesh, state, state_FE);
		// cout << "KR2_STEP2" << endl;
		RK2_step_2(mesh, state, state_FE);
		if (R_L_inf[i] < 1.0e-7)
		{
			cout << "iteration = " << i << ", |R|_L_inf = " << R_L_inf[i] << endl;
			auto_save(mesh, state, fname + "_state_final" + ".txt");
			save_history(R_L_inf, i, fname + "_history" + ".txt");
			break;
		}
		if (i % 1000 == 0)
		{
			cout << "iteration = " << i << ", |R|_L_inf = " << R_L_inf[i] << endl;
			auto_save(mesh, state, fname + "_state_iter_" + to_string(i) + ".txt");
		}
		
	}
}