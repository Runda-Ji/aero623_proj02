void refresh_R(MESH *mesh)
{
	//empty R for all cells
	int i, j;
	for (i = 0; i < mesh->nElem; i++)
		for (j = 0; j < 4; j++)
			mesh->Elems[i].R[j] = 0.0;
	EDGE edge;
	double n[2], delta_l, U_L[4], U_R[4], F[4], s;
	int t, t1, t2;
	double M_inf = 0.5;
	double T_t = 1 + ((gamma - 1) / 2) * M_inf * M_inf;
	double p_t = pow(T_t, gamma / (gamma - 1));
	double alpha = 0;
	double p_inf = 1.0;
	for (i = 0; i < mesh->nEdge; i++)
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
			for (j = 0; j < 4; j++)
			{
				U_L[j] = mesh->Elems[t1].state[j];
				U_R[j] = mesh->Elems[t2].state[j];
			}
			Roe_flux(U_L, U_R, n, F, &s); // F and s are return values
			for (j = 0; j < 4; j++)
			{
				mesh->Elems[t1].R[j] = mesh->Elems[t1].R[j] + F[j] * delta_l;
				mesh->Elems[t2].R[j] = mesh->Elems[t2].R[j] - F[j] * delta_l;
			}
			mesh->Edges[i].s = s;
		}
		if (edge.type == "Bottom" || edge.type == "Top")
		{
			// for wall
			t = edge.t1;
			for (j = 0; j < 4; j++)
				U_L[j] = mesh->Elems[t].state[j];
			wall_flux(U_L, n, F, &s);
			/* // ----- free stream test -----
			double U_inf[4];
			define_U(0.5, 0, gamma, U_inf);
			Roe_flux(U_L, U_inf, n, F, &s);
			// ----- free stream test ----- */
			for (j = 0; j < 4; j++)
				mesh->Elems[t].R[j] = mesh->Elems[t].R[j] + F[j] * delta_l;
			mesh->Edges[i].s = s;
		}
		if (edge.type == "Left")
		{
			// for inflow
			t = edge.t1;
			for (j = 0; j < 4; j++)
				U_L[j] = mesh->Elems[t].state[j];
			inlet_flux(U_L, n, T_t, p_t, alpha, F, &s);
			/* // ----- free stream test -----
			double U_inf[4];
			define_U(0.5, 0, gamma, U_inf);
			Roe_flux(U_L, U_inf, n, F, &s);
			// ----- free stream test ----- */
			for (j = 0; j < 4; j++)
				mesh->Elems[t].R[j] = mesh->Elems[t].R[j] + F[j] * delta_l;
			mesh->Edges[i].s = s;
		}
		if (edge.type == "Right")
		{
			// for outflow
			t = edge.t1;
			for (j = 0; j < 4; j++)
				U_L[j] = mesh->Elems[t].state[j];
			// flux evaluation
			outflow_flux(U_L, n, p_inf, F, &s);
			/*// ----- free stream test -----
			double U_inf[4];
			define_U(0.5, 0, gamma, U_inf);
			Roe_flux(U_L, U_inf, n, F, &s);
			// ----- free stream test -----*/
			for (j = 0; j < 4; j++)
				mesh->Elems[t].R[j] = mesh->Elems[t].R[j] + F[j] * delta_l;
			mesh->Edges[i].s = s;
		}
	}
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

double single_time_step(MESH *mesh)
{
	double max_R, dt, R[4], abs_R;
	int i, j;
	max_R = 0.0;
	for (i = 0; i < mesh->nElem; i++)
	{
		dt = mesh->Elems[i].dt;
		for (j = 0; j < 4; j++)
		{
			R[j] = mesh->Elems[i].R[j];
			mesh->Elems[i].state[j] = mesh->Elems[i].state[j] - dt*R[j];
			abs_R = abs(R[j]);
			if (abs_R > max_R)
				max_R = abs_R;
		}
	}
	return max_R;
}

void auto_save(MESH *mesh, string fname)
{
	ofstream output;
	output.open(fname);
	for (int i = 0; i < mesh->nElem; i++)
	{
		for (int j = 0; j < 4; j++)
			output << mesh->Elems[i].state[j] << " ";
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

void time_iteration(MESH *mesh, int nTime_step, string fname, double *R_L_inf)
{
	int i;
	double CFL = 0.85;
	for (i = 0;i < nTime_step;i++)
	{
		refresh_R(mesh);
		compute_time_step(mesh, CFL);
		R_L_inf[i] = single_time_step(mesh);
		if (R_L_inf[i] < 1.0e-7)
		{
			cout << "iteration = " << i << ", |R|_L_inf = " << R_L_inf[i] << endl;
			auto_save(mesh, fname + "_state_final" + ".txt");
			save_history(R_L_inf, i, fname + "_history" + ".txt");
			break;
		}
		if (i % 1000 == 0)
		{
			cout << "iteration = " << i << ", |R|_L_inf = " << R_L_inf[i] << endl;
			auto_save(mesh, fname + "_state_iter_" + to_string(i) + ".txt");
		}
	}
}