void Euler_flux(double U[4], PROP *prop)
{
	prop->rho = U[0];
	prop->v[0] = U[1] / U[0];
	prop->v[1] = U[2] / U[0];
	prop->E = U[3] / U[0];
	prop->p = (gamma - 1)*(prop->rho*prop->E - 0.5*prop->rho*(prop->v[0] *prop->v[0] + prop->v[1]*prop->v[1]));
	prop->c = sqrt(gamma*prop->p / prop->rho); // speed of sound
	prop->H = prop->E + prop->p / prop->rho;
	prop->F[0][0] = prop->rho * prop->v[0];
	prop->F[0][1] = prop->rho * prop->v[0] * prop->v[0] + prop->p;
	prop->F[0][2] = prop->rho * prop->v[0] * prop->v[1];
	prop->F[0][3] = prop->rho * prop->v[0] * prop->H;
	prop->F[1][0] = prop->rho * prop->v[1];
	prop->F[1][1] = prop->rho * prop->v[1] * prop->v[0];
	prop->F[1][2] = prop->rho * prop->v[1] * prop->v[1] + prop->p;
	prop->F[1][3] = prop->rho * prop->v[1] * prop->H;
}

void Roe_flux(double U_L[4], double U_R[4], double n[2], double F[4], double *s)
{
	PROP prop_L, prop_R;
	Euler_flux(U_L, &prop_L);
	Euler_flux(U_R, &prop_R);
	// find Roe avg velocity v
	double v[2];
	int i;
	for (i = 0;i < 2;i++)
		v[i] = (sqrt(U_L[0])*prop_L.v[i] + sqrt(U_R[0])*prop_R.v[i]) / (sqrt(U_L[0]) + sqrt(U_R[0]));
	// find Roe avg velocity magnitude q
	double q = sqrt(v[0] * v[0] + v[1] * v[1]);
	// find the normal velocity component u
	double u = v[0] * n[0] + v[1] * n[1];
	// find Roe avg enthalpy
	double H = (sqrt(U_L[0])*prop_L.H + sqrt(U_R[0])*prop_R.H) / (sqrt(U_L[0]) + sqrt(U_R[0]));
	double F_L[4], F_R[4];
	for (i = 0;i < 4;i++)
	{
		F_L[i] = prop_L.F[0][i] * n[0] + prop_L.F[1][i] * n[1];
		F_R[i] = prop_R.F[0][i] * n[0] + prop_R.F[1][i] * n[1];
	}
	double lambda[3], c, epsilon;
	c = sqrt((gamma - 1)*(H - 0.5*q*q));
	lambda[0] = abs(u + c);
	lambda[1] = abs(u - c);
	lambda[2] = abs(u);
	epsilon = 0.1*c;
	*s = 0;
	for (i = 0;i < 3;i++)
	{
		if (lambda[i] < epsilon) // entropy fix
			lambda[i] = (epsilon * epsilon + lambda[i] * lambda[i]) / (2 * epsilon);
		if (lambda[i] > *s) // find the max lambda, save for finding time step
			*s = lambda[i];
	}
	double C_1, C_2, G_1, G_2, s_1, s_2;
	s_1 = 0.5*(lambda[0] + lambda[1]);
	s_2 = 0.5*(lambda[0] - lambda[1]);
	G_1 = (gamma - 1)*(0.5*q*q*(U_R[0] - U_L[0])   \
		- (v[0] * (U_R[1] - U_L[1])  \
			+ v[1] * (U_R[2] - U_L[2])) \
		+ (U_R[3] - U_L[3]));
	G_2 = -u*(U_R[0] - U_L[0]) + ((U_R[1] - U_L[1])*n[0]
		+ (U_R[2] - U_L[2])*n[1]);
	C_1 = (G_1 / (c*c)) * (s_1 - lambda[2]) + (G_2 / c) * s_2;
	C_2 = (G_1 / c) * s_2 + (s_1 - lambda[2]) * G_2;
	/*-----debug-----debug-----debug-----debug-----*/
	/*
	cout << "-----Roe flux intermediate values-----" << endl;
	cout << "lambda" << endl;
	for (i = 0; i < 3; i++)
		cout << lambda[i] << " ";
	cout << endl;
	cout << "s_1 = " << s_1 << endl << "s_2 = " << s_2 << endl;
	cout << "G_1 = " << G_1 << endl << "G_2 = " << G_2 << endl;
	cout << "C_1 = " << C_1 << endl << "C_2 = " << C_2 << endl;
	cout << "-----Roe flux intermediate values-----" << endl;
	*/
	/*-----debug-----debug-----debug-----debug-----*/
	F[0] = 0.5*(F_L[0] + F_R[0]) - 0.5*(lambda[2] * (U_R[0] - U_L[0]) + C_1);
	F[1] = 0.5*(F_L[1] + F_R[1]) - 0.5*(lambda[2] * (U_R[1] - U_L[1]) + C_1*v[0] + C_2*n[0]);
	F[2] = 0.5*(F_L[2] + F_R[2]) - 0.5*(lambda[2] * (U_R[2] - U_L[2]) + C_1*v[1] + C_2*n[1]);
	F[3] = 0.5*(F_L[3] + F_R[3]) - 0.5*(lambda[2] * (U_R[3] - U_L[3]) + C_1*H + C_2*u);
}

void wall_flux(double U[4], double n[2], double F_b[4], double *s)
{
	PROP prop;
	Euler_flux(U, &prop);
	//-------------------------------------------------------------
	double v_b[2];
	for (int i=0;i<2;i++)
		v_b[i] = prop.v[i] - (prop.v[0] * n[0] + prop.v[1] * n[1]) * n[i];
	double p_b = (gamma - 1)*(prop.rho*prop.E - 0.5*prop.rho*(v_b[0]* v_b[0] + v_b[1]*v_b[1]));
	F_b[0] = 0;
	F_b[1] = p_b*n[0];
	F_b[2] = p_b*n[1];
	F_b[3] = 0;
	*s = prop.c; // alternatively we may use p_b ?
}

double root_finder(double a, double b, double c)
{
	double s[2];
	s[0] = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
	s[1] = (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
	if (s[1]<0)
		return s[0];
	else
		return s[1];
}

void inlet_flux(double U_plus[4], double n[2], double T_t, double p_t, double alpha, double F_b[4], double *s)
{
	PROP prop;
	Euler_flux(U_plus, &prop);
	double u_n_plus = prop.v[0] * n[0] + prop.v[1] * n[1];
	double J_plus = u_n_plus + 2 * prop.c / (gamma - 1);
	double d_n = cos(alpha)*n[0] + sin(alpha)*n[1];
	double A, B, C, M_b;
	A = gamma*R*T_t*d_n*d_n - 0.5*(gamma - 1)*J_plus*J_plus;
	B = (4 * gamma*R*T_t*d_n) / (gamma - 1);
	C = (4 * gamma*R*T_t) / ((gamma - 1)*(gamma - 1)) - J_plus*J_plus;
	M_b = root_finder(A, B, C);
	double T_b = T_t / (1 + 0.5*(gamma - 1)*M_b*M_b);
	double p_b = p_t * pow(T_b / T_t, gamma / (gamma - 1));
	double rho_b = p_b / (R*T_b);
	double c_b = sqrt(gamma*p_b / rho_b);
	double v_b[2];
	v_b[0] = M_b*c_b*cos(alpha);
	v_b[1] = M_b*c_b*sin(alpha);
	double rho_E_b = p_b / (gamma - 1) + 0.5*rho_b*(v_b[0] * v_b[0] + v_b[1] * v_b[1]);
	// note that this expression is inconsistant with rho E in AE523 Proj 2 ?
	double U_b[4];
	U_b[0] = rho_b;
	U_b[1] = rho_b*v_b[0];
	U_b[2] = rho_b*v_b[1];
	U_b[3] = rho_E_b;
	PROP prop_b;
	Euler_flux(U_b, &prop_b);
	for (int i = 0;i < 4;i++)
		F_b[i] = prop_b.F[0][i] * n[0] + prop_b.F[1][i] * n[1];
	*s = prop_b.c;
}

void outflow_flux(double U_plus[4], double n[2], double p_b, double F_b[4], double *s)
{
	PROP prop;
	Euler_flux(U_plus, &prop);
	double u_n_plus = prop.v[0] * n[0] + prop.v[1] * n[1];
	double J_plus = u_n_plus + 2 * prop.c / (gamma - 1);
	double S_plus = prop.p / (pow( prop.rho , gamma ));
	double rho_b = pow( p_b/S_plus , 1/gamma );
	double c_b = sqrt(gamma*p_b / rho_b);
	double u_n_b = J_plus - 2 * c_b / (gamma - 1);
	double v_b[2];
	for (int i = 0; i < 2; i++)
		v_b[i] = prop.v[i] - u_n_plus*n[i] + u_n_b*n[i];
	double rho_E_b = p_b / (gamma - 1) + 0.5*rho_b*(v_b[0] * v_b[0] + v_b[1] * v_b[1]);
	// note that this expression is inconsistant with rho E in AE523 Proj 2 ?
	double U_b[4];
	U_b[0] = rho_b;
	U_b[1] = rho_b*v_b[0];
	U_b[2] = rho_b*v_b[1];
	U_b[3] = rho_E_b;
	PROP prop_b;
	Euler_flux(U_b, &prop_b);
	for (int i = 0; i < 4; i++)
		F_b[i] = prop_b.F[0][i] * n[0] + prop_b.F[1][i] * n[1];
	*s = prop_b.c;
}
