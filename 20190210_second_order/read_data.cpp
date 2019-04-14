void read_data(string fname, MESH *mesh, TEMP *temp)
{
	ifstream input(fname);
	int i, j;
	int nNode, nElem, Dim;
	input >> nNode >> nElem >> Dim;
	double **node_pos;
	node_pos = new double *[nNode];
	/*-------------------------------------------------------------------------------*/
	for (i = 0; i < nNode; i++)
	{
		node_pos[i] = new double[2];
		input >> node_pos[i][0] >> node_pos[i][1];
	}
	/*-------------------------------------------------------------------------------*/
	int nBGroup;
	int n_0, n_1;
	input >> nBGroup;
	BOUNDARY *boundary;
	boundary = new BOUNDARY[nBGroup];
	// Use a temporary list "B_temp" to store the boundary nodes before the establishment of the global edge list
	// "B_temp" is a 3D list, nBGroup * nBFace * 2
	int ***B_temp;
	B_temp = new int**[nBGroup];
	for (i = 0; i < nBGroup; i++)
	{
		input >> boundary[i].nBFace >> boundary[i].nf >> boundary[i].Title;
		B_temp[i] = new int *[boundary[i].nBFace];
		for (j = 0; j < boundary[i].nBFace; j++)
		{
			input >> n_0 >> n_1;
			n_0 = n_0 - 1;
			n_1 = n_1 - 1;
			B_temp[i][j] = new int[2];
			B_temp[i][j][0] = n_0;
			B_temp[i][j][1] = n_1;
		}
	}
	/*-------------------------------------------------------------------------------*/
	int Order;
	int v_0, v_1, v_2;
	string Basis;
	CELL *E;
	E = new CELL[nElem];
	// Use a temporary list "edges_temp" to store the edges before the establishment of the global edge list
	int ***edges_temp;
	edges_temp = new int**[nElem];
	input >> nElem >> Order >> Basis;
	for (i = 0; i < nElem; i++)
	{
		input >> v_0 >> v_1 >> v_2;
		v_0 = v_0 - 1;
		v_1 = v_1 - 1;
		v_2 = v_2 - 1;
		E[i].tri = i;
		E[i].vertex[0] = v_0;
		E[i].vertex[1] = v_1;
		E[i].vertex[2] = v_2;
		edges_temp[i] = new int*[3];
		for (j = 0; j < 3; j++)
		{
			edges_temp[i][j] = new int[4];
			edges_temp[i][j][2] = i;
			edges_temp[i][j][3] = j;
		}
		edges_temp[i][0][0] = v_1; // e0 -> n0
		edges_temp[i][0][1] = v_2; // e0 -> n1
		edges_temp[i][1][0] = v_2; // e1 -> n0
		edges_temp[i][1][1] = v_0; // e1 -> n1
		edges_temp[i][2][0] = v_0; // e2 -> n0
		edges_temp[i][2][1] = v_1; // e2 -> n1
	}
	input.close();
	mesh->nNode = nNode;
	mesh->nElem = nElem;
	mesh->node_pos = node_pos;
	mesh->nBGroup = nBGroup;
	mesh->boundary = boundary;
	mesh->Elems = E;
	mesh->nEdge = 0;
	// mesh.Edges = no edges established;
	temp->B_temp = B_temp;
	temp->edges_temp = edges_temp;
}