# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 10:47:03 2019

@author: rundaji
"""
import numpy as np
import matplotlib.pyplot as plt

def read_gri_reduced(i):
    fname = 'mesh\\bump%d.gri'%i;
    f = open(fname, "r");
    #read some general info
    nNode,nElem,Dim = [int(string) for string in f.readline().split()];
    node_pos = [None]*nNode;
    #read the position of nodes
    for i in range(0,nNode):
        x,y = [float(string) for string in f.readline().split()];
        node_pos[i] = [x,y];
    node_pos =  np.asarray(node_pos);
    #read the number of boundary groups
    nBGroup = int(f.readline());
    #read the boundaries
    for i in range(0,nBGroup):
        nBFace,nf,Title = f.readline().split();
        nBFace = int(nBFace);
        for j in range(0,nBFace):
            n_0, n_1 = f.readline().split();
    #read cell info
    nElem,Order,Basis = f.readline().split();
    nElem = int(nElem);
    E = [None]*nElem;
    for i in range(0,nElem):
        v_0,v_1,v_2 = [int(string) for string in f.readline().split()];
        # the given index start from 1, we want the index start from 0
        v_0 = v_0-1;
        v_1 = v_1-1;
        v_2 = v_2-1;
        vertex = [v_0, v_1, v_2];
        # vertex, edge, tri, A, adj_cell, state, R, dt
        E[i] = vertex;
    E = np.asarray(E);
    mesh = {'nNode':nNode, 'nElem':nElem, 'node_pos':node_pos, 'Elems':E};
    f.close();
    return mesh;

#post_processing
def plot_mach(mesh, case_no):
    gamma = 1.4;
    x = mesh['node_pos'][:,0];
    y = mesh['node_pos'][:,1];
    tri = mesh['Elems'];
    M = np.zeros(mesh['nElem']);
    f = open('data\\task_2_c\\bump%d\\bump%d_state_final.txt' %(case_no,case_no), "r");
    for i in range(0,mesh['nElem']):
        state = [float(string) for string in f.readline().split()];
        rho = state[0];
        u = state[1]/state[0];
        v = state[2]/state[0];
        E = state[3]/state[0];
        p = (gamma-1)*(rho*E-0.5*rho*(u**2+v**2));
        a = np.sqrt(gamma*p/rho);
        M[i] = np.sqrt(u**2+v**2)/a;
    f1 = plt.figure(figsize=(15,5));
    plt.tripcolor(x, y, tri, M, cmap=plt.cm.jet);
    plt.axis('equal');
    plt.colorbar();
    plt.title('Mach number, bump%d' %case_no);
    plt.savefig('figure\\bump%d_mach.pdf' %case_no, dpi=150);
    plt.close(f1);
    return 0;

def plot_history(case_no):
    f = open('data\\task_2_c\\bump%d\\bump%d_history.txt' %(case_no,case_no), "r");
    nTime_step = int(f.readline());
    L_inf = [None]*nTime_step;
    for i in range(0,nTime_step):
        L_inf[i] = float(f.readline());
    f1 = plt.figure(figsize=(8,6));
    plt.semilogy(L_inf,'k-');
    plt.xlabel('Iteration',fontsize =16);
    plt.ylabel('$L_{\infty}$ error',fontsize =16);
    plt.grid();
    plt.title('Residual norm history, bump%d' %case_no);
    plt.savefig('figure\\bump%d_L_inf_error.pdf' %case_no, dpi=150);
    plt.close(f1);
    return 0;

def plot_cp(case_no):
    f = open('data\\task_2_c\\bump%d\\bump%d_cp_distribution.txt' %(case_no,case_no), "r");
    data = f.readlines();
    global cp;
    cp = [];
    for i in range(0,len(data)):
        cp.append([float(string) for string in data[i].split()]);
    cp = sorted(cp, key=lambda my_list: my_list[0]);
    cp = np.asarray(cp);
    x = cp[:,0];
    y = cp[:,1];
    f1 = plt.figure(figsize=(15,5));
    plt.plot(x,-y,'k-');
    plt.xlabel("x",fontsize =16);
    plt.ylabel('Pressure coefficient $-c_p$',fontsize =16);
    plt.grid();
    plt.title('Pressure coefficient distribution, bump%d' %case_no);
    plt.savefig('figure\\bump%d_cp_distribution.pdf' %case_no, dpi=150);
    plt.close(f1);
    return 0;

def main():
    global mesh;
    for i in range(0,5):
        mesh = read_gri_reduced(i);
        plot_mach(mesh, i);
        plot_history(i);
        plot_history(i);
        plot_cp(i);
    return 0;

if __name__=="__main__":
    main()