# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 14:37:40 2019

@author: runda
"""

import numpy as np
import matplotlib.pyplot as plt

cl_exact = 1.537095;
cd_exact = 2.94278 * 10 ** (-6);
Es_exact = 0.0;
    
def find_error(fname):
    global temp;
    f = open(fname, "r");
    f.readline();
    case_no = 4;
    DOF = np.zeros(case_no);
    cl_err = np.zeros(case_no);
    cd_err = np.zeros(case_no);
    Es_err = np.zeros(case_no);
    for i in range(0,case_no):
        title,DOF[i],cl,cd,Es = f.readline().split();
        DOF[i] = np.sqrt(int(DOF[i]));
        cl_err[i] = abs(float(cl) - cl_exact);
        cd_err[i] = abs(float(cd) - cd_exact);
        Es_err[i] = abs(float(Es) - Es_exact);
    err = {'DOF':DOF,'e_cl':cl_err,'e_cd':cd_err,'e_Es':Es_err};
    return err;

def plot_error(err_1,err_2):
    f1 = plt.figure(figsize=(10,8));
    plt.loglog(err_1['DOF'],err_1['e_cl'],'r-*', label='$1^{st}$ order, $c_l$');
    plt.loglog(err_2['DOF'],err_2['e_cl'],'r--x', label='$2^{nd}$ order, $c_l$');
    plt.loglog(err_1['DOF'],err_1['e_cd'],'g-*',label='$1^{st}$ order, $c_d$');
    plt.loglog(err_2['DOF'],err_2['e_cd'],'g--x',label='$2^{nd}$ order, $c_d$');
    plt.loglog(err_1['DOF'],err_1['e_Es'],'b-*',label='$1^{st}$ order, $E_s$');
    plt.loglog(err_2['DOF'],err_2['e_Es'],'b--x',label='$2^{nd}$ order, $E_s$');
    plt.legend(handlelength=5,loc='upper left',bbox_to_anchor=(1.0, 1.0));
    plt.xlabel('$\sqrt{dof}$',fontsize =20);
    plt.ylabel('Error',fontsize = 20);
    plt.grid();
    plt.tight_layout();
    plt.savefig('conv_study.pdf', dpi=150);
    plt.close(f1);
    return 0;

def func(x, a, b):
    return a + b*x;

def find_conv_rate(err):
    y = [None]*3;
    y[0] = err['e_cl'];
    y[1] = err['e_cd'];
    y[2] = err['e_Es'];
    for i in range(0,3):
        rate = np.abs(np.log2(y[i][3]/y[i][2]));
        print(rate);
    return 0;

def main():
    global err_1st, err_2nd;
    err_1st = find_error('cl_cd_Es_1ST.txt');
    err_2nd = find_error('cl_cd_Es_2ND.txt');
    plot_error(err_1st,err_2nd);
    find_conv_rate(err_1st);
    find_conv_rate(err_2nd);
    return 0;

if __name__=="__main__":
    main()