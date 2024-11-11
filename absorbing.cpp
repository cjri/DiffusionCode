#include "absorbing.h"
#include <iostream>
#include <string>
#include <cmath>

//Code for modelling 1D diffusion with absorbing boundaries.  NB Need to fix the diffusion code so that it does this also.

double expfunction (int n, double Pi, double A, double K, double lambda, double t) {
    double t1=pow(n*Pi/A,2);
    t1=(t1*K)+lambda/pow(A,2);
    t1=t1*t;
    double e=exp(-t1);
    return e;
}

double sinfunction (int n, double x, double L, double Pi) {
    double s=sin(n*Pi*x/L);
    return s;
}

double diffn1 (int n, double K, double Pi, double A, double lambda, double x, double x_0, double t) {
    double val=sinfunction(n,x,A,Pi)*expfunction(n,Pi,A,K,lambda,t)*sinfunction(n,x_0,A,Pi);
    return val;
}

double finddiff (double K, double Pi, double A, double lambda, double x, double x_0, double t) {
    double tot=0;
    int small=0;
    for (int n=1;n<=100000;n++) {
        double d1=diffn1(n,K,Pi,A,lambda,x,x_0,t);
        if (abs(d1)<pow(10,-30)) {
            if (tot>0&&abs(d1)<abs(tot*pow(10,-10))) { //New criterion
                small++;
            }
        } else {
            small=0;
            tot=tot+d1;
        }
        if (small==10) {
            //cout << "N was " << n << "\n";
            break;
        }
    }
    if (tot<pow(10,-25)) {  //What is happening here?  Change to relative measure?
        tot=0;
    }
    return tot;
}

double diffn2 (double K, double Pi, double A, double B, double lambda, double x, double y, double x_0, double y_0, double t) {
    double xval=finddiff(K,Pi,A,lambda,x,x_0,t);
    double yval=finddiff(K,Pi,B,lambda,y,y_0,t);
    double val=xval*yval;
    return val;
}


/*
void FindGrids (double d, double Pi, double L, double x_0, double t, vector<double>& xdat, vector<double>& cdat) {
    double tx=0;
    double max=L/0.02;
    //cout << "Max " << max << "\n";
    for (int i=0;i<max+0.01;i++) {
        xdat.push_back(tx);
        tx=tx+0.02;
    }
    double tot=0;
    for (int i=0;i<xdat.size();i++) {
        double f=finddiff(d,Pi,L,xdat[i],x_0,t);
        cdat.push_back(f);
        tot=tot+f;
        //cout << xdat[i] << " " << f << "\n";
    }
    for (int i=0;i<cdat.size();i++) {
        cdat[i]=cdat[i]/tot;
    }

}

double FindMean (vector<double>& xdat, vector<double>& cdat) {
    //Find mean
    double mn=0;
    for (int i=0;i<cdat.size();i++) {
        mn=mn+xdat[i]*cdat[i];
    }
    return mn;
}

double FindSD (double mn, vector<double>& xdat, vector<double>& cdat){
    double sd=0;
    for (int i=0;i<xdat.size();i++) {
        double s=cdat[i]*pow(xdat[i]-mn,2);
        sd=sd+s;
    }
    sd=sd/xdat.size();
    sd=sqrt(sd);
    return sd;
}
*/
