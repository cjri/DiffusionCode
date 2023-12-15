#include "absorbing.h"
#include <iostream>
#include <string>
#include <cmath>

//Code for modelling 1D diffusion with absorbing boundaries.  NB Need to fix the diffusion code so that it does this also.

double dref (double x, double x_0, double y, double y_0, double k, double lambda, double t) {
    double d=pow(x-x_0,2)+pow(y-y_0,2);
    d=d/(1+(4*k*t));
    d=d+(lambda*t);
    d=exp(-d);
    d=d/(1+(4*k*t));
    return d;
}

double sumreflections (double A, double B, double reps, double x, double x_0, double y, double y_0, double k, double lambda, double t) { //Here A and B are the dimensions of the room in the x and y directions
    //reps is the number of repititions of the sum
    
    //First normalise the location in case there is anything weird here.  First, within 0<=x<2A; 0<=y<2B
    while (x>(2*A)) {
        x=x-(2*A);
    }
    while (y>(2*B)) {
        y=y-2*B;
    }
    while (x<0) {
        x=x+(2*A);
    }
    while (y<0) {
        y=y+2*B;
    }
    
    //Now, within 0<=x<A; 0<=y<B
    if (x>A) {
        x=(2*A)-x;
    }
    if (y>B) {
        y=(2*B)-y;
    }
    
    double out=0;
    for (int i=-reps;i<=reps;i++) {
        for (int j=-reps;j<=reps;j++) {
            out=out+dref((2*i*A)+x,x_0,(2*j*B)+y,y_0,k,lambda,t);
            out=out+dref((2*i*A)-x,x_0,(2*j*B)+y,y_0,k,lambda,t);
            out=out+dref((2*i*A)+x,x_0,(2*j*B)-y,y_0,k,lambda,t);
            out=out+dref((2*i*A)-x,x_0,(2*j*B)-y,y_0,k,lambda,t);
        }
    }
    if (out<pow(10,-25)) {
	    out=0;
    }
    return out;
}
