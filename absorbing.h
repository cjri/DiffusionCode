
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>


using namespace std;

double expfunction (int n, double Pi, double A, double K, double lambda, double t);
double sinfunction (int n, double x, double L, double Pi);
double diffn1 (int n, double K, double Pi, double A, double lambda, double x, double x_0, double t);
double finddiff (double K, double Pi, double A, double lambda, double x, double x_0, double t);
double diffn2 (double K, double Pi, double A, double B, double lambda, double x, double y, double x_0, double y_0, double t);




void FindGrids (double d, double Pi, double L, double x_0, double t, vector<double>& xdat, vector<double>& cdat);
double FindMean (vector<double>& xdat, vector<double>& cdat);
double FindSD (double mn, vector<double>& xdat, vector<double>& cdat);





