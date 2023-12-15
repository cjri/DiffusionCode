#include <iostream>

void GetCMLOptions (options& o, int argc, const char **argv);
void ReadParameters(options& o, param& p);
void ReadLocations(options& o, param& p, vector<loc>& locations);
void ReadCough(options& o, param& p, vector<dist>& initial);
void ExportInitialState (options& o, vector< vector<double> >& initial_state);
void GetMethod (options& o, string& method);
void ExportExposures (options& o, vector<double>& times, vector< vector<double> >& exp_record);

void ExportTestExposures (options& o, double k, vector<dist>& exps);


//Analysis
void ReadExposures (options& o, const vector<loc>& locations, vector<double>& times, vector< vector< vector<double> > >& exposure_dat);
void SetupExposureData (options& o, const vector<loc>& locations, vector< vector< vector<double> > >& exposure_dat);
void ExportBottlenecks(options& o, vector<int>& bottlenecks);
