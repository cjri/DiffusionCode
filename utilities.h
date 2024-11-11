#include <iostream>


double FindK(double L, double gamma);
void FindEmissionMean (options& o, param& p, vector<dist> initial);
void FindLambda (options& o, param& p, double radius);
double EvaporationModel (double& radius);
void FindSedimentationRate (param& p, double r);
void FindEvacuationRate (param& p, double r);
void SetDegradationRate (param& p, options& o);
void OptimiseInitialTAbsorbing (options& o, param& p, double& t, double Pi, vector<dist>& initial, gsl_rng *rgen);
void OptimiseInitialTReflecting (options& o, param& p, double& t, double Pi, vector<dist>& initial, gsl_rng *rgen);
void CalculateTotalExposureAbsorbing(options& o, param& p, double Pi, double& total_exp_init);
void CalculateTotalExposureReflecting(options& o, param& p, double Pi, double& total_exp_init);
void MakeLocationGrids (options& o, param& p, vector<loc>& locations, vector< vector<loc> >& location_grids);
void CalcluateExposuresAbsorbing(options& o, param& p, double Pi, const double total_exp_init, const vector< vector<loc> >& location_grids, vector<double>& times, vector< vector<double> >& exp_record);
void CalcluateExposuresAbsorbingTest(options& o, param& p, double Pi, const double total_exp_init);
void CalcluateExposuresReflecting(options& o, param& p, double Pi, const double total_exp_init, const vector< vector<loc> >& location_grids, vector<double>& times, vector< vector<double> >& exp_record);
void CalcluateExposuresReflectingTest(options& o, param& p, double Pi, const double total_exp_init);


//Analysis
void CalculateVolumeDistribution (double Pi, vector<double>& sizedist);
void NormaliseExposureToVolume (options& o, const vector<double>& sizedist, vector< vector< vector<double> > >& exposure_dat);
void CalculateExposures (options& o, param& p, const vector<loc>& locations, const vector< vector< vector<double> > >& exposure_dat, vector<double>& exp_tot, vector< vector<double> >& exp_r_tot);
void CalculateExposuresMovement (options& o, param& p, const vector<loc>& locations, const vector< vector< vector<double> > >& exposure_dat, vector<double>& exp_tot, vector< vector<double> >& exp_r_tot);
void IndexShuffle(vector<int>& index);
void FindNP(options& o, double Pi, const vector< vector<double> >& exp_r_tot, vector< vector<double> >& np);
void FindEpp(options& o, double k, double Pi, vector<double>& epp);
void OptimiseViralLoad (options& o, double renv, const vector<double>& epp, const vector< vector<double> > np, double& k, gsl_rng *rgen);
double calcr (options& o, double k, const vector<double>& epp, const vector< vector<double> >& np);
void CalculateBottlenecks (options& o, vector< vector<double> >& np, vector<double>& epp, vector<loc>& locations, vector<int>& bottlenecks, gsl_rng *rgen);
void CalculateBottlenecksLimit (options& o, vector< vector<double> >& np, vector<double>& epp, vector<loc>& locations, vector<int>& bottlenecks, gsl_rng *rgen);

