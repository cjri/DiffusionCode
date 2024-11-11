
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>


using namespace std;


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>

struct options { //Options for running the code
    int seed;
    string method; //What we are doing here: Generate diffusion patterns, analyse data
    string env; //Name of environment e.g. Office
    string emm; //Emission type
    int absorbing; //Flag for absorbing boundary conditions versus reflective
    int opt_t0; //Flag to generate initial times t_0.  Alternative is to read them in
    double radius; //Radius of particle in microns;
    int verb; //Verbose output
    int maxr; //Maximum radius - depends on the emission type.
    double volume; //Volume of liquid in one emission event
    double viral_load; //Effective viral load
    int variable_emissions; //Flag for modelling a population who emit different amounts of virus
    string alpha; //Parameterisation for the gamma distribution of infectivity in a population.  As (a, 1/a)
    int move; //Moving around.  Option for nightclub
    double specify_k; //Option to specify the diffusion coefficient to equal a specific value
    int lowres; //Flag to use 20cm grid not 2cm
    int shorrt;
    double exp_per_hour; //Number of coughs or events per hour;
    double phi; //Value of phi_office.
    double vmult; //Multiplier for the volume when calcualting bottlenecks
    double r0; //Value of R_0.
    int test_output; //Code to generate outputs on a 2cm grid at one minute intervals.  Overrides standard calculation
    double deg_shift; //Modifier to rate of viral degradation in particles
};

struct param {
    double X;
    double Y;  //Note: Here Y is a horizontal coordinate while Z is the vertical direction
    double Z;
    double gamma; //Number of air changes per hour.  Contrast to use elsewhere?
    double x_emit; //Location of infected person
    double y_emit; //Location of infected person
    double x_0; //Centre of emission
    double y_0; //Centre of emission
    double L; //Characteristic length
    double K; //Diffusion constant
    double z_0; //Height at which people were sitting or standing
    double sed; //Sedimentation rate
    double evac; //Evacuation rate
    double deg; //Degradation rate
    double lambda; //Combination of -c parameters
    double t_0; //Correction for time: Begin at this value.
    double hours; //Number of hours of exposure.
    double inh; //Litres inhaled per minute.  Read this from Parameters
};

struct loc {
    double x;
    double y;
};

struct dist {
    double x;
    double y;
    double val;
};
