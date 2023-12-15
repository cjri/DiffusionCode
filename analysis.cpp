#include "diffusion.h"
#include "analysis.h"
#include "utilities.h"
#include "data.h"
#include "io.h"
#include <iostream>
#include <string>

//Analysis of pre-existing exposure values.  Read in and compile to generate total exposures

void RunAnalysis (options& o, param& p, double Pi, gsl_rng *rgen) {
    //Read in locations
    vector<loc> locations;
    ReadLocations(o,p,locations);
    //Read in exposures
    vector<double> times;
    vector< vector< vector<double> > > exposure_dat;
    ReadExposures (o,locations,times,exposure_dat);
    
    //We now want:
    
    //Find distribution of the proportion of each size of particle
    vector<double> sizedist;
    GetSizeDistribution (o,Pi,sizedist);
    if (o.verb==1) {
        cout << "Size dist\n";
        for (unsigned int i=0;i<sizedist.size();i++) {
            cout << i << " " << sizedist[i] << "\n";
        }
    }

    CalculateVolumeDistribution (Pi,sizedist);
    //2. Normalise to a total volume at emission
    //How we do this: The exposures are proportions of an initial volume.  We multiply them by an absolute volume.
    NormaliseExposureToVolume (o,sizedist,exposure_dat);

    //3.  Find exposures for each person based on time and data
    vector<double> exp_tot;  //Total exposures in litres
    vector< vector<double> > exp_r_tot; //Exposures by r in litres
    if (o.move==0) {
        CalculateExposures (o,p,locations,exposure_dat,exp_tot,exp_r_tot);
    } else {
        CalculateExposuresMovement (o,p,locations,exposure_dat,exp_tot,exp_r_tot);
    }
    
    //4.  Identify effective viral load
    //this requires an optimisation.
    //Information is in exp_tot vector.
    //NB this part of the code will be optional - we might also want to feed in a phi_office vlaue.
    
    //Optimisation of viral load if not specified
    if (o.viral_load==-1){
        //Setup here is default for office
        double renv=o.r0*(1./36.)*o.phi;
        
        cout << o.r0 << " " << o.phi << " " << o.r0*(1./36.)*o.phi << " Renv " << renv << "\n";

        
        //Find vector np.  Number of particles of given size
        vector< vector<double> > np;
        FindNP(o,Pi,exp_r_tot,np);

        //Find particle sizes epp
        vector<double> epp;
        FindEpp(o,1,Pi,epp);
        
        double k=pow(10,12);
        OptimiseViralLoad (o,renv,epp,np,k,rgen);
        cout << "Optimal k " << k << "\n";
        o.viral_load=k;
    }
    
    //Calculate samples from the random distribution with optimised k.
    cout << "Viral load " << o.viral_load << "\n";
    cout << "Locations size " << locations.size() << "\n";
    double k=o.viral_load;
    vector< vector<double> > np;

    FindNP(o,Pi,exp_r_tot,np);
    if (o.verb==1) {
        cout << "NP\n";
        for (unsigned int i=0;i<np.size();i++) {
            for (unsigned int j=0;j<np[i].size();j++) {
                cout << np[i][j] << " ";
            }
            cout << "\n";
        }
    }
    vector<double> epp;
    FindEpp(o,k,Pi,epp);
    if (o.verb==1) {
        cout << "Epp\n";
        for (unsigned int i=0;i<epp.size();i++) {
            cout << i << " " << epp[i] << "\n";
        }
    }
    //Construct set of random draws of each size of particle
    //Number of particles
    vector<int> bottlenecks;
    cout << "Calculating bottlenecks...\n";
    if (o.viral_load==1e+16){
	    cout << "Viral load limit\n";
        CalculateBottlenecksLimit (o,np,epp,locations,bottlenecks,rgen);
    } else {
        CalculateBottlenecks (o,np,epp,locations,bottlenecks,rgen);
    }
    //Export bottlenecks
    ExportBottlenecks(o,bottlenecks);
    
    
}

