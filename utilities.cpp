#include "diffusion.h"
#include "absorbing.h"
#include "reflecting.h"
#include "io.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <iterator>
#include <random>

double FindK(double L, double gamma) {
    double K=(0.52*gamma+0.31)*pow(L,2);
    return K;
}

void FindEmissionMean (options& o, param& p, vector<dist> initial) {
    double xm=0;
    double ym=0;
    for (unsigned int i=0;i<initial.size();i++) {
        xm=xm+initial[i].x*initial[i].val;
        ym=ym+initial[i].y*initial[i].val;
    }
    if (abs(xm)<1e-10) {
        xm=0;
    }
    p.x_0=xm;
    p.y_0=ym;
    if (o.verb==1) {
        cout << "Emission point " << p.x_emit << " " << p.y_emit << "\n";
        cout << "Central position " << p.x_0 << " " << p.y_0 << "\n";
    }
}

double EvaporationModel (double& radius) {
    double r=radius/4;
    return r;
}

void FindLambda (options& o, param& p, double radius) {
    //Sedimentation
    FindSedimentationRate(p,radius);
    FindEvacuationRate (p,radius);
    SetDegradationRate(p,o);
    
    p.lambda=p.sed+p.evac+p.deg;
    //Correct lambda with A and B
    if (o.absorbing==1) {
        p.lambda=p.lambda*pow(p.X,2)*pow(p.Y,2)/(pow(p.X,2)+pow(p.Y,2));
    }
    if (o.verb==1) {
        cout << "Rates " << p.sed << " " << p.evac << " " << p.deg << "\n";
        cout << "Lambda' " << p.lambda << "\n";
    }
}


void FindSedimentationRate (param& p, double r) {
    double t_max = (0.85*pow(10,-8))*(p.z_0/pow(r,2)); //0.85e-8 factor obtained from experimental data (UPDATE - CHECK OVERLEAF DOC)
    p.sed=3600/t_max;
}

void FindEvacuationRate (param& p, double r) {
    double mu_a=1.81*pow(10,-5); //dynamic viscosity of air [kg/(m s)]
    double g=9.81; //gravitational acceleration
    double rho_air=1.204; //density of air [kg/m^3]
    double rho_water=997.77; //density of water [kg/m^3]
    double delta_rho=rho_water-rho_air;
    double r_c=sqrt((9*p.gamma*p.L*mu_a)/(2*g*delta_rho));
    double r_star=max(r,r_c);
    p.evac=p.gamma*pow(r_c/r_star,2); //Rate per hour
}

void SetDegradationRate (param& p, options& o) {
    p.deg=0.630134; //Rate per hour
    p.deg=p.deg*o.deg_shift;
}

void OptimiseInitialTAbsorbing (options& o, param& p, double& t, double Pi, vector<dist>& initial, gsl_rng *rgen){
    t=0.0001;
    //Remove points outside the room
    vector<dist> new_initial;
    for (unsigned int i=0;i<initial.size();i++) {
        if (initial[i].x>0&&initial[i].y>0&&initial[i].x<=p.X&&initial[i].y<=p.Y) {
            new_initial.push_back(initial[i]);
        }
    }
    initial=new_initial;
    cout << "Size " << initial.size() << "\n";
    //Remake initial distribution to keep everything in the boundaries
    double delta=0.01;
    double t_best=0.01;
    double diff=1000;
    double diff_best=1000;
    double t_previous=0.01;
    int it_best=0;
    int first=1;
    int first_change=1;
    vector<dist> model=initial;
    for (unsigned int it=0;it<1000;it++) {
        if (first==0) {
            if (diff<diff_best) {
                if (first_change==0) {
                    t_previous=t_best;
                }
                it_best=it;
                diff_best=diff;
                t_best=t;
                if (o.verb==1) {
                    cout << "Iteration " << it << " New difference " << diff << " " << t << "\n";
                }
                first_change=0;
            } else {
                t=t_best;
            }
        }
        first=0;
        if (t_best<delta/10) {
            delta=t_best/10;
        }
        if (it-it_best==50) {
            double change=abs(t_best-t_previous);
            //cout << "Change delta " << delta << " " << change/2 << "\n";
            delta=change/2;
        }
        if (delta<pow(10,-10)) {
            break;
        }
        
        t=t+(gsl_rng_uniform(rgen)*delta)-(delta/2);
        if (t<0) {
            t=-t;
        }
        //Generate values on the initial emissions grid
        double tot=0;
        for (unsigned int i=0;i<model.size();i++) {
            model[i].val=diffn2(p.K,Pi,p.X,p.Y,p.lambda,model[i].x,model[i].y,p.x_0,p.y_0,t);
            tot=tot+model[i].val;
        }
        for (unsigned int i=0;i<model.size();i++) {
            model[i].val=model[i].val/tot;
        }
        diff=0;
        for (unsigned int i=0;i<model.size();i++) {
            diff=diff+sqrt(pow(model[i].val-initial[i].val,2));
        }
    }
    t=t_best;
    cout << "Best t " << t_best << " at " << diff_best << "\n";
}


void OptimiseInitialTReflecting (options& o, param& p, double& t, double Pi, vector<dist>& initial, gsl_rng *rgen){
    t=0.01;
    double min_t=-1/(4*p.K);
    //Remove points outside the room
    vector<dist> new_initial;
    for (unsigned int i=0;i<initial.size();i++) {
        if (initial[i].x>0&&initial[i].y>0&&initial[i].x<=p.X&&initial[i].y<=p.Y) {
            new_initial.push_back(initial[i]);
        }
    }
    initial=new_initial;
    cout << "Size " << initial.size() << "\n";
    double delta=0.01;
    double t_best=0.01;
    double diff=1000;
    double diff_best=1000;
    double t_previous=0.01;
    int it_best=0;
    int first=1;
    int first_change=1;
    vector<dist> model=initial;
    for (unsigned int it=0;it<100000;it++) {
        if (first==0) {
            if (diff<diff_best) {
                if (first_change==0) {
                    t_previous=t_best;
                }
                it_best=it;
                diff_best=diff;
                t_best=t;
                if (o.verb==1) {
                    cout << "Iteration " << it << " New difference " << diff << " " << t << "\n";
                }
                first_change=0;
            } else {
                t=t_best;
            }
        }
        first=0;
        if (abs(t_best)<delta/10) {
            delta=abs(t_best)/10;
        }
        if (it-it_best==50) {
            double change=abs(t_best-t_previous);
            //cout << "Change delta " << delta << " " << change/2 << "\n";
            delta=change/2;
        }
        if (delta<pow(10,-10)) {
            break;
        }
        if (delta<pow(10,-3)&&it%100==0) {
            delta=delta*10;
            t=t+(gsl_rng_uniform(rgen)*delta)-(delta/2);
            delta=delta/10;
        } else {
            t=t+(gsl_rng_uniform(rgen)*delta)-(delta/2);
        }
        if (t<min_t) {
            t=min_t-t;
        }
        //Generate values on the initial emissions grid
        double tot=0;
        for (unsigned int i=0;i<model.size();i++) {
            model[i].val=sumreflections(p.X,p.Y,5,model[i].x,p.x_0,model[i].y,p.y_0,p.K,p.lambda,t);
            tot=tot+model[i].val;
        }
        for (unsigned int i=0;i<model.size();i++) {
            model[i].val=model[i].val/tot;
        }
        diff=0;
        for (unsigned int i=0;i<model.size();i++) {
            diff=diff+sqrt(pow(model[i].val-initial[i].val,2));
        }
    }
    t=t_best;
    cout << "Best t " << t_best << " at " << diff_best << "\n";
}



void CalculateTotalExposureAbsorbing(options& o, param& p, double Pi, double& total_exp_init) {
    double max=0;
    vector< vector<double> > initial_state;
    for (double dx=0;dx<=p.X;dx=dx+0.02) {
        for (double dy=0;dy<=p.Y;dy=dy+0.02) {
            vector<double> d;
            double increment=diffn2(p.K,Pi,p.X,p.Y,p.lambda,dx,dy,p.x_0,p.y_0,p.t_0);
            d.push_back(dx);
            d.push_back(dy);
            d.push_back(increment);
            initial_state.push_back(d);
            total_exp_init=total_exp_init+increment;
            if (increment>max) {
                max=increment;
            }
        }
    }
    if (o.verb==1) {
        cout << "Total exposure at time t_0 " << total_exp_init << "\n";
        ExportInitialState(o,initial_state);
    }
}

void CalculateTotalExposureReflecting(options& o, param& p, double Pi, double& total_exp_init) {
    double max=0;
    vector< vector<double> > initial_state;
    for (double dx=0;dx<=p.X;dx=dx+0.02) {
        for (double dy=0;dy<=p.Y;dy=dy+0.02) {
            vector<double> d;
            double increment=sumreflections(p.X,p.Y,5,dx,p.x_0,dy,p.y_0,p.K,p.lambda,p.t_0);
            d.push_back(dx);
            d.push_back(dy);
            d.push_back(increment);
            initial_state.push_back(d);
            total_exp_init=total_exp_init+increment;
            if (increment>max) {
                max=increment;
            }
        }
    }
    if (o.verb==1) {
        cout << "Total exposure at time t_0 " << total_exp_init << "\n";
        ExportInitialState(o,initial_state);
    }
}


void MakeLocationGrids (options& o, param& p, vector<loc>& locations, vector< vector<loc> >& location_grids) {
    for (unsigned int i=0;i<locations.size();i++) {
        vector<loc> lg;
        loc l;
        for (double dx=-0.2;dx<=0.2;dx=dx+0.02) {
            for (double dy=-0.2;dy<=0.2;dy=dy+0.02) {
                l.x=locations[i].x+dx;
                l.y=locations[i].y+dy;
                lg.push_back(l);
            }
        }
        location_grids.push_back(lg);
    }
    if (o.verb==1) {
        cout << "Individual locations\n";
        for (unsigned int i=0;i<locations.size();i++) {
            cout << locations[i].x << " " << locations[i].y << " " << p.X << " " << p.Y << "\n";
        }
        cout << "\n";
    }
}

void CalcluateExposuresAbsorbing(param& p, double Pi, const double total_exp_init, const vector< vector<loc> >& location_grids, vector<double>& times, vector< vector<double> >& exp_record) {
    double delta_t=1./3600.;
    vector<int> is0;        
    for (unsigned int i=0;i<location_grids.size();i++) {
        is0.push_back(0);
    }
    for (double k=1;k<=3600;k++) {
        times.push_back(k);
        vector<double> exps;
        for (unsigned int i=0;i<location_grids.size();i++) {
            if (is0[i]==0) {
                double tot=0;
                double count=0;
                for (unsigned int j=0;j<location_grids[i].size();j++) {
                    tot=tot+diffn2(p.K,Pi,p.X,p.Y,p.lambda,location_grids[i][j].x,location_grids[i][j].y,p.x_0,p.y_0,p.t_0+(k*delta_t));
                    count++;
                }
                tot=tot/count;
                if (tot==0&&k>50) {
                    is0[i]=1;
                    //cout << "Found 0 at person " << i << "\n";
                }
                exps.push_back(tot/total_exp_init);
            } else {
                exps.push_back(0);
            }
        }
        exp_record.push_back(exps);
    }
}

void CalcluateExposuresAbsorbingTest(options& o, param& p, double Pi, const double total_exp_init) {
    double delta_t=1./3600.;
    for (double k=1;k<=3600;k=k+60) {
        vector<dist> exps;
        for (double x=0;x<p.X;x=x+0.02) {
            for (double y=0;y<p.Y;y=y+0.02) {
                dist d;
                d.x=x;
                d.y=y;
                d.val=diffn2(p.K,Pi,p.X,p.Y,p.lambda,x,y,p.x_0,p.y_0,p.t_0+(k*delta_t));
                d.val=d.val/total_exp_init;
                exps.push_back(d);
            }
        }
        //Export exps
        ExportTestExposures (o,k,exps);
    }
}


void CalcluateExposuresReflecting(options& o, param& p, double Pi, const double total_exp_init, const vector< vector<loc> >& location_grids, vector<double>& times, vector< vector<double> >& exp_record) {
    double delta_t=1./3600.;
    vector<int> is0;
    for (unsigned int i=0;i<location_grids.size();i++) {
        is0.push_back(0);
    }
    for (double k=1;k<=3600;k++) {
        times.push_back(k);
        vector<double> exps;
        for (unsigned int i=0;i<location_grids.size();i++) {
            if (is0[i]==0) {
                double tot=0;
                double count=0;
                for (unsigned int j=0;j<location_grids[i].size();j++) {
                    tot=tot+sumreflections(p.X,p.Y,5,location_grids[i][j].x,p.x_0,location_grids[i][j].y,p.y_0,p.K,p.lambda,p.t_0+(k*delta_t));
		    count++;
                }
                tot=tot/count;
                if (tot==0&&k>50) {
                    is0[i]=1;
                    //cout << "Found 0 at person " << i << "\n";
                }
		double rr=tot/total_exp_init;
                exps.push_back(rr);
		if (i==1&&o.verb==1) {
			cout << k << " " << i << " " << tot << " " << total_exp_init << "\n";
		}
            } else {
                exps.push_back(0);
            }
        }
        exp_record.push_back(exps);
    }
}

void CalcluateExposuresReflectingTest(options& o, param& p, double Pi, const double total_exp_init) {
    double delta_t=1./3600.;
    for (double k=1;k<=3600;k=k+60) {
        vector<dist> exps;
        for (double x=0;x<p.X;x=x+0.02) {
            for (double y=0;y<p.Y;y=y+0.02) {
                dist d;
                d.x=x;
                d.y=y;
                d.val=sumreflections(p.X,p.Y,5,x,p.x_0,y,p.y_0,p.K,p.lambda,p.t_0+(k*delta_t));
                d.val=d.val/total_exp_init;
                exps.push_back(d);
            }
        }
        //Export exps
        ExportTestExposures (o,k,exps);
    }
}


void CalculateVolumeDistribution (double Pi, vector<double>& sizedist) {
    double tot=0;
    for (unsigned int i=0;i<sizedist.size();i++) {
        sizedist[i]=sizedist[i]*(4/3)*Pi*pow(i+1,3);
        tot=tot+sizedist[i];
    }
    for (unsigned int i=0;i<sizedist.size();i++) {
        sizedist[i]=sizedist[i]/tot;
    }
}

void NormaliseExposureToVolume (options& o, const vector<double>& sizedist, vector< vector< vector<double> > >& exposure_dat) {
    for (int r=0;r<o.maxr;r++) {
        for (unsigned int i=0;i<exposure_dat[r].size();i++) {
            for (unsigned int j=0;j<exposure_dat[r][i].size();j++) {
                exposure_dat[r][i][j]=exposure_dat[r][i][j]*sizedist[r]*o.volume;
            }
        }
    }
}

void CalculateExposures (options& o, param& p, const vector<loc>& locations, const vector< vector< vector<double> > >& exposure_dat, vector<double>& exp_tot, vector< vector<double> >& exp_r_tot) {
    //Set up index vector
    for (unsigned int n=0;n<locations.size();n++) {
	if (o.verb==1) {
		cout << "Individual " << n << "\n";
	}
    	vector<double> er;
        for (int r=0;r<o.maxr;r++) {
            er.push_back(0);
        }
        //Calcualte exposures by radius for this individual
        //Have the total for a standard cough, measured over one hour
        for (int r=0;r<o.maxr;r++) {
            for (unsigned int t=0;t<exposure_dat[r][n].size();t++) {
                er[r]=er[r]+exposure_dat[r][n][t];
		if (exposure_dat[r][n][t]<0) {
			cout << "Negative at " << r << " " << n << " " << t << "\n";
		}
            }
        }
        for (int r=0;r<o.maxr;r++) {
            er[r]=er[r]*o.exp_per_hour*(p.hours-1);
        }
  	//Here we calculate the exposure in the final hour.  From zero to some end point.
        int step=3600/o.exp_per_hour;
      	for (int r=0;r<o.maxr;r++) {
            for (unsigned int i=0;i<o.exp_per_hour;i++) {
                //cout << "Step " << i << " limit " << exposure_dat[r][n].size()-(i*step) << "\n";
                for (unsigned int t=0;t<exposure_dat[r][n].size()-(i*step);t++) {
                    er[r]=er[r]+exposure_dat[r][n][t];
                }
            }
        }
        exp_r_tot.push_back(er);
        double e_tot=0;
        for (int r=0;r<o.maxr;r++) {
            e_tot=e_tot+er[r];
        }
        exp_tot.push_back(e_tot);
        if (o.verb==1) {
            cout << "Individual " << n << "\n";
            cout << "Total " << e_tot << "\n";
        }
    }
}

void CalculateExposuresMovement (options& o, param& p, const vector<loc>& locations, const vector< vector< vector<double> > >& exposure_dat, vector<double>& exp_tot, vector< vector<double> >& exp_r_tot) {
    //Set up index vector
    vector<int> index;
    for (unsigned int n=0;n<locations.size();n++) {
        index.push_back(n);
    }
    
    for (unsigned int n=0;n<locations.size();n++) {
        if (o.verb==1) {
            cout << "Individual " << n << "\n";
        }
        vector<double> er;
        for (int r=0;r<o.maxr;r++) {
            er.push_back(0);
        }
        
        //Calculate exposures by radius for this individual
        //Have the total for a standard cough, measured over one hour
        for (int r=0;r<o.maxr;r++) {//Radius
            for (int h=0;h<p.hours-1;h++) {//Explicit loop over multiple hours
                for (unsigned int t=0;t<exposure_dat[r][n].size();t++) {
                    if (t>0&&t%300==0) { //Shuffle index every five minutes
                        IndexShuffle(index);
                    }
                    er[r]=er[r]+exposure_dat[r][index[n]][t];
                }
            }
        }
        for (int r=0;r<o.maxr;r++) {
            er[r]=er[r]*o.exp_per_hour;
        }
      //Here we calculate the exposure in the final hour.  From zero to some end point.
        int step=3600/o.exp_per_hour;
          for (int r=0;r<o.maxr;r++) {
            for (unsigned int i=0;i<o.exp_per_hour;i++) {
                //cout << "Step " << i << " limit " << exposure_dat[r][n].size()-(i*step) << "\n";
                for (unsigned int t=0;t<exposure_dat[r][n].size()-(i*step);t++) {
                    er[r]=er[r]+exposure_dat[r][n][t];
                }
            }
        }
        exp_r_tot.push_back(er);
        double e_tot=0;
        for (int r=0;r<o.maxr;r++) {
            e_tot=e_tot+er[r];
        }
        exp_tot.push_back(e_tot);
        if (o.verb==1) {
            cout << "Individual " << n << "\n";
            cout << "Total " << e_tot << "\n";
        }
    }
}

void IndexShuffle(vector<int>& index) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(index.begin(), index.end(), g);
    /*std::copy(index.begin(), index.end(), std::ostream_iterator<int>(std::cout, " "));
    for (int i=0;i<index.size();i++) {
        cout << index[i] << " ";
    }
    cout << "\n";*/
}

void FindNP(options& o, double Pi, const vector< vector<double> >& exp_r_tot, vector< vector<double> >& np) {
    np=exp_r_tot;
    for (unsigned int i=0;i<np.size();i++) {
        for (int r=0;r<o.maxr;r++) {
            np[i][r]=np[i][r]/((4./3.)*Pi*pow((r+1)*pow(10,-5),3));
        }
    }
}

void FindEpp(options& o, double k, double Pi, vector<double>& epp) {
    for (int r=0;r<o.maxr;r++) {
        epp.push_back((4./3.)*Pi*k*pow((r+1)*pow(10,-5),3));
    }
}

void OptimiseViralLoad (options& o, double renv, const vector<double>& epp, const vector< vector<double> > np, double& k, gsl_rng *rgen) {
    double diff=100;
    double diff_best=100;
    int it_best=0;
    double k_best=pow(10,12);
    double k_previous=pow(10,12);
    double dk=k/50;
    int first=1;
    int first_change=1;
    for (unsigned int it=0;it<10000;it++) {
        if (first==0) {
            if (diff<diff_best) {
                if (first_change==0) {
                    k_previous=k_best;
                }
                diff_best=diff;
                k_best=k;
                first_change=0;
                it_best=it;
            } else {
                k=k_best;
            }
        }
        first=0;
        
        if (it-it_best==50) {
            double change=abs(k_best-k_previous);
            dk=change/2;
        }
        if (dk<1) {
            cout << "Check opt " << diff_best << "\n";
            break;
        }
        k=k+(gsl_rng_uniform(rgen)*dk)-(dk/2);
        
        //Find difference between formula and renv
        double tprodv=calcr(o,k,epp,np);
        diff=abs(tprodv-renv);
    }
}

double calcr (options& o, double k, const vector<double>& epp, const vector< vector<double> >& np) {
    double tprodv=0;
    for (unsigned int i=0;i<np.size();i++) {
        double prodv=1;
        for (int r=0;r<o.maxr;r++) {
            //cout << "r " << r << " " << np[i][r] << " " << epp[r] << "\n";
            double val=gsl_ran_poisson_pdf(0,np[i][r]);
            //cout << "Val " << val << "\n";
            for (unsigned int d=1;d<=10;d++) {
                double v=gsl_ran_poisson_pdf(d,np[i][r])*gsl_ran_poisson_pdf(0,d*k*epp[r]);
                val=val+v;
                //cout << "d " << d << " " << gsl_ran_poisson_pdf(d,np[i][r]) << " " << gsl_ran_poisson_pdf(0,d*k*epp[r]) << " " << v << " " << val << "\n";
            }
            //cout << "Previous " << prodv << "\n";
            prodv=prodv*val;
            //cout << "Here " << r << " " << val << " " << prodv << "\n";
        }
        prodv=1-prodv;
        tprodv=tprodv+prodv;
    }
    return tprodv;
}

void CalculateBottlenecks (options& o, vector< vector<double> >& np, vector<double>& epp, vector<loc>& locations, vector<int>& bottlenecks, gsl_rng *rgen) {
    //Pre-calculation of infectivity distribution
    vector<double> id;
    if (o.variable_emissions==1) {
        ifstream gam_file;
        gam_file.open("Gamma"+o.alpha+".dat");
        id.push_back(0);
        double idx;
        for (int i=1;i<1000;i++) {
            gam_file >> idx;
            id.push_back(idx);
            if (o.verb==1) {
                cout << "ID " << id[i] << "\n";
            }
        }
    } else {
        for (int i=0;i<1000;i++) {
            id.push_back(1);
        }
    }
    vector< vector<int> > draws;
    vector<int> d;
    for (int r=0;r<o.maxr;r++) {
        d.push_back(0);
    }
    for (unsigned int i=0;i<locations.size();i++) {
        draws.push_back(d);
    }
    for (unsigned int draw=1;draw<=1000000;draw++) {
        int idpos=draw%1000;
        for (unsigned int i=0;i<locations.size() ;i++) {
            for (int r=0;r<o.maxr;r++) {
                draws[i][r]=gsl_ran_poisson(rgen,id[idpos]*np[i][r]);
                //cout << draws[i][r] << " ";
            }
        }
        //Number of viruses per particle
        for (unsigned int i=0;i<locations.size();i++) {
            int size=0;
            for (int r=0;r<o.maxr;r++) {
                if (draws[i][r]>0) {
                    draws[i][r]=gsl_ran_poisson(rgen,draws[i][r]*epp[r]);
                    //cout << i << " " << r << " " << epp[r] << " " << draws[i][r] << "\n";
                    if (draws[i][r]>0) {
                        size=size+draws[i][r];
                    }
                }
            }
            if (size>0) {
                bottlenecks.push_back(size);
            }
        }
    }
    sort(bottlenecks.begin(),bottlenecks.end());
}

void CalculateBottlenecksLimit (options& o, vector< vector<double> >& np, vector<double>& epp, vector<loc>& locations, vector<int>& bottlenecks, gsl_rng *rgen) {
    //Version of the code which simulates a very high viral load: Exposure to a particle deterministically causes infection
    vector< vector<int> > draws;
    vector<int> d;
    for (int r=0;r<o.maxr;r++) {
        d.push_back(0);
    }
    for (unsigned int i=0;i<locations.size();i++) {
        draws.push_back(d);
    }
    for (unsigned int draw=1;draw<=1000000;draw++) {
        for (unsigned int i=0;i<locations.size() ;i++) {
            for (int r=0;r<o.maxr;r++) {
                draws[i][r]=gsl_ran_poisson(rgen,np[i][r]);
                //cout << draws[i][r] << " ";
            }
        }
        //Number of viruses per particle
        for (unsigned int i=0;i<locations.size();i++) {
            int size=0;
            for (int r=0;r<o.maxr;r++) {
                if (draws[i][r]>0) {
                    size++;
                }
            }
            if (size>0) {
                bottlenecks.push_back(size);
            }
        }
    }
    sort(bottlenecks.begin(),bottlenecks.end());
}

