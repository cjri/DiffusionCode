#include "diffusion.h"
#include "io.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <sstream>

void GetCMLOptions (options& o, int argc, const char **argv) {
    string p_switch;
    o.seed=(int) time(NULL);
    o.radius=pow(10,-5);
    o.opt_t0=0;
    o.absorbing=1;
    o.env="Office";
    o.emm="Cough";
    o.method="Generation";
    o.maxr=500;
    o.volume=38*pow(10,-12); //Liquid in cough
    o.exp_per_hour=10;
    o.viral_load=-1;
    o.phi=1;
    o.r0=2.5;
    o.vmult=1;
    o.specify_k=-1;
    o.test_output=0;
    o.verb=0;
    o.variable_emissions=0;
    o.alpha="";
    o.deg_shift=1;
    o.move=0; //Option for the nightclub
    int x=1;
    while (x < argc && (argv[x][0]=='-')) {
        p_switch=argv[x];
        if (p_switch.compare("--radius")==0) {
            x++;
            o.radius=atof(argv[x]);
            o.radius=o.radius*pow(10,-6); //Convert to microns
        } else if (p_switch.compare("--env")==0) {
            x++;
            o.env=argv[x];
        } else if (p_switch.compare("--emm")==0) {
            x++;
            o.emm=argv[x];
        } else if (p_switch.compare("--method")==0) {
            x++;
            o.method=argv[x];
        } else if (p_switch.compare("--opt_t0")==0) {
            x++;
            o.opt_t0=atoi(argv[x]);
        } else if (p_switch.compare("--absorbing")==0) {
            x++;
            o.absorbing=atoi(argv[x]);
        } else if (p_switch.compare("--reflecting")==0) {
            x++;
            o.absorbing=1-atoi(argv[x]);
        } else if (p_switch.compare("--deg_shift")==0) {
            x++;
            o.deg_shift=atof(argv[x]);
        } else if (p_switch.compare("--viral_load")==0) {
            x++;
            o.viral_load=atof(argv[x]);
        } else if (p_switch.compare("--specify_k")==0) {
            x++;
            o.specify_k=atof(argv[x]);
        } else if (p_switch.compare("--seed")==0) {
            x++;
            o.seed=atoi(argv[x]);
        } else if (p_switch.compare("--move")==0) {
            x++;
            o.move=atoi(argv[x]);
        } else if (p_switch.compare("--exp_per_hour")==0) {
            x++;
            o.exp_per_hour=atoi(argv[x]);
        } else if (p_switch.compare("--phi")==0) {
            x++;
            o.phi=atof(argv[x]);
        } else if (p_switch.compare("--alpha")==0) {
            x++;
            o.variable_emissions=1;
            o.alpha=argv[x];
        } else if (p_switch.compare("--r0")==0) {
            x++;
            o.r0=atof(argv[x]);
        } else if (p_switch.compare("--test")==0) {
            x++;
            o.test_output=atoi(argv[x]);
        } else if (p_switch.compare("--vmult")==0) {
                x++;
                o.vmult=atof(argv[x]);
        } else if (p_switch.compare("--verb")==0) {
            x++;
            o.verb=atoi(argv[x]);
        } else {
            cout << "Incorrect usage\n ";
            exit(1);
        }
        p_switch.clear();
        x++;
    }
    //Maximum radius
    if (o.emm.compare("Sneeze")==0||o.emm.compare("Sneeze_Rotated")==0) {
        o.maxr=1000;
        o.volume=0.007298332079915014;
        o.exp_per_hour=5;

    }
    if (o.emm.compare("Speaking")==0||o.emm.compare("Speaking_Rotated")==0) {
        o.maxr=500;
        o.volume=5.78726*pow(10,-13);
        o.exp_per_hour=3600;
    }
    o.volume=o.volume*o.vmult;
}

void ReadParameters(options& o, param& p) {
    //Parameters.dat is X, Y, Z, gamma, z_0, Hours, Inhalation (L)
    //N.B. Z is the height of the room in this framework.
    ifstream par_file;
    par_file.open("Data/"+o.env+"/Parameters.dat");
    double d;
    par_file >> d;
    p.X=d;
    par_file >> d;
    p.Y=d;
    par_file >> d;
    p.Z=d;
    par_file >> d;
    p.gamma=d;
    par_file >> d;
    p.z_0=d;
    par_file >> d;
    p.hours=d;
    par_file >> d;
    p.inh=d;
    //Characteristic length
    p.L=pow(p.X*p.Y*p.Z,0.333333333333);
    if (o.specify_k==-1) {
        p.K=FindK(p.L,p.gamma);
    } else {
        p.K=o.specify_k;
    }
    //Convert inhalation in litres per minute into metres cubed per second
    p.inh=p.inh/(60*1000);  //Inhalation per minute converted to m^3 per second
    //Further convert into fraction of box volume per second.  Box is 40cm^2 times height of room
    p.inh=p.inh/(p.Y*0.02*0.02);
    if (o.verb==1) {
        cout << "Params X " << p.X << " Y " << p.Y << " Z " << p.Z << " V " << p.X*p.Y*p.Z << "\n";
        cout << "Characteristic L " << p.L << "\n";
        if (o.specify_k!==-1) {
            cout << "Specified ";
        }
        cout << "K is " << p.K << "\n";
        cout << "Z_0 " << p.z_0 << "\n";
    }
}

void ReadLocations(options& o, param& p, vector<loc>& locations) {
    ifstream loc_file;
    loc_file.open("Data/"+o.env+"/Locations.dat");
    double x;
    double y;
    loc l;
    for (int i=0;i<1000;i++) {
        if (!(loc_file >> x)) break;
        if (!(loc_file >> y)) break;
        if (i==0) {
            p.x_emit=x;
            p.y_emit=y;
        } else {  //First point is the infected individual
            l.x=x;
            l.y=y;
            locations.push_back(l);
        }
    }
}

void ReadCough(options& o, param& p, vector<dist>& initial) {
    ifstream c_file;
    c_file.open("Data/"+o.env+"/Initial_distribution/"+o.emm+"_grid.dat");
    double x;
    double y;
    double v;
    double tot=0;
    dist d;
    for (int i=0;i<1000000;i++) {
        if (!(c_file >> x)) break;
        if (!(c_file >> y)) break;
        if (!(c_file >> v)) break;
        d.x=x+p.x_emit;
        d.y=y+p.y_emit;
        d.val=v;
        tot=tot+v;
        initial.push_back(d);
    }
    for (unsigned int i=0;i<initial.size();i++) {
        initial[i].val=initial[i].val/tot;
    }
}
        
void ExportInitialState (options& o, vector< vector<double> >& initial_state) {
    ofstream init_file;
    string method;
    GetMethod(o,method);
    init_file.open("Data/"+o.env+"/Initial_state_"+o.emm+"_"+method+".dat");
    for (unsigned int i=0;i<initial_state.size();i++) {
           init_file << initial_state[i][0] << " " << initial_state[i][1] << " " << initial_state[i][2] << "\n";
    }
}

void GetMethod (options& o, string& method) {
    if (o.absorbing==1) {
        method="Absorbing";
    } else {
        method="Reflecting";
    }
}


void ExportExposures (options& o, vector<double>& times, vector< vector<double> >& exp_record) {
    ofstream exp_file;
    double r=o.radius*pow(10,6);
    int rr=floor(r+0.001);
    stringstream ss;
    ss << rr;
    string str = ss.str();
    string method;
    GetMethod(o,method);
    exp_file.open("Data/"+o.env+"/Exposures_"+o.emm+"_"+method+"_"+str+".dat");
    for (int i=0;i<3600;i++) {
        exp_file << times[i] << " ";
        for (unsigned int j=0;j<exp_record[i].size();j++) {
            exp_file << exp_record[i][j] << " ";
        }
        exp_file << "\n";
    }

}

void ReadExposures (options& o, const vector<loc>& locations, vector<double>& times, vector< vector< vector<double> > >& exposure_dat){
    //Set up exposuree datastructure
    SetupExposureData (o,locations,exposure_dat);
	cout << "#to read: " << o.maxr << "\n";
    for (int r=1;r<=o.maxr;r++) {
        ifstream exp_file;
        stringstream ss;
        ss << r;
        string str = ss.str();
        string method;
        GetMethod(o,method);
        exp_file.open("Data/"+o.env+"/Exposures_"+o.emm+"_"+method+"_"+str+".dat");
        if (o.verb==1) {
            cout << "Read file " << "Data/"+o.env+"/Exposures_"+o.emm+"_"+method+"_"+str+".dat" << "\n";
        }
        double t;
        double e;
        for (int i=1;i<=3600;i++) {
            if (!(exp_file >> t)) break;
            if (r==0) {
                times.push_back(t);
            }
            //cout << "T " << t << " E ";
            for (unsigned int j=0;j<locations.size();j++) {
                if (!(exp_file >> e)) break;
              //  cout << e << " ";
                exposure_dat[r-1][j][i-1]=e;
                //Note exposure_dat[k] is radius k+1 microns
            }
            //cout << "\n";
        }
        exp_file.close();
    }
}

void SetupExposureData (options& o, const vector<loc>& locations, vector< vector< vector<double> > >& exposure_dat) {
    for (int r=1;r<=o.maxr;r++) {
        vector< vector<double> > exposures;
        for (unsigned int i=0;i<locations.size();i++) {
            vector<double> e;
            for (int j=0;j<3600;j++) {
                e.push_back(0);
            }
            exposures.push_back(e);
        }
        exposure_dat.push_back(exposures);
    }
}

void ExportBottlenecks(options& o, vector<int>& bottlenecks) {
    ofstream b_file;
    string method;
    stringstream ss;
    stringstream sss;
    ss << o.phi;
    string str = ss.str();
    sss << o.vmult;
    string vm = sss.str();
    GetMethod(o,method);
    b_file.open("Data/"+o.env+"/Bottlenecks_"+o.emm+"_"+method+"_phi"+str+"_vm"+vm+".dat");
    for (unsigned int i=0;i<bottlenecks.size();i++) {
        b_file << bottlenecks[i] << " ";
    }
    b_file << "\n";
}

void ExportTestExposures (options& o, double k, vector<dist>& exps) {
    stringstream ss;
    ss << floor(k+0.00001);
    string str = ss.str();
    string method;
    GetMethod(o,method);
    ofstream test_file;
    test_file.open("Data/"+o.env+"/Data"+o.emm+"_"+method+"_"+str+".dat");
    for (unsigned int i=0;i<exps.size();i++) {
        test_file << exps[i].x << " " << exps[i].y << " " << exps[i].val << "\n";
    }
}
