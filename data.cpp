#include "diffusion.h"
#include "data.h"
#include <iostream>
#include <string>

void GetPrecalcT0Absorbing (options& o, param& p) {
    if (o.env.compare("Office")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=3.11206*pow(10,-5);
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=3.3117*pow(10,-5);
        } else if (o.emm.compare("Speaking")==0) {
            p.t_0=9.96951*pow(10,-6);
        } else if (o.emm.compare("Speaking_Rotated")==0) {
            p.t_0=9.96951*pow(10,-6);
        } else if (o.emm.compare("Sneeze")==0) {
            p.t_0=2.59633*pow(10,-5);
        } else if (o.emm.compare("Sneeze_Rotated")==0) {
            p.t_0=2.59633*pow(10,-5);
        } else if (o.emm.compare("Cough0.4")==0) {
            p.t_0=5.01887*pow(10,-5);
        } else if (o.emm.compare("Cough0.6")==0) {
            p.t_0=6.90848*pow(10,-5);
        } else if (o.emm.compare("Cough0.8")==0) {
            p.t_0=8.3499*pow(10,-5);
        } else if (o.emm.compare("Cough1.0")==0) {
            p.t_0=9.89874*pow(10,-5);
        }
    }
    if (o.env.compare("Office_Death_x0.25")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=3.11206*pow(10,-5);
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=3.3117*pow(10,-5);
        }
    }
    if (o.env.compare("Office_Death_x0.5")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=3.11206*pow(10,-5);
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=3.3117*pow(10,-5);
        }
    }
    if (o.env.compare("Office_Death_x2")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=3.11206*pow(10,-5);
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=3.3117*pow(10,-5);
        }
    }
    if (o.env.compare("Office_Death_x4")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=3.11206*pow(10,-5);
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=3.3117*pow(10,-5);
        }
    }


    if (o.env.compare("Office_Gamma_x2")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=1.69666*pow(10,-5);
        }
    }
    if (o.env.compare("Office_Gamma_x4")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=8.88483*pow(10,-6);
        }
    }
    if (o.env.compare("Office_Gamma_x0.5")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=5.33905*pow(10,-5);
        }
    }
    if (o.env.compare("Office_Gamma_x0.25")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=8.31367*pow(10,-5);
        }
    }

    
    if (o.env.compare("Nightclub")==0||o.env.compare("Nightclub_Movement")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=2.10192*pow(10,-6);
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=3.3117*pow(10,-5);
        } else if (o.emm.compare("Speaking")==0) {
            p.t_0=6.73528*pow(10,-7);
        } else if (o.emm.compare("Speaking_Rotated")==0) {
            p.t_0=6.73376*pow(10,-7);
        } else if (o.emm.compare("Sneeze")==0) {
            p.t_0=1.75365*pow(10,-6);
        } else if (o.emm.compare("Sneeze_Rotated")==0) {
            p.t_0=1.75365*pow(10,-6);
        }
    }
    if (o.env.compare("Lounge")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=0.000217408;
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=0.000231354;
        } else if (o.emm.compare("Speaking")==0) {
            p.t_0=6.96465*pow(10,-5);
        } else if (o.emm.compare("Sneeze")==0) {
            p.t_0=0.000181378;
        }
    }
    if (o.env.compare("Lounge_Full")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=0.000217408;
        }
    }

    if (o.env.compare("Bus")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=2.81052*pow(10,-5);
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=2.99056*pow(10,-5);
        } else if (o.emm.compare("Speaking")==0) {
            p.t_0=9.00276*pow(10,-6);
        } else if (o.emm.compare("Speaking_Rotated")==0) {
            p.t_0=9.00279*pow(10,-6);
        } else if (o.emm.compare("Sneeze_Rotated")==0) {
            p.t_0=2.34456*pow(10,-5);
        }
    }
}

void GetPrecalcT0Reflecting (options& o, param& p) {
    if (o.env.compare("Office")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=-0.00295209;
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=-0.0029501;
        } else if (o.emm.compare("Speaking")==0) {
            p.t_0=-0.00297324;
        } else if (o.emm.compare("Speaking_Rotated")==0) {
            p.t_0=-0.00297324;
        } else if (o.emm.compare("Sneeze")==0) {
            p.t_0=-0.00295725;
        } else if (o.emm.compare("Sneeze_Rotated")==0) {
            p.t_0=-0.00295725;
        }
    }
    if (o.env.compare("Nightclub")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=-0.000199394;
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=-0.000199259;
        } else if (o.emm.compare("Speaking")==0) {
            p.t_0=-0.000200823;
        } else if (o.emm.compare("Sneeze")==0) {
            p.t_0=-0.000199743;
        }
    }
    if (o.env.compare("Lounge")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=-0.0206232;
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=-0.0206092;
        } else if (o.emm.compare("Speaking")==0) {
            p.t_0=-0.0207709;
        } else if (o.emm.compare("Sneeze")==0) {
            p.t_0=-0.0206592;
        }
    }
    if (o.env.compare("Bus")==0) {
        if (o.emm.compare("Cough")==0) {
            p.t_0=-0.00266582;
        } else if (o.emm.compare("Cough_Rotated")==0) {
            p.t_0=-0.00266402;
        } else if (o.emm.compare("Speaking")==0) {
            p.t_0=-0.00268492;
        } else if (o.emm.compare("Speaking_Rotated")==0) {
            p.t_0=-0.00268492;
        } else if (o.emm.compare("Sneeze_Rotated")==0) {
            p.t_0=-0.00267048;
        }
    }
}


void GetSizeDistribution (options& o, double Pi, vector<double>& sizedist) {
    if (o.emm.compare("Cough")==0||o.emm.compare("Cough_Rotated")==0) {
        GetLognormalCough(sizedist);
    } else if (o.emm.compare("Sneeze")==0||o.emm.compare("Sneeze_Rotated")==0) {
        GetBimodalSneeze(Pi,sizedist);
    } else if (o.emm.compare("Speaking")==0||o.emm.compare("Speaking_Rotated")==0) {
        GetLognormalSpeaking(sizedist);
    }
}

void GetLognormalCough (vector<double>& sizedist) {
    double p1=log(13.5);
    double p2=-log(0.5);
    double tot=0;
    for (int d=2;d<=1000;d=d+2){
        double s=gsl_cdf_lognormal_P(d,p1,p2)-gsl_cdf_lognormal_P(d-2,p1,p2);
        //cout << s << "\n";
        sizedist.push_back(s);
        tot=tot+s;
    }
    for (unsigned int i=0;i<sizedist.size();i++) {
        sizedist[i]=sizedist[i]/tot;
    }
}

void GetLognormalSpeaking (vector<double>& sizedist) {
    double p1=log(16.0);
    double p2=-log(0.55);
    double tot=0;
    for (int d=2;d<=1000;d=d+2){
        double s=gsl_cdf_lognormal_P(d,p1,p2)-gsl_cdf_lognormal_P(d-2,p1,p2);
        sizedist.push_back(s);
        tot=tot+s;
    }
    for (unsigned int i=0;i<sizedist.size();i++) {
        sizedist[i]=sizedist[i]/tot;
    }
}
    
double LNVar(double a, double s, double u, double d, double Pi) {
    double g1=sqrt(2*Pi);
    g1=a/(s*g1);
    double g2=-pow(log10(d)-u,2);
    g2=g2/(2*pow(s,2));
    g2=exp(g2);
    g1=g1*g2;
    return(g1);
}

void GetBimodalSneeze (double Pi, vector<double>& sizedist) {
    double a1=4.9359;
    double s1=0.1479;
    double u1=2.733;
    double a2=2.1649;
    double s2=0.1483;
    double u2=2.005;
    vector<double> pdf;
    double tot=0;
    for (int d=0;d<=2000;d++) {
        double g1=LNVar(a1,s1,u1,d,Pi);
        double g2=LNVar(a2,s2,u2,d,Pi);
        double g=max(g1,g2);
        pdf.push_back(g);
        tot=tot+g;
    }
    for (unsigned int i=0;i<pdf.size();i++) {
        pdf[i]=pdf[i]/tot;
    }
    vector<double> cdf=pdf;
    for (unsigned int i=1;i<cdf.size();i++) {
        cdf[i]=cdf[i]+cdf[i-1];
    }
    tot=0;
    for (unsigned int i=2;i<cdf.size();i=i+2) {
        double c=cdf[i]-cdf[i-2];
        sizedist.push_back(c);
        tot=tot+c;
    }
    for (unsigned int i=0;i<sizedist.size();i++) {
        sizedist[i]=sizedist[i]/tot;
    }
}



