
void GetPrecalcT0Absorbing (options& o, param& p);
void GetPrecalcT0Reflecting (options& o, param& p);
void GetSizeDistribution (options& o, double Pi, vector<double>& sizedist);
void GetLognormalCough (vector<double>& sizedist);
void GetLognormalSpeaking (vector<double>& sizedist);
double LNVar(double a, double s, double u, double d, double Pi);
void GetBimodalSneeze (double Pi, vector<double>& sizedist);
void CalculateVolumeDistribution (double Pi, vector<double>& sizedist);




