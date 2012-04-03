#ifndef CLEANING
#define CLEANING

enum BADLEVEL
{
	LooseMinusBad = 0,
	LooseBad = 1,
	MediumBad = 2,
	TightBad = 3
};

bool is_bad_jet(BADLEVEL criteria,
	double quality, double NegE,
	double emf,     double hecf,
	double time,    double fmax,
	double eta,     double chf ,
	double HecQ,    double LArQmean );

#endif
