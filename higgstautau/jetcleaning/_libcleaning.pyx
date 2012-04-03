from cpython cimport bool

cdef extern from "_cleaning.h":
        
        enum BADLEVEL:
            LooseMinusBad = 0,
            LooseBad = 1,
            MediumBad = 2,
            TightBad = 3

        bint is_bad_jet(BADLEVEL criteria,
            double quality, double NegE,
            double emf,     double hecf,
            double time,    double fmax,
            double eta,     double chf ,
            double HecQ,    double LArQmean )

def is_bad_jet_(criteria,
             quality, NegE,
             emf,     hecf,
             time,    fmax,
             eta,     chf ,
             HecQ,    LArQmean ):
    
    return is_bad_jet(criteria,
            quality, NegE,
            emf,     hecf,
            time,    fmax,
            eta,     chf ,
            HecQ,    LArQmean )
