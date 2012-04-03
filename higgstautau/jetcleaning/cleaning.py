from . import _libcleaning

LOOSER,\
LOOSE,\
MEDIUM,\
TIGHT = range(4)

def is_bad(level,
	quality, NegE,
	emf,     hecf,
	time,    fmax,
	eta,     chf ,
    HecQ,    LArQmean ):
    return _libcleaning.is_bad_jet_(level,
        quality, NegE,
        emf,     hecf,
        time,    fmax,
        eta,     chf ,
        HecQ,    LArQmean )
