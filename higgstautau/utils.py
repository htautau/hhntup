from .units import GeV
import math
import os
import fnmatch

sign_zero = lambda x: 1 if x > 0 else -1 if x < 0 else 0

sign = lambda x: 1 if x >= 0 else -1

dphi = lambda phi1, phi2 : abs(math.fmod((math.fmod(phi1, 2*math.pi) - math.fmod(phi2, 2*math.pi)) + 3*math.pi, 2*math.pi) - math.pi)

dR = lambda eta1, phi1, eta2, phi2: math.sqrt((eta1 - eta2)**2 + dphi(phi1, phi2)**2)

is_visible = lambda fourvect: (
    fourvect.Et() > 10 * GeV and abs(fourvect.Eta()) < 2.5)


def et2pt(et, eta, m):
    return math.sqrt(et**2 - (m**2)/(math.cosh(eta)**2))


def pt2et(pt, eta, m):
    return math.sqrt(pt**2 + (m**2)/(math.cosh(eta)**2))


def Mvis(et1, phi1, et2, phi2):
    return math.sqrt(2. * et1 * et2 * (1. - math.cos(dphi(phi1, phi2))))


def closest_reco_object(objects, thing, dR=0.2):
    closest_object = None
    closest_dR = 1111
    for other in objects:
        dr = utils.dR(other.eta, other.phi, thing.Eta(), thing.Phi())
        if dr < dR and dr < closest_dR:
            closest_object = other
            closest_dR = dr
    return closest_object, closest_dR


def all_files_matching(dir, pattern):
    matched = []
    for path, dirs, files in os.walk(dir):
        for file in files:
            if fnmatch.fnmatch(file, pattern):
                matched.append(os.path.join(path, file))
    return matched


def find_file(filename, search_path_var='PATH', include_working=True):
    """
    find filename in PATH with the option of including
    the current working directory in the search
    """
    if not os.environ.has_key(search_path_var):
        if os.path.exists(filename):
            return os.path.abspath(filename)
        return None
    search_path = os.environ[search_path_var]
    paths = search_path.split(os.pathsep)
    if include_working:
        paths = ['.'] + paths
    for path in paths:
        fullpath = os.path.join(path, filename)
        if os.path.exists(fullpath):
            return os.path.abspath(fullpath)
    return None
