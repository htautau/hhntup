#!/usr/bin/env python
import cPickle as pickle
import os

HERE = os.path.dirname(os.path.abspath(__file__))

SCALES_FILE = None
SCALES = {}
MODIFIED = False


def read_scales(name='background_scales.cache'):

    global SCALES_FILE
    global SCALES

    SCALES_FILE = os.path.join(HERE, name)
    print "reading background scale factors from %s" % SCALES_FILE

    if os.path.isfile(SCALES_FILE):
        with open(SCALES_FILE) as cache:
            SCALES = pickle.load(cache)


def get_scales(category, embedded, param, verbose=True):

    category = category.upper()
    param = param.upper()
    if has_category(category, embedded, param):
        qcd_scale, qcd_scale_error, \
        ztautau_scale, ztautau_scale_error = SCALES[category][embedded][param]
        if verbose:
            print "using the embedding scale factors: %s" % str(embedded)
            print "scale factors for %s category" % category
            print "fits were derived via %s parameters" % param
            print "    qcd scale: %.3f +/- %.4f" % (qcd_scale, qcd_scale_error)
            print "    ztautau scale: %.3f +/- %.4f" % (
                    ztautau_scale, ztautau_scale_error)
            print
        return qcd_scale, qcd_scale_error, ztautau_scale, ztautau_scale_error
    else:
        return None


def has_category(category, embedded, param):

    category = category.upper()
    param = param.upper()
    return (category in SCALES and
            embedded in SCALES[category] and
            param in SCALES[category][embedded])


def set_scales(category, embedded, param,
        qcd_scale, qcd_scale_error,
        ztautau_scale, ztautau_scale_error):

    global MODIFIED
    param = param.upper()
    category = category.upper()
    print "setting the embedding scale factors: %s" % str(embedded)
    print "setting scale factors for %s category" % category
    print "fits were derived via %s parameters" % param
    print "    qcd scale: %.3f +/- %.4f" % (qcd_scale, qcd_scale_error)
    print "    ztautau scale: %.3f +/- %.4f" % (ztautau_scale, ztautau_scale_error)
    print
    if has_category(category, embedded, param):
        qcd_scale_old, qcd_scale_error_old, \
        ztautau_scale_old, ztautau_scale_error_old = get_scales(
                category, embedded, param, verbose=False)
        print "scale factors were previously:"
        print "    qcd scale: %.3f +/- %.4f" % (
                qcd_scale_old,
                qcd_scale_error_old)
        print "    ztautau scale: %.3f +/- %.4f" % (
                ztautau_scale_old,
                ztautau_scale_error_old)
        print
    if category not in SCALES:
        SCALES[category] = {}
    if embedded not in SCALES[category]:
        SCALES[category][embedded] = {}
    SCALES[category][embedded][param] = (
            qcd_scale, qcd_scale_error,
            ztautau_scale, ztautau_scale_error)
    MODIFIED = True


if __name__ == '__main__':
    import sys
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--embedding', action='store_true', default=False)
    parser.add_argument('cache', default='background_scales.cache', nargs='?')
    args = parser.parse_args()

    read_scales(args.cache)
    for category in SCALES.keys():
        for embedding in SCALES[category].keys():
            if embedding != args.embedding:
                    continue
            for param, (qcd_scale, qcd_scale_error,
                    ztautau_scale, ztautau_scale_error) in \
                    SCALES[category][embedding].items():
                print "scale factors for embedding: %s" % str(embedding)
                print "scale factors for %s category" % category
                print "fits were derived via %s parameters" % param
                print "    qcd scale: %.3f +/- %.4f" % (qcd_scale, qcd_scale_error)
                print "    ztautau scale: %.3f +/- %.4f" % (ztautau_scale, ztautau_scale_error)
                print
else:
    import atexit

    read_scales()

    @atexit.register
    def write_scales():

        if MODIFIED:
            with open(SCALES_FILE, 'w') as cache:
                pickle.dump(SCALES, cache)
