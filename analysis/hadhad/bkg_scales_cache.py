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


def get_scales(category, embedded):

    if has_category(category, embedded):
        qcd_scale, qcd_scale_error, \
        ztautau_scale, ztautau_scale_error = SCALES[category][embedding]
        print "using the embedding scale factors: %s" % str(embedded)
        print "scale factors for %s category" % category
        print "    qcd scale: %.3f +/- %.4f" % (qcd_scale, qcd_scale_error)
        print "    ztautau scale: %.3f +/- %.4f" % (ztautau_scale, ztautau_scale_error)
        print
        return qcd_scale, qcd_scale_error, ztautau_scale, ztautau_scale_error
    else:
        return None


def has_category(category, embedded):

    return category in SCALES and embedded in SCALES[category]


def set_scales(category, embedding,
        qcd_scale, qcd_scale_error,
        ztautau_scale, ztautau_scale_error):

    global MODIFIED
    print "setting the embedding scale factors: %s" % str(embedding)
    print "setting scale factors for %s category" % category
    print "    qcd scale: %.3f +/- %.4f" % (qcd_scale, qcd_scale_error)
    print "    ztautau scale: %.3f +/- %.4f" % (ztautau_scale, ztautau_scale_error)
    print
    SCALES[category][embedding] = (
            qcd_scale, qcd_scale_error,
            ztautau_scale, ztautau_scale_error)
    MODIFIED = True


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        read_scales(sys.argv[1])
    else:
        read_scales()
    for category in SCALES.keys():
        for embedding in SCALES[category]:
            for (qcd_scale, qcd_scale_error,
                ztautau_scale, ztautau_scale_error) in SCALES[category][embedding]:
                print "scale factors for embedding: %s" % str(embedding)
                print "scale factors for %s category" % category
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
