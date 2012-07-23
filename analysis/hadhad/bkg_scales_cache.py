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


def get_scales(category):

    if category in SCALES:
        qcd_scale, ztautau_scale = SCALES[category]
        print "scale factors for %s category" % category
        print "    qcd scale: %.3f" % qcd_scale
        print "    ztautau scale: %.3f" % ztautau_scale
        print
        return qcd_scale, ztautau_scale
    else:
        return None


def has_category(category):

    return category in SCALES


def set_scales(category, qcd_scale, ztautau_scale):

    global MODIFIED
    print "setting scale factors for %s category" % category
    print "    qcd scale: %.3f" % qcd_scale
    print "    ztautau scale: %.3f" % ztautau_scale
    print
    SCALES[category] = (qcd_scale, ztautau_scale)
    MODIFIED = True


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        read_scales(sys.argv[1])
    else:
        read_scales()
    for category, (qcd_scale, ztautau_scale) in SCALES.items():
        print "scale factors for %s category" % category
        print "    qcd scale: %.5f" % qcd_scale
        print "    ztautau scale: %.5f" % ztautau_scale
        print
else:
    import atexit

    read_scales()

    @atexit.register
    def write_scales():

        if MODIFIED:
            with open(SCALES_FILE, 'w') as cache:
                pickle.dump(SCALES, cache)
