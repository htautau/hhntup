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
    parser.add_argument('--plot', action='store_true', default=False)
    parser.add_argument('cache', default='background_scales.cache', nargs='?')
    args = parser.parse_args()

    if args.plot:
        from matplotlib import pyplot as plt

    z_norms = {}
    qcd_norms = {}
    read_scales(args.cache)
    categories = sorted(SCALES.keys())
    for category in categories:
        for embedding in SCALES[category].keys():
            if embedding != args.embedding:
                    continue
            params = sorted(SCALES[category][embedding].keys())
            for param in params:
                qcd_scale, qcd_scale_error, \
                ztautau_scale, ztautau_scale_error = \
                SCALES[category][embedding][param]
                if param not in z_norms:
                    z_norms[param] = []
                    qcd_norms[param] = []
                z_norms[param].append((ztautau_scale, ztautau_scale_error))
                qcd_norms[param].append((qcd_scale, qcd_scale_error))
                print "scale factors for embedding: %s" % str(embedding)
                print "scale factors for %s category" % category
                print "fits were derived via %s parameters" % param
                print "    qcd scale: %.3f +/- %.4f" % (qcd_scale, qcd_scale_error)
                print "    ztautau scale: %.3f +/- %.4f" % (ztautau_scale, ztautau_scale_error)
                print
    if args.plot:

        plt.figure()
        fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)

        for param in params:
            z_n, z_e = zip(*z_norms[param])
            q_n, q_e = zip(*qcd_norms[param])
            x = range(len(z_n))

            ax = axs[0]
            ax.errorbar(x, z_n, z_e, fmt='o',
                    markersize=5, label='%s Fit' % param)
            ax = axs[1]
            ax.errorbar(x, q_n, q_e, fmt='o',
                    markersize=5, label='%s Fit' % param)

        axs[0].set_ylim(.5, 2)
        axs[0].set_ylabel('Z Scale Factor')

        axs[1].set_ylim(.5, 2.5)
        axs[1].set_xticklabels([''] + categories)
        axs[1].set_xlim(-0.5, len(z_norms[params[0]]) - .5)
        axs[1].set_ylabel('QCD Scale Factor')

        l1 = axs[0].legend(numpoints=1)
        l2 = axs[1].legend(numpoints=1)

        l1.get_frame().set_fill(False)
        l1.get_frame().set_linewidth(0)

        l2.get_frame().set_fill(False)
        l2.get_frame().set_linewidth(0)

        plt.savefig('bkg_norms.png')
        plt.savefig('bkg_norms.eps')

else:
    import atexit

    read_scales()

    @atexit.register
    def write_scales():

        if MODIFIED:
            with open(SCALES_FILE, 'w') as cache:
                pickle.dump(SCALES, cache)
