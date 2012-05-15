"""
This module handles retrieving the official Higgs cross sections [pb] and
branching ratios. See dat/README
"""
import os

__HERE = os.path.dirname(os.path.abspath(__file__))


def _read_xs_file(filename):

    xs = {}
    with open(filename) as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith('#'):
                continue
            mass, xs_mean, error_high, error_low = map(float, line[:4])
            xs[mass] = (xs_mean,
                        xs_mean * (1. + error_high),
                        xs_mean * (1. + error_low))
    return xs


def _read_br_file(filename):

    br = {}
    with open(filename) as f:
        for i, line in enumerate(f.readlines()):
            line = line.strip()
            if i == 0:
                # First line contains channel labels
                # Ignore first column which is the Higgs mass
                channels = line.split()[1:]
                for channel in channels:
                    br[channel] = {}
            else:
                line = map(float, line)
                for channel, value in zip(channels, line[1:]):
                    br[channel][line[0]] = value
    return br


MODES = [
    'ggf',
    'vbf',
    'wh',
    'zh',
    'tth'
]

CHANNEL_CATEGORIES = {
    '2fermions',
    '2gaugebosons',
    '4fermions_llll',
    '4fermions_llqq_qqqq'
]

__XS = {}
for mode in MODES:
    __XS[mode] = _read_xs_file(os.path.join(__HERE, 'dat', 'xs', '%s.txt' % mode))

__BR = {}
for channel_category in CHANNEL_CATEGORIES:
    __BR.update(_read_br_file(os.path.join(__HERE, 'dat', 'br', '%s.txt' %
                                           channel_category)))


def get_xs(mass, mode):
    """
    Return the production cross section [pb] in this mode in the form:
    (xs, xs_high, xs_low)
    """
    return __XS[mode][mass]


def get_br(mass, channel):
    """
    Return the branching ratio for this channel
    """
    return __BR[channel][mass]


def get_xsbr(mass, mode, channel):
    """
    Return the production cross section [pb] times branching ratio for this mode and
    channel in the form:
    (xsbr, xsbr_high, xsbr_low)
    """
    xs, xs_high, xs_low = get_xs(mass, mode)
    br = get_br(mass, channel)
    return (xs * br, xs_high * br, xs_low * br)
