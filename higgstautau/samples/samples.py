import os

__HERE = os.path.dirname(os.path.abspath(__file__))

import yaml
import yaml.constructor

try:
    # included in standard lib from Python 2.7
    from collections import OrderedDict
except ImportError:
    # try importing the backported drop-in replacement
    # it's available on PyPI
    from ordereddict import OrderedDict

import fnmatch
import itertools


class OrderedDictYAMLLoader(yaml.Loader):
    """
    A YAML loader that loads mappings into ordered dictionaries.
    https://gist.github.com/844388
    http://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts
    """
    def __init__(self, *args, **kwargs):

        yaml.Loader.__init__(self, *args, **kwargs)
        self.add_constructor(u'tag:yaml.org,2002:map', type(self).construct_yaml_map)
        self.add_constructor(u'tag:yaml.org,2002:omap', type(self).construct_yaml_map)

    def construct_yaml_map(self, node):

        data = OrderedDict()
        yield data
        value = self.construct_mapping(node)
        data.update(value)

    def construct_mapping(self, node, deep=False):

        if isinstance(node, yaml.MappingNode):
            self.flatten_mapping(node)
        else:
            raise yaml.constructor.ConstructorError(None, None,
                'expected a mapping node, but found %s' % node.id, node.start_mark)

        mapping = OrderedDict()
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            try:
                hash(key)
            except TypeError, exc:
                raise yaml.constructor.ConstructorError('while constructing a mapping',
                    node.start_mark, 'found unacceptable key (%s)' % exc, key_node.start_mark)
            value = self.construct_object(value_node, deep=deep)
            mapping[key] = value
        return mapping


SIGNALS_YML = {}
BACKGROUNDS_YML = {}
SAMPLES_YML = {}

SIGNALS = {}
BACKGROUNDS = {}
SAMPLES = {}

for channel in [o for o in os.listdir(__HERE) if
        os.path.isdir(os.path.join(__HERE, o))]:

    if channel not in SIGNALS_YML:
        SIGNALS_YML[channel] = {}
        BACKGROUNDS_YML[channel] = {}
        SAMPLES_YML[channel] = {}
        SIGNALS[channel] = {}
        BACKGROUNDS[channel] = {}
        SAMPLES[channel] = {}

    for year in [o
            for o in os.listdir(os.path.join(__HERE, channel)) if
            os.path.isdir(os.path.join(__HERE, channel, o))]:

        _year = int(year) % 1000
        SIGNALS_YML[channel][_year] = os.path.join(__HERE, channel, year,
                'signals.yml')
        BACKGROUNDS_YML[channel][_year] = os.path.join(__HERE, channel, year,
                'backgrounds.yml')
        SAMPLES_YML[channel][_year] = os.path.join(__HERE, channel, year,
                'samples.yml')

        SIGNALS[channel][_year] = yaml.load(
                open(SIGNALS_YML[channel][_year]),
                OrderedDictYAMLLoader)
        BACKGROUNDS[channel][_year] = yaml.load(
                open(BACKGROUNDS_YML[channel][_year]),
                OrderedDictYAMLLoader)
        SAMPLES[channel][_year] = yaml.load(
                open(SAMPLES_YML[channel][_year]))


def filter_with_patterns(samples, patterns):

    if not patterns:
        return samples
    if not isinstance(patterns, (list, tuple)):
        patterns = [patterns]
    matches = [fnmatch.filter(samples, pattern) for pattern in patterns]
    return list(set(itertools.chain(*matches)))


def iter_samples(channel, year, patterns=None, systematics=False):

    channel = channel.lower()
    year = year % 1000
    for sample_type, sample_info in SAMPLES[channel][year].items():
        samples = filter_with_patterns(sample_info['samples'][:], patterns)
        if systematics:
            if not samples:
                continue
            yield (samples,
                   [tuple(var.split(',')) for var in
                    sample_info['systematics']])
        else:
            if 'systematics_samples' in sample_info:
                for sample, sys_samp in sample_info['systematics_samples'].items():
                    samples += filter_with_patterns(sys_samp.keys(), patterns)
            if samples:
                yield samples


def samples(channel, year, patterns=None):

    return list(itertools.chain.from_iterable(
        iter_samples(channel, year, patterns, False)))


def get_systematics(channel, year, sample):

    channel = channel.lower()
    year = year % 1000
    terms = None
    sys_samples = None
    for sample_type, sample_info in SAMPLES[channel][year].items():
        if sample in sample_info['samples']:
            terms = [tuple(term.split(',')) for term in
                    sample_info['systematics']]
        else:
            continue
        if 'systematics_samples' in sample_info and sample in \
                sample_info['systematics_samples']:
            sys_samples = sample_info['systematics_samples'][sample]
        return terms, sys_samples
    raise ValueError("sample %s is not listed in samples.yml" % sample)


def get_sample(channel, year, sample_class, name):

    channel = channel.lower()
    year = year % 1000
    sample_class = sample_class.lower()
    if sample_class == 'signal':
        sample_class = SIGNALS[channel][year]
    elif sample_class == 'background':
        sample_class = BACKGROUNDS[channel][year]
    else:
        raise ValueError('sample class %s is not defined' % sample_class)
    if name in sample_class:
        return sample_class[name]
    raise ValueError('sample %s is not defined' % name)


if __name__ == '__main__':

    from higgstautau.datasets import Database
    from tabulartext import PrettyTable

    db = Database('datasets_hh')

    print "Signals"
    print "~~~~~~~"

    headers = ['dataset',
               'mean sigma [pb]', 'min sigma [pb]', 'max sigma [pb]',
               'sigma factor',
               'filter effic', 'K factor']

    for name, info in SIGNALS['hadhad'].items():
        print
        print ":math:`%s`" % info['math']
        print
        table = PrettyTable(headers)
        for sample in info['samples']:
            ds = db[sample]
            xsec, xsec_min, xsec_max, effic = ds.xsec_effic
            table.add_row([ds.ds, xsec*1E3, xsec_min*1E3, xsec_max*1E3,
                           ds.xsec_factor, effic, 1.])
        print table.get_string(hrules=1)

    print
    print "Backgrounds"
    print "~~~~~~~~~~~"


    for name, info in BACKGROUNDS['hadhad'].items():
        print
        print ":math:`%s`" % info['math']
        print
        table = PrettyTable(headers)
        for sample in info['samples']:
            ds = db[sample]
            xsec, xsec_min, xsec_max, effic = ds.xsec_effic
            table.add_row([ds.ds, xsec*1E3, xsec_min*1E3, xsec_max*1E3,
                           ds.xsec_factor, effic, 1.])
        print table.get_string(hrules=1)

