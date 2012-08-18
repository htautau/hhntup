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
#SAMPLES_YML = {}

SIGNALS = {}
BACKGROUNDS = {}
#SAMPLES = {}

for channel in [o for o in os.listdir(__HERE) if
        os.path.isdir(os.path.join(__HERE, o))]:

    SIGNALS_YML[channel] = os.path.join(__HERE, channel, 'signals.yml')
    BACKGROUNDS_YML[channel] = os.path.join(__HERE, channel, 'backgrounds.yml')
    #SAMPLES_YML[channel] = os.path.join(__HERE, channel, 'samples.yml')

    SIGNALS[channel] = yaml.load(open(SIGNALS_YML[channel]), OrderedDictYAMLLoader)
    BACKGROUNDS[channel] = yaml.load(open(BACKGROUNDS_YML[channel]), OrderedDictYAMLLoader)
    #SAMPLES[channel] = yaml.load(open(SAMPLES_YML[channel]), OrderedDictYAMLLoader)


def iter_samples(channel, type=None):

    channel = channel.lower()
    for sample_class in (SIGNALS[channel], BACKGROUNDS[channel]):
        for sample_type, sample_info in sample_class.items():
            if type is not None:
                if sample_type != type:
                    continue
            for sample_name in sample_info['samples']:
                yield sample_name


def samples(channel, type=None):

    return [name for name in iter_samples(channel, type)]


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

