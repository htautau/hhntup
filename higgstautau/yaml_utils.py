import inspect
import yaml
try:
    from collections import OrderedDict
except ImportError:
    # py 2.6
    from rootpy.extern.ordereddict import OrderedDict


class Serializable(yaml.YAMLObject):
    __metaclass__ = yaml.YAMLObjectMetaclass
    @property
    def _dict(self):
        dump_dict = OrderedDict()

        for var in inspect.getargspec(self.__init__).args[1:]:
            if getattr(self, var, None) is not None:
                dump_dict[var] = getattr(self, var)

        return dump_dict

    @classmethod
    def to_yaml(cls, dumper, data):
        return ordered_dump(dumper, '!{0}'.format(data.__class__.__name__),
                            data._dict)


    @classmethod
    def from_yaml(cls, loader, node):
        fields = loader.construct_mapping(node, deep=True)
        return cls(**fields)

def ordered_dump(dumper, tag, data):
    value = []
    node = yaml.nodes.MappingNode(tag, value)
    for key, item in data.iteritems():
        node_key = dumper.represent_data(key)
        node_value = dumper.represent_data(item)
        value.append((node_key, node_value))

    return node
