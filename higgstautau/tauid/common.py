import os


HERE = os.path.dirname(os.path.abspath(__file__))

LEVELS = {
    'loose': 1,
    'medium': 2,
    'tight': 3,
}

PRONGS = (1, 3)


def selection_file(year):

    if year == 2011:
        return os.path.join(HERE, 'selection', 'p851', 'bdt_selection.root')
    elif year == 2012:
        return os.path.join(HERE, 'selection', 'p1130', 'ParametrizedBDTSelection.root')
    else:
        raise ValueError("No BDT selection defined for year %d" % year)


def nprong(ntrack):

    if ntrack > 1:
        return 3
    return 1
