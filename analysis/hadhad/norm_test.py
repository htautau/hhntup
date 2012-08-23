

from utils import draw


import numpy as np
from rootpy.plotting import Hist
import ROOT
# Setting this to True (default in rootpy)
# changes how the histograms look in ROOT...
ROOT.TH1.SetDefaultSumw2(False)

# create normal distributions
mu1, mu2, sigma1, sigma2 = 100, 100, 10, 20
x1 = mu1 + sigma1 * np.random.randn(30000)
x2 = mu2 + sigma2 * np.random.randn(50000)
x1_obs = mu1 + sigma1 * np.random.randn(30000)
x2_obs = mu2 + sigma2 * np.random.randn(50000)


for a, name in ((1, 'good'), (1.2, 'high'), (0.8, 'low')):

    # create histograms
    h1 = Hist(50, 40, 200, title=r'$Z\rightarrow\tau\tau$')
    h2 = h1.Clone(title='QCD Multijets')
    h3 = h1.Clone(title='Data')

    # fill the histograms with our distributions
    map(h1.Fill, x1)
    map(h2.Fill, x2)
    map(h3.Fill, x1_obs)
    map(h3.Fill, x2_obs)

    h1 *= a
    b = (sum(h3) - sum(h1)) / float(sum(h2))
    h2 *= b
    # set visual attributes
    h1.fillstyle = 'solid'
    h1.fillcolor = 'green'
    h1.linecolor = 'green'
    h1.linewidth = 0

    h2.fillstyle = 'solid'
    h2.fillcolor = 'red'
    h2.linecolor = 'red'
    h2.linewidth = 0

    draw(data=h3, model=[h2, h1], name='some variable', output_name='norm_test',
            category='example_%s' % name, category_name='example %s' % name,
            show_ratio=True)
