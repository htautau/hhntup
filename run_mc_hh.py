#!/usr/bin/env python

import cluster


hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.noel.sfu.txt')

datasets = [
    # signal
    "PowHegPythia_VBFH*_tautauhh.mc11c",
    "PowHegPythia_ggH*_tautauhh.mc11c",
    "PythiaZH*_tautauhh.mc11c",
    "PythiaWH*_tautauhh.mc11c",
    # background
    "AlpgenJimmyZtautauNp[0-5]_pt20.mc11c",
    "AlpgenJimmyZtautauNp[0-5]_Mll10to40_pt20.mc11c",
    "AlpgenJimmyWtaunuNp[0-5]_pt20.mc11c",
    "AlpgenJimmyZmumuNp[0-5]_pt20.mc11c",
    "AlpgenJimmyZmumuNp[0-5]_Mll10to40_pt20.mc11c",
    "AlpgenJimmyWmunuNp[0-5]_pt20.mc11c",
    "AlpgenJimmyZeeNp[0-5]_pt20.mc11c",
    "AlpgenJimmyZeeNp[0-5]_Mll10to40_pt20.mc11c",
    "AlpgenJimmyWenuNp[0-5]_pt20.mc11c",
    "st_tchan_taunu_McAtNlo_Jimmy.mc11c",
    "st_schan_taunu_McAtNlo_Jimmy.mc11c",
    "st_Wt_McAtNlo_Jimmy.mc11c",
    "McAtNlo_JIMMY_WpWm_taunutaunu.mc11c",
    "gg2WW0240_JIMMY_WW_taunutaunu.mc11c",
    "McAtNlo_JIMMY_WpZ_taunutautau.mc11c",
    "McAtNlo_JIMMY_WmZ_taunutautau.mc11c",
    "McAtNlo_JIMMY_WpZ_qqtautau.mc11c",
    "McAtNlo_JIMMY_ZZ_4tau.mc11c",
    "McAtNlo_JIMMY_ZZ_tautaununu.mc11c",
    "McAtNlo_JIMMY_ZZ_tautauqq.mc11c",
    "T1_McAtNlo_Jimmy.mc11c",
    "TTbar_FullHad_McAtNlo_Jimmy.mc11c",
]

cluster.run('HHProcessor.py',
            db='datasets_hh',
            datasets=datasets,
            hosts=hosts,
            nproc=10,
            nice=10,
            setup=setup)
