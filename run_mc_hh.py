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
    # EW background
    "AlpgenJimmyZtautauNp[0-5]_pt20.mc11c",
    "AlpgenJimmyZtautauNp[0-5]_Mll10to40_pt20.mc11c",
    "AlpgenJimmyWtaunuNp[0-5]_pt20.mc11c",
    "AlpgenJimmyZmumuNp[0-5]_pt20.mc11c",
    "AlpgenJimmyZmumuNp[0-5]_Mll10to40_pt20.mc11c",
    "AlpgenJimmyWmunuNp[0-5]_pt20.mc11c",
    "AlpgenJimmyZeeNp[0-5]_pt20.mc11c",
    "AlpgenJimmyZeeNp[0-5]_Mll10to40_pt20.mc11c",
    "AlpgenJimmyWenuNp[0-5]_pt20.mc11c",
    # Top
    "st_*",
    "AcerMC_Wt.mc11c",
    "T1_McAtNlo_Jimmy.mc11c",
    "TTbar_FullHad_McAtNlo_Jimmy.mc11c",
    # Diboson
    "gg2WW0240_JIMMY_WW_*",
    "McAtNlo_JIMMY_WmZ_*",
    "McAtNlo_JIMMY_WpWm_*",
    "McAtNlo_JIMMY_WpZ_*",
    "McAtNlo_JIMMY_ZZ_",
]

cluster.run('HHProcessor.py',
            db='datasets_hh',
            datasets=datasets,
            hosts=hosts,
            nproc=10,
            nice=10,
            setup=setup)
