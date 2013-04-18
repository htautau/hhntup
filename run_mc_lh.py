#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--student', default='LHProcessorCN.py')
parser.add_argument('--systematics', default=None)
parser.add_argument('--nproc', type=int, default=5)
parser.add_argument('--queue', default='short')
parser.add_argument('--output-path', default='/cluster/data05/michel/Ntuples/lephad/2012')
parser.add_argument('--db', default='datasets_lh')
parser.add_argument('--nice', type=int, default=10)
parser.add_argument('--nominal-only', action='store_true', default=False)
parser.add_argument('--systematics-only', action='store_true', default=False)
parser.add_argument('--dry', action='store_true', default=False)
parser.add_argument('--use-ssh', dest='use_qsub', action='store_false', default=True)
parser.add_argument('samples', nargs='?', default=None)

args = parser.parse_args()

import sys
import cluster
from higgstautau import samples

hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.michel.sfu.txt')

datasets2011 = [
    'AcerMC_Wt.mc11c',
    'AlpgenJimmyWenuNp0_pt20.mc11c',
    'AlpgenJimmyWenuNp1_pt20.mc11c',
    'AlpgenJimmyWenuNp2_pt20.mc11c',
    'AlpgenJimmyWenuNp3_pt20.mc11c',
    'AlpgenJimmyWenuNp4_pt20.mc11c',
    'AlpgenJimmyWenuNp5_pt20.mc11c',
    'AlpgenJimmyWmunuNp0_pt20.mc11c',
    'AlpgenJimmyWmunuNp1_pt20.mc11c',
    'AlpgenJimmyWmunuNp2_pt20.mc11c',
    'AlpgenJimmyWmunuNp3_pt20.mc11c',
    'AlpgenJimmyWmunuNp4_pt20.mc11c',
    'AlpgenJimmyWmunuNp5_pt20.mc11c',
    'AlpgenJimmyWtaunuNp0_pt20.mc11c',
    'AlpgenJimmyWtaunuNp1_pt20.mc11c',
    'AlpgenJimmyWtaunuNp2_pt20.mc11c',
    'AlpgenJimmyWtaunuNp3_pt20.mc11c',
    'AlpgenJimmyWtaunuNp4_pt20.mc11c',
    'AlpgenJimmyWtaunuNp5_pt20.mc11c',
    'AlpgenJimmyZeeNp0Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZeeNp0VBFCut.mc11c',
    'AlpgenJimmyZeeNp0_Mll10to40_pt20.mc11c',
    #'AlpgenJimmyZeeNp0_pt20.mc11c',
    'AlpgenJimmyZeeNp1Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZeeNp1VBFCut.mc11c',
    'AlpgenJimmyZeeNp1_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZeeNp1_pt20.mc11c',
    'AlpgenJimmyZeeNp2Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZeeNp2VBFCut.mc11c',
    'AlpgenJimmyZeeNp2_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZeeNp2_pt20.mc11c',
    'AlpgenJimmyZeeNp3Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZeeNp3VBFCut.mc11c',
    'AlpgenJimmyZeeNp3_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZeeNp3_pt20.mc11c',
    'AlpgenJimmyZeeNp4Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZeeNp4VBFCut.mc11c',
    'AlpgenJimmyZeeNp4_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZeeNp4_pt20.mc11c',
    'AlpgenJimmyZeeNp5Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZeeNp5VBFCut.mc11c',
    'AlpgenJimmyZeeNp5_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZeeNp5_pt20.mc11c',
    'AlpgenJimmyZmumuNp0Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZmumuNp0VBFCut.mc11c',
    'AlpgenJimmyZmumuNp0_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZmumuNp0_pt20.mc11c',
    'AlpgenJimmyZmumuNp1Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZmumuNp1VBFCut.mc11c',
    'AlpgenJimmyZmumuNp1_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZmumuNp1_pt20.mc11c',
    'AlpgenJimmyZmumuNp2Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZmumuNp2VBFCut.mc11c',
    'AlpgenJimmyZmumuNp2_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZmumuNp2_pt20.mc11c',
    'AlpgenJimmyZmumuNp3Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZmumuNp3VBFCut.mc11c',
    'AlpgenJimmyZmumuNp3_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZmumuNp3_pt20.mc11c',
    'AlpgenJimmyZmumuNp4Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZmumuNp4VBFCut.mc11c',
    'AlpgenJimmyZmumuNp4_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZmumuNp4_pt20.mc11c',
    'AlpgenJimmyZmumuNp5Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZmumuNp5VBFCut.mc11c',
    'AlpgenJimmyZmumuNp5_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZmumuNp5_pt20.mc11c',
    'AlpgenJimmyZtautauNp0VBFCut.mc11c',
    'AlpgenJimmyZtautauNp0_Mll10to40_pt20.mc11c',
    #'AlpgenJimmyZtautauNp0_pt20.mc11c',
    'AlpgenJimmyZtautauNp1Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZtautauNp1VBFCut.mc11c',
    'AlpgenJimmyZtautauNp1_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZtautauNp1_pt20.mc11c',
    'AlpgenJimmyZtautauNp2Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZtautauNp2VBFCut.mc11c',
    'AlpgenJimmyZtautauNp2_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZtautauNp2_pt20.mc11c',
    'AlpgenJimmyZtautauNp3Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZtautauNp3VBFCut.mc11c',
    'AlpgenJimmyZtautauNp3_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZtautauNp3_pt20.mc11c',
    'AlpgenJimmyZtautauNp4Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZtautauNp4VBFCut.mc11c',
    'AlpgenJimmyZtautauNp4_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZtautauNp4_pt20.mc11c',
    'AlpgenJimmyZtautauNp5Mll7to40VBFCut.mc11c',
    'AlpgenJimmyZtautauNp5VBFCut.mc11c',
    'AlpgenJimmyZtautauNp5_Mll10to40_pt20.mc11c',
    'AlpgenJimmyZtautauNp5_pt20.mc11c',
    'McAtNlo_JIMMY_WmZ_lnull.mc11c',
    'McAtNlo_JIMMY_WmZ_lnuqq.mc11c',
    'McAtNlo_JIMMY_WmZ_lnutautau.mc11c',
    'McAtNlo_JIMMY_WmZ_qqll.mc11c',
    'McAtNlo_JIMMY_WmZ_qqtautau.mc11c',
    'McAtNlo_JIMMY_WmZ_taunull.mc11c',
    'McAtNlo_JIMMY_WmZ_taunutautau.mc11c',
    'McAtNlo_JIMMY_WpWm_enuenu.mc11c',
    'McAtNlo_JIMMY_WpWm_enumunu.mc11c',
    'McAtNlo_JIMMY_WpWm_enutaunu.mc11c',
    'McAtNlo_JIMMY_WpWm_munuenu.mc11c',
    'McAtNlo_JIMMY_WpWm_munumunu.mc11c',
    'McAtNlo_JIMMY_WpWm_munutaunu.mc11c',
    'McAtNlo_JIMMY_WpWm_taunuenu.mc11c',
    'McAtNlo_JIMMY_WpWm_taunumunu.mc11c',
    'McAtNlo_JIMMY_WpWm_taunutaunu.mc11c',
    'McAtNlo_JIMMY_WpZ_lnull.mc11c',
    'McAtNlo_JIMMY_WpZ_lnuqq.mc11c',
    'McAtNlo_JIMMY_WpZ_lnutautau.mc11c',
    'McAtNlo_JIMMY_WpZ_qqll.mc11c',
    'McAtNlo_JIMMY_WpZ_qqtautau.mc11c',
    'McAtNlo_JIMMY_WpZ_taunull.mc11c',
    'McAtNlo_JIMMY_WpZ_taunutautau.mc11c',
    'McAtNlo_JIMMY_ZZ_2l2tau.mc11c',
    'McAtNlo_JIMMY_ZZ_4tau.mc11c',
    'McAtNlo_JIMMY_ZZ_llll.mc11c',
    'McAtNlo_JIMMY_ZZ_llnunu.mc11c',
    'McAtNlo_JIMMY_ZZ_llqq.mc11c',
    'McAtNlo_JIMMY_ZZ_tautaununu.mc11c',
    'McAtNlo_JIMMY_ZZ_tautauqq.mc11c',
    'PowHegPythia_VBFH100_tautaulh.mc11c',
    'PowHegPythia_VBFH105_tautaulh.mc11c',
    'PowHegPythia_VBFH110_tautaulh.mc11c',
    'PowHegPythia_VBFH115_tautaulh.mc11c',
    'PowHegPythia_VBFH120_tautaulh.mc11c',
    'PowHegPythia_VBFH125_tautaulh.mc11c',
    'PowHegPythia_VBFH130_tautaulh.mc11c',
    'PowHegPythia_VBFH135_tautaulh.mc11c',
    'PowHegPythia_VBFH140_tautaulh.mc11c',
    'PowHegPythia_VBFH145_tautaulh.mc11c',
    'PowHegPythia_VBFH150_tautaulh.mc11c',
    'PowHegPythia_ggH100_tautaulh.mc11c',
    'PowHegPythia_ggH105_tautaulh.mc11c',
    'PowHegPythia_ggH110_tautaulh.mc11c',
    'PowHegPythia_ggH115_tautaulh.mc11c',
    'PowHegPythia_ggH120_tautaulh.mc11c',
    'PowHegPythia_ggH125_tautaulh.mc11c',
    'PowHegPythia_ggH130_tautaulh.mc11c',
    'PowHegPythia_ggH135_tautaulh.mc11c',
    'PowHegPythia_ggH140_tautaulh.mc11c',
    'PowHegPythia_ggH145_tautaulh.mc11c',
    'PowHegPythia_ggH150_tautaulh.mc11c',
    'PythiaWH100_tautaulh.mc11c',
    'PythiaWH105_tautaulh.mc11c',
    'PythiaWH110_tautaulh.mc11c',
    'PythiaWH115_tautaulh.mc11c',
    'PythiaWH120_tautaulh.mc11c',
    'PythiaWH125_tautaulh.mc11c',
    'PythiaWH130_tautaulh.mc11c',
    'PythiaWH135_tautaulh.mc11c',
    'PythiaWH140_tautaulh.mc11c',
    'PythiaWH145_tautaulh.mc11c',
    'PythiaWH150_tautaulh.mc11c',
    'PythiaZH100_tautaulh.mc11c',
    'PythiaZH105_tautaulh.mc11c',
    'PythiaZH110_tautaulh.mc11c',
    'PythiaZH115_tautaulh.mc11c',
    'PythiaZH120_tautaulh.mc11c',
    'PythiaZH125_tautaulh.mc11c',
    'PythiaZH130_tautaulh.mc11c',
    'PythiaZH135_tautaulh.mc11c',
    'PythiaZH140_tautaulh.mc11c',
    'PythiaZH145_tautaulh.mc11c',
    'PythiaZH150_tautaulh.mc11c',
    #'T1_McAtNlo_Jimmy.mc11c',
    'TTbar_FullHad_McAtNlo_Jimmy.mc11c',
    'gg2WW0240_JIMMY_WW_enuenu.mc11c',
    'gg2WW0240_JIMMY_WW_enutaunu.mc11c',
    'gg2WW0240_JIMMY_WW_munuenu.mc11c',
    'gg2WW0240_JIMMY_WW_munumunu.mc11c',
    'gg2WW0240_JIMMY_WW_munutaunu.mc11c',
    'gg2WW0240_JIMMY_WW_taunuenu.mc11c',
    'gg2WW0240_JIMMY_WW_taunutaunu.mc11c',
    'st_schan_enu_AcerMC.mc11c',
    'st_schan_munu_AcerMC.mc11c',
    'st_schan_taunu_AcerMC.mc11c',
    'st_tchan_enu_AcerMC.mc11c',
    'st_tchan_munu_AcerMC.mc11c',
    'st_tchan_taunu_AcerMC.mc11c',
    ]

datasets2012 = [
    'ggH100',
    'ggH105',
    'ggH110',
    'ggH115',
    'ggH120',
    'ggH125',
    'ggH130',
    'ggH135',
    'ggH140',
    'ggH145',
    'ggH150',

    'VBFH100',
    'VBFH105',
    'VBFH110',
    'VBFH115',
    'VBFH120',
    'VBFH125',
    'VBFH130',
    'VBFH135',
    'VBFH140',
    'VBFH145',
    'VBFH150',

    'WH100',
    'WH105',
    'WH110',
    'WH115',
    'WH120',
    'WH125',
    'WH130',
    'WH135',
    'WH140',
    'WH145',
    'WH150',

    'ZH100',
    'ZH105',
    'ZH110',
    'ZH115',
    'ZH120',
    'ZH125',
    'ZH130',
    'ZH135',
    'ZH140',
    'ZH145',
    'ZH150',

    'SingleTopSChanWenu',
    'SingleTopSChanWmunu',
    'SingleTopSChanWtaunu',
    'singletop_tchan_e',
    'singletop_tchan_mu',
    'singletop_tchan_tau',
    'SingleTopWtChanIncl',

    'ttbar',
    'ttbar_allhad',

    'WtaunuNp0',
    'WtaunuNp1',
    'WtaunuNp2',
    'WtaunuNp3',
    'WtaunuNp4',
    'WtaunuNp5',

    'WenuNp0',
    'WenuNp1',
    'WenuNp2',
    'WenuNp3',
    'WenuNp4',
    'WenuNp5',

    'WmunuNp0',
    'WmunuNp1',
    'WmunuNp2',
    'WmunuNp3',
    'WmunuNp4',
    'WmunuNp5',

    'WZ',
    'ZZ',
    'WWlnulnuNp0',
    'WWlnulnuNp1',
    'WWlnulnuNp2',
    'WWlnulnuNp3',
    'WWqqlnuNp0',
    'WWqqlnuNp1',
    'WWqqlnuNp2',
    'WWqqlnuNp3',

    'ZmumuNp0',
    'ZmumuNp1',

    'ZeeNp0',
    'ZeeNp1',

    'ZtautauNp0',
    'ZtautauNp1',
    ]

datasets2012_filter = [

    # 'ZmumuNp2',
    # 'ZmumuNp3',
    # 'ZmumuNp4',
    # 'ZmumuNp5',
    # 'ZeeNp2',
    # 'ZeeNp3',
    # 'ZeeNp4',
    # 'ZeeNp5',
    # 'ZtautauNp2',
    # 'ZtautauNp3',
    # 'ZtautauNp4',
    # 'ZtautauNp5',

    # 'VBF_ZmumuNp2',
    # 'VBF_ZmumuNp3',
    # 'VBF_ZmumuNp4',
    # 'VBF_ZmumuNp5',
    # 'VBF_ZeeNp2',
    # 'VBF_ZeeNp3',
    # 'VBF_ZeeNp4',
    # 'VBF_ZeeNp5',
    # 'VBF_ATau_ZtautauNp2',
    # 'VBF_ATau_ZtautauNp3',
    # 'VBF_ATau_ZtautauNp4',
    # 'VBF_ATau_ZtautauNp5',

    'TVBF_ZmumuNp2',
    'TVBF_ZmumuNp3',
    'TVBF_ZmumuNp4',
    'TVBF_ZmumuNp5',
    'TVBF_ZeeNp2',
    'TVBF_ZeeNp3',
    'TVBF_ZeeNp4',
    'TVBF_ZeeNp5',
    'TVBF_ATau_ZtautauNp2',
    'TVBF_ATau_ZtautauNp3',
    'TVBF_ATau_ZtautauNp4',
    'TVBF_ATau_ZtautauNp5',
]

datasets = datasets2012

if not args.systematics_only:
    # nominal values
    cluster.run(args.student,
                db=args.db,
                datasets=datasets,
                hosts=hosts,
                nproc=args.nproc,
                nice=args.nice,
                setup=setup,
                output_path=args.output_path,
                use_qsub=args.use_qsub,
                qsub_queue=args.queue,
                dry_run=args.dry)

if not args.nominal_only:
    if args.systematics is not None:
        args.systematics = [
                set(s.upper().split('+')) for s in
                args.systematics.split(',')]
    # systematics
    cluster.run_systematics('LEPHAD',
                args.student,
                db=args.db,
                systematics=args.systematics,
                datasets=datasets,
                hosts=hosts,
                nproc=args.nproc,
                nice=args.nice,
                setup=setup,
                output_path=args.output_path,
                use_qsub=args.use_qsub,
                qsub_queue=args.queue,
                dry_run=args.dry)
