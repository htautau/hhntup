
import samples
import features


s = samples.MC_VBF('2jet')

arr1 = s.array(features.hh_01jet_vars)

print arr1

arr2 = s.array(['mass_mmc_tau1_tau2'])

print arr2
