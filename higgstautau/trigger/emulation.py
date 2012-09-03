from rootpy.tree.filters import EventFilter

from externaltools import CoEPPTrigTool
from ROOT import CoEPP


class TauTriggerEmulation(EventFilter):
    """
    Tau trigger emulation (only apply on MC)
    """
    def __init__(self, year, outtree):

        self.outtree = outtree

        if year == 2011: # only can emulate 2011 currently...
            # initialize the trigger emulation tool
            self.trigger_tool_wrapper = CoEPP.OfficialWrapper()
            self.trigger_tool = CoEPP.TriggerTool()
            self.trigger_tool.setWrapper(trigger_tool_wrapper)
            trigger_config = CoEPPTrigTool.get_resource(
                    'config_EF_DiTau.xml')
            self.trigger_tool.setXMLFile(trigger_config)
            self.trigger_tool.initializeFromXML()

            self.trigger_A = trigger_tool.getTriggerChecked(
                    "EF_tau29_medium1_tau20_medium1_Hypo_00_02_42")
            self.trigger_B = trigger_tool.getTriggerChecked(
                    "EF_tau29_medium1_tau20_medium1_Hypo_00_03_02")
            self.trigger_C = trigger_tool.getTriggerChecked(
                    "EF_tau29T_medium1_tau20T_medium1_Hypo_00_03_02")

            self.trigger_run_dict = {
                180164: (self.trigger_A, 'EF_tau29_medium1_tau20_medium1'),
                183003: (self.trigger_B, 'EF_tau29_medium1_tau20_medium1'),
                186169: (self.trigger_B, 'EF_tau29_medium1_tau20_medium1'),
                189751: (self.trigger_C, 'EF_tau29T_medium1_tau20T_medium1'),
            }

            self.passes = self.passes_11
            self.finalize = self.finalize_11

    @staticmethod
    def update_trigger_trees(student, tool, name, file, tree):
        """
        This method must be called when each new tree is loaded in the chain
        """
        tool.trigger_tool_wrapper.loadMainTree(tree)
        tool.trigger_tool_wrapper.loadMetaTree(
                file.Get('%sMeta/TrigConfTree' % name))

    def passes_11(self, event):

        self.trigger_tool_wrapper.setEventNumber(event._entry.value)
        trigger, triggername = self.trigger_run_dict[event.RunNumber]
        trigger.switchOn()
        self.trigger_tool.executeTriggers()
        emulated_trigger_passed = trigger.passed()

        self.outtree.tau_trigger_match_index.clear()
        self.outtree.tau_trigger_match_thresh.clear()

        if emulated_trigger_passed:
            if triggername == 'EF_tau29_medium1_tau20_medium1':
                self.outtree.EF_tau29_medium1_tau20_medium1_EMULATED = True
            else:
                self.outtree.EF_tau29T_medium1_tau20T_medium1_EMULATED = True

            # trigger matching
            trig1 = trigger.getTrigger1() # EF_tau29(T)_medium1
            trig2 = trigger.getTrigger2() # EF_tau20(T)_medium1
            for tau in event.taus:
                thresh = 0
                idx = -1
                idx1 = trig1.matchIndex(tau.fourvect)
                idx2 = trig2.matchIndex(tau.fourvect)
                if idx1 == idx2 != -1:
                    idx = idx1
                    thresh = 29
                elif idx1 == -1 and idx2 > -1:
                    idx = idx2
                    thresh = 20
                elif idx2 == -1 and idx1 > -1:
                    idx = idx1
                    thresh = 29
                elif idx2 != idx1: # both >-1 and non-equal
                    # take index of closer one using dR
                    trigtau1TLV = self.trigger_tool.buildEFTauTLV(idx1)
                    trigtau2TLV = self.trigger_tool.buildEFTauTLV(idx2)
                    if trigtau1TLV.DeltaR(tau.fourvect) < trigtau2TLV.DeltaR(tau.fourvect):
                        idx = idx1
                        thresh = 29
                    else:
                        idx = idx2
                        thresh = 20

                self.outtree.tau_trigger_match_index.push_back(idx)
                self.outtree.tau_trigger_match_thresh.push_back(thresh)
        else:
            self.outtree.EF_tau29_medium1_tau20_medium1_EMULATED = False
            self.outtree.EF_tau29T_medium1_tau20T_medium1_EMULATED = False
            for tau in event.taus:
                self.outtree.tau_trigger_match_index.push_back(-1)
                self.outtree.tau_trigger_match_thresh.push_back(0)

        trigger.switchOff()
        return True

    def finalize_11(self):

        # turn on triggers so they show up as "active" in the report
        self.trigger_A.switchOn()
        self.trigger_B.switchOn()
        self.trigger_C.switchOn()

        # finalize the trigger_tool
        self.trigger_tool.finalize()
        self.trigger_tool.summary()
