from rootpy.tree.filtering import EventFilter


def update_trigger_trees(student, tool, name, file, tree):
    """
    This method must be called when each new tree is loaded in the chain
    """
    if not tool.passthrough and tool.year == 2011:
        tool.trigger_tool_wrapper.loadMainTree(tree)
        tool.trigger_tool_wrapper.loadMetaTree(
                file.Get('%sMeta/TrigConfTree' % name))


class TauTriggerEmulation(EventFilter):
    """
    Tau trigger emulation (only apply on MC)
    """
    def __init__(self, year, passthrough=False, **kwargs):

        if not passthrough:
            from externaltools import TauTriggerEmulation as TTE
            from ROOT import CoEPP

            self.year = year

            if year == 2011:
                # emulation not required in 2012 yet since the SFs are wrt
                # the default triggers
                # initialize the trigger emulation tool
                self.trigger_tool_wrapper = CoEPP.OfficialWrapper()
                self.trigger_tool = CoEPP.TriggerTool()
                self.trigger_tool.setWrapper(self.trigger_tool_wrapper)
                trigger_config = TTE.get_resource(
                        'config_EF_DiTau.xml')
                self.trigger_tool.setXMLFile(trigger_config)
                self.trigger_tool.initializeFromXML()

                self.trigger_A = self.trigger_tool.getTriggerChecked(
                        "EF_tau29_medium1_tau20_medium1_Hypo_00_02_42")
                self.trigger_B = self.trigger_tool.getTriggerChecked(
                        "EF_tau29_medium1_tau20_medium1_Hypo_00_03_02")
                self.trigger_C = self.trigger_tool.getTriggerChecked(
                        "EF_tau29T_medium1_tau20T_medium1_Hypo_00_03_02")

                self.trigger_run_dict = {
                    180164: (self.trigger_A, 'EF_tau29_medium1_tau20_medium1'),
                    183003: (self.trigger_B, 'EF_tau29_medium1_tau20_medium1'),
                    186169: (self.trigger_B, 'EF_tau29_medium1_tau20_medium1'),
                    189751: (self.trigger_C, 'EF_tau29T_medium1_tau20T_medium1'),
                }

                self.passes = self.passes_11
                self.finalize = self.finalize_11

        super(TauTriggerEmulation, self).__init__(
                passthrough=passthrough,
                **kwargs)

    def passes_11(self, event):

        self.trigger_tool_wrapper.setEventNumber(event._entry.value)
        trigger, triggername = self.trigger_run_dict[event.RunNumber]
        trigger.switchOn()
        self.trigger_tool.executeTriggers()

        if trigger.passed():
            if triggername == 'EF_tau29_medium1_tau20_medium1':
                event.EF_tau29_medium1_tau20_medium1 = True
            else:
                event.EF_tau29T_medium1_tau20T_medium1 = True

            # trigger matching
            trig1 = trigger.getTrigger1() # EF_tau29(T)_medium1
            trig2 = trigger.getTrigger2() # EF_tau20(T)_medium1
            for tau in event.taus:
                idx = -1
                idx1 = trig1.matchIndex(tau.fourvect)
                idx2 = trig2.matchIndex(tau.fourvect)
                if idx1 == idx2 != -1:
                    idx = idx1
                elif idx1 == -1 and idx2 > -1:
                    idx = idx2
                elif idx2 == -1 and idx1 > -1:
                    idx = idx1
                elif idx2 != idx1: # both >-1 and non-equal
                    # take index of closer one using dR
                    trigtau1TLV = self.trigger_tool.buildEFTauTLV(idx1)
                    trigtau2TLV = self.trigger_tool.buildEFTauTLV(idx2)
                    if trigtau1TLV.DeltaR(tau.fourvect) < trigtau2TLV.DeltaR(tau.fourvect):
                        idx = idx1
                    else:
                        idx = idx2
                tau.trigger_match_index = idx
        else:
            if triggername == 'EF_tau29_medium1_tau20_medium1':
                event.EF_tau29_medium1_tau20_medium1 = False
            else:
                event.EF_tau29T_medium1_tau20T_medium1 = False

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
