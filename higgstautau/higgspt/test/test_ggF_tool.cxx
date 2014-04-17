///#include "ggF_cross_section_uncertainty/ggF_XSecUncertTool.h"

// being lazy: to be able to compile in one line
#include "src/ggF_XSecUncertTool.cxx"

int main() {

  ggF_XSecUncertTool *ggF_uncert_jve  = new ggF_XSecUncertTool("JVE");
  ggF_XSecUncertTool *ggF_uncert_rist = new ggF_XSecUncertTool("RIST");
  ggF_XSecUncertTool *ggF_uncert_st   = new ggF_XSecUncertTool("ST");
  
  // Test
  printf("\n~~~~~~~~~~\nTest of the tool\n\n");
  printf("  For an event with 1 pT>25 GeV truth jet, the\n");
  printf("  event weight for uncertainty propagation is:\n");

  printf("    Yield JVE,  up: %.3f, down: %.3f\n",
	 ggF_uncert_jve->yieldShiftUp(1),ggF_uncert_jve->yieldShiftDown(1));
  printf("    Yield RIST, up: %.3f, down: %.3f\n",
	 ggF_uncert_rist->yieldShiftUp(1),ggF_uncert_rist->yieldShiftDown(1));
  printf("    Yield ST,   up: %.3f, down: %.3f\n",
	 ggF_uncert_st->yieldShiftUp(1),ggF_uncert_st->yieldShiftDown(1));
  
  printf("\n~~~~~~~~~~\n");
}
