#include "ggF_cross_section_uncertainty/ggF_XSecUncertTool.h"
#include <TROOT.h>

// constructor
ggF_XSecUncertTool::ggF_XSecUncertTool(Str method, Str configFile, bool verbose)
  : m_input(0), m_inputFN(configFile), m_method(method), m_verb(verbose) {
  init(configFile);
}

void ggF_XSecUncertTool::init(Str configFile) {

  printf("\n======================================\n");
  printf("========= ggF_XSecUncertTool =========\n");
  printf("  A tool for propagation of the ggF cross section uncertainty\n\n");
  
  // 1. open the input file
  m_input = new TEnv();
  int status = m_input->ReadFile(m_inputFN.Data(),EEnvLevel(0));
  if (status!=0) fatal("Cannot open "+m_inputFN);
  
  printf("  Reading input from %s\n\n",m_inputFN.Data());

  // Read in general information
  m_Njetbins = (int)getValue("NjetBins");
  m_pTcut = getValue("pTcut");
  m_jetR  = getValue("AntiKtDistanceParam");
  
  if (m_Njetbins!=3) fatal("Only support for 3 jet bins at the moment (NjetBins)");

  printf("  Jet definition:\n");
  printf("    anti-kt R = %.1f, pT > %.0f GeV\n",m_jetR,m_pTcut);
  printf("    built from truth particles, excluding jets originating\n");
  printf("    from the Higgs decay products (apply overlap removal)\n\n");

  // create an emtpy 2D-vector of Njetbins uncertainty components each with Njetbins uncertainty amplitudes
  // to be filled according to each method below
  for (int ijet=0;ijet<m_Njetbins;++ijet) m_relUncert.push_back(VecD(m_Njetbins,0.0));

  printf("  Initializing ggF uncertainties according to the %s\n\n",m_input->GetValue(m_method+".Description",""));
  if (m_method=="JVE") {
        
    // number of efficienies needed
    int Neffs = m_Njetbins-1;

    // Read in input parameters with accosiated relative uncertainites
    VecD eff, deff;
    double sigTot = getValue("JVE.SigmaTot"), dsigTot = getValue("JVE.SigmaTotRelUnc");
    for (int ijet=0;ijet<Neffs;++ijet) {
      add(eff,getValue(Form("JVE.Eff%d",ijet)));
      add(deff,getValue(Form("JVE.Eff%dRelUnc",ijet)));
    }
    
    printf("  %s\n",m_input->GetValue(m_method+".InputDetails",""));
    printf("    Input parameters:   sigTot  = %4.1f pb,",sigTot);
    for (int ieff=0;ieff<Neffs;++ieff) printf("  eff%d = %4.1f%%%s",ieff,eff[ieff]*100,ieff==Neffs-1?"\n":",");
    printf("    Rel. uncertainties: dsigTot =  %5.1f%%,",dsigTot*100);
    for (int ieff=0;ieff<Neffs;++ieff) printf(" deff%d = %4.1f%%%s",ieff,deff[ieff]*100,ieff==Neffs-1?"\n":",");
    
    // Calculating fractional uncertainty amplitudes according to the JVE method.
    // For formulas for the 3x3 case, see talk by David Hall, page 5:
    //   https://indico.cern.ch/conferenceDisplay.py?confId=296612
    //
    // Generalization to more jet bins is as follows:
    // the relative yield uncertainty is always dsigTot in all jet bins
    // the relative migration uncertainty between jet bins
    //   {i} -> {i+1}, i.e. =i jets -> i+1 or more jets
    // is dermined by eff{i} (P for =i jets given at least i jets)
    // its associated realtive uncertainty deff{i}, and become:
    //   0                           for <i jets (A)
    //   deff{i}                     for =i jets (B)
    //  -eff{i}/(1-eff{i})*deff{i}   for >i jets (C)
    //
    // Note the sign change across the {i}->{i+1} Njet boundary
    //
    // implemenation:
    for (int ijet=0;ijet<m_Njetbins;++ijet) {
      m_relUncert[0][ijet] = dsigTot;
      for (int ieff=0;ieff<Neffs;++ieff) {
	if (ijet<ieff) continue; // case (A) above
	int icomp = ieff+1;
	//printf("ieff = %d, ijet = %d\n",ieff,ijet);
	m_relUncert[icomp][ijet] =
	  ieff==ijet ? deff[ieff] : -eff[ieff]/(1.0-eff[ieff])*deff[ieff]; // cases (B) and (C) above
      }
    }

    // Calculate the corresponding cross exclusive and incluisve cross sections
    // (not directly needed for JVE uncertainies)
    for (int ijet=0;ijet<m_Njetbins;++ijet) {
      double sig=sigTot;

      // remove xsec for less than ijet jets
      for (int jjet=0;jjet<ijet;++jjet) sig*=(1.0-eff[jjet]);
      add(m_sigIncl,sig);
      
      // for all but last bin, remove xsec for >=ijet+1 jets
      if (ijet<m_Njetbins-1) sig*=eff[ijet];
      add(m_sigExcl,sig);
    }

  } else if (m_method=="RIST") {

    int Njb = m_Njetbins; // being lazy

    m_sigExcl = getValues("RIST.SigmaJetBins");
    if (m_sigExcl.size()!=m_Njetbins) fatal("Too few entries in RIST.SigmaJetBins");

    // sig_{\geq i} = sig_{i} + sig_{\geq i+1}
    m_sigIncl = m_sigExcl;
    for (int i=m_Njetbins-2;i>=0;--i) m_sigIncl[i] += m_sigIncl[i+1];

    VecD DsigMu = getValues("RIST.DSigmaMu");
    if (DsigMu.size()!=m_Njetbins) fatal("Too few entries in RIST.DSigmaMu");
    
    double D0cut = getValue("RIST.Delta0cut"),
      D1cut  = getValue("RIST.Delta1cut"),
      rho = getValue("RIST.Rho");
    
    // Dump input variables used to the screen
    printf("  %s\n",m_input->GetValue(m_method+".InputDetails",""));
    printf("    Cross sections:    ");
    for (size_t i=0;i<Njb;++i) printf("%8.2f pb%s",m_sigExcl[i],i==Njb-1?"\n":",");
    printf("    Yield uncert:      ");
    for (size_t i=0;i<Njb;++i) printf("%8.2f pb%s",DsigMu[i],i==Njb-1?"\n":",");
    printf("    Jet bin migration: %8.3f pb,%8.3f pb\n",D0cut,D1cut);
    printf("    Correlation param: %8.3f\n",rho);
    
    for (int ijet=0;ijet<Njb;++ijet) // yield part
      m_relUncert[0][ijet]=DsigMu[ijet]/m_sigExcl[ijet]; 

    // John's mig01 and mig12
    m_relUncert[1][0]=D0cut; m_relUncert[1][1]=-(1-rho)*D0cut; m_relUncert[1][2]=-rho*D0cut;
    m_relUncert[2][0]=0; m_relUncert[2][1]=D1cut; m_relUncert[2][2]=-D1cut;
    for (int ijet=0;ijet<Njb;++ijet) {
      m_relUncert[1][ijet] = m_relUncert[1][ijet]/m_sigExcl[ijet];
      m_relUncert[2][ijet] = m_relUncert[2][ijet]/m_sigExcl[ijet];
    }

    if (Njb>3) fatal("More than 3 jet bins not yet supported for RIST. Should be able to use ST for this though ...");
    
  } else if (m_method=="ST") {
    int Njb = m_Njetbins; // being lazy
    
    m_sigIncl = getValues("ST.SigmaIncl");
    if (m_sigIncl.size()!=m_Njetbins) fatal("Too few entries in ST.SigmaIncl");
    
    // sig_{\geq i} = sig_{i} + sig_{\geq i+1}
    m_sigExcl = m_sigIncl;
    for (int i=0;i<Njb-1;++i) m_sigExcl[i]=m_sigIncl[i]-m_sigIncl[i+1];
    VecD DsigIncl = getValues("ST.DSigmaIncl");
    
    // Dump input variables used to screen
    printf("  %s\n",m_input->GetValue(m_method+".InputDetails",""));
    printf("    Incl. jet-bin x-secs:");
    for (size_t i=0;i<Njb;++i) printf("%8.2f pb%s",m_sigIncl[i],i==Njb-1?"\n":",");
    printf("    Absolute uncert:     ");
    for (size_t i=0;i<Njb;++i) printf("%8.2f pb%s",DsigIncl[i],i==Njb-1?"\n":",");
    
    // Fixed-order ST
    m_relUncert[0][0] = DsigIncl[0]/m_sigExcl[0];
    for (int imig=1;imig<Njb;++imig) {
      double Dmig = DsigIncl[imig];
      m_relUncert[imig][imig-1]=  Dmig/m_sigExcl[imig-1];
      m_relUncert[imig][imig]  = -Dmig/m_sigExcl[imig];
    }
    
  } else
    fatal("Don't know how to init ggF uncertainty parametrization according to method "+m_method);
  
  // print the calculated uncertainty components and correlation details
  printNPs();
  if (m_verb) printCorrelationDetails();

  printf("\n\n======================================\n");
}

// access value from the TEnv database 
double ggF_XSecUncertTool::getValue(Str key) {
  double val = m_input->GetValue(key,-9999.0);
  if (val==-9999) fatal("Cannot access value for "+key+" in file: "+m_inputFN);
  return val;
}

// access vector of values from the TEnv database
VecD ggF_XSecUncertTool::getValues(Str key) {
  VecD vals; Str valStr = m_input->GetValue(key,"");
  if (valStr=="") fatal("Cannot access value for "+key+" in file: "+m_inputFN);
  // split string at " " and put each piece in the vector of doubles
  TObjArray *strings = valStr.Tokenize(" ");
  TIter istr(strings);
  while (TObjString* os=(TObjString*)istr()) add(vals,atof(os->GetString())); 
  delete strings;
  return vals;
}

void ggF_XSecUncertTool::printNP(Str uName, const VecD &relU) {
  printf("  %20s:",uName.Data());
  for (size_t bin=0;bin<relU.size();++bin) printf("%12.3f",relU[bin]);
  printf("\n");
}

void ggF_XSecUncertTool::printNPs() {
  printf("\n  Fractional unc. amplitude for each ggF uncertainty component (nuisance parameter):\n");
  printf("  %20s ","");
  for (int njet=0;njet<m_Njetbins;++njet) 
    if (njet+1==m_Njetbins) printf("%12s\n",Form(">=%d jets",njet));
    else printf("%12s",Form("=%d jets",njet));
  printNP("yield uncertainty",m_relUncert[0]);
  for (int ijet=1;ijet<m_Njetbins;++ijet) 
    printNP(Form("migration %d->%d",ijet-1,ijet),m_relUncert[ijet]);
  printf("  %26s-------------------------------\n","");
  printNP("total",totalXsecUncert());
}

double ggF_XSecUncertTool::relUncert(int iComp, int Ntruthjets) {
  if (iComp>=m_Njetbins||iComp<0) 
    fatal(Form("ggF_XSecUncertTool::relUncert Cannot access uncertainty component %d - there are only %d.",
	       iComp,m_Njetbins));
  if (Ntruthjets<0) 
    fatal(Form("ggF_XSecUncertTool::relUncert %d truth jets makes no sense. Don't feed me junk ;)",Ntruthjets));

  // return uncertainty of
  //  component iComp
  //  for jet bin: max(Ntruthjets,Njetbins), i.e. last bin is inclusive
  return m_relUncert[iComp][ Ntruthjets >= m_Njetbins ? m_Njetbins-1 : Ntruthjets ];
}

double ggF_XSecUncertTool::totalXsecUncert(int Ntruthjets) {
  double sum2=0;
  for (size_t icomp=0;icomp<m_Njetbins;++icomp)
    sum2+=pow(relUncert(icomp,Ntruthjets),2);
  return sqrt(sum2);
}

VecD ggF_XSecUncertTool::totalXsecUncert() {
  VecD unc;
  for (size_t ijet=0;ijet<m_Njetbins;++ijet)
    add(unc,totalXsecUncert(ijet));
  return unc;
}

void ggF_XSecUncertTool::printCorrelationDetails() {
  printf("\n  Covariance matrix (between jet bins):\n");
  for (int i=0;i<m_Njetbins;++i) 
    for (int j=0;j<m_Njetbins;++j)
      printf("%10.2f pb^2%s",cov(i,j),j+1==m_Njetbins?"\n":"");

  printf("\n  Correlation matrix:\n");
  for (int i=0;i<m_Njetbins;++i) 
    for (int j=0;j<m_Njetbins;++j)
      printf("%12.1f%%%s",corr(i,j)*100,j+1==m_Njetbins?"\n":"");


  double D0cut = migration01Uncert(0)*m_sigExcl[0],
    D1cut = migration12Uncert(1)*m_sigExcl[1],
    rho = -migration01Uncert(2)*m_sigExcl[2]/D0cut;
  printf("\n  Parametrization suggested by J. Walsh et al:\n");
  printf(  "    D0cut (0->1 migration):    %.3f pb\n",D0cut);
  printf(  "    D1cut (1->2 migration):    %.3f pb\n",D1cut);
  printf(  "    Rho   (0->2 migration):    %.3f\n",rho);
}

// For fun
double ggF_XSecUncertTool::cov(int ijet, int jjet) {
  int Nj=m_Njetbins; double cov_ij=0;
  for (int ic=0;ic<Nj;++ic)
    cov_ij+=m_relUncert[ic][ijet]*m_relUncert[ic][jjet]*m_sigExcl[ijet]*m_sigExcl[jjet];
  return cov_ij;
}
double ggF_XSecUncertTool::corr(int ijet, int jjet) {
  if (ijet==jjet) return 1.0;
  return cov(ijet,jjet)/sqrt(cov(ijet,ijet)*cov(jjet,jjet));
}
