def OverlapCheck(event):
    """
    Check for overlap between object, so that not a single object fires both
    the tau and the lepton requirements
    """

    # Overlap Check
    Objects = [] # Format : [0:object 4-vector, 1:isTau, 2:isElectron, 3:isMuon]

    print '##############################'
    print 'n. taus     ', len(event.taus)
    print 'n. electrons', len(event.electrons)
    print 'n. muons    ', len(event.muons)

    #collect all tau objects:
    for tau in event.taus:
        Objects.append([tau, 1, 0, 0])

    #collect non-overalpping electrons
    for e in event.electrons:
        for tau in event.taus:
            if utils.dR(e.eta, e.phi, tau.eta, tau.phi) > 0.2:
                Objects.append([e, 0, 1, 0])

    #collect non-overalpping muons
    for mu in event.muons:
        for tau in event.taus:
            if utils.dR(mu.eta, mu.phi, tau.eta, tau.phi) > 0.2:
                Objects.append([mu, 0, 0, 1])

    #Check if collected objects also have other types
    for obj in Objects:
                    
        #Check if object is also a tau, if not already a tau
        if obj[1] == 0:
            for tau in event.taus:
                if utils.dR(obj[0].eta, obj[0].phi, tau.eta, tau.phi) < 0.2:
                    obj[1] = 1

        #Check if object is also an electron, if not already an electron
        if obj[2] == 0:
            for e in event.electrons:
                if utils.dR(obj[0].eta, obj[0].phi, e.eta, e.phi) < 0.2:
                    obj[2] = 1

        #Check if object is also a muon, if not already a muon
        if obj[3] == 0:
            for mu in event.muons:
                if utils.dR(obj[0].eta, obj[0].phi, mu.eta, mu.phi) < 0.2:
                    obj[3] = 1
                            
                    
        print '=============================='
        print 'isTau isElectron isMuon'
        for obj in Objects:
            print obj[1], obj[2], obj[3]
