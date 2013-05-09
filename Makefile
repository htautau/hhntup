

HHSTUDENT ?= HHProcessor
HHNTUP ?= ntuples/hadhad/HHProcessor

default: clean lib

lib:
	cd higgstautau/jetcleaning && ./setup.py build_ext --inplace

kill-hung:
	-./pbs.py | grep HUNG | cut -d " " -f 1 | xargs qdel

hh-ntup-clean:
	rm -f $(HHNTUP)/$(HHSTUDENT).root
	rm -f $(HHNTUP)/$(HHSTUDENT).h5

hh-ntup-merge: hh-ntup-clean
	./ntup-merge -s $(HHSTUDENT) $(HHNTUP)
	(cd $(HHNTUP) && root2hd5 $(HHSTUDENT).root)

clean-pyc:                                                                      
	find higgstautau -name "*.pyc" | xargs rm -f
	rm -f *.pyc

clean: clean-pyc
	rm -f student*.root
	rm -f student*.profile
	rm -f TPileupReweighting.prw.root
	rm -f nohup.out
	rm -f supervisor*.log
	rm -f grid-setup.sh
	rm -f *.dot
