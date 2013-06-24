

HHSTUDENT ?= HHProcessor
HHNTUP ?= ntuples/hadhad/HHProcessor

.PHONY: dump

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
	find higgstautau -name "*.pyc" -exec rm {} \;
	rm -f *.pyc

clean: clean-pyc
	rm -f student*.root
	rm -f student*.profile
	rm -f TPileupReweighting.prw.root
	rm -f nohup.out
	rm -f supervisor*.log
	rm -f grid-setup.sh
	rm -f *.dot
	rm -f *.e[0-9]*
	rm -f *.o[0-9]*

dump:
	@./dump -t tau -s "hh_taus_pass && (RunNumber==207528)" --select-file etc/embed_select_ac.txt -o RunNumber,EventNumber hhskim.embed12_p1344_hadhad.root
