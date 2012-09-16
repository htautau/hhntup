
default: clean

clean:
	rm -f *.pyc
	rm -f cutflow.p
	rm -f student*.root
	rm -f student*.profile
	rm -f TPileupReweighting.prw.root
	rm -f nohup.out
	rm -f supervisor*.log
	rm -f grid-setup.sh
