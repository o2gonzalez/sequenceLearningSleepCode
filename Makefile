all: mindcurrent

mindcurrent: main.cpp CellSyn.cpp CellSyn.h currents.cpp currents.h io.cpp io.h network.h network.cpp bmalgo.h bmalgo.cpp
	g++ -g -O3 -lm -Wall -fopenmp  main.cpp CellSyn.cpp io.cpp currents.cpp network.cpp bmalgo.cpp -o mindcurrent
phi:
	icpc -mmic -O3 -Wall -fopenmp main.cpp CellSyn.cpp io.cpp currents.cpp network.cpp -o mindcurrent


#sanity check whether current version of code changes output (compared to previously stored test/.files)
check: mindcurrent
	./mindcurrent test/params.txt test connection_info2
	cd test; for f in *; do diff -u $$f .$$f; done

check-prepare: mindcurrent
	./mindcurrent test/params.txt test connection_info2
	cd test; for f in *; do cp -f  $$f .$$f; done

run: mindcurrent
	./mindcurrent params.txt out connection_info2

doxy:
	doxygen ./docs/Doxyfile

clean: 
	-rm mindcurrent 

network:
	g++ -O2 generate_network.cpp -o generate_network
	./generate_network > connection_info2

network_mri:
	g++ -O2 -g generate_network.cpp -o generate_network
	./generate_network /mb/home/psanda/DistanceFiles/DistanceMatrix_18-Dec-2013_LH.txt /mb/home/psanda/DistanceFiles/Map_1600_To_25k_18-Dec-2013_LH.txt > connection_info2
