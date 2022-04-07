ILS_NKland : aux_functions.o file_man.o global.o ILS.o perturbation.o selection.o statistics.o 
	g++ -Wall aux_functions.o file_man.o global.o ILS.o perturbation.o selection.o statistics.o -o ILS_NKland

aux_functions.o : aux_functions.cpp	
	g++ -Wall -o aux_functions.o -c aux_functions.cpp

file_man.o : file_man.cpp	
	g++ -Wall -o file_man.o -c file_man.cpp

global.o : global.cpp	
	g++ -Wall -o global.o -c global.cpp

ILS.o : ILS.cpp	
	g++ -Wall -o ILS.o -c ILS.cpp

perturbation.o : perturbation.cpp	
	g++ -Wall -o perturbation.o -c perturbation.cpp

selection.o : selection.cpp	
	g++ -Wall -o selection.o -c selection.cpp

statistics.o : statistics.cpp	
	g++ -Wall -o statistics.o -c statistics.cpp




