# Include CPLEX libraries.
CPX_INC = -I/opt/ibm/ILOG/CPLEX_Studio1210/cplex/include/ -I/opt/ibm/ILOG/CPLEX_Studio1210/cplex/include/ilcplex -I/opt/ibm/ILOG/CPLEX_Studio1210/concert/include/
CPX_LIB = -L/opt/ibm/ILOG/CPLEX_Studio1210/concert/lib/x86-64_linux/static_pic/ -L/opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/


.cpp.o:
	g++ -DIL_STD -O3 -c -I. $(CPX_INC) $<

# The following two lines are not necessary and can be removed.
#main.o Problem.o: Problem.h

#Common.o: Common.h

# Change the file names before .o below to the cpp file names you want to compile.
myexe: main.o GrandSubproblem2.o ThirdPhase.o Dijkstra.o GrandSubProbMasterPICEF.o GrandProblem.o GrandSubproblem.o Reading.o Class_Problem.o
	g++ -o 2StageROThisWork -DIL_STD $(CPX_LIB) -O3 $+ -lilocplex -lcplex -lconcert -lm -lpthread -ldl

# Finally rename this file to "Makefile", put this file to the same folder of the cpp files. Change to the directory of this file, and type in `make` in the Terminal. After that you should get an executable named "myexe". In the Terminal type in `.\myexe` to run the executable.

#-I/usr/local/c++/4.4.1/
