CXX=g++ -std=c++17 -I
OPT=-O3

SWTD: tools.o LPFunction.o graph.o TDH2H.o TDGTree.o DHL.o main.o -lboost_system -lpthread -lboost_thread
	$(CXX) -g -o SWTD *.o -lboost_system -lpthread -lboost_thread


tools.o:tools.cpp
	$(CXX) -g -c $(OPT) tools.cpp

LPFunction.o:LPFunction.cpp
	$(CXX) -g -c $(OPT) LPFunction.cpp

TDH2H.o:TDH2H.cpp
	$(CXX) -g -c $(OPT) TDH2H.cpp

TDGTree.o:TDGTree.cpp
	$(CXX) -g -c $(OPT) TDGTree.cpp

DHL.o:DHL.cpp
	$(CXX) -g -c $(OPT) DHL.cpp

main.o:main.cpp
	$(CXX) -g -c $(OPT) main.cpp

graph.o:graph.cpp
	$(CXX) -g -c $(OPT) graph.cpp


clean:
	rm *.o
	rm SWTD
