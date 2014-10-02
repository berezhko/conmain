all: compile run

compile:
	g++ -std=c++11 conmain.cpp -o conmain
run:
	./conmain
