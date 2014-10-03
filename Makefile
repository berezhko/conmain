all: compile run

compile:
	g++ gauss.c -c
	g++ -std=c++11 conmain.cpp -c
	g++ gauss.o conmain.o -o conmain -lm
run: compile
	./conmain
vim:
	vim conmain.cpp
