CC=g++
all: vim compile run

compile:
	${CC} gauss.cpp -c
	${CC} conmain.cpp -c -DDEBUG=0
	${CC} gauss.o conmain.o -o conmain -lm

run:
	./conmain

vim:
	vim conmain.cpp

speed:
	time ./conmain &>/dev/null

clean:
	rm *.o conmain
