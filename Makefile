CC=g++
all: vim compile run

compile:
	${CC} conmain.cpp -o conmain -lm -DDEBUG=0

run:
	./conmain

vim:
	vim conmain.cpp

speed:
	time ./conmain &>/dev/null

clean:
	rm *.o conmain
