CC = gcc
CPPC = g++
CFLAGS = -Wall -Werror -std=gnu99 -pedantic
CPPFLAGS = -Wall -Werror -pedantic -std=c++11
.PHONY: clean test run prof build
.DEFAULT_GOAL := run

all: run

prof: CFLAGS += -g -pg -O1
prof: CPPFLAGS += -g -pg -O1
prof: run

run: CFLAGS += -O3
run: CPPFLAGS += -O3
run: build

build: run.o polynomial.o polynomial.h io.o io.h fin_field.o fin_field.h fft.o fft.h
	$(CPPC) $(CPPFLAGS) run.o polynomial.o io.o fin_field.o fft.o -o run

test: test.o polynomial.o io.o fin_field.o fft.o
	$(CC) $(CFLAGS) test.o polynomial.o io.o fin_field.o fft.o -o test

fft.o: fft.c
	$(CC) $(CFLAGS) -c $<

fin_field.o: fin_field.c
	$(CC) $(CFLAGS) -c $<

polynomial.o: polynomial.c
	$(CC) $(CFLAGS) -c $<

run.o: run.cpp
	$(CPPC) $(CPPFLAGS) -c $<

test.o: test.c
	$(CC) $(CFLAGS) -c $<

io.o: io.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o test run
