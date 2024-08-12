CC = gcc
CPPC = g++
CFLAGS = -Wall -Werror -std=gnu99 -pedantic
CPPFLAGS = -Wall -Werror -pedantic -std=c++11
TEST = test
RUN = run
.PHONY: clean test run prof
.DEFAULT_GOAL := run

all: run

prof: CFLAGS += -g -pg
prof: CPPFLAGS += -g -pg
prof: run

run: run.o polynomial.o io.o fin_field.o fft.o
	$(CPPC) $(CPPFLAGS) run.o polynomial.o io.o fin_field.o fft.o -o $@

test: test.o polynomial.o io.o fin_field.o fft.o
	$(CC) $(CFLAGS) test.o polynomial.o io.o fin_field.o fft.o -o $@

fft.o: fft.c fft.h
	$(CC) $(CFLAGS) -c $<

fin_field.o: fin_field.c fin_field.h
	$(CC) $(CFLAGS) -c $<

polynomial.o: polynomial.c polynomial.h
	$(CC) $(CFLAGS) -c $<

run.o: run.cpp
	$(CPPC) $(CPPFLAGS) -c $<

test.o: test.c
	$(CC) $(CFLAGS) -c $<

io.o: io.c io.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(TEST) $(RUN)
