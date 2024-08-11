CC = gcc
CFLAGS = -Wall -Werror -std=gnu99 -pedantic
TEST = test
RUN = run
.PHONY: clean test run
.DEFAULT_GOAL := run

run: run.o polynomial.o io.o fin_field.o fft.o
	$(CC) $(CFLAGS) run.o polynomial.o io.o fin_field.o fft.o -o $@

test: test.o polynomial.o io.o fin_field.o fft.o
	$(CC) $(CFLAGS) test.o polynomial.o io.o fin_field.o fft.o -o $@

fft.o: fft.c fft.h
	$(CC) $(CFLAGS) -c $<

fin_field.o: fin_field.c fin_field.h
	$(CC) $(CFLAGS) -c $<

polynomial.o: polynomial.c polynomial.h
	$(CC) $(CFLAGS) -c $<

run.o: run.c
	$(CC) $(CFLAGS) -c $<

test.o: test.c
	$(CC) $(CFLAGS) -c $<

io.o: io.c io.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(TEST) $(RUN)
