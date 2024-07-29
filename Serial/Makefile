CC = gcc
CFLAGS = -Wall -Werror -std=gnu99 -pedantic
TEST = test
.PHONY: clean test
.DEFAULT_GOAL := test

test: test.o fft.o fin_field.o io.o polynomial.o
	$(CC) $(CFLAGS) test.o fft.o fin_field.o io.o polynomial.o -o $@

fft.o: fft.c fft.h
	$(CC) $(CFLAGS) -c $<

fin_field.o: fin_field.c fin_field.h
	$(CC) $(CFLAGS) -c $<

polynomial.o: polynomial.c polynomial.h
	$(CC) $(CFLAGS) -c $<

test.o: test.c
	$(CC) $(CFLAGS) -c $<

io.o: io.c io.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(TEST)
