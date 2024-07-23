CC = gcc
CFLAGS = -Wall -Werror -std=gnu99 -pedantic
TEST = test
.PHONY: test clean debug
.DEFAULT_GOAL := test

debug: CFLAGS += -g
debug: test

test: test.o fft.o fin_field.o
	$(CC) $(CFLAGS) test.o fft.o fin_field.o -o $@

fft.o: fft.c fft.h
	$(CC) $(CFLAGS) -c $<

fin_field.o: fin_field.c fin_field.h
	$(CC) $(CFLAGS) -c $<
