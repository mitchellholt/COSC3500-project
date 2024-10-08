#! /bin/bash

rm -f code.zip
zip code *.cpp *.h *.c Makefile goslurm*
scp code.zip rangpur:~/COSC3500/Milestone2
