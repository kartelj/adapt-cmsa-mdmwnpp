CC = g++
CFLAGS = -Wall -g -o4
 
# ****************************************************
# Targets needed to bring the executable up to date
 
main: VNS_GA.o Problem.o VNS.o GA.o
	$(CC) $(CFLAGS) -o VNS_GA VNS_GA.o Problem.o VNS.o GA.o
 
# The main.o target can be written more simply
 
main.o: VNS_GA.cpp Problem.h VNS.h GA.o
	$(CC) $(CFLAGS) -c VNS_GA.cpp
 
Problem.o: Problem.h
 
VNS.o: VNS.h Problem.h

GA.o: GA.h Problem.h
