#Makefile

# CXX = icc
# CXXFLAGS = -Wall -parallel -lm -fPIC

CC=gcc
all:rs

rs:graph.o bmp.o ford_fulkerson.o queue_stack.o rs_drv.o rangeswap.o -lm
	$(CC) graph.o bmp.o ford_fulkerson.o queue_stack.o rs_drv.o rangeswap.o -o rs -lm -g3
.c.o:
	$(CC) -c $< -Wall -l -g3
remove:
	rm -f *.o
clean:
	rm -f *.o *~ rs