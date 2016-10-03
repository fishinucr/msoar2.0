#**********************************************#
# Title: M-SOAR                                #
# Description: System for orthologs assignment #
# Author: Zheng Fu                             #
# Copyright: Copyright (c) 2005                # 
#**********************************************# 

CC = g++
FLAGS = -g -Wall 

MSOAR: main.o RelatedGenomes.o CompleteGraph.o PairmatchGraph.o
	$(CC) $(FLAGS) main.o RelatedGenomes.o CompleteGraph.o PairmatchGraph.o -o MSOAR

main.o: main.cc RelatedGenomes.h Genome.h Gene.h BlastHit.h Graph.h CompleteGraph.h PairmatchGraph.h VertexPair.h
	$(CC) $(FLAGS)-c main.cc

RelatedGenomes.o: RelatedGenomes.cc RelatedGenomes.h Genome.h Gene.h BlastHit.h Graph.h CompleteGraph.h PairmatchGraph.h VertexPair.h
	$(CC) $(FLAGS) -c RelatedGenomes.cc

CompleteGraph.o: CompleteGraph.cc CompleteGraph.h Graph.h Genome.h VertexPair.h
	$(CC) $(FLAGS) -c CompleteGraph.cc

PairmatchGraph.o: PairmatchGraph.cc PairmatchGraph.h Graph.h Genome.h VertexPair.h
	$(CC) $(FLAGS) -c PairmatchGraph.cc

clean: 
	rm -rf *.o *~ MSOAR
