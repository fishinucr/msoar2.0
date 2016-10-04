#***********************************************#
# Title: MSOAR 2.0                              #
# Description: To generate executable files     #
#              that are used in MSOAR2.0 script #
# Author: Guanqun Shi                           #
# Copyright: Copyright (c) 2009                 # 
#***********************************************# 

CC = g++
FLAGS = -O3

all: mapToPhylip familySeparator addOutgroup getTopHits postprocessing removeInparalogs normalizeScores tagGenerator treeTran

mapToPhylip:
	$(CC) $(FLAGS) MapToPhylip.cpp -o mapToPhylip

familySeparator:
	$(CC) $(FLAGS) FamilySeparator.cpp -o familySeparator

addOutgroup:
	$(CC) $(FLAGS) AddOutgroup.cpp -o addOutgroup

getTopHits:
	$(CC) $(FLAGS) GetTopHits.cpp -o getTopHits

treeTran:
	$(CC) $(FLAGS) TreeTran.cpp -o treeTran

postprocessing:
	$(CC) $(FLAGS) Postprocessing.cpp -o postprocessing

removeInparalogs:
	$(CC) $(FLAGS) RemoveInparalogs.cpp -o removeInparalogs

normalizeScores:
	$(CC) $(FLAGS) NormalizeScores.cpp -o normalizeScores

tagGenerator:
	$(CC) $(FLAGS) TagGenerator.cpp -o tagGenerator

.PHONY: clean
clean: 
	rm -f *~ mapToPhylip familySeparator addOutgroup getTopHits postprocessing removeInparalogs normalizeScores tagGenerator treeTran
