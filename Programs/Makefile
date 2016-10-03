#**********************************************#
# Title: MSOAR 2.0 									           #
# Description: To generate executable files    #
# 			 that are used in MSOAR2.0 script      #
# Author: Guanqun Shi                          #
# Copyright: Copyright (c) 2009                # 
#**********************************************# 

CC = g++

all: mapToPhylip familySeparator addOutgroup getTopHits postprocessing removeInparalogs normalizeScores tagGenerator treeTran

mapToPhylip:
	$(CC)  MapToPhylip.cpp -o mapToPhylip

familySeparator:
	$(CC)  FamilySeparator.cpp -o familySeparator

addOutgroup:
	$(CC)  AddOutgroup.cpp -o addOutgroup

getTopHits:
	$(CC)  GetTopHits.cpp -o getTopHits

treeTran:
	$(CC)  TreeTran.cpp -o treeTran

postprocessing:
	$(CC)  Postprocessing.cpp -o postprocessing

removeInparalogs:
	$(CC)  RemoveInparalogs.cpp -o removeInparalogs

normalizeScores:
	$(CC)  NormalizeScores.cpp -o normalizeScores

tagGenerator:
	$(CC)  TagGenerator.cpp -o tagGenerator

clean: 
	rm -f *~ mapToPhylip familySeparator addOutgroup getTopHits postprocessing removeInparalogs normalizeScores tagGenerator treeTran
