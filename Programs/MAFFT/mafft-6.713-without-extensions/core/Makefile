PREFIX = /usr/local
LIBDIR = $(PREFIX)/lib/mafft
BINDIR = $(PREFIX)/bin
MANDIR = $(PREFIX)/man/man1

#MNO_CYGWIN = -mno-cygwin

CC = gcc

#CFLAGS =  -O0  -pedantic -Wall -std=c99 -g -pg -DMALLOC_CHECK_=3
CFLAGS = -O3 $(MNO_CYGWIN)
#CFLAGS =  -mno-cygwin -O3 # for cygwin
#CFLAGS =  -fast -static      # for icc on Intel
#CFLAGS =  -O3 -ipo -static   # for icc on AMD
#CFLAGS =  -fast
#CFLAGS =  -g -O0 -pg 
#CFLAGS =   -ansi -pedantic
LIBS = -lm 

INSTALL = install

PROGS = dvtditr dndfast7 dndblast sextet5 mafft-distance pairlocalalign \
		pair2hat3s multi2hat3s rnatest pairash \
        splittbfast disttbfast tbfast mafft-profile f2cl mccaskillwrap contrafoldwrap countlen \
		seq2regtable regtable2seq score getlag dndpre dndpre2 setcore 
SCRIPTS = mafft mafft-homologs.rb
OBJREGTABLE2SEQ = mtxutl.o io.o regtable2seq.o defs.o mltaln9.o
OBJSEQ2REGTABLE = mtxutl.o io.o seq2regtable.o defs.o
OBJCOUNTLEN = mtxutl.o io.o countlen.o defs.o
OBJF2CL = mtxutl.o io.o f2cl.o constants.o defs.o
OBJMCCASKILLWRAP = mtxutl.o io.o mccaskillwrap.o constants.o defs.o mltaln9.o
OBJCONTRAFOLDWRAP = mtxutl.o io.o contrafoldwrap.o constants.o defs.o mltaln9.o
OBJMULTI2HAT3S = mtxutl.o io.o mltaln9.o tddis.o constants.o \
				 multi2hat3s.o defs.o fft.o fftFunctions.o
OBJPAIR2HAT3S = mtxutl.o io.o mltaln9.o tddis.o constants.o \
				pair2hat3s.o defs.o fft.o fftFunctions.o
OBJRNATEST = mtxutl.o io.o mltaln9.o tddis.o constants.o Lalignmm.o \
			 rna.o rnatest.o defs.o fft.o fftFunctions.o
OBJPAIRASH = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
			 Falign.o MSalignmm.o Galign11.o MSalign11.o suboptalign11.o genGalign11.o genalign11.o Lalign11.o SAalignmm.o \
			 Calignm1.o pairash.o defs.o fft.o fftFunctions.o
OBJPAIRLOCALALIGN = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		            Falign.o MSalignmm.o Galign11.o MSalign11.o suboptalign11.o genGalign11.o genalign11.o Lalign11.o SAalignmm.o \
					Calignm1.o pairlocalalign.o defs.o fft.o fftFunctions.o
OBJDUMMY = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o Calignm1.o MSalignmm.o \
			disttbfast_dummy.o dummy.o defs.o fft.o fftFunctions.o
OBJSPLITFROMALN = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o Calignm1.o MSalignmm.o \
			Lalign11.o splitfromaln.o defs.o fft.o fftFunctions.o
OBJSPLITTBFAST = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o Calignm1.o MSalignmm.o \
			Lalign11.o splittbfast.o defs.o fft.o fftFunctions.o
OBJSPLITTBFAST2 = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o Calignm1.o MSalignmm.o \
			Lalign11.o splittbfast2.o defs.o fft.o fftFunctions.o
OBJSPLITTBFASTP = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o Calignm1.o MSalignmm.o \
			Lalign11.o defs.o fft.o fftFunctions.o
OBJDISTTBFAST = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o Calignm1.o MSalignmm.o \
			disttbfast.o defs.o fft.o fftFunctions.o
OBJTBFAST = mtxutl.o io.o mltaln9.o tddis.o constants.o MSalignmm.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o Calignm1.o  \
			tbfast.o defs.o fft.o fftFunctions.o
OBJTBFAST2 = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o MSalignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o Calignm1.o \
			tbfast2.o defs.o fft.o fftFunctions.o
OBJSETCORE = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o Calignm1.o MSalignmm.o \
			setcore.o defs.o fft.o fftFunctions.o
OBJTDITR = mtxutl.o io.o mltaln9.o tddis.o constants.o nj.o  partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		   Calignm1.o Falign.o Falign_localhom.o Galign11.o fftFunctions.o fft.o  \
		   tditeration.o tditr.o defs.o SAalignmm.o treeOperation.o
OBJDVTDITR = mtxutl.o io.o mltaln9.o tddis.o constants.o nj.o  partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
			Calignm1.o Falign.o Falign_localhom.o Galign11.o MSalignmm.o fftFunctions.o fft.o  \
			tditeration.o dvtditr.o defs.o SAalignmm.o treeOperation.o
OBJGETLAG = mtxutl.o io.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
		    Falign.o Falign_localhom.o Galign11.o SAalignmm.o Calignm1.o MSalignmm.o \
			getlag.o defs.o fft.o fftFunctions.o
OBJGAPFILL = mtxutl.o io.o constants.o gapfill.o defs.o
OBJDNDFAST5 = dndfast5.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJDNDBLAST = dndblast.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJDNDFAST7 = dndfast7.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJDNDFAST6 = dndfast6.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJDNDFAST4 = dndfast4.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJDNDFAST6 = dndfast6.o io.o constants.o mtxutl.o mltaln9.o tddis.o defs.o
OBJSEXTET5 = io.o constants.o mtxutl.o mltaln9.o tddis.o sextet5.o defs.o
OBJDISTANCE = io.o constants.o mtxutl.o mltaln9.o tddis.o mafft-distance.o defs.o
OBJTRIPLET6 = io.o constants.o mtxutl.o mltaln9.o tddis.o triplet6.o defs.o
OBJTRIPLET5 = io.o constants.o mtxutl.o mltaln9.o tddis.o triplet5.o defs.o
OBJOCTET4 = io.o constants.o mtxutl.o mltaln9.o tddis.o octet4.o defs.o
OBJDNDPRE = dndpre.o io.o constants.o mtxutl.o mltaln9.o defs.o
OBJDNDPRESCORE = dndpre2.o io.o constants.o mtxutl.o mltaln9.o defs.o
OBJGALN = io.o mtxutl.o mltaln9.o tddis.o constants.o partQalignmm.o partSalignmm.o MSalignmm.o Lalignmm.o rna.o Salignmm.o Ralignmm.o Qalignmm.o Halignmm.o \
          Calignm1.o SAalignmm.o Galign11.o Falign.o Falign_localhom.o fftFunctions.o fft.o mafft-profile.o defs.o
OBJSCORE = io.o mtxutl.o mltaln9.o score.o constants.o defs.o

HEADER = mltaln.h mtxutl.h
DPHEADER = dp.h
FFTHEADER = fft.h

MANPAGES = mafft.1 mafft-homologs.1



all : $(PROGS) $(SCRIPTS)
	cp $(SCRIPTS) ../scripts
	chmod 755 ../scripts/*
	cp $(PROGS) $(MANPAGES) ../binaries
	chmod 755 ../binaries/*
	@echo done.

univscript: univscript.tmpl Makefile
	sed "s:_PROGS:$(PROGS):" univscript.tmpl  > univscript

mafft: mafft.tmpl
	sed "s:_LIBDIR:$(LIBDIR):" mafft.tmpl  > mafft

mafft-homologs.rb: mafft-homologs.tmpl
#	cp mafft-homologs.tmpl mafft-homologs.rb
	sed "s:_BINDIR:$(BINDIR):" mafft-homologs.tmpl  > mafft-homologs.rb

mltaln.h : functions.h
	touch mltaln.h

ord : $(OBJORD)
	$(CC) -o $@ $(OBJORD) $(CFLAGS) $(LIBS)

tbfast : $(OBJTBFAST)
	$(CC) -o $@ $(OBJTBFAST) $(CFLAGS) $(LIBS)

tbfast2 : $(OBJTBFAST2)
	$(CC) -o $@ $(OBJTBFAST2) $(CFLAGS) $(LIBS)

disttbfast : $(OBJDISTTBFAST)
	$(CC) -o $@ $(OBJDISTTBFAST) $(CFLAGS) $(LIBS)

splittbfast : $(OBJSPLITTBFAST)
	$(CC) -o $@ $(OBJSPLITTBFAST) $(CFLAGS) $(LIBS)

splitfromaln : $(OBJSPLITFROMALN)
	$(CC) -o $@ $(OBJSPLITFROMALN) $(CFLAGS) $(LIBS)

splittbfast2 : $(OBJSPLITTBFAST2)
	$(CC) -o $@ $(OBJSPLITTBFAST2) $(CFLAGS) $(LIBS)

dummy : $(OBJDUMMY)
	$(CC) -o $@ $(OBJDUMMY) $(CFLAGS) $(LIBS)

setcore : $(OBJSETCORE)
	$(CC) -o $@ $(OBJSETCORE) $(CFLAGS) $(LIBS)

countlen : $(OBJCOUNTLEN)
	$(CC) -o $@ $(OBJCOUNTLEN) $(CFLAGS) $(LIBS)

seq2regtable : $(OBJSEQ2REGTABLE)
	$(CC) -o $@ $(OBJSEQ2REGTABLE) $(CFLAGS) $(LIBS)

regtable2seq : $(OBJREGTABLE2SEQ)
	$(CC) -o $@ $(OBJREGTABLE2SEQ) $(CFLAGS) $(LIBS)

f2cl : $(OBJF2CL)
	$(CC) -o $@ $(OBJF2CL) $(CFLAGS) $(LIBS)

mccaskillwrap : $(OBJMCCASKILLWRAP)
	$(CC) -o $@ $(OBJMCCASKILLWRAP) $(CFLAGS) $(LIBS)

contrafoldwrap : $(OBJCONTRAFOLDWRAP)
	$(CC) -o $@ $(OBJCONTRAFOLDWRAP) $(CFLAGS) $(LIBS)

pairlocalalign : $(OBJPAIRLOCALALIGN)
	$(CC) -o $@ $(OBJPAIRLOCALALIGN) $(CFLAGS) $(LIBS)

pairash : $(OBJPAIRASH)
	$(CC) -o $@ $(OBJPAIRASH) $(CFLAGS) $(LIBS)

rnatest : $(OBJRNATEST)
	$(CC) -o $@ $(OBJRNATEST) $(CFLAGS) $(LIBS)

pair2hat3s : $(OBJPAIR2HAT3S)
	$(CC) -o $@ $(OBJPAIR2HAT3S) $(CFLAGS) $(LIBS)

multi2hat3s : $(OBJMULTI2HAT3S)
	$(CC) -o $@ $(OBJMULTI2HAT3S) $(CFLAGS) $(LIBS)

getlag : $(OBJGETLAG)
	$(CC) -o $@ $(OBJGETLAG) $(CFLAGS) $(LIBS)

tditr : $(OBJTDITR)
	$(CC) -o $@ $(OBJTDITR) $(CFLAGS) $(LIBS)

dvtditr : $(OBJDVTDITR)
	$(CC) -o $@ $(OBJDVTDITR) $(CFLAGS) $(LIBS)

mafft-profile : $(OBJGALN)
	$(CC) -o $@ $(OBJGALN) $(CFLAGS) $(LIBS)

gapfill : $(OBJGAPFILL)
	$(CC) -o $@ $(OBJGAPFILL) $(CFLAGS) $(LIBS)

dndfast4 : $(OBJDNDFAST4)
	$(CC) -o $@ $(OBJDNDFAST4) $(CFLAGS) $(LIBS)

dndfast5 : $(OBJDNDFAST5)
	$(CC) -o $@ $(OBJDNDFAST5) $(CFLAGS) $(LIBS)

dndfast6 : $(OBJDNDFAST6)
	$(CC) -o $@ $(OBJDNDFAST6) $(CFLAGS) $(LIBS)

dndfast7 : $(OBJDNDFAST7)
	$(CC) -o $@ $(OBJDNDFAST7) $(CFLAGS) $(LIBS)

dndblast : $(OBJDNDBLAST)
	$(CC) -o $@ $(OBJDNDBLAST) $(CFLAGS) $(LIBS)

dndfast3 : $(OBJDNDFAST3)
	$(CC) -o $@ $(OBJDNDFAST3) $(CFLAGS) $(LIBS)

triplet : $(OBJTRIPLET)
	$(CC) -o $@ $(OBJTRIPLET) $(CFLAGS) $(LIBS)

triplet3 : $(OBJTRIPLET3)
	$(CC) -o $@ $(OBJTRIPLET3) $(CFLAGS) $(LIBS)

sextet3 : $(OBJSEXTET3)
	$(CC) -o $@ $(OBJSEXTET3) $(CFLAGS) $(LIBS)

sextet4 : $(OBJSEXTET4)
	$(CC) -o $@ $(OBJSEXTET4) $(CFLAGS) $(LIBS)

sextet5 : $(OBJSEXTET5)
	$(CC) -o $@ $(OBJSEXTET5) $(CFLAGS) $(LIBS)

mafft-distance : $(OBJDISTANCE)
	$(CC) -o $@ $(OBJDISTANCE) $(CFLAGS) $(LIBS)

triplet5 : $(OBJTRIPLET5)
	$(CC) -o $@ $(OBJTRIPLET5) $(CFLAGS) $(LIBS)

triplet6 : $(OBJTRIPLET6)
	$(CC) -o $@ $(OBJTRIPLET6) $(CFLAGS) $(LIBS)

octet4 : $(OBJOCTET4)
	$(CC) -o $@ $(OBJOCTET4) $(CFLAGS) $(LIBS)

dndpre : $(OBJDNDPRE)
	$(CC) -o $@ $(OBJDNDPRE) $(CFLAGS) $(LIBS)

dndpre2 : $(OBJDNDPRESCORE)
	$(CC) -o $@ $(OBJDNDPRESCORE) $(CFLAGS) $(LIBS)

score : $(OBJSCORE)
	$(CC) -o $@ $(OBJSCORE) $(CFLAGS) $(LIBS)

genMtx : $(OBJGENMTX)
	$(CC) -o $@ $(OBJGENMTX) $(CFLAGS) $(LIBS)

gapfill.o : gapfill.c $(HEADER)
	$(CC) $(CFLAGS) -c gapfill.c

mltaln9.o : mltaln9.c $(HEADER)
	$(CC) $(CFLAGS) -c mltaln9.c

tddis.o : tddis.c $(HEADER)
	$(CC) $(CFLAGS) -c tddis.c

constants.o : constants.c miyata.h miyata5.h blosum.c DNA.h JTT.c $(HEADER)
	$(CC) $(CFLAGS) -c constants.c

defs.o : defs.c 
	$(CC) $(CFLAGS) -c defs.c

Calignm1.o : Calignm1.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c Calignm1.c

#A+++alignmm.o : SA+++alignmm.c $(HEADER) $(DPHEADER)
#	$(CC) $(CFLAGS) -c SA+++alignmm.c -o A+++alignmm.o

Salignmm.o : Salignmm.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c Salignmm.c 

Halignmm.o : Halignmm.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c Halignmm.c 

Ralignmm.o : Ralignmm.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c Ralignmm.c 

Qalignmm.o : Qalignmm.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c Qalignmm.c 

MSalignmm.o : MSalignmm.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c MSalignmm.c 

partSalignmm.o : partSalignmm.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c partSalignmm.c  

partQalignmm.o : partQalignmm.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c partQalignmm.c  

Lalign11.o : Lalign11.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c Lalign11.c 

genalign11.o : genalign11.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c genalign11.c 

genGalign11.o : genGalign11.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c genGalign11.c 

suboptalign11.o : suboptalign11.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c suboptalign11.c 

Galign11.o : Galign11.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c Galign11.c 

MSalign11.o : MSalign11.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c MSalign11.c 

SAalignmm.o : SAalignmm.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c SAalignmm.c -o SAalignmm.o

Lalignmm.o : Lalignmm.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c Lalignmm.c 

rna.o : rna.c $(HEADER) $(DPHEADER)
	$(CC) $(CFLAGS) -c rna.c 

ord.o : ord.c $(HEADER)
	$(CC) $(CFLAGS) -c ord.c

1itr.o : 1itr.c $(HEADER)
	$(CC) $(CFLAGS) -c 1itr.c

tb.o : tb.c $(HEADER)
	$(CC) $(CFLAGS) -c tb.c

disttbfast.o : disttbfast.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c disttbfast.c

splitfromaln.o : splitfromaln.c $(HEADER) $(FFTHEADER) 
	$(CC) $(CFLAGS) -c splitfromaln.c

splittbfast.o : splittbfast.c $(HEADER) $(FFTHEADER) 
	$(CC) $(CFLAGS) -c splittbfast.c

splittbfast2.o : splittbfast2.c $(HEADER) $(FFTHEADER) 
	$(CC) $(CFLAGS) -c splittbfast2.c

disttbfast_dummy.o : disttbfast_dummy.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c disttbfast_dummy.c

dummy.o : dummy.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c dummy.c

tbfast.o : tbfast.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c tbfast.c

tbfast2.o : tbfast2.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c tbfast2.c

setcore.o : setcore.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c setcore.c

getlag.o : getlag.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c getlag.c

tditr.o : tditr.c $(HEADER)
	 $(CC) $(CFLAGS) -c tditr.c

dvtditr.o : dvtditr.c $(HEADER)
	 $(CC) $(CFLAGS) -c dvtditr.c

tditeration.o : tditeration.c $(HEADER)
	 $(CC) $(CFLAGS) -c tditeration.c

mafft-profile.o : mafft-profile.c $(HEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c mafft-profile.c

dndfast4.o : dndfast4.c $(HEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c dndfast4.c

dndfast5.o : dndfast5.c $(HEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c dndfast5.c

dndfast6.o : dndfast6.c $(HEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c dndfast6.c

dndfast7.o : dndfast7.c $(HEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c dndfast7.c

dndblast.o : dndblast.c $(HEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c dndblast.c

dndfast3.o : dndfast3.c $(HEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c dndfast3.c

dndpre.o : dndpre.c $(HEADER)
	$(CC) $(CFLAGS) -c dndpre.c

countlen.o : countlen.c $(HEADER)
	$(CC) $(CFLAGS) -c countlen.c

seq2regtable.o : seq2regtable.c $(HEADER)
	$(CC) $(CFLAGS) -c seq2regtable.c

regtable2seq.o : regtable2seq.c $(HEADER)
	$(CC) $(CFLAGS) -c regtable2seq.c

f2cl.o : f2cl.c $(HEADER)
	$(CC) $(CFLAGS) -c f2cl.c

mccaskillwrap.o : mccaskillwrap.c $(HEADER)
	$(CC) $(CFLAGS) -c mccaskillwrap.c

contrafoldwrap.o : contrafoldwrap.c $(HEADER)
	$(CC) $(CFLAGS) -c contrafoldwrap.c

pairlocalalign.o : pairlocalalign.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c pairlocalalign.c

pairash.o : pairash.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c pairash.c

rnatest.o : rnatest.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c rnatest.c

multi2hat3s.o : multi2hat3s.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c multi2hat3s.c

pair2hat3s.o : pair2hat3s.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c pair2hat3s.c

dndpre2.o : dndpre2.c $(HEADER)
	$(CC) $(CFLAGS) -c dndpre2.c

io.o : io.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c io.c

nj.o : nj.c $(HEADER)
	$(CC) $(CFLAGS) -c nj.c

treeOperation.o : treeOperation.c $(HEADER)
	$(CC) $(CFLAGS) -c treeOperation.c

sextet5.o : sextet5.c $(HEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c sextet5.c

mafft-distance.o : mafft-distance.c $(HEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c mafft-distance.c

triplet5.o : triplet5.c $(HEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c triplet5.c

triplet6.o : triplet6.c $(HEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c triplet6.c

fft.o : fft.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c fft.c 

fftFunctions.o : fftFunctions.c $(HEADER) $(FFTHEADER)
	$(CC) $(CFLAGS) -c fftFunctions.c

Falign.o : Falign.c $(HEADER) $(FFTHEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c Falign.c

Falign_localhom.o : Falign_localhom.c $(HEADER) $(FFTHEADER) $(MTXHEADER)
	$(CC) $(CFLAGS) -c Falign_localhom.c

mtxutl.o : mtxutl.c 
	$(CC) $(CFLAGS) -c mtxutl.c

score.o : score.c $(HEADER)
	$(CC) $(CFLAGS) -c score.c

clean :
	rm -f *.o *.a *.exe *~ $(PROGS)  $(SCRIPTS)
#	rm -f ../binaries/* ../scripts/*

install : all
	mkdir -p $(LIBDIR)
	chmod 755 $(LIBDIR)
	mkdir -p $(BINDIR)
	chmod 755 $(BINDIR)
	chmod 755 $(SCRIPTS)
	$(INSTALL) $(SCRIPTS)  $(BINDIR)
	chmod 755 $(PROGS)
	$(INSTALL) $(PROGS) $(MANPAGES) $(LIBDIR)
	( cd $(BINDIR); \
rm -f linsi ginsi einsi fftns fftnsi nwns nwnsi xinsi qinsi; \
rm -f mafft-linsi mafft-ginsi mafft-einsi mafft-fftns mafft-fftnsi mafft-nwns mafft-nwnsi mafft-xinsi mafft-qinsi; \
ln -s mafft linsi; ln -s mafft ginsi; ln -s mafft fftns; \
ln -s mafft fftnsi; ln -s mafft nwns; ln -s mafft nwnsi; \
ln -s mafft einsi; \
ln -s mafft mafft-linsi; ln -s mafft mafft-ginsi; ln -s mafft mafft-fftns; \
ln -s mafft mafft-fftnsi; ln -s mafft mafft-nwns; ln -s mafft mafft-nwnsi; \
ln -s mafft mafft-einsi; ln -s mafft mafft-xinsi; ln -s mafft mafft-qinsi;\
rm -f mafft-profile mafft-profile.exe; ln -s ../lib/mafft/mafft-profile .; \
rm -f mafft-distance mafft-distance.exe; ln -s ../lib/mafft/mafft-distance . )

	mkdir -p $(MANDIR)
	chmod 755 $(MANDIR)
	$(INSTALL) $(MANPAGES) $(MANDIR)
