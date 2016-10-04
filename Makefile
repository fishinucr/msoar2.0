BLAST_VERSION  = 2.2.13
PHYLIP_VERSION = 3.67

BLAST_URL     = ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/$(BLAST_VERSION)/blast-$(BLAST_VERSION)-x64-linux.tar.gz
BLAST_DIR     = ./Programs/blast-$(BLAST_VERSION)
BLAST         = $(BLAST_DIR)/bin/blastall \
				$(BLAST_DIR)/bin/formatdb

MAFFT_DIR = ./Programs/MAFFT/mafft-6.713-without-extensions
MAFFT     = $(MAFFT_DIR)/binaries/mafft

PAL2NAL = ./Programs/PAL2NAL/pal2nal.v12/pal2nal.pl

MCL_DIR = ./Programs/mcl-06-058
MCL     = $(MCL_DIR)/bin/mcxdeblast $(MCL_DIR)/bin/mcl

PHYLIP_DIR    = ./Programs/phylip-$(PHYLIP_VERSION)
PHYLIP_SRC    = $(PHYLIP_DIR)/src
PHYLIP_EXE    = $(PHYLIP_DIR)/exe
DNADIST       = $(PHYLIP_EXE)/dnadist

MSOAR_DIR = Programs/MSOAR

all: $(MAFFT) $(PAL2NAL) $(DNADIST)
	[[ ! -z `type -P blastall` || -f $(BLAST_DIR)/bin/blastall ]] || ${MAKE} blast
	${MAKE} mcl
	cd Programs && ${MAKE}
	cd $(MSOAR_DIR) && ${MAKE}

$(PAL2NAL):
	echo "PAL2NAL not found" > /dev/stderr
	exit 1

$(DNADIST):
	echo "Building PHYLIP" > /dev/stderr
	cd $(PHYLIP_SRC) && ${MAKE} clean && ${MAKE} install

$(MAFFT):
	echo "Building MAFFT" > /dev/stderr
	cd $(MAFFT_DIR)/core && ${MAKE} clean && ${MAKE}

.PHONY: blast
blast:
	echo "Building blast" > /dev/stderr
	rm -rf $(BLAST_DIR)
	echo "Retrieving BLAST binary for x64-linux" > /dev/stderr
	wget -O Programs/blast.tar.gz $(BLAST_URL) && \
		cd Programs                            && \
		tar -xzf blast.tar.gz                  && \
		rm blast.tar.gz

.PHONY: mcl
mcl:
	echo "Building MCL" > /dev/stderr
	cd $(MCL_DIR) && ${MAKE}

.PHONY: runclean
runclean:
	rm -f *.codon outfile *.map *.blastp *.cluster *.pep
	rm -rf Families

.PHONY: distclean
clean:
	rm -rf $(BLAST_DIR)
	cd Programs          && ${MAKE} clean
	cd $(MAFFT_DIR)/core && ${MAKE} clean
	cd $(MSOAR_DIR)      && ${MAKE} clean
	cd $(MCL_DIR)        && ${MAKE} clean
	cd $(PHYLIP_SRC)     && ${MAKE} clean
