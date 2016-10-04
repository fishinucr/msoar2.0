#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0

#define END_OF_VEC -1

static int treein;
static int topin;
static int treeout;
static int noalign;
static int distout;
static float lenfaca, lenfacb, lenfacc, lenfacd;
#if 0
#define PLENFACA 0.0123
#define PLENFACB 10252
#define PLENFACC 10822
#define PLENFACD 0.5
#define DLENFACA 0.01
#define DLENFACB 2445
#define DLENFACC 2412
#define DLENFACD 0.1
#else
#define PLENFACA 0.01
#define PLENFACB 10000
#define PLENFACC 10000
#define PLENFACD 0.1
#define DLENFACA 0.01
#define DLENFACB 2500
#define DLENFACC 2500
#define DLENFACD 0.1
#endif


void arguments( int argc, char *argv[] )
{
    int c;

	topin = 0;
	treein = 0;
	treeout = 0;
	distout = 0;
	noalign = 0;
	nevermemsave = 0;
	inputfile = NULL;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0;
	force_fft = 0;
	fftscore = 1;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
    weight = 3;
    utree = 1;
	tbutree = 1;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 1;
    alg = 'A';
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
	treemethod = 'X';
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	dorp = NOTSPECIFIED;
	ppenalty = -1530;
	ppenalty_ex = NOTSPECIFIED;
	poffset = -123;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	TMorJTT = JTT;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = atoi( *++argv );
					fprintf( stderr, "kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = atoi( *++argv );
					scoremtx = 1;
//					fprintf( stderr, "blosum %d / kimura 200 \n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = atoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					fprintf( stderr, "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = atoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					fprintf( stderr, "tm %d\n", pamN );
					--argc;
					goto nextoption;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'y':
					distout = 1;
					break;
				case 't':
					treeout = 1;
					break;
				case 'T':
					noalign = 1;
					break;
				case 'r':
					fmodel = -1;
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'e':
					fftscore = 0;
					break;
				case 'O':
					fftNoAnchStop = 1;
					break;
#if 0
				case 'R':
					fftRepeatStop = 1;
					break;
#endif
				case 's':
					treemethod = 's';
					break;
				case 'X':
					treemethod = 'X'; // mix
					break;
				case 'E':
					treemethod = 'E'; // upg (average)
					break;
				case 'q':
					treemethod = 'q'; // minimum
					break;
#if 0
				case 'a':
					alg = 'a';
					break;
#endif
				case 'R':
					alg = 'R';
					break;
				case 'Q':
					alg = 'Q';
					break;
				case 'H':
					alg = 'H';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'N':
					nevermemsave = 1;
					break;
				case 'M':
					alg = 'M';
					break;
				case 'B':
					break;
				case 'S':
					alg = 'S';
					break;
				case 'C':
					alg = 'C';
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'G':
					use_fft = 1;
					force_fft = 1;
					break;
				case 'V':
					topin = 1;
					break;
				case 'U':
					treein = 1;
					break;
				case 'u':
					weight = 0;
					tbrweight = 0;
					break;
				case 'v':
					tbrweight = 3;
					break;
				case 'd':
					disp = 1;
					break;
				case 'o':
					outgap = 0;
					break;
				case 'J':
					tbutree = 0;
					break;
				case 'z':
					fftThreshold = atoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = atoi( *++argv );
					--argc;
					goto nextoption;
				case 'Z':
					checkC = 1;
					break;
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}
    if( argc == 1 )
    {
        cut = atof( (*argv) );
        argc--;
    }
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : o, m or u\n" );
		exit( 1 );
	}
	if( alg == 'C' && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : C, o\n" );
		exit( 1 );
	}
}

static int maxl;
static int tsize;

void seq_grp_nuc( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 4 )
			*grp++ = tmp;
		else
			fprintf( stderr, "WARNING : Unknown character %c\r", *(seq-1) );
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
		fprintf( stderr, "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void seq_grp( int *grp, char *seq )
{
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
			fprintf( stderr, "WARNING : Unknown character %c\r", *(seq-1) );
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
		fprintf( stderr, "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void makecompositiontable_p( short *table, int *pointt )
{
	int point;

	while( ( point = *pointt++ ) != END_OF_VEC )
		table[point]++;
}

int commonsextet_p( short *table, int *pointt )
{
	int value = 0;
	short tmp;
	int point;
	static short *memo = NULL;
	static int *ct = NULL;
	static int *cp;

	if( *pointt == -1 )
		return( 0 );

	if( !memo )
	{
		memo = (short *)calloc( tsize, sizeof( short ) );
		if( !memo ) ErrorExit( "Cannot allocate memo\n" );
		ct = (int *)calloc( MIN( maxl, tsize )+1, sizeof( int ) );
		if( !ct ) ErrorExit( "Cannot allocate memo\n" );
	}

	cp = ct;
	while( ( point = *pointt++ ) != END_OF_VEC )
	{
		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
		if( tmp == 0 ) *cp++ = point;
	}
	*cp = END_OF_VEC;
	
	cp =  ct;
	while( *cp != END_OF_VEC )
		memo[*cp++] = 0;

	return( value );
}

void makepointtable_nuc( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 1024;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  7776;
	point += *n++ *  1296;
	point += *n++ *   216;
	point += *n++ *    36;
	point += *n++ *     6;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 7776;
		point *= 6;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void treebase( int *nlen, char **aseq, char **mseq1, char **mseq2, int ***topol, double *effarr, int *alloclen )
{
	int l, len1, len2;
	int clus1, clus2;
	float pscore, tscore;
	static char *indication1, *indication2;
	static double *effarr1 = NULL;
	static double *effarr2 = NULL;
	static int *fftlog; // fixed at 2006/07/26
	float dumfl = 0.0;
	int ffttry;
	int m1, m2;
#if 0
	int i, j;
#endif


	if( effarr1 == NULL ) 
	{
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		fftlog = AllocateIntVec( njob );
	}
	for( l=0; l<njob; l++ ) fftlog[l] = 1;

#if 0
	fprintf( stderr, "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	tscore = 0.0;
	for( l=0; l<njob-1; l++ )
	{
		m1 = topol[l][0][0];
		m2 = topol[l][1][0];
		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen < len1 + len2 )
		{
			fprintf( stderr, "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
			fprintf( stderr, "done. *alloclen = %d\n", *alloclen );
		}

#if 1 // CHUUI@@@@
		clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1 );
		clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2 );
#else
		clus1 = fastconjuction_noweight( topol[l][0], aseq, mseq1, effarr1,  indication1 );
		clus2 = fastconjuction_noweight( topol[l][1], aseq, mseq2, effarr2,  indication2 );
#endif


#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            fprintf( stderr, "i = %d / %d\n", i, clus1 ); 
            fprintf( stderr, "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            fprintf( stderr, "j = %d / %d\n", j, clus2 ); 
            fprintf( stderr, "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
#endif

		free( topol[l][0] );
		free( topol[l][1] );

#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            fprintf( stderr, "i = %d / %d\n", i, clus1 ); 
            fprintf( stderr, "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            fprintf( stderr, "j = %d / %d\n", j, clus2 ); 
            fprintf( stderr, "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
#endif


//		fprintf( trap_g, "\nSTEP-%d\n", l );
//		fprintf( trap_g, "group1 = %s\n", indication1 );
//		fprintf( trap_g, "group2 = %s\n", indication2 );

//		fprintf( stderr, "\rSTEP % 5d / %d %d-%d", l+1, njob-1, clus1, clus2 );
		fprintf( stderr, "\rSTEP % 5d / %d ", l+1, njob-1 );

#if 0
		fprintf( stderr, "STEP %d /%d\n", l+1, njob-1 );
		fprintf( stderr, "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, "\n" );
		fprintf( stderr, "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, "\n" );
#endif

/*
		fprintf( stderr, "before align all\n" );
		display( aseq, njob );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/


		if( !nevermemsave && ( alg != 'M' && ( len1 > 10000 || len2 > 10000  ) ) )
		{
			fprintf( stderr, "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000); // v6.708
//		fprintf( stderr, "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (float)len1/fftlog[m1], clus1, (float)len2/fftlog[m2], clus2 );

		if( force_fft || ( use_fft && ffttry ) )
		{
			fprintf( stderr, "f" );
			if( alg == 'M' )
			{
				fprintf( stderr, "m" );
//				fprintf( stderr, "%d-%d", clus1, clus2 );
				pscore = Falign_udpari_long( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
//				fprintf( stderr, "%d-%d", clus1, clus2 );
				pscore = Falign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, fftlog+m1 );
			}
		}
		else
		{
			fprintf( stderr, "d" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					fprintf( stderr, "m" );
//					fprintf( stderr, "%d-%d", clus1, clus2 );
					pscore = MSalignmm( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL );
					break;
				case( 'Q' ):
					if( clus1 == 1 && clus2 == 1 && 0 )
					{
//							fprintf( stderr, "%d-%d", clus1, clus2 );
							pscore = G__align11( mseq1, mseq2, *alloclen );
					}
					else
					{
						pscore = Q__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					}
					break;
				case( 'R' ):
					pscore = R__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				case( 'H' ):
					pscore = H__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						fprintf( stderr, "%d-%d", clus1, clus2 );
						pscore = G__align11( mseq1, mseq2, *alloclen );
					}
					else
					{
//						fprintf( stderr, "%d-%d", clus1, clus2 );
						pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					}
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
#if SCOREOUT
		fprintf( stderr, "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

//		writePre( njob, name, nlen, aseq, 0 );

		if( disp ) display( aseq, njob );
//		fprintf( stderr, "\n" );
	}
#if SCOREOUT
	fprintf( stderr, "totalscore = %10.2f\n\n", tscore );
#endif
}

static void WriteOptions( FILE *fp )
{

	if( dorp == 'd' ) fprintf( fp, "DNA\n" );
	else
	{
		if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
		else if( scoremtx ==  1 ) fprintf( fp, "BLOSUM %d\n", nblosum );
		else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	}
    fprintf( stderr, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
    if( use_fft ) fprintf( fp, "FFT on\n" );

	fprintf( fp, "tree-base method\n" );
	if( tbrweight == 0 ) fprintf( fp, "unweighted\n" );
	else if( tbrweight == 3 ) fprintf( fp, "clustalw-like weighting\n" );
	if( tbitr || tbweight ) 
	{
		fprintf( fp, "iterate at each step\n" );
		if( tbitr && tbrweight == 0 ) fprintf( fp, "  unweighted\n" ); 
		if( tbitr && tbrweight == 3 ) fprintf( fp, "  reversely weighted\n" ); 
		if( tbweight ) fprintf( fp, "  weighted\n" ); 
		fprintf( fp, "\n" );
	}

   	 fprintf( fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );

	if( alg == 'a' )
		fprintf( fp, "Algorithm A\n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Algorithm A+\n" );
	else if( alg == 'S' ) 
		fprintf( fp, "Apgorithm S\n" );
	else if( alg == 'C' ) 
		fprintf( fp, "Apgorithm A+/C\n" );
	else
		fprintf( fp, "Unknown algorithm\n" );

	if( treemethod == 'X' )
		fprintf( fp, "Tree = UPGMA (mix).\n" );
	else if( treemethod == 'E' )
		fprintf( fp, "Tree = UPGMA (average).\n" );
	else if( treemethod == 'q' )
		fprintf( fp, "Tree = Minimum linkage.\n" );
	else
		fprintf( fp, "Unknown tree.\n" );

    if( use_fft )
    {
        fprintf( fp, "FFT on\n" );
        if( dorp == 'd' )
            fprintf( fp, "Basis : 4 nucleotides\n" );
        else
        {
            if( fftscore )
                fprintf( fp, "Basis : Polarity and Volume\n" );
            else
                fprintf( fp, "Basis : 20 amino acids\n" );
        }
        fprintf( fp, "Threshold   of anchors = %d%%\n", fftThreshold );
        fprintf( fp, "window size of anchors = %dsites\n", fftWinSize );
    }
	else
        fprintf( fp, "FFT off\n" );
	fflush( fp );
}
	 

int main( int argc, char *argv[] )
{
	static int  *nlen;	
	static int  *nogaplen;	
	static char **name, **seq;
	static char **mseq1, **mseq2;
	static char **bseq;
	static double *eff;
	int i, j;
	static int ***topol;
	static float **len;
	FILE *infp;
	char c;
	int alloclen;
	float longer, shorter;
	float lenfac;
	float bunbo;

	FILE *orderfp, *hat2p;
	int *grpseq;
	char *tmpseq;
	int  **pointt;
	float **mtx = NULL; // by D. Mathog
	static short *table1;
	char b[B];
	int ien;



	arguments( argc, argv );

	if( inputfile )
	{
		infp = fopen( inputfile, "r" );
		if( !infp )
		{
			fprintf( stderr, "Cannot open %s\n", inputfile );
			exit( 1 );
		}
	}
	else
		infp = stdin;

	getnumlen( infp );
	rewind( infp );

	if( njob > 20000 )
	{
		fprintf( stderr, "The number of sequences must be < %d\n", 20000 );
		fprintf( stderr, "Please try the --parttree option for such large data.\n" );
		exit( 1 );
	}

	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	seq = AllocateCharMtx( njob, nlenmax*1+1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );

	topol = AllocateIntCub( njob, 2, 0 );
	len = AllocateFloatMtx( njob, 2 );
	eff = AllocateDoubleVec( njob );


#if 0
	Read( name, nlen, seq );
	readData( infp, name, nlen, seq );
#else
    name = AllocateCharMtx( njob, B+1 );
    nlen = AllocateIntVec( njob ); 
    nogaplen = AllocateIntVec( njob ); 
	readData_pointer( infp, name, nlen, seq );
#endif
	fclose( infp );

	constants( njob, seq );


#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();

	WriteOptions( trap_g );

	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}

	fprintf( stderr, "\n" );

	if( !treein )
	{
		fprintf( stderr, "\n\nMaking a distance matrix ..\n" );

	    tmpseq = AllocateCharVec( nlenmax+1 );
		grpseq = AllocateIntVec( nlenmax+1 );
		pointt = AllocateIntMtx( njob, nlenmax+1 );
    	mtx = AllocateFloatHalfMtx( njob ); 
		if( dorp == 'd' ) tsize = (int)pow( 4, 6 );
		else              tsize = (int)pow( 6, 6 );

		if( dorp == 'd' )
		{
			lenfaca = DLENFACA;
			lenfacb = DLENFACB;
			lenfacc = DLENFACC;
			lenfacd = DLENFACD;
		}
		else    
		{
			lenfaca = PLENFACA;
			lenfacb = PLENFACB;
			lenfacc = PLENFACC;
			lenfacd = PLENFACD;
		}

		maxl = 0;
		for( i=0; i<njob; i++ ) 
		{
			gappick0( tmpseq, seq[i] );
			nogaplen[i] = strlen( tmpseq );
			if( nogaplen[i] < 6 )
			{
//				fprintf( stderr, "Seq %d, too short, %d characters\n", i+1, nogaplen[i] );
//				fprintf( stderr, "Please use mafft-ginsi, mafft-linsi or mafft-ginsi\n\n\n" );
//				exit( 1 );
			}
			if( nogaplen[i] > maxl ) maxl = nogaplen[i];
			if( dorp == 'd' ) /* nuc */
			{
				seq_grp_nuc( grpseq, tmpseq );
				makepointtable_nuc( pointt[i], grpseq );
			}
			else                 /* amino */
			{
				seq_grp( grpseq, tmpseq );
				makepointtable( pointt[i], grpseq );
			}
		}
		for( i=0; i<njob; i++ )
		{
			table1 = (short *)calloc( tsize, sizeof( short ) );
			if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
			if( i % 10 == 0 )
			{
				fprintf( stderr, "\r% 5d / %d", i+1, njob );
			}
			makecompositiontable_p( table1, pointt[i] );
	
			for( j=i; j<njob; j++ ) 
			{
				mtx[i][j-i] = (float)commonsextet_p( table1, pointt[j] );
			} 
			free( table1 );
		}
		fprintf( stderr, "\ndone.\n\n" );
		ien = njob-1;

		for( i=0; i<ien; i++ )
		{
			for( j=i+1; j<njob; j++ ) 
			{
				if( nogaplen[i] > nogaplen[j] )
				{
					longer=(float)nogaplen[i];
					shorter=(float)nogaplen[j];
				}
				else
				{
					longer=(float)nogaplen[j];
					shorter=(float)nogaplen[i];
				}
//				lenfac = 3.0 / ( LENFACA + LENFACB / ( longer + LENFACC ) + shorter / longer * LENFACD );
				lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
//				lenfac = 1.0;
//				fprintf( stderr, "lenfac = %f (%.0f,%.0f)\n", lenfac, longer, shorter );
				bunbo = MIN( mtx[i][0], mtx[j][0] );
				if( bunbo == 0.0 )
					mtx[i][j-i] = 1.0;
				else
					mtx[i][j-i] = ( 1.0 - mtx[i][j-i] / bunbo ) * lenfac;
//				fprintf( stdout, "##### mtx = %f, mtx[i][0]=%f, mtx[j][0]=%f, bunbo=%f\n", mtx[i][j-i], mtx[i][0], mtx[j][0], bunbo );
			}
		}
		if( disopt )
		{
			for( i=0; i<njob; i++ ) 
			{
				sprintf( b, "=lgth = %04d", nogaplen[i] );
				strins( b, name[i] );
			}
		}
		free( grpseq );
		free( tmpseq );
		FreeIntMtx( pointt );

#if 1 // writehat2 wo kakinaosu
		if( distout )
		{
			hat2p = fopen( "hat2", "w" );
			WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, mtx );
			fclose( hat2p );
		}
#endif

	}
	else {
#if 0 // readhat2 wo kakinaosu
		fprintf( stderr, "Loading 'hat2' ... " );
		prep = fopen( "hat2", "r" );
		if( prep == NULL ) ErrorExit( "Make hat2." );
		readhat2_float( prep, njob, name, mtx ); // name chuui
		fclose( prep );
		fprintf( stderr, "done.\n" );
#endif
	}

	if( treein )
	{
		fprintf( stderr, "Loading a tree ... " );
		loadtree( njob, topol, len, name, nogaplen );
	}
	else if( topin )
	{
		fprintf( stderr, "Loading a topology ... " );
		loadtop( njob, mtx, topol, len );
		FreeFloatHalfMtx( mtx, njob );
	}
	else if( treeout )
	{
		fprintf( stderr, "Constructing a UPGMA tree ... " );

		fixed_musclesupg_float_realloc_nobk_halfmtx_treeout( njob, mtx, topol, len, name, nogaplen );
//		veryfastsupg_float_realloc_nobk_halfmtx_treeout( njob, mtx, topol, len, name, nogaplen );


		FreeFloatHalfMtx( mtx, njob );
	}
	else
	{
		fprintf( stderr, "Constructing a UPGMA tree ... " );
		fixed_musclesupg_float_realloc_nobk_halfmtx( njob, mtx, topol, len );
		FreeFloatHalfMtx( mtx, njob );
	}
//	else 
//		ErrorExit( "Incorrect tree\n" );
	fprintf( stderr, "\ndone.\n\n" );

	orderfp = fopen( "order", "w" );
	if( !orderfp )
	{
		fprintf( stderr, "Cannot open 'order'\n" );
		exit( 1 );
	}
	for( i=0; (j=topol[njob-2][0][i])!=-1; i++ )
	{
		fprintf( orderfp, "%d\n", j );
	}
	for( i=0; (j=topol[njob-2][1][i])!=-1; i++ )
	{
		fprintf( orderfp, "%d\n", j );
	}
	fclose( orderfp );

	if( ( treeout || distout )  && noalign ) 
	{
		writeData_pointer( stdout, njob, name, nlen, seq );
		fprintf( stderr, "\n" );
		SHOWVERSION;
		return( 0 );
	}
	

	if( tbrweight )
	{
		weight = 3; 
#if 0
		utree = 0; counteff( njob, topol, len, eff ); utree = 1;
#else
		counteff_simple_float( njob, topol, len, eff );
#endif
	}
	else
	{
		for( i=0; i<njob; i++ ) eff[i] = 1.0;
	}

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stdout, "eff[%d] = %20.16f\n", i, eff[i] );
	exit( 1 );
#endif


	FreeFloatMtx( len );

	bseq = AllocateCharMtx( njob, nlenmax*2+1 );
	alloclen = nlenmax*2;
	for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );

	fprintf( stderr, "Progressive alignment ... \n" );
	treebase( nlen, bseq, mseq1, mseq2, topol, eff, &alloclen );
	fprintf( stderr, "\ndone.\n\n" );
#if DEBUG
	fprintf( stderr, "closing trap_g\n" );
#endif
	fclose( trap_g );

//	writePre( njob, name, nlen, aseq, !contin );
#if DEBUG
	fprintf( stderr, "writing alignment to stdout\n" );
#endif
	writeData_pointer( stdout, njob, name, nlen, bseq );
#if 0
	writeData( stdout, njob, name, nlen, bseq );
#endif
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif
	SHOWVERSION;
	return( 0 );
}
