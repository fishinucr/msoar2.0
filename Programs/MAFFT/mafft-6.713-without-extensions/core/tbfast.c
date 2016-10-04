#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0

static int treein;
static int topin;
static int treeout;
static int distout;
static int noalign;

void arguments( int argc, char *argv[] )
{
    int c;

	treein = 0;
	topin = 0;
	rnaprediction = 'm';
	rnakozo = 0;
	nevermemsave = 0;
	inputfile = NULL;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0; // chuui
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
	ppenalty = NOTSPECIFIED;
	ppenalty_ex = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	RNAppenalty = NOTSPECIFIED;
	RNAppenalty_ex = NOTSPECIFIED;
	RNApthr = NOTSPECIFIED;
	TMorJTT = JTT;
	consweight_multi = 1.0;
	consweight_rna = 0.0;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( ( c = *++argv[0] ) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
                    goto nextoption;
				case 'e':
					RNApthr = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'o':
					RNAppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
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
					fprintf( stderr, "blosum %d / kimura 200\n", nblosum );
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
				case 'l':
					fastathreshold = atof( *++argv );
					constraint = 2;
					--argc;
					goto nextoption;
				case 'r':
					consweight_rna = atof( *++argv );
					rnakozo = 1;
					--argc;
					goto nextoption;
				case 'c':
					consweight_multi = atof( *++argv );
					--argc;
					goto nextoption;
				case 'R':
					rnaprediction = 'r';
					break;
				case 's':
					RNAscoremtx = 'r';
					break;
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
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'O':
					fftNoAnchStop = 1;
					break;
#if 0
				case 'e':
					fftscore = 0;
					break;
				case 'r':
					fmodel = -1;
					break;
				case 'R':
					fftRepeatStop = 1;
					break;
				case 's':
					treemethod = 's';
					break;
#endif
				case 'X':
					treemethod = 'X';
					break;
				case 'E':
					treemethod = 'E';
					break;
				case 'q':
					treemethod = 'q';
					break;
#if 0
				case 'a':
					alg = 'a';
					break;
#endif
				case 'Q':
					alg = 'Q';
					break;
				case 'H':
					alg = 'H';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'S':
					alg = 'S';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'N':
					nevermemsave = 1;
					break;
				case 'B':
					break;
				case 'C':
					alg = 'C';
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'G':
					force_fft = 1;
					use_fft = 1;
					break;
				case 'U':
					treein = 1;
					break;
				case 'V':
					topin = 1;
					break;
				case 'u':
					tbrweight = 0;
					weight = 0;
					break;
				case 'v':
					tbrweight = 3;
					break;
				case 'd':
					disp = 1;
					break;
#if 0
				case 'o':
					outgap = 0;
					break;
#endif
/* Modified 01/08/27, default: user tree */
				case 'J':
					tbutree = 0;
					break;
/* modification end. */
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


void treebase( int nlen[M], char **aseq, char **mseq1, char **mseq2, int ***topol, double *effarr, int *alloclen, LocalHom **localhomtable, RNApair ***singlerna, double *effarr_kozo )
{
	int i, l;
	int len1, len2;
	int clus1, clus2;
	float pscore, tscore;
	static char *indication1, *indication2;
	static double *effarr1 = NULL;
	static double *effarr2 = NULL;
	static double *effarr1_kozo = NULL;
	static double *effarr2_kozo = NULL;
	static LocalHom ***localhomshrink = NULL;
	static int *fftlog;
	int m1, m2;
	float dumfl = 0.0;
	int ffttry;
	RNApair ***grouprna1, ***grouprna2;

	if( rnakozo && rnaprediction == 'm' )
	{
		grouprna1 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		grouprna2 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
	}
	else
	{
		grouprna1 = grouprna2 = NULL;
	}

	if( effarr1 == NULL ) 
	{
		fftlog = AllocateIntVec( njob );
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
#if 0
#else
		if( constraint )
		{
			localhomshrink = (LocalHom ***)calloc( njob, sizeof( LocalHom ** ) );
			for( i=0; i<njob; i++)
				localhomshrink[i] = (LocalHom **)calloc( njob, sizeof( LocalHom *) );
		}
#endif
		effarr1_kozo = AllocateDoubleVec( njob ); //tsuneni allocate sareru.
		effarr2_kozo = AllocateDoubleVec( njob ); //tsuneni allocate sareru.
		for( i=0; i<njob; i++ ) effarr1_kozo[i] = 0.0;
		for( i=0; i<njob; i++ ) effarr2_kozo[i] = 0.0;
	}

	for( l=0; l<njob; l++ ) fftlog[l] = 1;

#if 0
	fprintf( stderr, "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif


	if( constraint )
		calcimportance( njob, effarr, aseq, localhomtable );


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

		if( effarr_kozo )
		{
			clus1 = fastconjuction_noname_kozo( topol[l][0], aseq, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
			clus2 = fastconjuction_noname_kozo( topol[l][1], aseq, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
		}
		else
		{
			clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1 );
			clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2 );
		}


		fprintf( trap_g, "\nSTEP-%d\n", l );
		fprintf( trap_g, "group1 = %s\n", indication1 );
		fprintf( trap_g, "group2 = %s\n", indication2 );

#if 1
		fprintf( stderr, "\rSTEP % 5d /%d ", l+1, njob-1 );
#else
		fprintf( stdout, "STEP %d /%d\n", l+1, njob-1 );
		fprintf( stderr, "STEP %d /%d\n", l+1, njob-1 );
		fprintf( stderr, "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, "\n" );
		fprintf( stderr, "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, "\n" );
#endif



//		for( i=0; i<clus1; i++ ) fprintf( stderr, "## STEP%d-eff for mseq1-%d %f\n", l+1, i, effarr1[i] );

		if( constraint )
		{
			fastshrinklocalhom( topol[l][0], topol[l][1], localhomtable, localhomshrink );
//			msfastshrinklocalhom( topol[l][0], topol[l][1], localhomtable, localhomshrink );
//			fprintf( stderr, "localhomshrink =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
//			weightimportance4( clus1, clus2, effarr1, effarr2, localhomshrink );
//			fprintf( stderr, "after weight =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
		}
		if( rnakozo && rnaprediction == 'm' )
		{
			makegrouprna( grouprna1, singlerna, topol[l][0] );
			makegrouprna( grouprna2, singlerna, topol[l][1] );
		}

		free( topol[l][0] );
		free( topol[l][1] );

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

		if( !nevermemsave && ( constraint != 2  && alg != 'M'  && ( len1 > 10000 || len2 > 10000 ) ) )
		{
			fprintf( stderr, "\nlen1=%d, len2=%d, Switching to the memsave mode.\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}
		

//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000 );
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000 ); // v6.708
//		fprintf( stderr, "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (float)len1/fftlog[m1], clus1, (float)len2/fftlog[m2], clus2 );
//		fprintf( stderr, "f=%d, clus1=%d, fftlog[m1]=%d, clus2=%d, fftlog[m2]=%d\n", ffttry, clus1, fftlog[m1], clus2, fftlog[m2] );
		if( constraint == 2 )
		{
			if( alg == 'M' )
			{
				fprintf( stderr, "\n\nMemory saving mode is not supported.\n\n" );
				exit( 1 );
			}
			fprintf( stderr, "c" );
			if( alg == 'A' )
			{
				imp_match_init_strict( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
				if( rnakozo ) imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, NULL, NULL, NULL );
				pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL );
			}
			else if( alg == 'H' )
			{
				imp_match_init_strictH( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
				pscore = H__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL );
			}
			else if( alg == 'Q' )
			{
				imp_match_init_strictQ( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
				if( rnakozo ) imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, NULL, NULL, NULL );
				pscore = Q__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL );
			}
			else if( alg == 'R' )
			{
				imp_match_init_strictR( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
				pscore = R__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumfl, NULL, NULL, NULL, NULL );
			}
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			fprintf( stderr, "f" );
			if( alg == 'M' )
			{
				fprintf( stderr, "m" );
				pscore = Falign_udpari_long( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
				pscore = Falign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, fftlog+m1 );
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
					pscore = MSalignmm( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL );
					break;
				case( 'A' ):
					pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				case( 'Q' ):
					pscore = Q__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				case( 'R' ):
					pscore = R__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				case( 'H' ):
					pscore = H__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					break;
				case ( 'C' ):
					if( outgap && ( ( clus1 == 1 && clus2 != 1 ) || ( clus1 != 1 && clus2 == 1 ) ) )
					{
						pscore = translate_and_Calign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					}
					else
					{
						pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
					}
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}

		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

#if SCOREOUT
		fprintf( stderr, "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
/*
		fprintf( stderr, "after align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "after align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/

//		writePre( njob, name, nlen, aseq, 0 );

		if( disp ) display( aseq, njob );
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
	static float *selfscore;
	int nogaplen;
	static char **name, **seq;
	static char **mseq1, **mseq2;
	static char **bseq;
	static float **iscore, **iscore_kozo;
	static double *eff, *eff_kozo, *eff_kozo_mapped = NULL;
	int i, j, ien, ik, jk;
	static int ***topol, ***topol_kozo;
	static float **len, **len_kozo;
	FILE *prep;
	FILE *infp;
	FILE *orderfp;
	FILE *hat2p;
	
	char c;
	int alloclen;
	LocalHom **localhomtable;
	RNApair ***singlerna;
	float ssi, ssj, bunbo;
	static char *kozoarivec;
	int nkozo;

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

	nkozo = 0;

	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	seq = AllocateCharMtx( njob, nlenmax+1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );

	name = AllocateCharMtx( njob, B+1 );
	nlen = AllocateIntVec( njob );
	selfscore = AllocateFloatVec( njob );

	topol = AllocateIntCub( njob, 2, 0 );
	len = AllocateFloatMtx( njob, 2 );
	iscore = AllocateFloatHalfMtx( njob );
	eff = AllocateDoubleVec( njob );
	kozoarivec = AllocateCharVec( njob );
	if( constraint )
	{
		localhomtable = (LocalHom **)calloc( njob, sizeof( LocalHom *) );
		for( i=0; i<njob; i++)
		{
			localhomtable[i] = (LocalHom *)calloc( njob, sizeof( LocalHom ) );
			for( j=0; j<njob; j++)
			{
				localhomtable[i][j].start1 = -1;
				localhomtable[i][j].end1 = -1;
				localhomtable[i][j].start2 = -1;
				localhomtable[i][j].end2 = -1;
				localhomtable[i][j].overlapaa = -1.0;
				localhomtable[i][j].opt = -1.0;
				localhomtable[i][j].importance = -1.0;
				localhomtable[i][j].next = NULL;
				localhomtable[i][j].korh = 'h';
			}
		}

		fprintf( stderr, "Loading 'hat3' ... " );
		prep = fopen( "hat3", "r" );
		if( prep == NULL ) ErrorExit( "Make hat3." );
		readlocalhomtable( prep, njob, localhomtable, kozoarivec );
		fclose( prep );
		fprintf( stderr, "\ndone.\n" );


		nkozo = 0;
		for( i=0; i<njob; i++ ) 
		{
//			fprintf( stderr, "kozoarivec[%d] = %d\n", i, kozoarivec[i] );
			if( kozoarivec[i] ) nkozo++;
		}
		if( nkozo )
		{
			topol_kozo = AllocateIntCub( nkozo, 2, 0 );
			len_kozo = AllocateFloatMtx( nkozo, 2 );
			iscore_kozo = AllocateFloatHalfMtx( nkozo );
			eff_kozo = AllocateDoubleVec( nkozo );
			eff_kozo_mapped = AllocateDoubleVec( njob );
		}


//		outlocalhom( localhomtable, njob );

#if 0
		fprintf( stderr, "Extending localhom ... " );
		extendlocalhom2( njob, localhomtable );
		fprintf( stderr, "done.\n" );
#endif
	}

#if 0
	readData( infp, name, nlen, seq );
#else
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

//	writePre( njob, name, nlen, seq, 0 );

	if( treein )
	{
#if 0
		if( nkozo )
		{
			fprintf( stderr, "Both structure and user tree have been given. Not yet supported!\n" );
			exit( 1 );
		}
#endif
		fprintf( stderr, "Loading a tree ... " );
		loadtree( njob, topol, len, name, nlen );
		fprintf( stderr, "\ndone.\n\n" );
	}
	else
	{
		if( tbutree == 0 )
		{
			for( i=1; i<njob; i++ ) 
			{
				if( nlen[i] != nlen[0] ) 
				{
					fprintf( stderr, "Input pre-aligned seqences or make hat2.\n" );
					exit( 1 );
				}
			}
	
			fprintf( stderr, "Making a distance matrix .. \n" );
			ien = njob-1;
#if 0 //060424
			for( i=0; i<njob; i++ ) 
			{
				fprintf( stderr, "\r% 5d / %d", i, ien );
				for( j=i; j<njob; j++ ) 
				{
					iscore[i][j-i] = score_calcp( seq[i], seq[j], nlen[0] );
					fprintf( stderr, "i=%d,j=%d### %f\n", i, j, iscore[i][j-i] );
				}
			}
			fprintf( stderr, "\ndone.\n\n" );
   	    	for( i=0; i<ien; i++ )
	   	    {
   		        for( j=i+1; j<njob; j++ ) 
   		        {
   		            iscore[i][j-i] = ( 1.0 - iscore[i][j-i] / MIN( iscore[i][0], iscore[j][0] ) ) * 3;
					fprintf( stderr, "i=%d,j=%d### %f\n", i, j, iscore[i][j-i] );
   		        }
   		    }
#else
			for( i=0; i<njob; i++ ) 
			{
				selfscore[i] = naivepairscore11( seq[i], seq[i], penalty );
			}
			for( i=0; i<ien; i++ ) 
			{
				if( i % 10 == 0 ) fprintf( stderr, "\r% 5d / %d", i, ien );
				ssi = selfscore[i];
				for( j=i+1; j<njob; j++ ) 
				{
					ssj = selfscore[j];
					bunbo = MIN( ssi, ssj );
					if( bunbo == 0.0 )
						iscore[i][j-i] = 1.0;
					else
						iscore[i][j-i] = 1.0 - naivepairscore11( seq[i], seq[j], penalty ) / MIN( selfscore[i], selfscore[j] );

#if 0
					fprintf( stderr, "### ssj = %f\n", ssj );
					fprintf( stderr, "### selfscore[i] = %f\n", selfscore[i] );
					fprintf( stderr, "### selfscore[j] = %f\n", selfscore[j] );
					fprintf( stderr, "### rawscore = %f\n", naivepairscore11( seq[i], seq[j], penalty ) );
#endif
				}
			}
			fprintf( stderr, "\ndone.\n\n" );
#endif
	
		}
		else
		{
			fprintf( stderr, "Loading 'hat2' ... " );
			prep = fopen( "hat2", "r" );
			if( prep == NULL ) ErrorExit( "Make hat2." );
			readhat2_floathalf_pointer( prep, njob, name, iscore );
			fclose( prep );
			fprintf( stderr, "done.\n" );
		}
#if 1
		if( distout )
		{
			hat2p = fopen( "hat2", "w" );
			WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, iscore );
			fclose( hat2p );
		}
#endif
		if( nkozo )
		{
			ien = njob-1;
			ik = 0;
			for( i=0; i<ien; i++ )
			{
				jk = ik+1;
				for( j=i+1; j<njob; j++ ) 
				{
					if( kozoarivec[i] && kozoarivec[j] )
					{
						iscore_kozo[ik][jk-ik] = iscore[i][j-i];
					}
					if( kozoarivec[j] ) jk++;
				}
				if( kozoarivec[i] ) ik++;
			}
		}

		fprintf( stderr, "Constructing a UPGMA tree ... " );
		if( topin )
		{
			fprintf( stderr, "Loading a topology ... " );
			loadtop( njob, iscore, topol, len );
			fprintf( stderr, "\ndone.\n\n" );
		}
		else if( treeout )
		{
			fixed_musclesupg_float_realloc_nobk_halfmtx_treeout( njob, iscore, topol, len, name, nlen );
		}
		else
		{
			fixed_musclesupg_float_realloc_nobk_halfmtx( njob, iscore, topol, len );
		}
//		else 
//			ErrorExit( "Incorrect tree\n" );

		if( nkozo )
		{
//			for( i=0; i<nkozo-1; i++ )
//				for( j=i+1; j<nkozo; j++ )
//					fprintf( stderr, "iscore_kozo[%d][%d] =~ %f\n", i, j, iscore_kozo[i][j-i] );
			fixed_musclesupg_float_realloc_nobk_halfmtx( nkozo, iscore_kozo, topol_kozo, len_kozo );
		}
		fprintf( stderr, "\ndone.\n\n" );
	}


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

	if( treeout && noalign ) 
	{
		writeData_pointer( prep_g, njob, name, nlen, seq );
		fprintf( stderr, "\n" ); 
		SHOWVERSION;
		return( 0 );
	}

//	countnode( njob, topol, node0 );
	if( tbrweight )
	{
		weight = 3; 
#if 0
		utree = 0; counteff( njob, topol, len, eff ); utree = 1;
#else
		counteff_simple_float( njob, topol, len, eff );

		if( nkozo )
		{
//			counteff_simple_float( nkozo, topol_kozo, len_kozo, eff_kozo ); // single weight nanode iranai
			for( i=0,j=0; i<njob; i++ )
			{
				if( kozoarivec[i] )
				{
//					eff_kozo_mapped[i] = eff_kozo[j]; //
					eff_kozo_mapped[i] = eff[i];      // single weight
					j++;
				}
				else
					eff_kozo_mapped[i] = 0.0;
//				fprintf( stderr, "eff_kozo_mapped[%d] = %f\n", i, eff_kozo_mapped[i] );
//				fprintf( stderr, "            eff[%d] = %f\n", i, eff[i] );
			}
		}


#endif
	}
	else
	{
		for( i=0; i<njob; i++ ) eff[i] = 1.0;
		if( nkozo ) 
		{
			for( i=0; i<njob; i++ ) 
			{
				if( kozoarivec[i] ) 
					eff_kozo_mapped[i] = 1.0;
				else
					eff_kozo_mapped[i] = 0.0;
			}
		}
	}

	FreeFloatHalfMtx( iscore, njob );
	FreeFloatMtx( len );

	bseq = AllocateCharMtx( njob, nlenmax*2+1 );
	alloclen = nlenmax*2;

	for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );

	if( rnakozo && rnaprediction == 'm' )
	{
		singlerna = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		prep = fopen( "hat4", "r" );
		if( prep == NULL ) ErrorExit( "Make hat4 using mccaskill." );
		fprintf( stderr, "Loading 'hat4' ... " );
		for( i=0; i<njob; i++ )
		{
			nogaplen = strlen( bseq[i] );
			singlerna[i] = (RNApair **)calloc( nogaplen, sizeof( RNApair * ) );
			for( j=0; j<nogaplen; j++ )
			{
				singlerna[i][j] = (RNApair *)calloc( 1, sizeof( RNApair ) );
				singlerna[i][j][0].bestpos = -1;
				singlerna[i][j][0].bestscore = -1.0;
			}
			readmccaskill( prep, singlerna[i], nogaplen );
		}
		fclose( prep );
		fprintf( stderr, "\ndone.\n" );
	}
	else
		singlerna = NULL;


	fprintf( stderr, "Progressive alignment ... \n" );
	treebase( nlen, bseq, mseq1, mseq2, topol, eff, &alloclen, localhomtable, singlerna, eff_kozo_mapped );
	fprintf( stderr, "\ndone.\n" );

#if 0
	if( constraint )
	{
		LocalHom *tmppt1, *tmppt2;
		for( i=0; i<njob; i++)
		{
			for( j=0; j<njob; j++)
			{
				tmppt1 = localhomtable[i]+j;
				while( tmppt2 = tmppt1->next )
				{
					free( (void *)tmppt1 );
					tmppt1 = tmppt2;
				}
				free( (void *)tmppt1 );
			}
			free( (void *)(localhomtable[i]+j) );
		}
		free( (void *)localhomtable );
	}
#endif

	fprintf( trap_g, "done.\n" );
	fclose( trap_g );

	writeData_pointer( prep_g, njob, name, nlen, bseq );
#if 0
	writeData( stdout, njob, name, nlen, bseq );
	writePre( njob, name, nlen, bseq, !contin );
	writeData_pointer( prep_g, njob, name, nlen, aseq );
#endif
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif

	if( constraint ) FreeLocalHomTable( localhomtable, njob );

	SHOWVERSION;
	return( 0 );
}
