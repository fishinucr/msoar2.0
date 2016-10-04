#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0


void arguments( int argc, char *argv[] )
{
    int c;

	thrinter = 1.0;
	inputfile = NULL;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0;
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
	treemethod = 'x';
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
	TMorJTT = JTT;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( c = *++argv[0] )
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
					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
					fprintf( stderr, "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = atoi( *++argv );
					fprintf( stderr, "kimuraR = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = atoi( *++argv );
					scoremtx = 1;
					fprintf( stderr, "blosum %d\n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = atoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					fprintf( stderr, "jtt %d\n", pamN );
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
					fprintf( stderr, "weighti = %f\n", fastathreshold );
					--argc;
					goto nextoption;
				case 'c':
					thrinter = atof( *++argv );
					fprintf( stderr, "thrinter %f\n", thrinter );
					--argc;
					goto nextoption;
#if 0
				case 'm':
					fmodel = 1;
					break;
#endif
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
				case 'R':
					fftRepeatStop = 1;
					break;
				case 'Q':
					calledByXced = 1;
					break;
				case 's':
					treemethod = 's';
					break;
				case 'x':
					treemethod = 'x';
					break;
				case 'p':
					treemethod = 'p';
					break;
				case 'a':
					alg = 'a';
					break;
				case 'A':
					alg = 'A';
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
				case 'v':
					tbrweight = 3;
					break;
				case 'd':
					disp = 1;
					break;
				case 'o':
					outgap = 0;
					break;
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


void treebase( char name[M][B], int nlen[M], char **seq, char **aseq, char **mseq1, char **mseq2, int ***topol, double *effarr, int alloclen, LocalHom **localhomtable )
{
	int i, j, l;
	int clus1, clus2;
	int s1, s2, r1, r2;
	float pscore, tscore;
	static char *indication1, *indication2;
	FILE *trap;
	static char name1[M][B], name2[M][B];
	static double *effarr1 = NULL;
	static double *effarr2 = NULL;
	static LocalHom ***localhomshrink = NULL;
	float dumfl = 0.0;
	double total;


	if( effarr1 == NULL ) 
	{
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
	}

#if 0
	fprintf( stderr, "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );

	if( constraint )
		calcimportance( njob, effarr, aseq, localhomtable );


	writePre( njob, name, nlen, aseq, 0 );

	tscore = 0.0;
	for( l=0; l<njob-1; l++ )
	{

		clus1 = fastconjuction( topol[l][0], aseq, mseq1, effarr1, effarr, name, name1, indication1 );
		clus2 = fastconjuction( topol[l][1], aseq, mseq2, effarr2, effarr, name, name2, indication2 );


		fprintf( trap_g, "\nSTEP-%d\n", l );
		fprintf( trap_g, "group1 = %s\n", indication1 );
		fprintf( trap_g, "group2 = %s\n", indication2 );

#if 1
		fprintf( stderr, "STEP %d /%d\r", l+1, njob-1 );
#else
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
//			fprintf( stderr, "localhomshrink =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
//			weightimportance4( clus1, clus2, effarr1, effarr2, localhomshrink );
//			fprintf( stderr, "after weight =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
		}

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

		if( constraint == 2 )
		{
			imp_match_init_strict( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
			pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, localhomshrink, &dumfl );
		}
		else if( use_fft )
		{
			pscore = Falign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );  // _noudp ??
		}
		else
		{
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
					break;
				case( 'A' ):
					pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, NULL, &dumfl );
					break;
				case ( 'C' ):
					if( outgap && ( clus1 == 1 && clus2 != 1 || clus1 != 1 && clus2 == 1 ) )
					{
						pscore = translate_and_Calign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
					}
					else
					{
						pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, NULL, &dumfl );
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
/*
		fprintf( stderr, "after align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "after align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/

		writePre( njob, name, nlen, aseq, 0 );

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

	if( treemethod == 'x' )
		fprintf( fp, "Tree = UPGMA (3).\n" );
	else if( treemethod == 's' )
		fprintf( fp, "Tree = UPGMA (2).\n" );
	else if( treemethod == 'p' )
		fprintf( fp, "Tree = UPGMA (1).\n" );
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
	static int  nlen[M];	
	static char name[M][B], **seq;
	static char **mseq1, **mseq2;
	static char **aseq;
	static char **bseq;
	static double **iscore;
	static double *eff;
	int i, j;
	int identity;
	static int ***topol;
	static double **len;
	FILE *prep;
	FILE *infp;
	FILE *hat3p;
	FILE *orderfp;
	char c;
	int alloclen;
	double total;
	LocalHom **localhomtable;
	LocalHom *tmpptr;

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

	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	seq = AllocateCharMtx( njob, nlenmax*9+1 );
	aseq = AllocateCharMtx( njob, nlenmax*9+1 );
	bseq = AllocateCharMtx( njob, nlenmax*9+1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	alloclen = nlenmax*9;

	topol = AllocateIntCub( njob, 2, njob );
	len = AllocateDoubleMtx( njob, 2 );
	iscore = AllocateDoubleMtx( njob, njob );
	eff = AllocateDoubleVec( njob );
	if( constraint )
	{
		localhomtable = (LocalHom **)calloc( njob, sizeof( LocalHom *) );
		for( i=0; i<njob; i++)
		{
			localhomtable[i] = (LocalHom *)calloc( njob, sizeof( LocalHom ) );
			for( j=0; j<njob; j++)
			{
				localhomtable[i][j].nokori = 0;
				localhomtable[i][j].start1 = -1;
				localhomtable[i][j].end1 = -1;
				localhomtable[i][j].start2 = -1;
				localhomtable[i][j].end2 = -1;
				localhomtable[i][j].overlapaa = -1.0;
				localhomtable[i][j].opt = -1.0;
				localhomtable[i][j].importance = -1.0;
				localhomtable[i][j].next = NULL;
				localhomtable[i][j].extended = -1;
				localhomtable[i][j].last = localhomtable[i]+j;
			}
		}

		fprintf( stderr, "Loading 'hat3' ... " );
		prep = fopen( "hat3", "r" );
		if( prep == NULL ) ErrorExit( "Make hat3." );
		readlocalhomtable2( prep, njob, localhomtable );
		fclose( prep );
		fprintf( stderr, "done.\n" );

#if 0
		for( i=0; i<njob; i++ )
			for( j=0; j<njob; j++ )
				fprintf( stderr, "[%d][%d].nokori = %d\n", i, j, localhomtable[i][j].nokori );
#endif

//		outlocalhom( localhomtable, njob );

	}




#if 0
	Read( name, nlen, seq );
#else
	readData( infp, name, nlen, seq );
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
		fprintf( stderr, "Illeagal character %c\n", c );
		exit( 1 );
	}

	writePre( njob, name, nlen, seq, 0 );

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
		for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) 
		{
		/*
			pscore[i][j] = (double)score_calc1( seq[i], seq[j] );
		*/
			iscore[i][j] = (double)( substitution_hosei( seq[i], seq[j] ) );
		}
	}
	else
	{
		fprintf( stderr, "Loading 'hat2' ... " );
		prep = fopen( "hat2", "r" );
		if( prep == NULL ) ErrorExit( "Make hat2." );
		readhat2( prep, njob, name, iscore );
		fclose( prep );
		fprintf( stderr, "done.\n" );

#if 0
		prep = fopen( "hat2_check", "w" );
		WriteHat2( prep, njob, name, pscore );
		fclose( prep );
#endif

	}

#if 0
	for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ )
		iscore[j][i] = iscore[i][j];
#endif

	fprintf( stderr, "Constructing dendrogram ... " );
	if( treemethod == 'x' )
		veryfastsupg( njob, iscore, topol, len );
	else if( treemethod == 's' )
		spg( njob, iscore, topol, len );
	else if( treemethod == 'p' )
		upg2( njob, iscore, topol, len );
	else 
		ErrorExit( "Incorrect tree\n" );
	fprintf( stderr, "done.\n" );

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


//	countnode( njob, topol, node0 );
	if( tbrweight )
	{
		weight = 3; 
#if 0
		utree = 0; counteff( njob, topol, len, eff ); utree = 1;
#else
		counteff_simple( njob, topol, len, eff );
#endif
	}
	else
	{
		for( i=0; i<njob; i++ ) eff[i] = 1.0;
	}

	if( constraint )
	{
#if 1
		if( thrinter > 0.0 )
		{
			fprintf( stderr, "Extending localhom ... " );
			extendlocalhom2( njob, localhomtable, iscore );
			fprintf( stderr, "done.\n" );
// job de watteoku
	    	for( i=0; i<njob-1; i++ ) 
			{
				for( j=i+1; j<njob; j++ )
				{
					if( localhomtable[i][j].opt == -1.0 ) continue;
					for( tmpptr=localhomtable[i]+j; tmpptr; tmpptr=tmpptr->next )
						tmpptr->opt /= localhomtable[i][j].nokori;
				}
			}
		}
#endif
//		resetlocalhom( njob, localhomtable ); // chuui !!!
#if 1
	    fprintf( stderr, "##### writing hat4, nlenmax = %d\n", nlenmax );
    	hat3p = fopen( "hat4", "w" );
    	if( !hat3p ) ErrorExit( "Cannot open hat4." );
    	for( i=0; i<njob-1; i++ ) 
    	{
        	for( j=i+1; j<njob; j++ )
        	{
				if( i == j ) continue;
            	for( tmpptr=localhomtable[i]+j; tmpptr; tmpptr=tmpptr->next )
            	{
                	fprintf( hat3p, "%d %d %d %6.3f %d %d %d %d %p nokori=%d ext=%d\n", i, j, tmpptr->overlapaa, tmpptr->opt * 5.8 / 600, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->next, tmpptr->nokori, tmpptr->extended ); 
					fflush( hat3p );
                	if( tmpptr->opt == -1.0 ) continue;
            	}
        	}
    	}
    	fclose( hat3p );
	    fprintf( stderr, "done.\n" );
#endif
	}

	FreeDoubleMtx( len );
	FreeDoubleMtx( iscore );

	for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );

	treebase( name, nlen, bseq, aseq, mseq1, mseq2, topol, eff, alloclen, localhomtable );


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

	writePre( njob, name, nlen, aseq, !contin );
#if 0
	writeData( stdout, njob, name, nlen, aseq );
#endif
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif

	if( constraint ) FreeLocalHomTable( localhomtable, njob );

	SHOWVERSION;
	return( 0 );
}
