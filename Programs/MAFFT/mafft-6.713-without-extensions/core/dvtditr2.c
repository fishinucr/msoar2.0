 /* Tree-dependent-iteration */
 /* Devide to segments       */ 

#include "mltaln.h"

extern char **seq_g;
extern char **res_g;

void arguments( int argc, char *argv[] )
{
	int c;

	inputfile = NULL;
	score_check = 1;
	fftkeika = 1;
	constraint = 0;
	fmodel = 0;
	kobetsubunkatsu = 1;
	bunkatsu = 1;
	nblosum = 80;
	niter = 100;
	calledByXced = 0;
	devide = 1;
	divWinSize = 20; /* 70 */
	divThreshold = 65;
	fftscore = 1;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
    scmtd = 5;
	cooling = 0;
    weight = 4;
    utree = 1;
    refine = 1;
    check = 1;
    cut = 0.0;
	disp = 0;
	outgap = 1;
	use_fft = 0;
	alg = 'A';  /* chuui */
	mix = 0;
	checkC = 0;
	tbitr = 0;
	treemethod = 'x';
	scoremtx = 0;
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
				case 'I':
					niter = atoi( *++argv );
					fprintf( stderr, "niter = %d\n", niter );
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
					fprintf( stderr, "fastathreshold %f\n", fastathreshold );
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
				case 'Q':
					calledByXced = 1;
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
				case 'T':
					kobetsubunkatsu = 0;
					break;
				case 'B':
					bunkatsu = 0;
					break;
				case 'c':
					cooling = 1;
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
				case 't':
					weight = 4;
					break;
				case 'u':
					weight = 0;
					break;
				case 'J':
					utree = 0;
					break;
				case 'd':
					disp = 1;
					break;
				case 'Z':
					score_check = 0;
					break;
				case 'Y':
					score_check = 2;
					break;
				case 'n' :
					treemethod = 'n';
					break;
				case 's' :
					treemethod = 's';
					break;
				case 'x' :
					treemethod = 'x';
					break;
				case 'p' :
					treemethod = 'p';
					break;
				case 'z':
					fftThreshold = atoi( *++argv );
					--argc;
					goto nextoption;
				case 'w':
					fftWinSize = atoi( *++argv );
					--argc;
					goto nextoption;
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
		fprintf( stderr, "options : Check source file!\n" );
		exit( 1 );
	}
#if 0
	if( alg == 'A' && weight == 0 ) 
		ErrorExit( "ERROR : Algorithm A+ and un-weighted\n" ); 
#endif
}


int main( int argc, char *argv[] )
{
    int identity;
	static int nlen[M];
	static char name[M][B], **seq, **aseq, **bseq;
	static Segment *segment = NULL;
	static int anchors[MAXSEG];
	int i, j;
	int iseg, nseg;
	int ***topol;
	double **len;
	double **eff;
	FILE *prep;
	FILE *infp;
	int alloclen;
	int returnvalue;
	char c;
	int ocut;
	char **seq_g_bk;
	LocalHom **localhomtable;

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

#if 0
	PreRead( stdin, &njob, &nlenmax );
#else
	getnumlen( infp );
#endif
	rewind( infp );

	if( njob < 2 ) 
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	ocut = cut;

	segment = (Segment *)calloc( MAXSEG, sizeof( Segment ) );
	topol = AllocateIntCub( njob, 2, njob );
	len = AllocateDoubleMtx( njob, 2 );
	eff = AllocateDoubleMtx( njob, njob );
	seq = AllocateCharMtx( njob, nlenmax*5+1 );
	seq_g = AllocateCharMtx( njob, nlenmax*5+1 );
	res_g = AllocateCharMtx( njob, nlenmax*5+1 );
	aseq = AllocateCharMtx( njob, nlenmax*5+1 );
	bseq = AllocateCharMtx( njob, nlenmax*5+1 );
	alloclen = nlenmax * 5;
	seq_g_bk = AllocateCharMtx( njob, 0 );
	for( i=0; i<njob; i++ ) seq_g_bk[i] = seq_g[i];

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
			}
		}
		fprintf( stderr, "Loading 'hat3' ... " );
		prep = fopen( "hat3", "r" );
		if( prep == NULL ) ErrorExit( "Make hat3." );
		readlocalhomtable( prep, njob, localhomtable );
		fclose( prep ); 
//		for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ )
//			fprintf( stdout, "%d %d %d %d %d %d %d\n", i, j, localhomtable[i][j].opt, localhomtable[i][j].start1, localhomtable[i][j].end1, localhomtable[i][j].start2, localhomtable[i][j].end2 );
		fprintf( stderr, "done.\n" );
#if 0
		fprintf( stderr, "Extending localhom ... " );
		extendlocalhom( njob, localhomtable );
		fprintf( stderr, "done.\n" );
#endif
	}

#if 0
	Read( name, nlen, seq_g );
#else
	readData( infp, name, nlen, seq_g );
#endif
	fclose( infp );


	for( i=0; i<njob; i++ )
	{
		res_g[i][0] = 0;
	}

	identity = 1;
	for( i=0; i<njob; i++ ) 
	{
		identity *= ( nlen[i] == nlen[0] );
	}
	if( !identity ) 
	{
		fprintf( stderr, "Input pre-aligned data\n" );
		exit( 1 );
	}
	constants( njob, seq_g );

#if 0
	fprintf( stderr, "penalty = %d\n", penalty ); 
	fprintf( stderr, "penalty_ex = %d\n", penalty_ex ); 
	fprintf( stderr, "offset = %d\n", offset ); 
#endif

	initSignalSM();

	initFiles();

#if 1
	if( njob == 2 )
	{
		writePre( njob, name, nlen, seq_g, 1 );
		SHOWVERSION;
		return( 0 );
	}
#endif

	c = seqcheck( seq_g );
	if( c )
	{
		fprintf( stderr, "Illeagal character %c\n", c );
		exit( 1 );
	}
	commongappick( njob, seq_g );



	if( utree )
	{
		prep = fopen( "hat2", "r" );
		if( !prep ) ErrorExit( "Make hat2." );
		readhat2( prep, njob, name, eff );
		fclose( prep );
#if DEBUG
		for( i=0; i<njob-1; i++ ) 
		{
			for( j=i+1; j<njob; j++ ) 
			{
				printf( " %f", eff[i][j] );
			}
			printf( "\n" );
		}
#endif
		if     ( treemethod == 'x' ) 
			veryfastsupg( njob, eff, topol, len );
		else if( treemethod == 'n' ) 
			nj( njob, eff, topol, len );
		else if( treemethod == 's' )
			spg( njob, eff, topol, len );
		else if( treemethod == 'p' )
			upg2( njob, eff, topol, len );
		else ErrorExit( "Incorrect treemethod.\n" );
	}
#if DEBUG
	printf( "utree = %d\n", utree );
#endif

	if( ( !utree && kobetsubunkatsu ) || constraint || !bunkatsu )
	{
		nseg = 0;
		anchors[0] =0;
		anchors[1] =strlen( seq_g[0] );
		nseg += 2;
	}
	else
	{
		nseg = searchAnchors( njob, seq_g, segment );
#if 0
		fprintf( stderr, "nseg = %d\n", nseg );
		fprintf( stderr, "seq_g[0] = %s\n", seq_g[0] );
		fprintf( stderr, "nlenmax = %d\n", nlenmax );
		fprintf( stderr, "len = %d\n", strlen( seq_g[0] ) );
#endif

		anchors[0] = 0;
		for( i=0; i<nseg; i++ ) anchors[i+1] = segment[i].center;
		anchors[nseg+1] = strlen( seq_g[0] );
		nseg +=2;

#if 0
		for( i=0; i<nseg; i++ )
			fprintf( stderr, "anchor[%d] = %d\n", i, anchors[i] );
#endif
	}

	for( i=0; i<njob; i++ ) res_g[i][0] = 0;

	for( iseg=0; iseg<nseg-1; iseg++ )
	{
		int tmplen = anchors[iseg+1]-anchors[iseg];
		int pos = strlen( res_g[0] );
		for( j=0; j<njob; j++ )
		{
			strncpy( seq[j], seq_g[j], tmplen );
			seq[j][tmplen]= 0;
			seq_g[j] += tmplen;	

		}
		fprintf( stderr, "Segment %3d/%3d %4d-%4d\n", iseg+1, nseg-1, pos+1, pos+1+tmplen );
		fprintf( trap_g, "Segment %3d/%3d %4d-%4d\n", iseg+1, nseg-1, pos+1, pos+1+tmplen );
	
		cut = ocut;
		returnvalue = TreeDependentIteration( njob, name, nlen, seq, bseq, topol, len, alloclen, localhomtable );

		for( i=0; i<njob; i++ )
			strcat( res_g[i], bseq[i] );
	}
	FreeCharMtx( seq_g_bk );
	FreeIntCub( topol );
	FreeDoubleMtx( len );
	FreeDoubleMtx( eff );
	fprintf( stderr, "constraint = %d\n", constraint );
	if( constraint ) FreeLocalHomTable( localhomtable, njob );

#if 0
	Write( stdout, njob, name, nlen, bseq );
#endif

	fprintf( stderr, "done\n" );
	fprintf( trap_g, "done\n" );
	fclose( trap_g );


	devide = 0; 
	writePre( njob, name, nlen, res_g, 1 );
#if 0
	writeData( stdout, njob, name, nlen, res_g, 1 );
#endif


	SHOWVERSION;
	return( 0 );
}

#if 0
signed int main( int argc, char *argv[] )
{
	int i, nlen[M];
	char b[B];
	char a[] = "=";
	int value;

	gets( b ); njob = atoi( b );

/*
	scoremtx = 0;
	if( strstr( b, "ayhoff" ) ) scoremtx = 1;
	else if( strstr( b, "dna" ) || strstr( b, "DNA" ) ) scoremtx = -1;
	else if( strstr( b, "M-Y" ) || strstr( b, "iyata" ) ) scoremtx = 2;
	else scoremtx = 0;
*/
	if( strstr( b, "constraint" ) ) cnst = 1;

	nlenmax = 0;
	i = 0;
	while( i<njob )
	{
		gets( b );
		if( !strncmp( b, a, 1 ) ) 
		{
			gets( b ); nlen[i] = atoi( b );
			if( nlen[i] > nlenmax ) nlenmax = nlen[i];
			i++;
		}
	}
	if( nlenmax > N || njob > M ) 
	{
		fprintf( stderr, "ERROR in main\n" );
		exit( 1 );
	}
	/*
	nlenmax = Na;
	*/
	rewind( stdin );
	value = main1( nlen, argc, argv );
	exit( 0 );
}
#endif
