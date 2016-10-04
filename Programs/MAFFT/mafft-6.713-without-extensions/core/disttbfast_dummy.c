#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0

#define END_OF_VEC -1


void arguments()
{
    int c;

	inputfile = NULL;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 1;
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
    alg = 'M';
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
	treemethod = 'x';
	contin = 0;
	kobetsubunkatsu = 0;
	makedistmtx = 1;
	ppenalty = -1530;
	ppenalty_ex = NOTSPECIFIED;
	poffset = -123;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	TMorJTT = JTT;
}

static int maxl;
static int tsize;

void seq_grp_nuc( int *grp, char *seq )
{
	int tmp;
	while( *seq )
	{
		tmp = amino_grp[*seq++];
		if( tmp < 4 )
			*grp++ = tmp;
		else
			fprintf( stderr, "WARNING : Unknown character %c\r", *(seq-1) );
	}
	*grp = END_OF_VEC;
}

void seq_grp( int *grp, char *seq )
{
	int tmp;
	while( *seq )
	{
		tmp = amino_grp[*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
			fprintf( stderr, "WARNING : Unknown character %c\r", *(seq-1) );
	}
	*grp = END_OF_VEC;
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

	if( !memo )
	{
		memo = calloc( tsize, sizeof( short ) );
		if( !memo ) ErrorExit( "Cannot allocate memo\n" );
		ct = calloc( MIN( maxl, tsize )+1, sizeof( int ) );
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

void treebase( char name[M][B], int nlen[M], char **seq, char **aseq, char **mseq1, char **mseq2, int ***topol, double *effarr, int alloclen )
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
	float dumfl = 0.0;
	double total;
	int intdum;


	if( effarr1 == NULL ) 
	{
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
	}

#if 0
	fprintf( stderr, "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );



//	writePre( njob, name, nlen, aseq, 0 );

	tscore = 0.0;
	for( l=0; l<njob-1; l++ )
	{

		clus1 = fastconjuction( topol[l][0], aseq, mseq1, effarr1, effarr, name, name1, indication1 );
		clus2 = fastconjuction( topol[l][1], aseq, mseq2, effarr2, effarr, name, name2, indication2 );


		fprintf( stderr, "STEP %d /%d\r", l+1, njob-1 );
//		fprintf( stderr, "STEP %d /%d\n", l+1, njob-1 );
//		fprintf( stderr, "group1 = %.66s", indication1 );
//		if( strlen( indication1 ) > 66 ) fprintf( stderr, "..." );
//		fprintf( stderr, "\n" );
//		fprintf( stderr, "group2 = %.66s", indication2 );
//		if( strlen( indication2 ) > 66 ) fprintf( stderr, "..." );
//		fprintf( stderr, "\n" );
//		for( i=0; i<clus1; i++ ) fprintf( stderr, "## STEP%d-eff for mseq1-%d %f\n", l+1, i, effarr1[i] );

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

		if( use_fft )
		{
			if( alg == 'M' )
				pscore = Falign_noudp( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, &intdum );
			else
				pscore = Falign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, &intdum );
		}
		else
		{
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
					break;
				case( 'M' ):
					pscore = MSalignmm( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, NULL, NULL, NULL, NULL );
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
						pscore = G__align11( mseq1, mseq2, alloclen );
					}
					else
					{
						pscore = A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL );
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
	 

int disttbfast( char **in, int nlen[M], char name[M][B] )
{
	static char **seq;
	static char **mseq1, **mseq2;
	static char **aseq;
	static char **bseq;
	static int **pscore;
	static double *eff;
	int i, j;
	int identity;
	static int ***topol;
	static double **len;
	FILE *prep;
	FILE *infp;
	char c;
	int alloclen;
	double total;

	FILE *fp;
	FILE *orderfp;
	int *grpseq;
	char *tmpseq;
	int  **pointt;
	double **mtx;
	double **mtx2; 
	int score0; 
	static short *table1;
	char b[B];

	arguments();

	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	seq = in;
	aseq = AllocateCharMtx( njob, nlenmax*5+1 );
	bseq = AllocateCharMtx( njob, nlenmax*5+1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	alloclen = nlenmax*5;

	topol = AllocateIntCub( njob, 2, njob );
	len = AllocateDoubleMtx( njob, 2 );
	pscore = AllocateIntMtx( njob, njob );
	eff = AllocateDoubleVec( njob );



	constants( njob, seq );


#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();


	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illeagal character %c\n", c );
		exit( 1 );
	}

	fprintf( stderr, "\n" );

//	writePre( njob, name, nlen, seq, 0 );
	if( makedistmtx )
	{
		fprintf( stderr, "Making a distance matrix ..\n" );

	    tmpseq = AllocateCharVec( nlenmax+1 );
		grpseq = AllocateIntVec( nlenmax+1 );
		pointt = AllocateIntMtx( njob, nlenmax+1 );
    	mtx = AllocateDoubleMtx( njob, njob ); 
		if( dorp == 'd' ) tsize = (int)pow( 4, 6 );
		else              tsize = (int)pow( 6, 6 );

		maxl = 0;
		for( i=0; i<njob; i++ ) 
		{
			gappick0( tmpseq, seq[i] );
			nlen[i] = strlen( tmpseq );
			if( nlen[i] < 6 )
			{
				fprintf( stderr, "Seq %d, too short, %d characters\n", i+1, nlen[i] );
				exit( 1 );
			}
			if( nlen[i] > maxl ) maxl = nlen[i];
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
			table1 = calloc( tsize, sizeof( short ) );
			if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
			if( i % 10 == 0 )
			{
				fprintf( stderr, "%#4d / %#4d\r", i+1, njob );
			}
			makecompositiontable_p( table1, pointt[i] );
	
			for( j=i; j<njob; j++ ) 
			{
				mtx[i][j] = commonsextet_p( table1, pointt[j] );
			} 
			free( table1 );
		}
		for( i=0; i<njob; i++ )
		{
			score0 = mtx[i][i];
			for( j=0; j<njob; j++ ) 
				pscore[i][j] = (int)( ( score0 - mtx[MIN(i,j)][MAX(i,j)] ) / score0 * 3 * INTMTXSCALE + 0.5 );
		}
		for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) 
		{
			pscore[i][j] = MIN( pscore[i][j], pscore[j][i] );
		}
	
		if( disopt )
		{
			for( i=0; i<njob; i++ ) 
			{
				sprintf( b, "=lgth = %#04d", nlen[i] );
				strins( b, name[i] );
			}
		}
		FreeDoubleMtx( mtx );
		free( tmpseq );
		free( grpseq );
		FreeIntMtx( pointt );
		fprintf( stderr, "\ndone.\n\n" );
	}
	else if( tbutree == 0 )
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
			pscore[i][j] = (int)( substitution_hosei( seq[i], seq[j] ) * INTMTXSCALE + 0.5 );
		}
	}
	else
	{
#if 1
		fprintf( stderr, "Loading 'hat2' ... " );
		prep = fopen( "hat2", "r" );
		if( prep == NULL ) ErrorExit( "Make hat2." );
		readhat2_int( prep, njob, name, pscore );
		fclose( prep );
		fprintf( stderr, "done.\n" );
#endif
	}
#if 0
		prep = fopen( "hat2_check", "w" );
		WriteHat2_int( prep, njob, name, pscore );
		fclose( prep );
#endif

	fprintf( stderr, "Constructing dendrogram ... " );
	if( treemethod == 'x' )
	{
		veryfastsupg_int( njob, pscore, topol, len );
	}
	else 
		ErrorExit( "Incorrect tree\n" );
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
	

	for( i=0; i<njob; i++ )
	{
		len[i][0] /= INTMTXSCALE;
		len[i][1] /= INTMTXSCALE;
	}

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


	for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );

	FreeIntMtx( pscore );
	FreeDoubleMtx( len );
	fprintf( stderr, "Progressive alignment ... \n" );
	treebase( name, nlen, bseq, aseq, mseq1, mseq2, topol, eff, alloclen );
	fprintf( stderr, "\ndone.\n\n" );

#if DEBUG
	fprintf( stderr, "writing alignment to stdout\n" );
#endif


	for( i=0; i<njob; i++ ) strcpy( seq[i], aseq[i] );

	fprintf( stderr, "mafft v%s\n", VERSION );
	return( 0 );
}
