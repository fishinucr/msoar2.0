#include "mltaln.h"
#include "mtxutl.h"

#define DEBUG 0
#define TEST  0

#define END_OF_VEC -1

static int maxl;
static int tsize;
static float lenfaca, lenfacb, lenfacc, lenfacd;
#define PLENFACA 0.01
#define PLENFACB 10000
#define PLENFACC 10000
#define PLENFACD 0.1
#define DLENFACA 0.01
#define DLENFACB 2500
#define DLENFACC 2500
#define DLENFACD 0.1

void arguments( int argc, char *argv[] )
{
	int c;

	inputfile = NULL;
	scoremtx = 1;
	nblosum = 62;
	dorp = NOTSPECIFIED;
	alg = 'X';

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
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
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
	if( inputfile == NULL )
	{
		argc--;
		inputfile = *argv;
		fprintf( stderr, "inputfile = %s\n", inputfile );
	}
    if( argc != 0 )
    {
        fprintf( stderr, "Usage: mafft-distance [-PD] [-i inputfile] inputfile > outputfile\n" );
        exit( 1 );
    }
}

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
		ct = (int *)calloc( MIN( maxl, tsize)+1, sizeof( int ) );
		if( !ct ) ErrorExit( "Cannot allocate memo\n" );
	}

	cp = ct;
	while( ( point = *pointt++ ) != END_OF_VEC )
	{
		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
		if( tmp == 0 ) *cp++ = point;
//		fprintf( stderr, "cp - ct = %d (tsize = %d)\n", cp - ct, tsize );
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

int main( int argc, char **argv )
{
	int i, j;
	FILE *infp;
	char **seq;
	int *grpseq;
	char *tmpseq;
	int  **pointt;
	static char **name;
	static int *nlen;
	double *mtxself;
	float score;
	static short *table1;
	float longer, shorter;
	float lenfac;
	float bunbo;

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

	tmpseq = AllocateCharVec( nlenmax+1 );
	seq = AllocateCharMtx( njob, nlenmax+1 );
	grpseq = AllocateIntVec( nlenmax+1 );
	pointt = AllocateIntMtx( njob, nlenmax+1 );
	mtxself = AllocateDoubleVec( njob );
	pamN = NOTSPECIFIED;
	name = AllocateCharMtx( njob, B );
	nlen = AllocateIntVec( njob );

#if 0
	FRead( infp, name, nlen, seq );
#else
	readData_pointer( infp, name, nlen, seq );
#endif

	fclose( infp );

	constants( njob, seq );

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
		nlen[i] = strlen( tmpseq );
//		if( nlen[i] < 6 )
//		{
//			fprintf( stderr, "Seq %d, too short, %d characters\n", i+1, nlen[i] );
//			exit( 1 );
//		}
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
	fprintf( stderr, "\nCalculating i-i scores ... " );
	for( i=0; i<njob; i++ )
	{
		table1 = (short *)calloc( tsize, sizeof( short ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		makecompositiontable_p( table1, pointt[i] );

		score = commonsextet_p( table1, pointt[i] );
		mtxself[i] = score;
		free( table1 );
	}

	fprintf( stderr, "done.\n" );
	fprintf( stderr, "\nCalculating i-j scores ... \n" );
	for( i=0; i<njob; i++ )
	{
		table1 = (short *)calloc( tsize, sizeof( short ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 10 == 0 )
		{
			fprintf( stderr, "%4d / %4d\r", i+1, njob );
		}
		makecompositiontable_p( table1, pointt[i] );


        for( j=i+1; j<njob; j++ ) 
        {
			if( nlen[i] > nlen[j] )
			{
				longer=(float)nlen[i];
				shorter=(float)nlen[j];
			}
			else
			{
				longer=(float)nlen[j];
				shorter=(float)nlen[i];
			}
//			lenfac = 3.0 / ( LENFACA + LENFACB / ( longer + LENFACC ) + shorter / longer * LENFACD );
			lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
//			lenfac = 1.0;
//			fprintf( stderr, "lenfac = %f (%.0f,%.0f)\n", lenfac, longer, shorter );
			score = commonsextet_p( table1, pointt[j] );
			bunbo = MIN( mtxself[i], mtxself[j] );
			if( bunbo == 0.0 )
				fprintf( stdout, "%d-%d d=%4.2f l=%d,%d\n", i+1, j+1, 1.0, nlen[i], nlen[j] );
			else
				fprintf( stdout, "%d-%d d=%4.2f l=%d,%d\n", i+1, j+1, ( 1.0 - score / bunbo ) * lenfac, nlen[i], nlen[j] );
//			fprintf( stderr, "##### mtx = %f, mtx[i][0]=%f, mtx[j][0]=%f, bunbo=%f\n", mtx[i][j-i], mtx[i][0], mtx[j][0], bunbo );
//          score = (double)commonsextet_p( table1, pointt[j] );
//			fprintf( stdout, "%d-%d d=%4.2f l=%d,%d\n", i+1, j+1, ( 1.0 - score / MIN( mtxself[i], mtxself[j] ) ) * 3, nlen[i], nlen[j] );


		} 
		free( table1 );
	}
		
	fprintf( stderr, "\n" );
	SHOWVERSION;
	exit( 0 );
}
