#include "mltaln.h"

#define DEBUG 0


static char *comment;
static char *orderfile;
static int format;

static void fillspace( char *seq, int lenmax )
{
	int len = strlen( seq );
	seq += len;
	lenmax -= len;
	while( lenmax-- ) *seq++ = ' ';
	*seq = 0;
}

void setmark_clustal( int nlen, int nseq, char **seq, char *mark )
{
	int i, j, k;

	char *strong[] = { 
					"STA",
					"NEQK",
					"NHQK",
					"NDEQ",
					"QHRK",
					"MILV",
					"MILF",
					"HY",
					"FYW",
				  };
	int nstrong = 9;
	char *weaker[] = { 
					"CSA",
					"ATV",
					"SAG",
					"STNK",
					"STPA",
					"SGND",
					"SNDEQK",
					"NDEQHK",
					"NEQHRK",
					"FVLIM",
					"HFY",
				  };
	int nweaker = 11;

	for( i=0; i<nlen; i++ )
	{
		mark[i] = ' ';
		for( j=0; j<nseq; j++ )
			if( '-' == seq[j][i] ) break;
		if( j != nseq ) 
		{
			continue;
		}
		for( j=0; j<nseq; j++ )
			if( seq[0][i] != seq[j][i] ) break;
		if( j == nseq ) 
		{
			mark[i] = '*';
			continue;
		}
		for( k=0; k<nstrong; k++ )
		{
			for( j=0; j<nseq; j++ )
			{
				if( !strchr( strong[k], seq[j][i] ) ) break;
			}
			if( j == nseq ) break;
		}
		if( k < nstrong )
		{
			mark[i] = ':';
			continue;
		}
		for( k=0; k<nweaker; k++ )
		{
			for( j=0; j<nseq; j++ )
			{
				if( !strchr( weaker[k], seq[j][i] ) ) break;
			}
			if( j == nseq ) break;
		}
		if( k < nweaker )
		{
			mark[i] = '.';
			continue;
		}
	}
	mark[nlen] = 0;
}

void setmark( int nlen, int nseq, char **seq, char *mark )
{
	int i, j;

	for( i=0; i<nlen; i++ )
	{
		mark[i] = ' ';
		for( j=0; j<nseq; j++ )
			if( '-' == seq[j][i] ) break;
		if( j != nseq ) 
		{
			continue;
		}
		for( j=0; j<nseq; j++ )
			if( seq[0][i] != seq[j][i] ) break;
		if( j == nseq ) 
		{
			mark[i] = '*';
			continue;
		}
		for( j=0; j<nseq; j++ )
			if( amino_grp[(int)seq[0][i]] != amino_grp[(int)seq[j][i]] ) break;
		if( j == nseq ) 
		{
			mark[i] = '.';
			continue;
		}
	}
	mark[nlen] = 0;
}

void arguments( int argc, char *argv[] )
{
    int c;
	scoremtx = 1;
	nblosum = 62;
	dorp = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	inputfile = NULL;
	comment = NULL;
	orderfile = NULL;
	format = 'c';

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
				case 'c':
					comment = *++argv;
					fprintf( stderr, "comment = %s\n", comment );
					--argc;
					goto nextoption;
				case 'r':
					orderfile = *++argv;
					fprintf( stderr, "orderfile = %s\n", orderfile );
					--argc;
					goto nextoption;
				case 'f':
					format = 'f';
					break;
				case 'm':
					format = 'm';
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
    if( argc != 0 ) 
    {
        fprintf( stderr, "options: Check source file !\n" );
        exit( 1 );
    }
}


int main( int argc, char *argv[] )
{
	static int  *nlen;	
	static char **name, **seq, *mark;
	static int *order;
	int i;
	FILE *infp;
	FILE *orderfp;
	char gett[B];

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

	seq = AllocateCharMtx( njob, nlenmax*2+1 );
	mark = AllocateCharVec( nlenmax*2+1 );
	order = AllocateIntVec( njob );
	name = AllocateCharMtx( njob, B+1 );
    nlen = AllocateIntVec( njob );


	if( orderfile )
	{
		orderfp = fopen( orderfile, "r" );
		if( !orderfile )
		{
			fprintf( stderr, "Cannot open %s\n", orderfile );
			exit( 1 );
		}
		for( i=0; i<njob; i++ )
		{
			fgets( gett, B-1, orderfp );
			order[i] = atoi( gett );
		}
		fclose( orderfp );
	}
	else
	{
		for( i=0; i<njob; i++ ) order[i] = i;
	}

	readData_pointer( infp, name, nlen, seq );
	fclose( infp );

	if( format == 'c' ) for( i=0; i<njob; i++ ) fillspace( seq[i], nlenmax );
	constants( njob, seq );

//	initSignalSM();

//	initFiles();



//	setmark( nlenmax, njob, seq, mark );
	setmark_clustal( nlenmax, njob, seq, mark );

	if( format == 'f' )
		writeData_reorder_pointer( stdout, njob, name, nlen, seq, order );
	else if( format == 'c' )
		clustalout_pointer( stdout, njob, nlenmax, seq, name, mark, comment, order );
	else if( format == 'm' )
		miyataout_reorder_pointer( stdout, njob, nlenmax, name, nlen, seq, order );
	else
		fprintf( stderr, "Unknown format\n" );

//	SHOWVERSION;
	return( 0 );
}
