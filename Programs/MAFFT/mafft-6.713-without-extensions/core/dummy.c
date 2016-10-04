#include "mltaln.h"

#define ALNLEN 1e5

int main()
{
	int res;
	static int nlen[M];
	static char name[M][B];
	static char **seq;

	seq = AllocateCharMtx( 2, ALNLEN );

	njob = 2;
	strcpy( name[0], ">name1" );
	strcpy( name[1], ">name2" );
	strcpy( seq[0], "atgcctaatgccgta" );
	strcpy( seq[1], "atgccgtaatgccgta" );

	nlen[0] = strlen( seq[0] );
	nlen[1] = strlen( seq[1] );
	nlenmax = MAX( nlen[0], nlen[1] );

	dorp = 'd';
	scoremtx = -1;

	res = disttbfast( seq, nlen, name );

	fprintf( stdout, "%s\n", name[0] );
	fprintf( stdout, "%s\n", seq[0] );
	fprintf( stdout, "%s\n", name[1] );
	fprintf( stdout, "%s\n", seq[1] );
	return( res );
}
