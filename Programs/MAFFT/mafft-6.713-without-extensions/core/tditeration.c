
/* 
	tree-dependent iteration   
    algorithm A+ when group-to-group, C when group-to-singleSeqence 
	                 OR
    algorithm A+
*/

#include "mltaln.h"

#define DEBUG 0
#define RECORD 0

extern char **seq_g;
extern char **res_g;

static void Writeoption2( FILE *fp, int cycle, double cut )
{
	fprintf( fp, "%dth cycle\n", cycle );
    fprintf( fp, "marginal score to search : current score * (100-%d) / 100\n", (int)cut );
}

static void Writeoptions( FILE *fp )
{
	fprintf( fp, "Tree-dependent-iteration\n" );
    if( scoremtx == 1 )
        fprintf( fp, "Blosum %d\n", nblosum );
    else if( scoremtx == -1 )
        fprintf( fp, "DNA\n" );
    else if( scoremtx == 2 )
        fprintf( fp, "Miyata-Yasunaga\n" );
	else
		fprintf( fp, "JTT %dPAM\n", pamN );

	if( scoremtx == 0 || scoremtx == 1 )
    	fprintf( fp, "Gap Penalty = %+5.3f, %5.2f, %+5.3f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
	else
		fprintf( fp, "Gap Penalty = %+5.3f\n", (double)penalty/1000 );
		

    if( scmtd == 3 )
        fprintf( fp, "score of rnd or sco\n" );
    else if( scmtd == 4 )
        fprintf( fp, "score = sigma( score for a pair of homologous amino acids ) / ( number of amino acids pairs )\n" );
    else if( scmtd == 5 )
        fprintf( fp, "score : SP\n" );
    if( mix )
        fprintf( fp, "?\n" );
    else
    { 
        if( weight == 2 )
            fprintf( fp, "weighted rationale-1,  geta2 = %f\n", geta2 );
        else if( weight == 3 )
            fprintf( fp, "weighted like ClustalW," );
        else if( weight == 4 )
            fprintf( fp, "weighted rationale-2,  geta2 = %f\n", geta2 );
        else
            fprintf( fp, "unweighted\n" );
    }
    if( weight && utree )
        fprintf( fp, "using tree defined by the file hat2.\n" );
    if( weight && !utree )
        fprintf( fp, "using temporary tree.\n" );

	if( treemethod == 'n' )
		fprintf( fp, "Tree is calculated with Neighbor-Joining Method.\n" );
	else if( treemethod == 'q' )
		fprintf( fp, "Tree is calculated with Minimum linkage.\n" );
	else if( treemethod == 'X' )
		fprintf( fp, "Tree is calculated with simplified UPG Method and UPG Method.\n" );
	else if( treemethod == 'E' )
		fprintf( fp, "Tree is calculated with UPG Method.\n" );
	else
		fprintf( fp, "Tree is calculated with unknown Method.\n" );
		
	if( alg == 'C' ) 
		fprintf( fp, "Algorithm A+ / C\n" );
	else if( alg == 'S' ) 
		fprintf( fp, "Algorithm S \n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Algorithm A+ \n" );
	else if( alg == 'a' ) 
		fprintf( fp, "Algorithm A \n" );
	else 
		fprintf( fp, "Algorithm ? \n" );

	if( use_fft )
	{
		if( scoremtx == -1 )
		{
			fprintf( fp, "Basis : 4 nucleotides\n" );
		}
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
}


int TreeDependentIteration( int locnjob, char name[M][B], int nlen[M], 
                             char **aseq, char **bseq, int ***topol, double **len, 
                             int alloclen, LocalHom **localhomtable, 
							 RNApair ***singlerna,
							 int nkozo, char *kozoarivec )
{
	int i, j, k, l, iterate, ii, iu, ju;
	int lin, ldf, length; 
	int clus1, clus2;
	int s1, s2;
	static double **imanoten;
	static Node *stopol;
	static double *effarrforlocalhom = NULL;
	static double *effarr = NULL;
	static double *effarr1 = NULL;
	static double *effarr2 = NULL;
	static double *effarr_kozo = NULL;
	static double *effarr1_kozo = NULL;
	static double *effarr2_kozo = NULL;
	static double **mtx = NULL;
	static int **node = NULL;
	static int *branchnode = NULL;
	static double **branchWeight = NULL;
	static char **mseq1, **mseq2;
	static float ***history;
	FILE *trap;
	double tscore, mscore;
	int identity;
	int converged;
	int oscillating;
	float naivescore0 = 0.0; // by D.Mathog, a guess
	float naivescore1;
#if 0
	char pair[njob][njob];
#else
	static char **pair;
#endif
#if DEBUG + RECORD
	double score_for_check0, score_for_check1;
	static double **effmtx = NULL;
	extern double score_calc0();
#endif
	static char *indication1, *indication2;
	static LocalHom ***localhomshrink = NULL;
	float impmatch, oimpmatch;
	static int *gapmap1;
	static int *gapmap2;
	double tmpdouble;
	int intdum;
	static RNApair *rnapairboth;
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

	Writeoptions( trap_g );
	fflush( trap_g );

	if( effarr == NULL ) /* locnjob == njob ni kagiru */
	{
		indication1 = AllocateCharVec( njob*3+50 );
		indication2 = AllocateCharVec( njob*3+50 );
		effarr = AllocateDoubleVec( locnjob );
		effarrforlocalhom = AllocateDoubleVec( locnjob );
		effarr1 = AllocateDoubleVec( locnjob );
		effarr2 = AllocateDoubleVec( locnjob );
		mseq1 = AllocateCharMtx( locnjob, 0 );
		mseq2 = AllocateCharMtx( locnjob, 0 );
		mtx = AllocateDoubleMtx( locnjob, locnjob );
		node = AllocateIntMtx( locnjob, locnjob );
		branchnode = AllocateIntVec( locnjob );
		branchWeight = AllocateDoubleMtx( locnjob, 2 );
		history = AllocateFloatCub( niter, locnjob, 2 );
		stopol = (Node *)calloc( locnjob * 2, sizeof( Node ) );
		gapmap1 = AllocateIntVec( alloclen );
		gapmap2 = AllocateIntVec( alloclen );
		if( score_check == 2 ) imanoten = AllocateDoubleMtx( njob, njob );

		effarr1_kozo = AllocateDoubleVec( locnjob ); // tsuneni allocate suru.
		effarr2_kozo = AllocateDoubleVec( locnjob ); // tsuneni allocate suru.
		effarr_kozo = AllocateDoubleVec( locnjob ); 
		for( i=0; i<locnjob; i++ )
			effarr_kozo[i] = effarr1_kozo[i] = effarr2_kozo[i] = 0.0;

#if 0
#else
		pair = AllocateCharMtx( locnjob, locnjob );
		if( rnakozo ) rnapairboth = (RNApair *)calloc( alloclen, sizeof( RNApair ) );

		if( constraint )
		{
			localhomshrink = (LocalHom ***)calloc( njob, sizeof( LocalHom ** ) );
			for( i=0; i<njob; i++)
			{
				localhomshrink[i] = (LocalHom **)calloc( njob, sizeof( LocalHom * ) );
			}
		}
#endif
	}
#if DEBUG + RECORD
	if( !effmtx ) effmtx = AllocateDoubleMtx( locnjob, locnjob );
	for( i=0; i<locnjob; i++ ) for( j=0; j<locnjob; j++ ) effmtx[i][j] = 1.0;
#endif

	for( i=0; i<locnjob; i++ ) strcpy( bseq[i], aseq[i] );

	writePre( locnjob, name, nlen, aseq, 0 );

	if( utree )
	{
		if( constraint )
		{
			counteff_simple( locnjob, topol, len, effarrforlocalhom );
			calcimportance( locnjob, effarrforlocalhom, aseq, localhomtable );
		}

		if( weight == 2 ) 
		{
			countnode_int( locnjob, topol, node );
			if( nkozo )
			{
				fprintf( stderr, "Not supported, weight=%d nkozo=%d.\n", weight, nkozo );
			}
		}
		else if( weight == 4 )
		{
			treeCnv( stopol, locnjob, topol, len, branchWeight );
			calcBranchWeight( branchWeight, locnjob, stopol, topol, len );
		}
	}

	converged = 0;
	if( cooling ) cut *= 2.0;
	for( iterate = 0; iterate<niter; iterate++ ) 
	{
		if( cooling ) cut *= 0.5; /* ... */

		fprintf( trap_g, "\n" );
		Writeoption2( trap_g, iterate, cut );
		fprintf( trap_g, "\n" );

		if( utree == 0 )
		{
			if( nkozo )
			{
				fprintf( stderr, "Dynamic tree and kozo is not supported.\n" );
				exit( 1 );
			}
			if( devide )
			{
				static char *buff1 = NULL;
				static char *buff2 = NULL;
				if( !buff1 )
				{
					buff1 = AllocateCharVec( alloclen );
					buff2 = AllocateCharVec( alloclen );
				}	

				for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) 	
				{
					buff1[0] = buff2[0] = 0;
					strcat( buff1, res_g[i] ); strcat( buff2, res_g[j] );
					strcat( buff1,  bseq[i] ); strcat( buff2,  bseq[j] );
					strcat( buff1, seq_g[i] ); strcat( buff2, seq_g[j] );

					mtx[i][j] = (double)substitution_hosei( buff1, buff2 );
				}
			}
			else
			{
				for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) 	
					mtx[i][j] = (double)substitution_hosei( bseq[i], bseq[j] );
			}

			if     ( treemethod == 'n' )
				nj( locnjob, mtx, topol, len );
			else if( treemethod == 's' )
				spg( locnjob, mtx, topol, len );
			else if( treemethod == 'X' )
				supg( locnjob, mtx, topol, len );
			else if( treemethod == 'p' )
				upg2( locnjob, mtx, topol, len );
			/* veryfastsupgは、今のところ使えません。*/
			/* 順番の問題があるので                  */

			if( weight == 2 )
				countnode_int( locnjob, topol, node );
			else if( weight == 4 )
			{
				treeCnv( stopol, locnjob, topol, len, branchWeight );
				calcBranchWeight( branchWeight, locnjob, stopol, topol, len );
			}
			trap = fopen( "hat2", "w" );
			if( !trap ) ErrorExit( "Cannot open hat2." );
			WriteHat2( trap, locnjob, name, mtx );
			fclose( trap );
			if( constraint )
			{
				counteff_simple( locnjob, topol, len, effarrforlocalhom );
				calcimportance( locnjob, effarrforlocalhom, aseq, localhomtable );
			}
		}

		if( iterate % 2 == 0 ) 
		{
			lin = 0; ldf = +1;
		}
		else
		{
			lin = locnjob - 2; ldf = -1;
		}	

		if( score_check == 2 )
		{
			effarr1[0] = 1.0;
			effarr2[0] = 1.0;
			length = strlen( bseq[0] );
			for( i=0; i<locnjob-1; i++ )
				for( j=i+1; j<locnjob; j++ )
					intergroup_score( bseq+i, bseq+j, effarr1, effarr2, 1, 1, length, imanoten[i]+j );
		}

		for( l=lin; l < locnjob-1 && l >= 0 ; l+=ldf )
		{
			for( k=0; k<2; k++ ) 
			{
				if( l == locnjob-2 ) k = 1;
#if 1
				fprintf( stderr, "STEP %03d-%03d-%d ", iterate+1, l+1, k );
#else
				fprintf( stderr, "STEP %03d-%03d-%d %s", iterate+1, l+1, k, use_fft?"\n":"\n" );
#endif
				for( i=0; i<locnjob; i++ ) for( j=0; j<locnjob; j++ ) pair[i][j] = 0;

				OneClusterAndTheOther( locnjob, pair, &s1, &s2, topol, l, k );
#if 0
				fprintf( stderr, "STEP%d-%d\n", l, k );
				for( i=0; i<locnjob; i++ ) 
				{
					for( j=0; j<locnjob; j++ ) 
					{
						fprintf( stderr, "%#3d", pair[i][j] );
					}
					fprintf( stderr, "\n" );
				}
#endif
				if( !weight )
				{
					for( i=0; i<locnjob; i++ ) effarr[i] = 1.0;
					if( nkozo )
					{
						for( i=0; i<locnjob; i++ ) 
						{
							if( kozoarivec[i] )
								effarr_kozo[i] = 1.0;
							else
								effarr_kozo[i] = 0.0;
						}
					}
				}
				else if( weight == 2 ) 
				{
					nodeFromABranch( locnjob, branchnode, node, topol, len, l, k );
					node_eff( locnjob, effarr, branchnode );
				}
				else if( weight == 4 )
				{
					weightFromABranch( locnjob, effarr, stopol, topol, l, k );
#if 0
					if( nkozo )
					{
						assignstrweight( locnjob, effarr_kozo, stopol, topol, l, k, kozoarivec, effarr );
					}

#else
					if( nkozo ) // hitomadu single weight.
						for( i=0; i<locnjob; i++ ) 
						{
							if( kozoarivec[i] ) effarr_kozo[i] = effarr[i]; 
							else effarr_kozo[i] = 0.0;
						}
#endif
#if 0
					fprintf( stderr, "\n" );
					fprintf( stderr, "effarr_kozo = \n" );
					for( i=0; i<locnjob; i++ ) fprintf( stderr, "%5.3f ", effarr_kozo[i] );
					fprintf( stderr, "\n" );
					fprintf( stderr, "effarr = \n" );
					for( i=0; i<locnjob; i++ ) fprintf( stderr, "%5.3f ", effarr[i] );
					fprintf( stderr, "\n\n" );
#endif
				}

				for( i=0; i<locnjob; i++ ) strcpy( aseq[i], bseq[i] );
				length = strlen( aseq[0] );

				if( nkozo )
				{
#if 1
					double tmptmptmp;
					tmptmptmp = 0.0;
					clus1 = conjuctionfortbfast_kozo( &tmptmptmp, pair, s1, aseq, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
					for( i=0; i<clus1; i++ ) effarr1_kozo[i] *= 1.0; // 0.5 ga sairyo ?
					tmptmptmp = 0.0;
					clus2 = conjuctionfortbfast_kozo( &tmptmptmp, pair, s2, aseq, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
					for( i=0; i<clus2; i++ ) effarr2_kozo[i] *= 1.0; // 0.5 ga sairyo ?

#if 0
					fprintf( stderr, "\ngroup1 = %s\n", indication1 );
					for( i=0; i<clus1; i++ ) fprintf( stderr, "effarr1_kozo[%d], effarr1[]  = %f, %f\n", i, effarr1_kozo[i], effarr1[i] );
					fprintf( stderr, "\ngroup2 = %s\n", indication2 );
					for( i=0; i<clus2; i++ ) fprintf( stderr, "effarr2_kozo[%d], effarr2[]  = %f, %f\n", i, effarr2_kozo[i], effarr2[i] );
#endif





#else
					clus1 = conjuctionfortbfast_kozo_BUG( pair, s1, aseq, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
					clus2 = conjuctionfortbfast_kozo_BUG( pair, s2, aseq, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
#endif
				}
				else
				{
					clus1 = conjuctionfortbfast( pair, s1, aseq, mseq1, effarr1, effarr, indication1 );
					clus2 = conjuctionfortbfast( pair, s2, aseq, mseq2, effarr2, effarr, indication2 );
				}



		        if( rnakozo && rnaprediction == 'm' )
		        {       
					makegrouprnait( grouprna1, singlerna, pair, s1 );
					makegrouprnait( grouprna2, singlerna, pair, s2 );
		        }       

				if( score_check == 2 )
				{
					if( constraint )
					{
//						msshrinklocalhom( pair, s1, s2, localhomtable, localhomshrink );
						shrinklocalhom( pair, s1, s2, localhomtable, localhomshrink );
						oimpmatch = 0.0;
						if( use_fft )
						{
							if( alg == 'Q' )
							{
								part_imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
								if( rnakozo ) part_imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
								for(  i=length-1; i>=0; i-- ) oimpmatch += part_imp_match_out_scQ( i, i );
							}
							else
							{
								part_imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
								if( rnakozo ) part_imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
								for(  i=length-1; i>=0; i-- ) oimpmatch += part_imp_match_out_sc( i, i );
							}
						}
						else
						{
							if( alg == 'Q' )
							{
								imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
								if( rnakozo ) imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
								for(  i=length-1; i>=0; i-- ) oimpmatch += imp_match_out_scQ( i, i );
							}
							else
							{
								imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
								fprintf( stderr, "not supported\n" );
								exit( 1 );
							}
						}
//						fprintf( stderr, "### oimpmatch = %f\n", oimpmatch );
					}
					else
					{
						oimpmatch = 0.0;
					}
					tmpdouble = 0.0;
					iu=0; 
					for( i=s1; i<locnjob; i++ ) 
					{
						if( !pair[s1][i] ) continue;
						ju=0;
						for( j=s2; j<locnjob; j++ )
						{
							if( !pair[s2][j] ) continue;
//							fprintf( stderr, "i = %d, j = %d, effarr1=%f, effarr2=%f\n", i, j, effarr1[iu], effarr2[ju] );
							tmpdouble += effarr1[iu] * effarr2[ju] * imanoten[MIN(i,j)][MAX(i,j)];
							ju++;
						}
						iu++;
					}
					mscore = oimpmatch + tmpdouble;
				}
				else if( score_check )
				{
#if 1
					if( RNAscoremtx == 'r' )
						intergroup_score_gapnomi( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame
					else
						intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame
#else
					intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble ); // gappick mae denaito dame
#endif

					if( constraint )
					{
						shrinklocalhom( pair, s1, s2, localhomtable, localhomshrink );
	//					weightimportance4( clus1, clus2,  effarr1, effarr2, localhomshrink ); // >>>
						oimpmatch = 0.0;
						if( use_fft )
						{
							if( alg == 'Q' )
							{
								part_imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
								if( rnakozo ) part_imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
								for(  i=length-1; i>=0; i-- )
								{
									oimpmatch += part_imp_match_out_scQ( i, i );
//									fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
								}
							}
							else
							{
								part_imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
								if( rnakozo ) part_imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
								for(  i=length-1; i>=0; i-- )
								{
									oimpmatch += part_imp_match_out_sc( i, i );
	//								fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
								}
							}
	//						fprintf( stderr, "otmpmatch = %f\n", oimpmatch );
						}
						else
						{
							if( alg == 'Q' )
							{
								imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
								if( rnakozo ) imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );

								for(  i=length-1; i>=0; i-- )
								{
									oimpmatch += imp_match_out_scQ( i, i );
//									fprintf( stderr, "#### i=%d, initial impmatch = %f\n", i, oimpmatch );
								}
							}
							else
							{
								imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );

								fprintf( stderr, "not supported\n" );
								exit( 1 );

								for(  i=length-1; i>=0; i-- )
								{
									oimpmatch += imp_match_out_sc( i, i );
	//								fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
								}
							}
	//						fprintf( stderr, "otmpmatch = %f\n", oimpmatch );
						}
	//					fprintf( stderr, "#### initial impmatch = %f\n", oimpmatch );
					}
					else
					{
						oimpmatch = 0.0;
					}


//					fprintf( stderr, "#### tmpdouble = %f\n", tmpdouble );
					mscore = (double)oimpmatch + tmpdouble;
				}
				else
				{
//					fprintf( stderr, "score_check = %d\n", score_check );
/* atode kousokuka */
					intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
					mscore = tmpdouble;
/* atode kousokuka */

					if( constraint )
					{
						oimpmatch = 0.0;
						shrinklocalhom( pair, s1, s2, localhomtable, localhomshrink );
						if( use_fft )
						{
							if( alg == 'Q' )
							{
								part_imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
								if( rnakozo ) part_imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
							}
							else
							{
								part_imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
								if( rnakozo ) part_imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
							}
						}
						else
						{
							if( alg == 'Q' )
							{
								imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
								if( rnakozo ) imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, gapmap1, gapmap2, NULL );
							}
							else
							{
								imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
								fprintf( stderr, "Not supported\n" );
								exit( 1 );
							}
						}
					}
				}

//				oimpmatch = 0.0;
				if( constraint )
				{
#if 0 // iranai
					if( alg == 'Q' )
					{
						imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
						for(  i=length-1; i>=0; i-- )
						{
							oimpmatch += imp_match_out_scQ( i, i );
//							fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
						}
					}
					else
					{
						imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
						for(  i=length-1; i>=0; i-- )
						{
							oimpmatch += imp_match_out_sc( i, i );
//							fprintf( stderr, "#### i=%d, initial impmatch = %f seq1 = %c, seq2 = %c\n", i, oimpmatch, mseq1[0][i], mseq2[0][i] );
						}
					}
#endif
				}
#if 0
				if( alg == 'H' )
					naivescore0 = naivepairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty ) + oimpmatch;
				else if( alg == 'Q' )
					naivescore0 = naiveQpairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty ) + oimpmatch;
				else if( alg == 'R' )
					naivescore0 = naiveRpairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty ) + oimpmatch;
#endif

//				if( rnakozo ) foldalignedrna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, rnapairboth );

//				if( !use_fft && !rnakozo )
				if( !use_fft )
				{
					commongappick_record( clus1, mseq1, gapmap1 );
					commongappick_record( clus2, mseq2, gapmap2 );
				}

#if 0
				fprintf( stderr, "##### mscore = %f\n", mscore );
#endif
	
#if DEBUG
				if( !devide )
				{
		       		fprintf( trap_g, "\nSTEP%d-%d-%d\n", iterate+1, l+1, k );
       		    	fprintf( trap_g, "group1 = %s\n", indication1 );
   	   		    	fprintf( trap_g, "group2 = %s\n", indication2 );
					fflush( trap_g );
				}
	
#endif
#if 0
				printf( "STEP %d-%d-%d\n", iterate, l, k );
				for( i=0; i<clus2; i++ ) printf( "%f ", effarr2[i] );
				printf( "\n" );
#endif
				if( constraint == 2 )
				{
					if( use_fft )
					{
//						if( alg == 'Q' )
//							part_imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 );
//						else
//							part_imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
						Falign_localhom( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, localhomshrink, &impmatch, gapmap1, gapmap2 );
//						fprintf( stderr, "##### impmatch = %f\n", impmatch );
					}
					else
					{
						if( alg == 'Q' )
						{
							float wm;
//							imp_match_init_strictQ( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 ); // Ichijiteki, gapmap  wo tsukuttakara iranai.
//							if( rnakozo ) imp_rnaQ( clus1, clus2, mseq1, mseq2, effarr1, effarr2, gapmap1, gapmap2, rnapairboth );

							wm = Q__align_gapmap( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, localhomshrink, &impmatch, NULL, NULL, NULL, NULL, gapmap1, gapmap2 );
							fprintf( stderr, "wm = %f\n", wm );
#if 0
							fprintf( stderr, "##### impmatch = %f->%f\n", oimpmatch, impmatch );
							naivescore1 = naiveQpairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty ) + impmatch;
							fprintf( stderr, "##### naivscore1 = %f\n", naivescore1 );

							if( naivescore1 > naivescore0 )
								fprintf( stderr, "%d-%d, ns: UP!\n", clus1, clus2 );
							else if( naivescore1 < naivescore0 )
								fprintf( stderr, "%d-%d, ns: DOWN! %f->%f\n", clus1, clus2, naivescore0, naivescore1 );
							else
								fprintf( stderr, "%d-%d, ns: IDENTICAL\n", clus1, clus2 );
#if 0 // chuui
							if( abs( wm - naivescore1 ) > 100 )
							{
//								fprintf( stderr, "WARNING, wm=%f but naivescore=%f\n", wm, naivescore1 );
//								rewind( stdout );
//								for( i=0; i<clus1; i++ ) printf( ">\n%s\n", mseq1[i] );
//								for( i=0; i<clus2; i++ ) printf( ">\n%s\n", mseq2[i] );
//								exit( 1 );
							}
#endif
#endif
						}
						else if( alg == 'R' )
						{
							float wm;
							imp_match_init_strictR( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 ); // Ichijiteki, gapmap ha mada
							wm = R__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, localhomshrink, &impmatch, NULL, NULL, NULL, NULL );
//							fprintf( stderr, "##### impmatch = %f->%f\n", oimpmatch, impmatch );
							naivescore1 = naiveRpairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty ) + impmatch;
//							fprintf( stderr, "##### naivscore1 = %f\n", naivescore1 );

							if( naivescore1 > naivescore0 )
								fprintf( stderr, "%d-%d, ns: %f->%f UP!\n", clus1, clus2, naivescore0, naivescore1 );
							else if( naivescore1 < naivescore0 )
								fprintf( stderr, "%d-%d, ns: DOWN! %f->%f\n", clus1, clus2, naivescore0, naivescore1 );
							else
								fprintf( stderr, "%d-%d, ns: IDENTICAL\n", clus1, clus2 );
#if 0 // chuui
							if( abs( wm - naivescore1 ) > 100 )
							{
//								fprintf( stderr, "WARNING, wm=%f but naivescore=%f\n", wm, naivescore1 );
								rewind( stdout );
								for( i=0; i<clus1; i++ ) printf( ">\n%s\n", mseq1[i] );
								for( i=0; i<clus2; i++ ) printf( ">\n%s\n", mseq2[i] );
								exit( 1 );
							}
#endif
						}
						else if( alg == 'H' )
						{
							float wm;
							imp_match_init_strictH( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, localhomshrink, 1 ); // Ichijiteki, gapmap ha mada
							wm = H__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, localhomshrink, &impmatch, NULL, NULL, NULL, NULL );
							fprintf( stderr, "##### impmatch = %f->%f\n", oimpmatch, impmatch );
							naivescore1 = naivepairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty ) + impmatch;
							fprintf( stderr, "##### naivscore1 = %f\n", naivescore1 );

							if( naivescore1 > naivescore0 )
								fprintf( stderr, "%d-%d, ns: UP!\n", clus1, clus2 );
							else if( naivescore1 < naivescore0 )
								fprintf( stderr, "%d-%d, ns: DOWN! %f->%f\n", clus1, clus2, naivescore0, naivescore1 );
							else
								fprintf( stderr, "%d-%d, ns: IDENTICAL\n", clus1, clus2 );
#if 0 // chuui
							if( abs( wm - naivescore1 ) > 100 )
							{
//								fprintf( stderr, "WARNING, totalwm=%f but naivescore=%f\n", totalwm, naivescore1 );
//								for( i=0; i<clus1; i++ ) printf( ">\n%s\n", mseq1[i] );
//								for( i=0; i<clus2; i++ ) printf( ">\n%s\n", mseq2[i] );
//								exit( 1 );
							}
#endif
						}
						else
						{
//							imp_match_init_strict( NULL, clus1, clus2, length, length, mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, 1 );
							A__align_gapmap( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, localhomshrink, &impmatch, gapmap1, gapmap2 );
//							fprintf( stderr, "##### impmatch = %f\n", impmatch );
						}
					}
				}
				else if( use_fft )
				{
					float totalwm;
					totalwm = Falign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, &intdum );

//					fprintf( stderr, "totalwm = %f\n", totalwm );
#if 0
					if( alg == 'Q' )
					{
						fprintf( stderr, "totalwm = %f\n", totalwm );
						naivescore1 = naiveQpairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty );
	
						if( naivescore1 > naivescore0 )
							fprintf( stderr, "%d-%d, ns: UP!\n", clus1, clus2 );
						else if( naivescore1 < naivescore0 )
							fprintf( stderr, "%d-%d, ns: DOWN! %f->%f\n", clus1, clus2, naivescore0, naivescore1 );
						else
							fprintf( stderr, "%d-%d, ns: IDENTICAL\n", clus1, clus2 );
#if 1 // chuui
						if( totalwm != 0.0 && abs( totalwm - naivescore1 ) > 100 )
						{
//							fprintf( stderr, "WARNING, totalwm=%f but naivescore=%f\n", totalwm, naivescore1 );
//							for( i=0; i<clus1; i++ ) printf( ">\n%s\n", mseq1[i] );
//							for( i=0; i<clus2; i++ ) printf( ">\n%s\n", mseq2[i] );
//							exit( 1 );
						}
#endif
					}
#endif
					if( alg == 'R' )
					{
						fprintf( stderr, "totalwm = %f\n", totalwm );
						naivescore1 = naiveRpairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty );
	
						if( naivescore1 > naivescore0 )
							fprintf( stderr, "%d-%d, ns: UP!\n", clus1, clus2 );
						else if( naivescore1 < naivescore0 )
							fprintf( stderr, "%d-%d, ns: DOWN! %f->%f\n", clus1, clus2, naivescore0, naivescore1 );
						else
							fprintf( stderr, "%d-%d, ns: IDENTICAL\n", clus1, clus2 );
#if 1 // chuui
						if( totalwm != 0.0 && abs( totalwm - naivescore1 ) > 100 )
						{
//							fprintf( stderr, "WARNING, totalwm=%f but naivescore=%f\n", totalwm, naivescore1 );
//							for( i=0; i<clus1; i++ ) printf( ">\n%s\n", mseq1[i] );
//							for( i=0; i<clus2; i++ ) printf( ">\n%s\n", mseq2[i] );
//							exit( 1 );
						}
					}
#endif
				}
				else
				{
					if( alg == 'M' )
					{
						MSalignmm( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, NULL, NULL, NULL, NULL );
					}
					else if( alg == 'A' )
					{
						A__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, NULL, &impmatch, NULL, NULL, NULL, NULL );
					}
					else if( alg == 'Q' )
					{
						float wm;
						wm = Q__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, NULL, &impmatch, NULL, NULL, NULL, NULL );
						fprintf( stderr, "wm = %f\n", wm );
						fprintf( stderr, "impmatch = %f\n", impmatch );
						naivescore1 = naiveQpairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty );

						if( naivescore1 > naivescore0 )
							fprintf( stderr, "%d-%d, ns: UP!\n", clus1, clus2 );
						else if( naivescore1 < naivescore0 )
							fprintf( stderr, "%d-%d, ns: DOWN!\n", clus1, clus2 );
						else
							fprintf( stderr, "%d-%d, ns: IDENTICAL\n", clus1, clus2 );
#if 1 // chuui
						if( abs( wm - naivescore1 ) > 100 )
						{
//							fprintf( stderr, "WARNING, wm=%f but naivescore=%f\n", wm, naivescore1 );
//							rewind( stderr );
//							rewind( stdout );
//							for( i=0; i<clus1; i++ ) printf( ">\n%s\n", mseq1[i] );
//							for( i=0; i<clus2; i++ ) printf( ">\n%s\n", mseq2[i] );
//							exit( 1 );
						}
#endif
					}
					else if( alg == 'R' )
					{
						float wm;
						wm = R__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, NULL, &impmatch, NULL, NULL, NULL, NULL );
						naivescore1 = naiveRpairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty );

						if( naivescore1 > naivescore0 )
							fprintf( stderr, "%d-%d, ns: UP!\n", clus1, clus2 );
						else if( naivescore1 < naivescore0 )
							fprintf( stderr, "%d-%d, ns: DOWN! %f->%f\n", clus1, clus2, naivescore0, naivescore1 );
						else
							fprintf( stderr, "%d-%d, ns: IDENTICAL\n", clus1, clus2 );
#if 1 // chuui
						if( abs( wm - naivescore1 ) > 100 )
						{
//							fprintf( stderr, "WARNING, wm=%f but naivescore=%f\n", wm, naivescore1 );
//							rewind( stderr );
//							rewind( stdout );
//							for( i=0; i<clus1; i++ ) printf( ">\n%s\n", mseq1[i] );
//							for( i=0; i<clus2; i++ ) printf( ">\n%s\n", mseq2[i] );
//							exit( 1 );
						}
#endif
					}
					else if( alg == 'H' )
					{
						float wm;
						wm = H__align( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen, NULL, &impmatch, NULL, NULL, NULL, NULL );
						naivescore1 = naivepairscore( clus1, clus2, mseq1, mseq2, effarr1, effarr2, penalty );

						if( naivescore1 > naivescore0 )
							fprintf( stderr, "%d-%d, ns: UP!\n", clus1, clus2 );
						else if( naivescore1 < naivescore0 )
						{
							fprintf( stderr, "%d-%d, ns: DOWN!\n", clus1, clus2 );
						}
						else
							fprintf( stderr, "%d-%d, ns: IDENTICAL\n", clus1, clus2 );

#if 0 // chuui
						if( abs( wm - naivescore1 ) > 100 )
						{
							rewind( stdout );
							for( i=0; i<clus1; i++ ) printf( ">\n%s\n", mseq1[i] );
							for( i=0; i<clus2; i++ ) printf( ">\n%s\n", mseq2[i] );
							exit( 1 );
						}
#endif
					}
					else if( alg == 'a' ) 
					{
						Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, alloclen );
					}
					else ErrorExit( "Sorry!" );
				}
//				fprintf( stderr, "## impmatch = %f\n", impmatch );
							
					if( checkC )
					{
						extern double DSPscore();
						extern double SSPscore();
						static double cur;
						static double pre;
	
/*
						pre = SSPscore( locnjob, bseq );
						cur = SSPscore( locnjob, aseq );
*/
						pre = DSPscore( locnjob, bseq );
						cur = DSPscore( locnjob, aseq );
	
						fprintf( stderr, "Previous Sscore = %f\n", pre );
						fprintf( stderr, "Currnet  Sscore = %f\n\n", cur );
					}
					
//				fprintf( stderr, "## impmatch = %f\n", impmatch );
				identity = !strcmp( aseq[s1], bseq[s1] );
				identity *= !strcmp( aseq[s2], bseq[s2] );


/* Bug?  : idnetitcal but score change when scoreing mtx != JTT  */

				length = strlen( mseq1[0] );
	
				if( identity )
				{
					tscore = mscore;
					if( !devide ) fprintf( trap_g, "tscore =  %f   identical.\n", tscore );
					fprintf( stderr, " identical." );
					converged++;
				}
				else
				{
					if( score_check )
					{
						if( constraint == 2 )
						{
#if 1
							if( RNAscoremtx == 'r' )
								intergroup_score_gapnomi( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
							else
								intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
#else
							intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
#endif

							tscore = impmatch + tmpdouble;

//							fprintf( stderr, "tmpdouble=%f, impmatch = %f -> %f, tscore = %f\n", tmpdouble, oimpmatch, impmatch, tscore );
						}
						else
						{
							intergroup_score( mseq1, mseq2, effarr1, effarr2, clus1, clus2, length, &tmpdouble );
							tscore = tmpdouble;
						}
//						fprintf( stderr, "#######ii=%d, iterate=%d score = %f -> %f \n", ii, iterate , mscore, tscore );
	#if 0
						for( i=0; i<1; i++ )
							fprintf( stderr, "%s\n", mseq1[i] );
						fprintf( stderr, "+++++++\n" );
						for( i=0; i<1; i++ )
							fprintf( stderr, "%s\n", mseq2[i] );
	#endif
	
					}
					else
					{
						tscore = mscore + 1.0;
//						tscore = 0.0;
//						fprintf( stderr, "in line 705, tscore=%f\n", tscore );
//						for( i=0; i<length; i++ )
//							tscore = tscore + (double)mseq1[0][i];
//						mscore = tscore - 1.0;
					}

					if( isnan( mscore ) )
					{
						fprintf( stderr, "\n\nmscore became NaN\n" );
						exit( 1 );
					}
					if( isnan( tscore ) )
					{
						fprintf( stderr, "\n\ntscore became NaN\n" );
						exit( 1 );
					}



//					fprintf( stderr, "@@@@@ mscore,tscore = %f,%f\n", mscore, tscore );

					if( tscore > mscore - cut/100.0*mscore ) 
					{
						writePre( locnjob, name, nlen, aseq, 0 );
						for( i=0; i<locnjob; i++ ) strcpy( bseq[i], aseq[i] );
						if( score_check == 2 )
						{
							effarr1[0] = 1.0;
							effarr2[0] = 1.0;
							for( i=0; i<locnjob-1; i++ )
								for( j=i+1; j<locnjob; j++ )
									intergroup_score( bseq+i, bseq+j, effarr1, effarr2, 1, 1, length, imanoten[i]+j );
						}
	
#if 0
						fprintf( stderr, "tscore =  %f   mscore = %f  accepted.\n", tscore, mscore );
#endif
						fprintf( stderr, " accepted." );
						converged = 0;
	
					}
					else 
					{
#if 0
						fprintf( stderr, "tscore =  %f   mscore = %f  rejected.\n", tscore, mscore );
#endif
						fprintf( stderr, " rejected." );
						tscore = mscore;
						converged++;
					}
				}
				fprintf( stderr, "\r" );


				history[iterate][l][k] = (float)tscore;

//				fprintf( stderr, "tscore = %f\n", tscore );
	
				if( converged >= locnjob * 2 )
				{
					fprintf( trap_g, "Converged.\n\n" );
					fprintf( stderr, "\nConverged.\n\n" );
					return( 0 );
				}
				if( iterate >= 1 )
				{
	/*   oscillation check    */
					oscillating = 0;
					for( ii=iterate-2; ii>=0; ii-=2 ) 
					{
						if( (float)tscore == history[ii][l][k] )
						{
							oscillating = 1;
							break;
						}
					}
					if( ( oscillating && !cooling ) || ( oscillating && cut < 0.001 && cooling ) )
					{
						fprintf( trap_g, "Oscillating.\n" );
						fprintf( stderr, "\nOscillating.\n\n" );
#if 1 /* hujuubun */
						return( -1 );
#endif
					}
				}      /* if( iterate ) */
			}          /* for( k ) */
		}              /* for( l ) */
	}                  /* for( iterate ) */
	return( 2 );
}                  	   /* int Tree... */
