/*********************************************************************/
/*Object: RelatedGenomes 					     */
/*This is the core object in MSOAR, which contain two related genomes*/
/*S and T and all useful informations. ("related genomes" simply     */
/*refer to two genomes under compairison)                            */
/*********************************************************************/
#ifndef _RELATEDGENOMES_H_
#define _RELATEDGENOMES_H_

#include <vector>
#include <string>
#include <fstream>
#include "Genome.h"
#include "Gene.h"
#include "CompleteGraph.h"
#include "PairmatchGraph.h"

using namespace std;

#define WEIGHT_PARALOG 2
#define EVALUE_DIFF_CUTOFF 0.2

class RelatedGenomes
{
 private:
 

 public:
   typedef struct { int A;int B; } EndPair;
   typedef vector<EndPair> GeneFamily;
   Genome S;
   Genome T;
   vector<int>matchesS; //matchesS[i] = j means the i-th gene in genome S is matched to th j-th gene in T
   vector<int>matchesT; //matchesT[j] = i means the j-th gene in genome T is matched to th i-th gene in S
   vector<int>modified; //store the ids of matches modified by newest actions.
   CompleteGraph G;     //Complete graph for S and T
   PairmatchGraph P;    //Pairmatch Graph for S and T
   vector<GeneFamily>GFs;	//Gene Families after assigning main orthologs pairs
   vector<EndPair>Paralogs;	//Paralogs vector
   
   RelatedGenomes( );
   RelatedGenomes( vector<int>A, vector<int>B, Genome s, Genome t );
   RelatedGenomes( string infile1, string infile2 );
   void RemoveFP( );
   void SuboptimalRules( );
   void ApproximateMCP( );
   bool ApproximateMGD( );
   bool ApproximateMGD(int start );
   void UpdateMatches_MGD( );
   void UpdateMatches_MCP( );
   void UpdateCompleteGraph( );
   void SuboptimalRule1( );
   void SuboptimalRule2( );
   void SuboptimalRule3( );
   void selectBlastHit( int id, int obj );
   bool tobeMatchedT( int i );
   void RandomMatch( );
   void RandomMatch( int i );
   void print (  );
   void printresult( string outfile );
   bool isSame( const char* s1, const char* s2 );
   void CheckBlastHits( );
   void Update( );
   void GeneFamiliesConstruction( vector<Gene>A, vector<Gene>B);
   void ParalogsDetection( );
   void ParalogsPrint( );
   void OrthologCandidate( );
};

#endif
