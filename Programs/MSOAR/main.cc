/************************************************/
/* Title: M-SOAR                                */
/* Description: System for orthologs assignment */
/* Author: Zheng Fu                             */
/* Copyright: Copyright (c) 2005                */ 
/************************************************/ 

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "RelatedGenomes.h"
#include "Genome.h"
#include "Gene.h"
#include "BlastHit.h"
#include "Graph.h"
#include "CompleteGraph.h"
#include "VertexPair.h"
#include "PairmatchGraph.h"

using namespace std;
#define ROUNDCUTOFF 1000

int main(int argc, char *argv[])
{
   cout.flush( );
  if( argc != 5 ){
     /* 
      Inputfile1: the blastp results for Genome1 vs. Genome2 including gene location, e.g., hm_msoar
      Inputfile2: the blastp results for Genome2 vs. Genome1 including gene location, e.g., mh_msoar
      Resultfile1: the ortholog pairs assigned before noise gene detection 
      Resultfile2: the ortholog pairs assigned after noise gene detection
     */
     cout << "Usage: "<< argv[0] << " Inputfile1" << " Inputfile2" << " Resultfile1" <<" Resultfile2"<<endl;
     return 0;
  }

  cout << "Load input data ......" << endl;
  RelatedGenomes relatedgenomes( argv[1], argv[2] );
	cout << "Here5 " << endl;

  RelatedGenomes relatedgenomes_copy( argv[1], argv[2] );

	cout << "Here2 " << endl;
  relatedgenomes.RemoveFP( );

	cout << "Here3 " << endl;
  relatedgenomes_copy.RemoveFP( );

	cout << "Here4 " << endl;
  relatedgenomes.OrthologCandidate( );  
  cout << "Approximate MCP ......" << endl;
  relatedgenomes.ApproximateMCP( ); // Minimum Common Partitions

  cout << "Apply three sub-optimal rules ......" << endl;
  relatedgenomes.SuboptimalRules( ); //Apply three sub-optimal rules

  cout << "Maximum Graph Decomposition"<<endl;
  relatedgenomes.G.ConstructGraph( relatedgenomes.S, relatedgenomes.T );
  // relatedgenomes.print(  ); //print out the results
  
  int start = 0;
  while(1){
         
     if( relatedgenomes.ApproximateMGD( start++  ) == true) break ; //Maximum Cycle Decomposition 
     if( start >= (int)relatedgenomes.T.genes.size( )  )
        break;     
  }

  cout <<"Random Match ......" << endl;
  relatedgenomes.RandomMatch( );

  cout << "Print out the results before noise gene pairs detection......" <<endl;
  relatedgenomes.printresult( argv[3] ); //print out the results
  
  // cout<<"relatedgenomes.print( )"<<endl;
  // relatedgenomes.print( );
  
  cout << "Noise gene pairs detetion......"<<endl;
  cout << "\t Related genomes update......"<<endl;
  relatedgenomes_copy.matchesS = relatedgenomes.matchesS ;
  relatedgenomes_copy.matchesT = relatedgenomes.matchesT;
  // cout<<"\t relatedgenomes_copy.print( )"<<endl;
  // relatedgenomes_copy.print( );
  relatedgenomes_copy.Update( );
  cout<<"\t relatedgenomes_copy.print( ) after update"<<endl;
  relatedgenomes_copy.print( );
  cout << "\t Detect paralogs......"<<endl;
  relatedgenomes_copy.ParalogsDetection( );  
  
  cout << "Print out the results after noise gene pairs detection......" <<endl;
  relatedgenomes_copy.ParalogsPrint( );
  relatedgenomes_copy.printresult( argv[4] ); //print out the results

}
