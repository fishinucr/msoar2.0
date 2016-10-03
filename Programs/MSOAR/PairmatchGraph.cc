#include <list>
#include <string>
#include <vector>
#include "Gene.h"
#include "PairmatchGraph.h"

using namespace std;

PairmatchGraph::PairmatchGraph( ){
}

/*Constructing pairmatch graph for a pair of genome SS and TT*/
void PairmatchGraph::ConstructGraph( Genome SS, Genome TT ){
  S = SS;
  T = TT;
  vertices.clear( );
  vertexPairS.clear( );
  vertexPairT.clear( );
  vmPairs.clear( );
  cout << "CreateVertexPairs"<<endl;
  CreateVertexPairs( );
  cout << "CreateVertexMatchPairs"<<endl;
  CreateVertexMatchPairs( );
  
  cout<< "insert vertices into the pairmatch graph "<<endl;
  for( int i = 0; i < (int)vmPairs.size( ); i++ )
    insert_vertex( i );
  
  cout<<"insert edges into the pairmatch graph"<<endl;
  for( int i = 0; i < (int)vmPairs.size( ) - 1; i++ ){
    for( int j = i + 1; j < (int)vmPairs.size( ); j++ ){
      if( ConflictingEdges( vmPairs[i], vmPairs[j] ) ) //insert edges between confliting vertex match pairs
	insert_edge( i, j );
    }
  }
  // Print( );
}

/*for every adjacent gene pair (i, i+1) in genome S, search for the gene pairs */
/*in genome t (j,j+1) matchs (i,i+1) or its reversal. If such gene pair (j,j+1)*/
/*exists, push (i,i+1) to "vertexPairS" and push(j,j+1) to "vertexPairT". And  */
/* keep the match infromation between this two vertex pairs.                   */
void PairmatchGraph::CreateVertexPairs( ){
  for( int i = 0; i < (int)S.genes.size( ) - 1; i ++ ){
    //the two genes of a VertexPair have to be on the same chromosome
    if( S.genes[i+1].ChromId != S.genes[i].ChromId)
      continue;
    //check every VertexPair in T that might match to (i,i+1)
    for( int j = 0; j < (int)S.genes[i].BlastHits.size( ); j ++ ){
      int k = S.genes[i].BlastHits[j].ObjectIndex;
      //same direction
      if( S.genes[i].Sign == T.genes[k].Sign ){
	int ii = i+1;
	int kk = k+1;
	//the two genes of a VertexPair have to be on the same chromosome
	if( T.genes[kk].ChromId != T.genes[k].ChromId)
	  continue;
	if( ii >= (int)S.genes.size( ) || kk >= (int)T.genes.size( ))
	  continue;	
	if( !S.genes[ii].isBlastHit( kk ) )
	  continue;
	if( S.genes[ii].Sign != T.genes[kk].Sign )
	  continue;
      }     
      //different direction
      else{
	int ii = i+1;
	int kk = k-1;
	if( ii >= (int)S.genes.size( ) || kk < 0 )
	  continue;
	if( !S.genes[ii].isBlastHit( kk ) )
	  continue;
	if( S.genes[ii].Sign == T.genes[kk].Sign )
	  continue;
      }
      //push (i,i+1) to "vertexPairS" and assign (k,k+1) to be one of its matches
      if( S.genes[i].Sign != T.genes[k].Sign )
	k --;
      if( ((int)vertexPairS.size( ) > 0 && vertexPairS[(int)vertexPairS.size( )-1].index != i ) ||(int)vertexPairS.size( ) == 0 )
	vertexPairS.push_back( VertexPair( i ) );      
      if( !vertexPairS[(int)vertexPairS.size( )-1].isExistMatch(k) )
	vertexPairS[(int)vertexPairS.size( )-1].MATCHES.push_back( k );
	
      //push (k,k+1) to "vertexPairT" and assign (i,i+1) to be one of its matches
      getVertexPairT(k, i);     
     
      
    }
  }
}

/*push (k,k+1) to "vertexPairT" and assign (i,i+1) to be one of its matches*/
void PairmatchGraph::getVertexPairT( int k ,int i ){
   for( int j = 0; j <(int)vertexPairT.size( ); j ++ ){
    if( vertexPairT[j].index == k ){
      if( !vertexPairT[j].isExistMatch(i) )
	vertexPairT[j].MATCHES.push_back( i );
      return;
    }
  
  }
  VertexPair t( k );
  vertexPairT.push_back( t );
  vertexPairT[(int)vertexPairT.size( )-1].MATCHES.push_back( i );
}

/*Based on "vertexPairS" and "vertexPairT", assign value to "vmp"*/
void PairmatchGraph::CreateVertexMatchPairs( ){
   cout << "(int)vertexPairS.size( )"<<(int)vertexPairS.size( )<<endl;
   for( int i = 0; i < (int)vertexPairS.size( ); i++ ){
     // cout <<"i="<<i<<endl;
     for( int j = 0; j < (int)vertexPairS[i].MATCHES.size( ); j++ ){
        // cout<<"j="<<j<<endl;
        VertexMatchPair vmp;
        vmp.vertexPairS = vertexPairS[i].index;
        vmp.vertexPairT = vertexPairS[i].MATCHES[j];
        vmp.orient = -1;
        if(  S.genes[vmp.vertexPairS].isBlastHit( vmp.vertexPairT ) && S.genes[vmp.vertexPairS + 1].isBlastHit( vmp.vertexPairT + 1 ) && S.genes[vmp.vertexPairS].Sign == T.genes[vmp.vertexPairT].Sign && S.genes[vmp.vertexPairS+1].Sign == T.genes[vmp.vertexPairT+1].Sign )
           vmp.orient = 0;
        if(  S.genes[vmp.vertexPairS].isBlastHit( vmp.vertexPairT + 1 ) && S.genes[vmp.vertexPairS + 1].isBlastHit( vmp.vertexPairT ) && S.genes[vmp.vertexPairS].Sign != T.genes[vmp.vertexPairT + 1].Sign && S.genes[vmp.vertexPairS+1].Sign != T.genes[vmp.vertexPairT].Sign ){
           if( vmp.orient == -1 )
              vmp.orient = 1;
           else
              vmp.orient = 2; // both ways
        }          
        vmPairs.push_back( vmp );
    } 
  }
} 

/*Check if two VertexmatchPairs are conflict with each other*/
bool PairmatchGraph::ConflictingEdges( VertexMatchPair p, VertexMatchPair q ){
 
  int u = p.vertexPairS;
  int v = q.vertexPairS;
  if( u == v )
    return true;
  int w = p.vertexPairT;
  int x = q.vertexPairT;
  if( w == x )
    return true;
  
  //dealing with special cases
  if( v == u + 1 )
    return ConflictingEdgesS( u, v, w, x );
  if( v == u - 1 )
    return ConflictingEdgesS( v, u, x, w );
  if( w == x + 1 || w == x - 1 )
    return true;
 
  return false;
}

bool PairmatchGraph::ConflictingEdgesS( int u, int v, int w, int x ){
  int i = u;
  int j = w;
  
  if( j+2 < (int)T.genes.size( ) && i+2 < (int)S.genes.size( ) && S.genes[i].isBlastHit( j ) && S.genes[i].Sign == T.genes[j].Sign && S.genes[i+1].isBlastHit( j+1 ) && S.genes[i+1].Sign == T.genes[j+1].Sign && S.genes[i+2].isBlastHit( j+2 ) && S.genes[i+2].Sign == T.genes[j+2].Sign){
    if( x == j+1 )
      return false;
 
  }
  else if( j+1 < (int)T.genes.size( ) && j-1 >= 0 && i-1 >= 0 && i+1 < (int)S.genes.size( ) && S.genes[i].isBlastHit( j+1 ) && S.genes[i].Sign != T.genes[j+1].Sign && S.genes[i+1].isBlastHit( j ) && S.genes[i+1].Sign != T.genes[j].Sign && S.genes[i+2].isBlastHit( j-1 ) && S.genes[i+2].Sign != T.genes[j-1].Sign){
    if( x == j-1 )
      return false;
  }

  return true;
}

void PairmatchGraph::Print( ){
  //print vertexPairS
  cout << "(vertexPairS)"<<endl;
  for( int i = 0; i< (int)vertexPairS.size( ); i++ )
    vertexPairS[i].print( );
  //print vertexPairT
  cout << endl<<"(vertexPairT)"<<endl;
  for( int i = 0; i< (int)vertexPairT.size( ); i++ )
      vertexPairT[i].print( );
  //print vmPairs
  cout << endl<<"(vmPairs)"<<endl;
  for( int i = 0; i < (int)vmPairs.size( ); i ++ ){
    cout<< vmPairs[i].vertexPairS<<"---";
    cout<< vmPairs[i].vertexPairT<<endl;
  }
  //print graph
  cout << endl<<"(PairmatchGraph)"<<endl;
  print( );
}
