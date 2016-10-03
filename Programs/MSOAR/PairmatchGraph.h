/*********************************************************************/
/*Object: PairmatchGraph 					     */
/*Derived from class "Graph".					     */
/*Constructing a pairmatch graph for a pair of genome S and T to     */
/*approximate Minimum Common Partition problem. As for the definition*/
/*of the pairmatch graph, please refer to the paper		     */
/*********************************************************************/

#ifndef PAIRMATCHGRAPH_H
#define PAIRMATCHGRAPH_H

#include <list>
#include <string>
#include <vector>
#include "Graph.h"
#include "Genome.h"
#include "VertexPair.h"

using namespace std;

class PairmatchGraph: public Graph
{
private:
  
public:
  typedef struct {
     int vertexPairS;
     int vertexPairT;
     int orient;
  }VertexMatchPair;
  
  //a vector of all the adjacent vertex pairs in genome S with at least a match in another genome
  vector<VertexPair>vertexPairS; 
  //a vector of all the adjacent vertex pairs in genome T with at least a match in another genome
  vector<VertexPair>vertexPairT; 
  //a vector of all the vertex match pairs between genome S and T 
  vector<VertexMatchPair>vmPairs;
  
  Genome S;
  Genome T;

  PairmatchGraph( );
  
  void ConstructGraph( Genome SS, Genome TT );
  void CreatePairmatches(  );
  void CreateVertexPairs(  );
  void getVertexPairT( int k, int i );
  void CreateVertexMatchPairs( );
  bool ConflictingEdges( VertexMatchPair p, VertexMatchPair q );
  bool ConflictingEdgesS( int u, int v, int w, int x );
  void Print( );
 
  
  PairmatchGraph& operator = ( const PairmatchGraph &g ){
    vertices = g.vertices;
    return  *this;
  }

  
};// end of class CompleteGraph


#endif

//Local Variables:
//mode: c++
//End:
