/*********************************************************************/
/*Object: CompleteGraph 					     */
/*Derived from class "Graph".					     */
/*Each gene is represented by two vertices and there is a cap verte  */
/*at each end of the Chromosome. The blastp hit is represented by    */
/*cross-genome edges. And there is an edge between two vertices of   */
/*adjacent genes.                                                    */
/*This Object could read in a input file in certain format and       */
/*transform it to a graph in adjacent-list structure                 */
/*********************************************************************/

#ifndef COMPLETEGRAPH_H
#define COMPLETEGRAPH_H

#include <list>
#include <string>
#include <vector>
#include "Graph.h"
#include "Genome.h"
#include "VertexPair.h"
using namespace std;
#define STEPCUTOFF 200

class CompleteGraph: public Graph
{
private:
  
public:
   int endgeneS; 	//id of the last gene vertex of genome S 
   int endchromS;	//id of the last chromosome vertex of genome S
   int endgeneT;	//id of the last gene vertex of genome T
   int endchromT;	//id of the last chromosome vertex of genome T
   CompleteGraph( );
   void ConstructGraph( Genome S, Genome T );
   void ClearGraph( );	
   bool BFSShortest( int i );
   bool ShortestCycle( int i );
   bool isidlevertex( Vertices V, int i );
   void remove_idle_vertices( Vertices& V, int i );
   bool ShortestPath( int i );
   void remove_incompatible_edge(Vertices& V, vertex i, vertex j );
   int ApproximateDistance( );
   int ApproximateDistance( int i, int j );
   void deletePair( int i, int j );
   CompleteGraph& operator = ( const CompleteGraph &g ){
     endgeneS = g.endgeneS;
     endchromS = g.endchromS;
     endgeneT = g.endgeneT;
     endchromT = g.endchromT;
     vertices = g.vertices;
     return  *this;
   }  
  
};// end of class CompleteGraph


#endif

//Local Variables:
//mode: c++
//End
