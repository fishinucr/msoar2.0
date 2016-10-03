/*********************************************************************/
/*Object: BlastHit for each gene				     */
/*Each gene has a vector of BlastHits, which saves all the qulified  */
/*Blast Hits of this gene. Each BlastHit object keeps the information*/
/*of Object gene index, E-value and BitScore.			     */
/*********************************************************************/
 
#ifndef _BLASTHIT_H_
#define _BLASTHIT_H_

#include <vector>
#include <iostream>

using namespace std;

class BlastHit 
{
 private:

 public:
  int ObjectIndex;//Object Gene Index
  double Evalue; 
  double BitScore;
  BlastHit( ) { ObjectIndex = -1; Evalue = 0.0; BitScore = 0.0; }
  BlastHit( int o, double e, double b ) { ObjectIndex = o; Evalue = e; BitScore = b; }
 
  void print( ){
    cout << ObjectIndex << "\t" << Evalue << "\t" << BitScore << "\n";
  }
};

#endif
