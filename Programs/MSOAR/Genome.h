/*********************************************************************/
/*Object: Genome						     */
/*Each Genome contains a vertor of genes and another vector          */
/*indicating the size of each chromosome			     */
/*********************************************************************/

#ifndef _GENOME_H_
#define _GENOME_H_

#include <vector>
#include <fstream>
#include "Gene.h"


using namespace std;

class Genome
{
 private:

 public:
  vector<Gene>genes;
  vector<int>chromsizes;

  Genome( ){ };
  
  /*calculate the chromosome sizes*/
  void getChromSizes( ){
    vector<int>tmp(40,0);
    for( int i = 0; i < (int)genes.size( ); i ++ )
      tmp[genes[i].ChromId - 1]++;
    
    for( int i = 0; i < (int)tmp.size( ); i ++ ){
      if( tmp[i] != 0)
         chromsizes.push_back( tmp[i] );
    }
  }


   const Genome& operator = ( const Genome& g ){
      genes = g.genes;
      chromsizes = g.chromsizes;
      return *this;
   }
 

  void print( ) {
     for( int i = 0; i < (int)chromsizes.size( ); i ++ )
        cout << chromsizes[i] <<endl;
    cout<<"**************"<<endl;
    for( int i = 0; i < (int)genes.size( ); i ++ ){
       cout <<"id=" <<i<<":  "<<endl;
       genes[i].print( );
    }
  }
  
  

};

#endif
