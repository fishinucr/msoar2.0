/*********************************************************************/
/*Object: VertexPair 					     	    */
/*This object is used in Parimatch Graph			    */
/********************************************************************/

#ifndef VERTEXPAIR_H
#define VERTEXPAIR_H

#include <list>
#include <string>
#include <vector>


using namespace std;

class VertexPair 
{
private:
  
public:
  int index; //The beginning position of a neighboring pair in this genome
  vector<int>MATCHES; //all the possible matches (VertexPairs) in another genome
  int selected; // id of the VertexPair in vector "MATCHES" that match to this VetexPair

  VertexPair( ){
    index = -1;
    selected = -1;
  }

  VertexPair( int i ){
    index =i;
    selected = -1;
  }
  
  bool isSelected( ){
    if( selected >= 0 ) return true;
    return false;
  }  
  
 

  bool isExistMatch( int v ){
    for( int i = 0; i < (int)MATCHES.size( ); i++ ){
      if( MATCHES[i] == v )
	return true;
    }
    return false;
  }
  
  void print( ){
    cout <<"index = "<<index<<": ";
    for( int i = 0 ; i < (int)MATCHES.size( ); i ++ )
      cout << MATCHES[i] <<" ";
    cout<<endl;
  }
  
  
};

#endif

