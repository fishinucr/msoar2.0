/*********************************************************************/
/*Object: Graph							     */
/*This is a base class for classes CompleteGrah and PairmatchGrah    */
/*********************************************************************/

#ifndef GRAPH_H
#define GRAPH_H

#include <list>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

class Graph
{
public:
   
   typedef int vertex;
   typedef int edge;
 
   //vertex structure
   struct Vertex
   {
     //for special use
      string vtype; 

      //value
      int value;
      
      //direction
      int direction; //0: forward direction; 1: backword direction
     
      //storing edge number
      list<edge> adj_list;


      //construct a vertex with empty adjacent list
      Vertex( string _t, int _v, int _d ): vtype(_t), value(_v), direction(_d) {}

      //insert a edge to the adjacent list and keep it in order
      bool insert_edge( edge e){
        list<int>::iterator it = adj_list.begin( );
        for(; it != adj_list.end( ); ++it ){
           if( *it == e )
	    return false;
           if( *it >= e ){
	      adj_list.insert( it, e );
	      return true;
           } 
	}
        adj_list.push_back( e );
        return true;
      }
     
      void delete_edge( int e ) 
      { 
         adj_list.remove( e );
      }
	  
  };// end of vertex

   //edge structure
   struct Edge
    {
     //left vertex
     vertex left;

     //right vertex
     vertex right;  

     Edge( vertex l, vertex r)
       : left( l ), right( r ) { }
    
   };//end of edge

   ///list of vertices
   typedef vector<Vertex> Vertices;
   Vertices vertices;

  ///list of edges
  //vector<Edge> edges;

 public:

  Graph( ){ }
  
  Graph( int n ){
    vertices.reserve( n );
  }

 
 
  //insert a new vertex with a specific type and value
  bool insert_vertex( string t, int v, int d ){
    vertices.push_back( Vertex( t, v, d ) );
    return true;
  }

  bool insert_vertex( int v ){
    vertices.push_back( Vertex( "", v, 1 ) );
    return true;
  }  
 
  //insert a new edge iff it is not already present in the graph
  bool insert_edge( vertex l, vertex r ){
    if( vertices[l].insert_edge( r ) && vertices[r].insert_edge( l ) )
      return true;
    return false;
  }
   

  Graph& operator = ( const Graph &g ){
    vertices = g.vertices;
    return  *this;
  }
   

  /* A greedy algorithm to approximate independent set problem */
  /* with approximation ratio 2*/
  void ApproximateIndependentSet( ){
    cout << "VertexCover"<<endl;
    int degree_max, s, t;
    list<int>::iterator it;
    vector<Vertex> vertices_bak = vertices;
    
    while( 1 ){     
      degree_max = 0;
      for( int i = 0; i < (int)vertices.size( ); i ++ ){
	if( vertices[i].value == -1 )
	  continue;
	it = vertices[i].adj_list.begin( );
        for(; it != vertices[i].adj_list.end( ); ++it ){
	  if( vertices[*it].value == -1 )
	    continue;

	  int degree = vertices[i].adj_list.size( )+vertices[*it].adj_list.size( ) - 1; 
	  if( degree > degree_max ){
	    s = i;
	    t = *it;
	    degree_max = degree;
	  }
	}
      }
      if( degree_max == 0)
	break;
      
      it = vertices[s].adj_list.begin( );
      for(; it != vertices[s].adj_list.end( ); ++it )
	  vertices[*it].delete_edge( s );

      it = vertices[t].adj_list.begin( );
      for(; it != vertices[t].adj_list.end( ); ++it )
	  vertices[*it].delete_edge( t );
      
      vertices[s].value = -1;
      vertices[t].value = -1;

      for( int i = 0; i < (int)vertices.size( ); i ++ ){
	if( vertices[i].value != -1 )
	  continue;
	int conflictflag = 0;
	it = vertices_bak[i].adj_list.begin( );
	for(; it != vertices_bak[i].adj_list.end( ); ++it ){
	  if( vertices[*it].value != -1 ){
	    conflictflag = 1;
	    break;
	  }
	}
	if( conflictflag == 0 )
	  vertices[i] = vertices_bak[i];
	
      }
            
       cout<<"s="<<s<<" t="<<t<<endl;
      print();
    }
  }

   void Clear( ){
      vertices.clear( );
   }
   
  void print( ){
    for( int i = 0; i < (int)vertices.size( ); i++ )
    {
      cout << vertices[i].vtype <<" "<< i << "("<<vertices[i].value<<"): ";
      list<int> tmp = vertices[i].adj_list;
      int t = (int)tmp.size( );
      for( int j = 0; j < t; j ++ )
	{
	  cout << tmp.front( ) <<" ";
	  tmp.pop_front( );
	}
      cout <<endl;
    }
    cout <<"******************************"<<endl;
 }
 
};//endof class Graph  
  
#endif

//Local Variables:
//mode: c++Com
//End:
