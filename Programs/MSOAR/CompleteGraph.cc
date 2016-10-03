#include <list>
#include <string>
#include <vector>
#include <queue>
#include "CompleteGraph.h"
#include "Genome.h"

using namespace std;

CompleteGraph::CompleteGraph( ){
}

void CompleteGraph::ClearGraph( ){
   vertices.clear( );
   endgeneS = 0;
   endchromS = 0;
   endgeneT = 0;
   endchromT = 0;
}

void CompleteGraph::ConstructGraph( Genome S, Genome T ){
  
  // Add empty chromosomes to make two genomes have same num of chroms 
  if( (int)S.chromsizes.size( ) >= (int)T.chromsizes.size( ) ){
    int gap = (int)S.chromsizes.size( ) - (int)T.chromsizes.size( );
    for( int i = 0; i < gap; i ++ )
      T.chromsizes.push_back( 0 );
  }
  else{
    cout<<"Error: Genome S should have more chomosomes than T."<<endl;
  }
  
  // Insert gene vertices for Genome S 
  for( int i = 0; i < (int)S.genes.size( ); i++ ){
    if( S.genes[i].Sign == 1 ) { // gene strand is +
      insert_vertex( "S_gen", 2*i, S.genes[i].Sign);
      insert_vertex( "S_gen", 2*i + 1, S.genes[i].Sign );
    }
    else { // gene strand is -
      insert_vertex( "S_gen", 2*i + 1, S.genes[i].Sign);
      insert_vertex( "S_gen", 2*i, S.genes[i].Sign);
    }
  }
  
  // Insert chromosome vertices for genome S
  char tmp[10];
  string vt ;
  for( int i = 0; i < (int)S.chromsizes.size( ); i++ ){
    sprintf( tmp, "%d", i );
    vt = string( "S_chr_" ) + string( tmp );
    insert_vertex( vt, 2*i, 1 );
    insert_vertex( vt, 2*i + 1, 1 );
   
  }
 
 // Insert gene vertices for Genome T 
 for( int i = 0; i < (int)T.genes.size( ); i++ ){
    if( T.genes[i].Sign == 1 ) { // gene strand is +
      insert_vertex( "T_gen", 2*i, T.genes[i].Sign );
      insert_vertex( "T_gen", 2*i + 1, T.genes[i].Sign );
    }
    else { // gene strand is -
      insert_vertex( "T_gen", 2*i + 1, T.genes[i].Sign );
      insert_vertex( "T_gen", 2*i, T.genes[i].Sign );
    }
  }
  
  // Insert chromosome vertices for genome S
  for( int i = 0; i < (int)T.chromsizes.size( ); i++ ){
    sprintf( tmp, "%d", i );
    vt = string( "T_chr_" ) + string( tmp );
    insert_vertex( vt, 2*i, 1 );
    insert_vertex( vt, 2*i + 1, 1 );
  }
  
  // Insert Edges within Genome S
  int genenodeindex = 0;
  int chromnodeindex = (int)S.genes.size( ) * 2;
  for( int i = 0; i < (int)S.chromsizes.size( ); i++ ){
    insert_edge( chromnodeindex++ ,genenodeindex++ );
    for( int j = 1; j < S.chromsizes[i]; j++, genenodeindex+=2 ){
      insert_edge( genenodeindex, genenodeindex + 1 );
    }    
    insert_edge( genenodeindex++, chromnodeindex++ );
  }
  
  // Insert Edges within Genome T
  genenodeindex = chromnodeindex;   
  chromnodeindex = genenodeindex + (int)T.genes.size( ) * 2;
  for( int i = 0; i < (int)T.chromsizes.size( ); i++ ){
   //deal with empty chromosome
    if( T.chromsizes[i] == 0 )
      continue;
    vertices[genenodeindex++].vtype = "T_gam";
    for( int j = 1; j < T.chromsizes[i]; j++, genenodeindex+=2 ){
      insert_edge( genenodeindex, genenodeindex + 1 );
    }
    vertices[genenodeindex++].vtype = "T_gam";
  }
  
  // Insert Edges between Genomes 
  int Schromstart = (int)S.genes.size( ) * 2;
  int Tgenestart = ((int)S.genes.size( )+ (int)S.chromsizes.size( )) * 2;
  int Tchromstart = Tgenestart + (int)T.genes.size( ) * 2;
  endgeneS = Schromstart - 1;
  endchromS = Tgenestart - 1;
  endgeneT = Tchromstart -1;
  endchromT = (int)vertices.size( ) - 1;
  // insert edges between blastp hits
  for( int i = 0; i < (int)S.genes.size( ); i++ ){
    for( int j = 0; j < (int)S.genes[i].BlastHits.size( ); j ++ ){
      int obj = S.genes[i].BlastHits[j].ObjectIndex;
      if( S.genes[i].Sign == T.genes[obj].Sign ){
	insert_edge( 2*i, Tgenestart+ obj * 2 );
	insert_edge( 2*i + 1, Tgenestart+ obj * 2 + 1 );
      }
      else{
      	insert_edge( 2*i, Tgenestart+ obj * 2 + 1 );
	insert_edge( 2*i + 1, Tgenestart+ obj * 2 );
      }
    }   
  }  
  // insert edges between chromosom nodes
  for( int i = 0; i < (int)S.chromsizes.size( ); i++ ){
    insert_edge( Schromstart+2*i, Tchromstart+2*i );
    insert_edge( Schromstart+2*i+1, Tchromstart+2*i+1 );
  } 
  //  print( );
}

/*Find shortest cycle/path in a BFS fashion*/
bool CompleteGraph::BFSShortest( int i ){
  int leftstart = i;
  int rightstart = i + 1;
  string vt;
  bool flag1, flag2;

  cout<<"find the shortest cycle/path starting with the left start vertex"<<endl;
  vt = vertices[leftstart].vtype.substr(0,5);
  if( vt == "T_gen")
    flag1 = ShortestCycle( leftstart );
  else
    flag1 = ShortestPath( leftstart ); 
  
  vt = vertices[rightstart].vtype.substr(0,5);
  if( vt == "T_gen")
    flag2 = ShortestCycle( rightstart );
  else
    flag2 = ShortestPath( rightstart );

  return flag1||flag2;
} 

/* Search for shortest cycle starting at vertices i in a BFS fashion*/
/*   At first, push i into the Queue Q. In order to satisfy the     */	
/*   pairing consition, for every element pop out from Q, the       */
/*   algorithm always traverse its black edge first then push all   */
/*   possible grey edges in a Queue Q. Then goes to the next round. */
/*   For every vertex pushed into Q, there is a corresponding       */
/*   complete graph pushed into another Q to keep the consistency   */
bool CompleteGraph::ShortestCycle( int i ){
   
   cout<<"find shortest cycle starting with" <<i<<endl;
   queue<Vertices>Q_graph; //graph Queue
   queue<int>Q_id;//id Queue
   int step = 0;
   
   Vertices CurrentGraph = vertices; 
   int next = CurrentGraph[i].adj_list.back( );//the vertex connect i by a black edge
   if( isidlevertex( CurrentGraph,next ) ) //if next is a idlee vertex, remove it and modify CurrentGraph
      remove_idle_vertices( CurrentGraph, next );
   next = CurrentGraph[i].adj_list.back( );
   //Q_id.push( i );
   Q_id.push( next );
   //cout<<"\t push:"<<i<<endl;
   cout<<"\t push:"<<next<<endl;
   Q_graph.push( CurrentGraph );  
  
   while( !Q_id.empty( )){     
      //Q_id.pop( );
      int u = Q_id.front( );
      cout<<"\t pop:"<<u<<endl;
      Q_id.pop( );
      int utype = vertices[u].vtype[0];
      Vertices uG0 = Q_graph.front( );
      Q_graph.pop( );
    
      //go through every possible grey edges
      list<int> U_adj = uG0[u].adj_list;
      list<int>::iterator it = U_adj.begin( ) ; 
      for(; it != U_adj.end( ); it++ ){
         step ++;
         // terminate the process if the number of steps exceeds cutoff
	 if( step > STEPCUTOFF ){
            cout <<"TERMINATE"<<endl;
            return false;
         } 
         //make a copy of current graph
	 Vertices uG = uG0;
         //exclude the black edge
         if( vertices[*it].vtype[0] == utype )
            continue;
         cout<<"*it="<<*it<<endl;
         //check if it goes back to the start node to make a cycle
         if( *it == i ){
            remove_incompatible_edge( uG, u, *it );            
            vertices = uG; //update complete graph
            cout <<"find cycle *it="<<*it<<endl;
            return true;
         }
         else{
            string tmpstr = uG[*it].vtype.substr(0,5);
	    
	    //if *it is a T_gen or in S chromosome, it definitely has a black edge
	    //so it is necessary to check the if its black edge connecting an idle gene                       
            if( tmpstr == "T_gen" ){
               next = uG[*it].adj_list.back( );
               if( isidlevertex( uG,next ) )
                  remove_idle_vertices( uG, next );        
            }
            else if( tmpstr[0] == 'S' ){
               next = uG[*it].adj_list.front( );
               if( isidlevertex( uG,next ) )
              remove_idle_vertices( uG, next );
            }            
            
	    //push next non idle gene into Queue
            tmpstr = uG[*it].vtype.substr(0,5);
            if( tmpstr == "T_gen" ){
               next = uG[*it].adj_list.back( );
               remove_incompatible_edge( uG, u, *it );          
               //Q_id.push( *it );
               Q_id.push( next );
               //cout<<"\t push:"<<*it<<endl;
               cout<<"\t push:"<<next<<endl;
               Q_graph.push( uG );
            }
            else if( tmpstr[0] == 'S' ){
               next = uG[*it].adj_list.front( );
               remove_incompatible_edge( uG, u, *it );               
               //Q_id.push( *it );
               Q_id.push( next );
               //cout<<"\t push:"<<*it<<endl;
               cout<<"\t push:"<<next<<endl;
               Q_graph.push( uG );	  
            }
         }                 
      }
   }
   // if no cycle found when Queue gets empty, return false
   cout << "Not find cycle" <<endl;
   return false;
}

/* If a gene vertex doesn't have any grey edge connecting it to a gene */
/* vertex in another genome, this gene vertex is called a idle vertex. */
bool CompleteGraph::isidlevertex( Vertices V, int i ){
  if( V[i].vtype[0] == 'S' && (int)V[i].adj_list.size( ) <= 1)
    return true;
  if( V[i].vtype.substr(0,5) == "T_gen" && (int)V[i].adj_list.size( ) <= 1)
    return true;
  if( V[i].vtype.substr(0,5) == "T_gam" && (int)V[i].adj_list.size( ) == 0 )
    return true;
  return false;
}

/* For a idle gene g, its two vertices i and i+1 are idle vertices.  */
/* Remove idle gene g by removing i and i+1 from complete graph and  */
/* connect the vertices i-1 and i+2 by a black edge. If i-1 or i+2   */
/* is idle gene, keep doing the same process until there is no idle  */
/* vertices left in this block. */
void CompleteGraph::remove_idle_vertices( Vertices& V, int i ){
  int left, right, leftmost, rightmost;

  // find the partner idle vertex of i
  // vertex "left" is the left vertex of this idle gene while
  // "right" is the right vertex of this idle gene.
  if( i % 2 == 0 ){
    left = i;
    right = i + 1;    
  }
  else{
    left = i - 1;
    right = i;
  }

  // If the idle vertices is in genome S, then try to find the vertiex 
  // "leftmost" and "rightmost". "Leftmost" is the first non-idle gene
  // vertex to the left of vertex "left"; "rightmost" is the first non-idle
  // gene vertex to the right of vertex "right".
  // Delete all the edges incident on the vertices between "leftmost" and
  // "righatmost" and connect "leftmost" and "rightmost" by a black edge
  if( V[left].vtype.substr( 0,1 ) == "S" ){
    leftmost =  V[left].adj_list.front( );
    while( (int)V[leftmost].adj_list.size( ) == 1){
      leftmost = V[leftmost-1].adj_list.front( ); 
    }
    rightmost =  V[right].adj_list.front( );
    while( (int)V[rightmost].adj_list.size( ) == 1){
      rightmost = V[rightmost+1].adj_list.front( ) ; 
    }
   
    for( int t= V[leftmost].adj_list.front( ); t <= V[rightmost].adj_list.front( ); t ++ ){
       V[t].adj_list.clear( );
    }

    V[leftmost].adj_list.pop_front( );
    V[leftmost].insert_edge( rightmost );
    V[rightmost].adj_list.pop_front( );
    V[rightmost].insert_edge( leftmost );
    
  }
  
  // If the idle vertices is in genome T, the vertex "leftmost" either
  // be the first non-idle gene vertex to the left of vertex "left"
  // or the first "T_gam" gene vertex to the left of the vertex "left";
  // "rightmost" is defined the similar way. 
  
  else{
    leftmost = left;
    while( (int)V[leftmost].adj_list.size( ) == 1 && V[leftmost].vtype.substr( 2, 3 ) == "gen")
      leftmost = V[leftmost].adj_list.back( ) - 1; 
    
    rightmost = right;
    while( (int)V[rightmost].adj_list.size( ) == 1 && V[rightmost].vtype.substr( 2, 3 ) == "gen")
       rightmost = V[rightmost].adj_list.back( ) + 1; 

    if( V[leftmost].vtype.substr( 2, 3 ) == "gam" && (int)V[leftmost].adj_list.size( ) == 0 ){
       leftmost++;
       rightmost--;
//       for( int t= leftmost; t <= rightmost; t ++ ){
	for( int t= leftmost; t <= V[rightmost].adj_list.back( ); t ++ ){
          V[t].adj_list.clear( );
       }
       V[rightmost].vtype.replace(2, 3 ,"gam");
       V[rightmost].adj_list.pop_back( );	
    }
    else if( V[rightmost].vtype.substr( 2, 3 ) == "gam" &&(int)V[rightmost].adj_list.size( ) == 0 ){
       leftmost++;
       rightmost--;
//       for( int t= leftmost; t <= rightmost; t ++ ){
	for( int t= rightmost; t >=V[leftmost].adj_list.back( ); t -- ){
          V[t].adj_list.clear( );
       }
       V[leftmost].vtype.replace( 2, 3 , "gam");
       V[leftmost].adj_list.pop_back( );
    }
    else {
       leftmost++;
       rightmost--;

       for( int t= V[leftmost].adj_list.back( ); t <= V[rightmost].adj_list.back( ); t ++ ){
          V[t].adj_list.clear( );
       }   
       V[leftmost].adj_list.pop_back( );
       V[leftmost].insert_edge( rightmost );
       V[rightmost].adj_list.pop_back( );
       V[rightmost].insert_edge( leftmost ); 
    }
  }
} 

/* Search for shortest path starting at vertices i in a BFS fashion    */
/*   At first, push i into the Queue Q_id. In order to satisfy the     */	
/*   pairing consition, for every element pop out from Q_id, the       */
/*   algorithm always traverse its black edge first then push all      */
/*   possible grey edges in a Queue Q_id. Then goes to the next round. */
/*   For every vertex pushed into Q_id, there is a corresponding       */
/*   complete graph pushed into another Q_id to keep the consistency   */
bool CompleteGraph::ShortestPath( int i ){
  cout<<"find shortest path starting with" <<i<<endl;
  queue<Vertices>Q_graph; //graph Queue
  queue<int>Q_id;//id Queue
  int step = 0;
  
  //Q_id.push( i );
  Q_id.push( i );
  Q_graph.push( vertices );
  
  while( !Q_id.empty( )){
    //Q_id.pop( );
    int u = Q_id.front( );
    cout<<"\t pop:"<<u<<endl;
    Q_id.pop( );
    int utype = vertices[u].vtype[0];
    Vertices uG0 = Q_graph.front( );
    Q_graph.pop( );
    
    list<int> U_adj = uG0[u].adj_list;
    list<int>::iterator it = U_adj.begin( ) ; 
    for(; it != U_adj.end( ); it++ ){
       step ++;
       if( step > STEPCUTOFF ){
          cout <<"TERMINATE"<<endl;
          return false;
       } 
      Vertices uG = uG0;
      string tmpstr = vertices[*it].vtype;
      //exclude the black edge
      if( tmpstr[0] == utype )
	continue;
      cout<<"*it="<<*it<<endl;
      //check if it goes to a dead end and form a path
      if( tmpstr.substr( 0,5)=="T_chr" ||tmpstr.substr( 0,5)=="T_gam" ){
	remove_incompatible_edge( uG, u, *it );
	vertices = uG;

	cout<<"find path *it="<<*it<<endl;

	return true;
      }
      else{
         int next;
	string tmpstr = uG[*it].vtype.substr(0,5);
        if( tmpstr == "T_gen" )
           next = uG[*it].adj_list.back( );
        else if( tmpstr[0] == 'S' )
           next = uG[*it].adj_list.front( );

        if( isidlevertex( uG,next ) )
           remove_idle_vertices( uG, next );
        
        tmpstr = uG[*it].vtype.substr(0,5);
	if( tmpstr == "T_gen" ){
           next = uG[*it].adj_list.back( );
           remove_incompatible_edge( uG, u, *it );

	 // Q_id.push( *it );
	  Q_id.push( next );
          //cout<<"\t push:"<<*it<<endl;
	  cout<<"\t push:"<<next<<endl;
	  Q_graph.push( uG );
	}
	else if( tmpstr[0] == 'S' ){
	  next = uG[*it].adj_list.front( );	 
	  remove_incompatible_edge( uG, u, *it );
	  //Q_id.push( *it );
	  Q_id.push( next );
          //cout<<"\t push:"<<*it<<endl;
	  cout<<"\t push:"<<next<<endl;
	  Q_graph.push( uG );	  
	}
      }
    }
  }
  cout << "Not find path" <<endl;
  return false;
}

/*If a gene g of genome S is matched to gene h in genome T, then*/
/*the two vertices of g "up" and "up2" will be matched to the two vertices*/
/*of h "down" and "down2" respectively. In this case, we need to update*/
/*the complete graph by removing incompativle grey edges*/
void CompleteGraph::remove_incompatible_edge(Vertices& V, int i, int j ){
  int up,down,up2, down2;
  if( i < j ){
    up = i; 
    down = j;
  }
  else{
    up = j;
    down = i;
  }
  
  //if they are chromosome vertices, do nothing and return
  if( up > endgeneS )
    return;

  //find two related vertices of up and down
  up2 = ( up % 2 == 0 )? up + 1: up - 1;
  down2 = ( down % 2 == 0 )? down + 1: down - 1;
 
  //remove incompatible grey edges
  list<int>::iterator it = V[up2].adj_list.begin( );
  for(; it != V[up2].adj_list.end( ); ++it ){
    if( *it != down2 && *it > endchromS ){	    
      V[up2].delete_edge( *it );
      V[*it].delete_edge( up2 );
      it --;
    }
  }
  it = V[down2].adj_list.begin( );
  for(; it != V[down2].adj_list.end( ); ++it ){
    if( *it != up2 && *it <= endchromS ){
      V[down2].delete_edge( *it );
      V[*it].delete_edge( down2 );
      it --;
    }
  }

  it = V[up].adj_list.begin( );
  for(; it != V[up].adj_list.end( ); ++it ){
    if( *it != down && *it > endchromS ){
      V[up].delete_edge( *it );
      V[*it].delete_edge( up );
      it --;
    }
  }
 
  it = V[down].adj_list.begin( );
  for(; it != V[down].adj_list.end( ); ++it ){
    if( *it != up && *it <= endchromS ){
      V[down].delete_edge( *it );
      V[*it].delete_edge( down );
      it --;
    }
  }
}

/*Calculate the approximate distance between two genomes */
/*d = b - c - p1 + p2*/
/*b: the number of black edges*/
/*c: the number of all the cycles*/
/*p1: the number of all the paths*/
/*p2: the number of all the Pai-Pai paths*/
int CompleteGraph::ApproximateDistance( ){

   //Calculate b
   int b = 0;
   for( int i = 0; i < endgeneS; i+=2 ){
      if( (int)vertices[i].adj_list.size( ) == 2 )
         b++;
   }
   b = b + ( endchromS - endgeneS ) / 2;

   vector<int>flag( (int)vertices.size( ), 0 );

   /*Calculate p1 and p2 */
   int p1 = 0;
   int p2 = 0; 
   for( int i = endchromS+1; i < (int)vertices.size( ); i++ ){
      if( flag[i] != 0 )
         continue;
      if( (int)vertices[i].adj_list.size( ) == 0 )
         continue;      
      string tmpstr = vertices[i].vtype.substr(0,5);
      int start, end;
      if( tmpstr == "T_gam" )
         start = 0;
      else if( tmpstr == "T_chr" )
         start = 1;
      else
         continue;
      queue<int> Q;
      Q.push( i );
      flag[i] = 1;
      while( !Q.empty( ) ){
         int u = Q.front( );
         Q.pop( );
         if( !Q.empty( ) ){
            cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE1"<<endl;
            exit( 1 );
         }                
         list<int> U_adj = vertices[u].adj_list;
         list<int>::iterator it = U_adj.begin( ) ; 
         for(; it != U_adj.end( ); it++ ){
            if( flag[*it] != 0 )
               continue;
            flag[*it] = 1;
            tmpstr = vertices[*it].vtype.substr(0,5);
            if( tmpstr == "T_gam" ){
               end = 0;
               break;
            }
            else if( tmpstr == "T_chr" ){
               end = 1;
               break;
            }
            Q.push( *it );
           
         }
       
      }
      if( !Q.empty( ) ){
         cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE1.1"<<endl;
         exit( 1 );
      }  
      p1 ++;
      if( start == end )
         p2 ++;
   }

   /*Calculate c */
   int c = 0;
   for( int i = 0; i < (int)endchromS; i++ ){
      if( flag[i] != 0 )
         continue;
      if( (int)vertices[i].adj_list.size( ) == 0 )
         continue;
      string tmpstr = vertices[i].vtype.substr(0,5);
      if( tmpstr == "T_gam" || tmpstr == "T_chr" ){
         cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE2"<<endl;
         exit( 1 );
      }
      queue<int> Q;
      int next = vertices[i].adj_list.front( );//the vertex connect i by a
      Q.push(next);
      flag[i] = 1;
      flag[next] = 1;
      
      while( !Q.empty( ) ){                    
         int u = Q.front( );
         Q.pop( );
         
         if( !Q.empty( ) ){            
            cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE3"<<endl;
            exit( 1 );
         }
         
         list<int> U_adj = vertices[u].adj_list;
         list<int>::iterator it = U_adj.begin( ) ; 
         for(; it != U_adj.end( ); it++ ){
            if(u == i )              
               break;   
            if( flag[*it] != 0 )
               continue;
            Q.push( *it );
            flag[*it] = 1;
         }
      }
     
      c++;
   }
   cout<<"b="<<b<<"\t c="<<c<<"\t p1="<<p1<<"\t p2="<<p2<<endl;
   return b - c - p1 + p2;     
}
   
 
/*delete S.genes[i] and T.genes[j] then compute the distance*/
int CompleteGraph::ApproximateDistance( int i, int j ){
   
   int vidil = i * 2;
   int vidir = vidil + 1;
   int vidjl = endchromS + 1 + j * 2;
   int vidjr = vidjl + 1;

   if( vertices[vidil].adj_list.back( ) != vidjl && vertices[vidil].adj_list.back( ) != vidjr ){
      cout << vidil << "\t" <<vidir <<"\t"<<vidjl<<"\t"<<vidjr<<endl;
      cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE4"<<endl;
      exit( 1 );
   }
   if( vertices[vidir].adj_list.back( ) != vidjl && vertices[vidir].adj_list.back( ) != vidjr ){
      cout << vidil << "\t" <<vidir <<"\t"<<vidjl<<"\t"<<vidjr<<endl;
      cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE5"<<endl;
      exit( 1 );
   }
   if( vertices[vidjl].adj_list.front( ) != vidil && vertices[vidjl].adj_list.front( ) != vidir ){
      cout << vidil << "\t" <<vidir <<"\t"<<vidjl<<"\t"<<vidjr<<endl;
      cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE6"<<endl;
      exit( 1 );
   }
   if( vertices[vidjr].adj_list.front( ) != vidil && vertices[vidjr].adj_list.front( ) != vidir ){
      cout << vidil << "\t" <<vidir <<"\t"<<vidjl<<"\t"<<vidjr<<endl;
      cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE7"<<endl;
      exit( 1 );
   }

   int left = vertices[vidil].adj_list.front( );
   int right = vertices[vidir].adj_list.front( );
   if( vertices[left].adj_list.front( ) != vidil || vertices[right].adj_list.front( ) != vidir ){
      cout << left << "\t" <<vidil <<"\t"<<vidir<<"\t"<<right<<endl;
      cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE8"<<endl;
      exit( 1 );
   }
   
   vertices[left].adj_list.pop_front( );
   vertices[left].insert_edge( right );
   vertices[right].adj_list.pop_front( );
   vertices[right].insert_edge( left );
   vertices[vidil].adj_list.clear( );
   vertices[vidir].adj_list.clear( );

   string tmpstrl = vertices[vidjl].vtype.substr(0,5);
   string tmpstrr = vertices[vidjr].vtype.substr(0,5);
   if( tmpstrl != "T_gam" && tmpstrr != "T_gam"){
      left = vertices[vidjl].adj_list.back( );
      right = vertices[vidjr].adj_list.back( );
      if( vertices[left].adj_list.back( ) != vidjl || vertices[right].adj_list.back( ) != vidjr ){
         cout << left << "\t" <<vidjl <<"\t"<<vidjr<<"\t"<<right<<endl;
         cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE9"<<endl;
         exit( 1 );
      }
      
      vertices[left].adj_list.pop_back( );
      vertices[left].insert_edge( right );
      vertices[right].adj_list.pop_back( );
      vertices[right].insert_edge( left );  
   }
   else if( tmpstrl == "T_gam" && tmpstrr == "T_gam" ){
       vertices[vidil].adj_list.clear( );
       vertices[vidir].adj_list.clear( );
   }
   else if( tmpstrl == "T_gam" && tmpstrr != "T_gam" ){
      right = vertices[vidjr].adj_list.back( );
      if(  vertices[right].adj_list.back( ) != vidjr ){
         cout << left << "\t" <<vidjl <<"\t"<<vidjr<<"\t"<<right<<endl;
         cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE10"<<endl;
         exit( 1 );
      }
      vertices[right].adj_list.pop_back( );
      vertices[right].vtype.replace(0,5,"T_gam");
      
   }
   else if( tmpstrr == "T_gam" && tmpstrl != "T_gam"){
      left = vertices[vidjl].adj_list.back( );
      if(  vertices[left].adj_list.back( ) != vidjl ){
         cout << left << "\t" <<vidjl <<"\t"<<vidjr<<"\t"<<right<<endl;
         cout<<"ERROR@COMPLETEGRAPH@APPROXIMATEDISTANCE11"<<endl;
         exit( 1 );
      }
      vertices[left].adj_list.pop_back( );
      vertices[left].vtype.replace(0,5,"T_gam");      
      
   }
   else{
      cout<<"STRANGE!"<<endl;
      exit( 1 );
   }
            
   vertices[vidjl].adj_list.clear( );
   vertices[vidjr].adj_list.clear( );
   return  ApproximateDistance( );
   
}

/* Delete a gene pair from the complete graph: */
/* gene i from genome S and gene j from genome T*/
void CompleteGraph::deletePair( int i, int j ){
   int vidil = i * 2;
   int vidir = vidil + 1;
   int vidjl = endchromS + 1 + j * 2;
   int vidjr = vidjl + 1;

   
   int left = vertices[vidil].adj_list.front( );
   int right = vertices[vidir].adj_list.front( );
    
   vertices[left].adj_list.pop_front( );
   vertices[left].insert_edge( right );
   vertices[right].adj_list.pop_front( );
   vertices[right].insert_edge( left );
   vertices[vidil].adj_list.clear( );
   vertices[vidir].adj_list.clear( );

   string tmpstrl = vertices[vidjl].vtype.substr(0,5);
   string tmpstrr = vertices[vidjr].vtype.substr(0,5);
   if( tmpstrl != "T_gam" && tmpstrr != "T_gam"){
      left = vertices[vidjl].adj_list.back( );
      right = vertices[vidjr].adj_list.back( );
      vertices[left].adj_list.pop_back( );
      vertices[left].insert_edge( right );
      vertices[right].adj_list.pop_back( );
      vertices[right].insert_edge( left );
    
   }
   else if( tmpstrl == "T_gam" ){
      right = vertices[vidjr].adj_list.back( );
      vertices[right].adj_list.pop_back( );
      vertices[right].vtype.replace(0,5, "T_gam");                                                           
   }
   else if( tmpstrr == "T_gam" ){
      left = vertices[vidjl].adj_list.back( );
      vertices[left].adj_list.pop_back( );
      vertices[left].vtype.replace(0,5,"T_gam");      
      
   }
   else{
      cout<<"STRANGE!"<<endl;
      exit( 1 );
   }
   vertices[vidjl].adj_list.clear( );
   vertices[vidjr].adj_list.clear( );
   
}
