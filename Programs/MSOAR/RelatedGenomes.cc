#include <vector>
#include <fstream>
#include <string>
#include <queue>
#include <time.h>
#include <math.h>
#include <algorithm>
#include "Gene.h"
#include "Genome.h"
#include "CompleteGraph.h"
#include "PairmatchGraph.h"
#include "RelatedGenomes.h"

using namespace std;

RelatedGenomes::RelatedGenomes( vector<int>A, vector<int>B, Genome s, Genome t ){
   matchesS = A;
   matchesT = B;
   S = s;
   T = t;
}

/*Construct RelatedGenomes by reading data from infile1 and infile2*/
RelatedGenomes::RelatedGenomes( string infile1, string infile2 ){
   
   //Read in the gene lists for two genomes
   cout.flush( );
   vector<Gene>genes;
   vector<BlastHit>blasthits;
   string query_index;
   string query_accession;
   string query_symbol;
   int query_chromid;
   int query_sign;
   int query_posonchrom;
   int obj_index;
   double evalue;
   double bitscore;
   string tmp;
   string pre_index = "0";
   string previous_chrom = "chr1";
   int previous_chromid = 1;
   
   ifstream fin1( infile1.c_str( ) );
   ifstream fin2( infile2.c_str( ) );
   ifstream fin;
   
   //Load two input files respectively
   for( int i = 0; i < 2; i ++ ){
      pre_index = "0";
      previous_chrom = "chr0";
      previous_chromid = 0;
      tmp = ( i == 0 )? infile1:infile2;
      ifstream fin( tmp.c_str( ) );
      
      genes.clear( );
      blasthits.clear( );
			int indexnum=0;
      while( getline( fin, query_index,'\t' ) ){
				indexnum++;
         
         if( query_index != pre_index && pre_index != "0" ){
            Gene currentgene ( query_accession, query_symbol, query_chromid, query_posonchrom, query_sign, -1,blasthits );
            genes.push_back( currentgene );
            blasthits.clear( );
            pre_index = query_index;   	
      }
         
         if( pre_index == "0" )
            pre_index = query_index;
         
         getline( fin, query_accession,'\t' );
         
         getline( fin, query_symbol,'\t' );

         getline( fin, tmp,'\t' );

         if ( tmp != previous_chrom ){
            previous_chrom = tmp;           
            previous_chromid ++;
         }
         query_chromid = previous_chromid;
        
         getline( fin, tmp,'\t' );
         query_sign = (tmp =="+")? 1:0;
         
         getline( fin, tmp,'\t' );
         query_posonchrom = atoi( tmp.c_str( ) );
         
         getline( fin, tmp,'\t' );
         if( tmp != "" ){
            obj_index = atoi( tmp.c_str( ) ) - 1;
            
            getline( fin, tmp,'\t' );
            evalue = atof( tmp.c_str( ) );
            
            getline( fin, tmp,'\n' );
            bitscore = atof( tmp.c_str( ) );
            
            BlastHit currentblasthit( obj_index, evalue, bitscore);
            blasthits.push_back( currentblasthit );
         }
         else{         
            getline( fin, tmp,'\t' );
            getline( fin, tmp,'\n' );         
            
         }
      }
      
      Gene currentgene ( query_accession, query_symbol, query_chromid, query_posonchrom, query_sign, -1,blasthits );
      genes.push_back( currentgene );
      
      if( i == 0 ){     
         S.genes = genes;
         S.getChromSizes( );
      }
      else {
         T.genes = genes;
         T.getChromSizes( );
      }
   }
//S.print( );
//cout<<"&&&"<<endl;
//T.print( );
   
  
   
   //only keep the bi-dirctional hits  
   for( int i = 0; i < (int)S.genes.size( ); i++ ){
      for( int j = 0; j < (int)S.genes[i].BlastHits.size( ); j++ ){
         
         int obj = S.genes[i].BlastHits[j].ObjectIndex;
         
				 cout << i << " " << j << " " << obj << endl;

         if( !T.genes[obj].isBlastHit( i ) ){
            S.genes[i].deleteBlastHit( obj );
            j--;
         }      

      }
   }
   
	// cout << "I am here " << endl;
   for( int i = 0; i < (int)T.genes.size( ); i++ ){
      for( int j = 0; j < (int)T.genes[i].BlastHits.size( ); j++ ){
         int obj = T.genes[i].BlastHits[j].ObjectIndex;
         
         if( !S.genes[obj].isBlastHit( i ) ){
            T.genes[i].deleteBlastHit( obj );
            j--;
         }
      }
			//cout << i << endl;
   } 

   
    //S.print( );
   //cout <<"&&&&&&&&&&&&&"<<endl;
    //T.print( );
   
   //Initialize the vector of matches
   for( int i = 0; i < (int)T.genes.size( ); i++ )
      matchesT.push_back( -1 );
  
   for( int i = 0; i < (int)S.genes.size( ); i++ ){
      if( (int)S.genes[i].BlastHits.size( ) != 1 ){
         matchesS.push_back( -1 );
         continue;
      }
      int j = S.genes[i].BlastHits[0].ObjectIndex;
      if( (int)T.genes[j].BlastHits.size( ) == 1  ){
         matchesS.push_back( j );
         matchesT[j] =i;
         S.genes[i].Match = j;
         T.genes[j].Match = i;
      }
      else
         matchesS.push_back( -1 );    
   }
   
   //Initialize the vector of modified
   modified.clear( );
   
   //Construct pairmatch graph P
   //P.ConstructGraph( S,T );
   
   //Construct complete graph G
   // G.ConstructGraph( S,T );
   
}

/*Used for debugging*/
void RelatedGenomes::CheckBlastHits( ){
   cout <<"One direction"<<endl;
   for( int i = 0; i < (int)S.genes.size( ); i++ ){
      for( int j = 0; j < (int)S.genes[i].BlastHits.size( ); j++ ){         
         int obj = S.genes[i].BlastHits[j].ObjectIndex;         
         if( !T.genes[obj].isBlastHit( i ) )
            cout<<"ERROR@CHECKBLASTHITS: "<< i<<"~"<<obj<<endl;      
      }
   }
   cout<<"Another direction"<<endl;
   for( int i = 0; i < (int)T.genes.size( ); i++ ){
      for( int j = 0; j < (int)T.genes[i].BlastHits.size( ); j++ ){         
         int obj = T.genes[i].BlastHits[j].ObjectIndex;         
         if( !S.genes[obj].isBlastHit( i ) )
            cout<<"ERROR@CHECKBLASTHITS: "<< obj<<"~"<<i<<endl;      
      }
   }
}

/* Elimenate the obvious false positive*/
void RelatedGenomes::RemoveFP( ){
   for( int i = 0; i < (int)S.genes.size( ); i++ ){
      double min = S.genes[i].getBestEvalue( );
      for( int j = 0; j < (int)S.genes[i].BlastHits.size( ); j++ ){
 	double cutoff;
        if( min == 0.0) 
	    cutoff = 0.0;
        else
	    cutoff = pow(10,(log10( min ) * (1 - EVALUE_DIFF_CUTOFF)));
        if( S.genes[i].BlastHits[j].Evalue > cutoff){
            int obj = S.genes[i].BlastHits[j].ObjectIndex;
            S.genes[i].deleteBlastHit( obj );
            T.genes[obj].deleteBlastHit( i );
            j--;
         }
      }
   }
   
   for( int i = 0; i < (int)T.genes.size( ); i++ ){
      double min = T.genes[i].getBestEvalue( );
      for( int j = 0; j < (int)T.genes[i].BlastHits.size( ); j++ ){
         double cutoff;
         if( min == 0.0)
            cutoff = 0.0;
         else
            cutoff = pow(10,(log10( min ) * (1 - EVALUE_DIFF_CUTOFF)));
         if( T.genes[i].BlastHits[j].Evalue > cutoff){
            int obj = T.genes[i].BlastHits[j].ObjectIndex;
            T.genes[i].deleteBlastHit( obj );
            S.genes[obj].deleteBlastHit( i );
            j--;
         }
      }
   }

   //only keep the bi-dirctional hits  
   for( int i = 0; i < (int)S.genes.size( ); i++ ){
      for( int j = 0; j < (int)S.genes[i].BlastHits.size( ); j++ ){
         
         int obj = S.genes[i].BlastHits[j].ObjectIndex;
         
         if( !T.genes[obj].isBlastHit( i ) ){
            S.genes[i].deleteBlastHit( obj );
            j--;
         }      
      }
   }
   
   for( int i = 0; i < (int)T.genes.size( ); i++ ){
      for( int j = 0; j < (int)T.genes[i].BlastHits.size( ); j++ ){
         int obj = T.genes[i].BlastHits[j].ObjectIndex;
         
         if( !S.genes[obj].isBlastHit( i ) ){
            T.genes[i].deleteBlastHit( obj );
            j--;
         }
      }
   } 
}

/*Three suboptimal Rules*/
void RelatedGenomes::SuboptimalRules( ){
  cout<<"Sub optimal rule 1"<<endl;
  SuboptimalRule1( );
  cout<<"Sub optimal rule 2"<<endl;
  SuboptimalRule2( );  
  cout<<"Sub optimal rule 3"<<endl;
  SuboptimalRule3( );
  
}  

/*Approximate minimum common partition problem*/
void RelatedGenomes::ApproximateMCP( ){
  P.ConstructGraph( S, T );
  P.ApproximateIndependentSet( );
  UpdateMatches_MCP( );
}

/*Approximate maximum graph decomposition*/
bool RelatedGenomes::ApproximateMGD( ){
  modified.clear( );
  int i;
  cout<<"debug:check if any gene unmatched left"<<endl;
  for( i = 0; i < (int)T.genes.size( ); i++ ){
    if( tobeMatchedT( i ))
      break;
  }
  if( i == (int)T.genes.size( ))
    return true; // no genes left for matching

  cout<<"debug:compute the left vertex id for this gene in Complete Graph"<<endl;
  int gid = G.endchromS + 1 + i * 2; 
  cout<<"debug:find shortest path/cycle"<<endl;
  
  if( G.BFSShortest( gid ) ){
     cout<<"debug:update mathces"<<endl;
     UpdateMatches_MGD( );
  }
  else{
     RandomMatch( i ); 
  }
  return false;
}

bool RelatedGenomes::ApproximateMGD( int start ){
  modified.clear( );
  int i;
  cout<<"debug:check if any gene unmatched left"<<endl;
  for( i = start; i < (int)T.genes.size( ); i++ ){
    if( tobeMatchedT( i ))
      break;
  }
  if( i == (int)T.genes.size( ))
    return true; // no genes left for matching

  cout<<"debug:compute the left vertex id for this gene in Complete Graph"<<endl;
  int gid = G.endchromS + 1 + i * 2; 
  cout<<"debug:find shortest path/cycle"<<endl;
 
  if( G.BFSShortest( gid ) ){
     cout<<"debug:update mathces"<<endl;
     UpdateMatches_MGD( );
  }
  else{
     RandomMatch( i );
     UpdateCompleteGraph( );
  }
  return false;
}

/*Randomly match genes if unmatched genes are still available after MGD*/
void RelatedGenomes::RandomMatch( ){
  srand( (unsigned)time( NULL ) );
   for( int i = 0; i < (int)matchesT.size( ); i ++ ){
      if( !tobeMatchedT( i ) )
         continue;
      int id = rand( )%(int)T.genes[i].BlastHits.size( );
      int j = T.genes[i].BlastHits[id].ObjectIndex;
      matchesS[j] = i;
      matchesT[i] = j;
      S.genes[j].Match = i;
      T.genes[i].Match = j;
      selectBlastHit( j, i );
   }
}

/*Randomly match gene "i" */
void RelatedGenomes::RandomMatch( int i ){
  srand( (unsigned)time( NULL ) );
  if( !tobeMatchedT( i ) )
        return;
  double maxscore = 0;
  int id = 0;
  for( int j = 0; j < (int)T.genes[i].BlastHits.size( ); j++ ){
     if( T.genes[i].BlastHits[j].BitScore > maxscore ){
        id = j;
        maxscore = T.genes[i].BlastHits[j].BitScore;
     }   
  }
 
  int j = T.genes[i].BlastHits[id].ObjectIndex;
  matchesS[j] = i;
  matchesT[i] = j;
  S.genes[j].Match = i;
  T.genes[i].Match =j;          
  selectBlastHit( j, i );
  modified.push_back( j );
  
}

/*Check if gene "i" is available for matching*/
bool RelatedGenomes::tobeMatchedT( int i ){
  if( matchesT[i] != -1 )
    return false; //already matched
  if( (int)T.genes[i].BlastHits.size( ) == 0)
    return false; //no one to match
  return true;
}

/*Update two vectors "matchesS" and "matchesT" after MGD*/
void RelatedGenomes::UpdateMatches_MGD( ){
   for( int i = 0; i < (int)matchesS.size( ); i ++ ){
      if( matchesS[i] != -1 )
         continue;
      int gid = i * 2; 
      if( (int)G.vertices[gid].adj_list.size( ) == 2 ){
         int obj = G.vertices[gid].adj_list.back( );
         if(( ((int)G.vertices[obj].adj_list.size( )== 2 && G.vertices[obj].vtype.substr(0,5) == "T_gen") ||( (int)G.vertices[obj].adj_list.size( )== 1 && G.vertices[obj].vtype.substr(0,5) == "T_gam") ) && G.vertices[obj].adj_list.front( ) == gid ){
            int j = (obj -G.endchromS - 1)/2;
            matchesS[i] = j;
            matchesT[j] = i;
            S.genes[i].Match = j;
            T.genes[j].Match = i;
            selectBlastHit( i, j );
         }
      }
   }
}

/*Update two vectors "matchesS" and "matchesT" after MCP*/
void RelatedGenomes::UpdateMatches_MCP( ){
  //update matches modified genome S and T  
   for( int i = 0; i < (int)P.vertices.size( ); i++ ){
      if( P.vertices[i].value == -1 )
         continue;
      
      int j = P.vertices[i].value;
      int ps = P.vmPairs[j].vertexPairS;
      int pt= P.vmPairs[j].vertexPairT;
      int po = P.vmPairs[j].orient;
    
      if( matchesS[ps] != -1 && matchesS[ps+1] != -1 )
         continue;
    
      if( po == 0 ){
         matchesS[ps] = pt;
         matchesT[pt] = ps;
         S.genes[ps].Match = pt;
         T.genes[pt].Match = ps;
         matchesS[ps+1] = pt+1;
         matchesT[pt+1] = ps+1;
         S.genes[ps+1].Match = pt+1;
         T.genes[pt+1].Match = ps+1;
         cout<< "ps="<<ps<<"\t pt="<<pt<<endl;   
         selectBlastHit( ps, pt );
         cout<< "ps+1="<<ps+1<<"\t pt+1="<<pt+1<<endl;  
         selectBlastHit( ps+1, pt+1 );
         if( (int)modified.size( )== 0 || ps != modified[(int)modified.size( )-1] )
            modified.push_back(ps);
         modified.push_back(ps+1);
         P.vertices[i].value = -1;
      }
      else if ( po == 1 ){
         matchesS[ps] = pt+1;
         matchesT[pt+1] = ps;
         S.genes[ps].Match = pt+1;
         T.genes[pt+1].Match = ps;
         matchesS[ps+1] = pt;
         matchesT[pt] = ps+1;
         S.genes[ps+1].Match = pt;
         T.genes[pt].Match = ps+1;
         cout<< "ps="<<ps<<"\t pt+1="<<pt<<endl;  
         selectBlastHit( ps, pt+1 );
         cout<< "ps+1="<<ps<<"\t pt="<<pt<<endl;  
         selectBlastHit( ps+1, pt );
         if( (int)modified.size( )== 0|| ps != modified[(int)modified.size( ) - 1] )
            modified.push_back(ps);
         modified.push_back(ps+1);
         P.vertices[i].value = -1;
      }
      else if ( po == 2 )
         continue;
      else{
         cout<<"ERROR@RELATEDGENOMES@UPDATEMATCHES_MCP1"<<endl;
         exit( 1 );
      }
  }

   //dealing with po == 2
   for( int i = 0; i < (int)P.vertices.size( ); i++ ){
      if( P.vertices[i].value == -1 )
         continue;
      
      int j = P.vertices[i].value;
      int ps = P.vmPairs[j].vertexPairS;
      int pt= P.vmPairs[j].vertexPairT;
      int po = P.vmPairs[j].orient;
    
      if( matchesS[ps] != -1 || matchesS[ps+1] != -1 )
         continue;
    
      if( po == 2 ){
         if( S.genes[ps].isBlastHit( pt ) && S.genes[ps+1].isBlastHit( pt+1 ) ){
            matchesS[ps] = pt;
            matchesT[pt] = ps;
            S.genes[ps].Match = pt;
            T.genes[pt].Match = ps;
            matchesS[ps+1] = pt+1;
            matchesT[pt+1] = ps+1;
            S.genes[ps+1].Match = pt+1;
            T.genes[pt+1].Match = ps+1;
            cout<< "ps="<<ps<<"\t pt="<<pt<<endl;   
            selectBlastHit( ps, pt );
            cout<< "ps+1="<<ps+1<<"\t pt+1="<<pt+1<<endl;  
            selectBlastHit( ps+1, pt+1 );
            if( (int)modified.size( )== 0 || ps != modified[(int)modified.size( )-1] )
               modified.push_back(ps);
            modified.push_back(ps+1);
            P.vertices[i].value = -1;
         }
         else if( S.genes[ps].isBlastHit( pt+1 ) && S.genes[ps+1].isBlastHit( pt ) ){
            matchesS[ps] = pt+1;
            matchesT[pt+1] = ps;
            S.genes[ps].Match = pt+1;
            T.genes[pt+1].Match = ps;
            matchesS[ps+1] = pt;
            matchesT[pt] = ps+1;
            S.genes[ps+1].Match = pt;
            T.genes[pt].Match = ps+1;
            cout<< "ps="<<ps<<"\t pt+1="<<pt<<endl;  
            selectBlastHit( ps, pt+1 );
            cout<< "ps+1="<<ps<<"\t pt="<<pt<<endl;  
            selectBlastHit( ps+1, pt );
            if( (int)modified.size( )== 0|| ps != modified[(int)modified.size( ) - 1] )
               modified.push_back(ps);
            modified.push_back(ps+1);
            P.vertices[i].value = -1;            
         }
         else{
            cout<< "ps="<<ps<<"\t pt="<<pt<<"\t po="<<po<<endl;  
            cout<<"ERROR@RELATEDGENOMES@UPDATEMATCHES_MCP2"<<endl;
            exit( 1 );
         }
      }
      else{
         cout<<"ERROR@RELATEDGENOMES@UPDATEMATCHES_MCP3"<<endl;
         exit( 1 );
      }
   }
}


/*Use vector "modified" to modify complete graph "G"*/
void RelatedGenomes::UpdateCompleteGraph( ){
  //Use matches and modified to update CompleteGraph G
   for( int i = 0; i < (int)modified.size( ); i ++ ){

      int geneS = modified[i];
      int geneT = matchesS[geneS]; 
      int leftS = geneS*2;
      int rightS = leftS + 1;
      int leftT, rightT;
      int Tgenestart = G.endchromS + 1;

      if( S.genes[geneS].Sign == T.genes[geneT].Sign ){
         leftT = Tgenestart + geneT * 2;
         rightT = leftT + 1;
      }
      else{
         rightT = Tgenestart + geneT * 2;
         leftT = rightT + 1;
      }
      
      list<int>::iterator it;

      it = G.vertices[leftS].adj_list.begin( );
      for(; it != G.vertices[leftS].adj_list.end( ); ++it ){
         string vt = G.vertices[*it].vtype.substr(0,1);
         if( vt == "T" && *it !=leftT){
            G.vertices[leftS].delete_edge( *it );
            G.vertices[*it].delete_edge( leftS );
            it --;
         }
      }
      
      it = G.vertices[rightS].adj_list.begin( );
      for(; it != G.vertices[rightS].adj_list.end( ); ++it ){
         string vt = G.vertices[*it].vtype.substr(0,1);
         if( vt == "T" && *it !=rightT){
            G.vertices[rightS].delete_edge( *it );
            G.vertices[*it].delete_edge( rightS );
            it --;
         }
      }
      
      it = G.vertices[leftT].adj_list.begin( );
      for(; it != G.vertices[leftT].adj_list.end( ); ++it ){
         string vt = G.vertices[*it].vtype.substr(0,1);
         if( vt == "S" && *it !=leftS){
            G.vertices[leftT].delete_edge( *it );
            G.vertices[*it].delete_edge( leftT );
            it --;
         }
      }

      it = G.vertices[rightT].adj_list.begin( );
      for(; it != G.vertices[rightT].adj_list.end( ); ++it ){
         string vt = G.vertices[*it].vtype.substr(0,1);
         if( vt == "S" && *it !=rightS){
            G.vertices[rightT].delete_edge( *it );
            G.vertices[*it].delete_edge( rightT );
            it --;
         }
      }
   }
  
   modified.clear();
}


/*Select gene "obj" of genome T to be the main ortholog to gene "id" of genome S*/
void RelatedGenomes::selectBlastHit( int id, int obj ){
  for( int i = 0; i < (int)S.genes[id].BlastHits.size( ); i++ ){
    int j = S.genes[id].BlastHits[i].ObjectIndex;
    if( j != obj )
      T.genes[j].deleteBlastHit( id );
  }
  
  S.genes[id].selectBlastHit( obj );

 for( int i = 0; i < (int)T.genes[obj].BlastHits.size( ); i++ ){
    int j = T.genes[obj].BlastHits[i].ObjectIndex;
    if( j != id )
      S.genes[j].deleteBlastHit( obj );
  }
  
  T.genes[obj].selectBlastHit( id );

}

/*First sub-optimal rule*/
void RelatedGenomes::SuboptimalRule1( ){
  for( int i = 0; i < (int)S.genes.size( ); i++ ){
    //Find the left singleton pair
    if( (int)S.genes[i].BlastHits.size( ) != 1 )
      continue;
   
    int leftS = i;
    int leftT = S.genes[i].BlastHits[0].ObjectIndex;
    int rightS = 0;
    if( (int)T.genes[leftT].BlastHits.size( ) != 1 )
      continue;
    int ChromidS = S.genes[leftS].ChromId;
    int ChromidT = T.genes[leftT].ChromId;

    if( S.genes[leftS].Sign == T.genes[leftT].Sign ){ 
      int dS = leftS+1;
      int dT = leftT+1;
      if( dS >= (int)S.genes.size( ) || dT >= (int)T.genes.size( ) )
         continue;
      int flag = 0;
      
      if( (int)S.genes[dS].BlastHits.size( ) <= 1 || (int)S.genes[dS].BlastHits.size( ) <= 1 )
	continue;
      while(  dS < (int)S.genes.size( ) && dT<(int)T.genes.size( ) && ( (int)S.genes[dS].BlastHits.size( ) > 1 || (int)T.genes[dT].BlastHits.size( ) > 1) ){ 
	if( S.genes[dS].isBlastHit(dT) && S.genes[dS].Sign == T.genes[dT].Sign && S.genes[dS].ChromId == ChromidS && T.genes[dT].ChromId == ChromidT){
	  dS++;
	  dT++;	  
	}
	else{
	  flag = 1;
	  break;
	}
      }
      if( flag == 1 )
	continue;
      if( (int)S.genes[dS].BlastHits.size( ) != 1 || (int)T.genes[dT].BlastHits.size( ) != 1)
	continue;
      if( !S.genes[dS].isBlastHit(dT))
	continue;
      if( S.genes[dS].ChromId != ChromidS || T.genes[dT].ChromId != ChromidT )
	continue;
      if(  S.genes[dS].Sign != T.genes[dT].Sign )
         continue;
      //update matches and modified
      rightS = dS;
      //int rightT = dT;
      
      for(dS = leftS, dT = leftT; dS<= rightS; dS++, dT++ ){
         
	if( matchesS[dS] == dT )
	  continue;
	if( matchesS[dS]!=-1 ){
	  cout<<"ERROR@RELATEDGENOMES@SUBOPTIMALRULE1"<<endl;
	  exit(1);
	}
        cout <<dS<<"---"<<dT<<endl;
	matchesS[dS] = dT;
	matchesT[dT] = dS;
        S.genes[dS].Match = dT;
        T.genes[dT].Match = dS;
	modified.push_back( dS );
	selectBlastHit( dS, dT );
      }

    }
    else{
      int dS = leftS+1;
      int dT = leftT-1;
      if( dS >= (int)S.genes.size( ) || dT < 0 )
         continue;      
      int flag = 0;
      if( (int)S.genes[dS].BlastHits.size( ) <= 1 || (int)S.genes[dS].BlastHits.size( ) <= 1 )
	continue;
      while( dS < (int)S.genes.size( ) && dT >= 0 && ((int)S.genes[dS].BlastHits.size( ) > 1 || (int)T.genes[dT].BlastHits.size( ) > 1 )){ 
	if( S.genes[dS].isBlastHit(dT) && S.genes[dS].Sign != T.genes[dT].Sign && S.genes[dS].ChromId == ChromidS && T.genes[dT].ChromId == ChromidT){
	  dS++;
	  dT--;	  
	}
	else{
	  flag = 1;
	  break;
	}
      }
      if( flag == 1 )
	continue;
      if( (int)S.genes[dS].BlastHits.size( ) != 1 || (int)T.genes[dT].BlastHits.size( ) != 1)
	continue;
      if( !S.genes[dS].isBlastHit(dT))
	continue;
      if( S.genes[dS].ChromId != ChromidS || T.genes[dT].ChromId != ChromidT )
	continue;
      if(  S.genes[dS].Sign == T.genes[dT].Sign )
         continue;
      //update matches and modified
      rightS = dS;
      //int rightT = dT;
      for(dS = leftS, dT = leftT; dS<= rightS; dS++, dT-- ){
	if( matchesS[dS] == dT )
	  continue;
	if( matchesS[dS]!=-1 ){
	  cout<<"ERROR@RELATEDGENOMES@SUBOPTIMALRULE1"<<endl;
	  exit(1);
	}
        cout <<dS<<"---"<<dT<<endl;
	matchesS[dS] = dT;
	matchesT[dT] = dS;
        S.genes[dS].Match = dT;
        T.genes[dT].Match = dS;
	modified.push_back( dS );
	selectBlastHit( dS, dT );	
      }
    }
    i = rightS - 1;
  } 
}

/*Second sub-optimal rule*/
void RelatedGenomes::SuboptimalRule2( ){
   
   list<int> endA;
   list<int> endB;
   EndPair tmppair;
   list<EndPair> endpairs;
   vector<int> duplicatesA1;
   vector<int> duplicatesB1;
   vector<int> duplicatesA2;
   vector<int> duplicatesB2;
   //list<int>::iterator it;
   int flag ;
     
   for( int i = 0; i < (int)S.genes.size( ); i++ ){

      endA.clear();
      endB.clear();
      duplicatesA1.clear();
      duplicatesA2.clear();
      duplicatesB1.clear();
      duplicatesB2.clear();
      endpairs.clear();
      if( (int)S.genes[i].BlastHits.size( ) != 1 )
         continue;
      int leftS = i;
      int leftT = S.genes[i].BlastHits[0].ObjectIndex;
      if( (int)T.genes[leftT].BlastHits.size( ) != 1 )
         continue;
      int ChromidS = S.genes[leftS].ChromId;
      int ChromidT = T.genes[leftT].ChromId;     

      if( S.genes[leftS].Sign == T.genes[leftT].Sign ){
         int dS = leftS + 1;
         int dT = leftT + 1;
         if( dS >= (int)S.genes.size( ) || dT >= (int)T.genes.size( ) )
            continue;
         int gap = dS - dT;
         if( S.genes[dS].ChromId != ChromidS || T.genes[dT].ChromId != ChromidT)
            continue;
         if(( (int)S.genes[dS].BlastHits.size( ) > 1 || (int)T.genes[dT].BlastHits.size( ) > 1 )&& S.genes[dS].ChromId == ChromidS && T.genes[dT].ChromId == ChromidT){
            for( int ii = 0; ii < (int)S.genes[dS].BlastHits.size( ); ii ++ ){
               int j = S.genes[dS].BlastHits[ii].ObjectIndex;
               if( S.genes[dS].Sign != T.genes[j].Sign && T.genes[j].ChromId == ChromidT )
                  endB.push_back( j );
            }
            for( int ii = 0; ii < (int)T.genes[dT].BlastHits.size( ); ii ++ ){
               int j = T.genes[dT].BlastHits[ii].ObjectIndex;
               if( T.genes[dT].Sign != S.genes[j].Sign && S.genes[j].ChromId == ChromidS )
                  endA.push_back( j );
            }            
         }
         // sorting two lists to ascending order
         endA.sort( );
         endB.sort( );
         // find endpair candidates
         list<int>::iterator itA;
         list<int>::iterator itB;
         for( itA = endA.begin( ); itA != endA.end( ); ++itA ){
            flag = 0;
            for( itB = endB.begin( ); itB != endB.end( ); ++itB ){
               if( *itA - *itB == gap ) {
                  tmppair.A = *itA;
                  tmppair.B = *itB;
                  endpairs.push_front( tmppair );
                  flag = 1;
                  break;
               }            
            }
         }
         
         //find the matching duplicates block
         list<EndPair>::iterator itEP;
         for( itEP = endpairs.begin( ); itEP != endpairs.end( ); ++itEP ) {
            int NA = dS;
            int NB = (*itEP).B;
            
            duplicatesA1.push_back( NA );
            duplicatesB1.push_back( NB );
            
            while( NB > dT  ) {
               flag = 0;
               NA ++;
               NB --;
               if( !S.genes[NA].isBlastHit(NB) || ((int)S.genes[NA].BlastHits.size( ) <= 1 && (int)T.genes[NB].BlastHits.size( ) <= 1)  || S.genes[NA].Sign == T.genes[NB].Sign){
                  duplicatesA1.clear( );
                  duplicatesB1.clear( );
                  break;
               }
               duplicatesA1.push_back( NA );
               duplicatesB1.push_back( NB );
            }
            if( (int)duplicatesA1.size( ) > 0 )
               break;
         }

         for( int ii = 0; ii < (int)duplicatesA1.size( ); ii ++){
            
            int m = duplicatesA1[ii];
            int n = duplicatesB1[ii];
            if( matchesS[m] == n )
               continue;
            if( matchesS[m]!=-1 ){
               cout<<"ERROR@RELATEDGENOMES@SUBOPTIMALRULE2"<<endl;
               exit(1);
            } 
            cout<<m<<"-----"<<n<<endl;
            matchesS[m] = n;
            matchesT[n] = m;
            S.genes[m].Match = n;
            T.genes[n].Match = m;
            selectBlastHit( m,n );
            modified.push_back( m );
         }
         endA.clear( );
         endB.clear( );
         endpairs.clear( );          

         // check backword direction
         dS = leftS - 1;
         dT = leftT - 1;
         if( dS < 0 || dT < 0 )
            continue;
         gap = dS - dT;
         if( S.genes[dS].ChromId != ChromidS || T.genes[dT].ChromId != ChromidT)
            continue;
         if( ((int)S.genes[dS].BlastHits.size( ) > 1 || (int)T.genes[dT].BlastHits.size( ) > 1) && S.genes[dS].ChromId == ChromidS && T.genes[dT].ChromId == ChromidT){
            for( int ii = 0; ii < (int)S.genes[dS].BlastHits.size( ); ii ++ ){
               int j = S.genes[dS].BlastHits[ii].ObjectIndex;
               if( S.genes[dS].Sign != T.genes[j].Sign && T.genes[j].ChromId == ChromidT )
                  endB.push_back( j );
            }
            for( int ii = 0; ii < (int)T.genes[dT].BlastHits.size( ); ii ++ ){
               int j = T.genes[dT].BlastHits[ii].ObjectIndex;
               if( T.genes[dT].Sign != S.genes[j].Sign && S.genes[j].ChromId == ChromidS )
                  endA.push_back( j );
            }            
         }

         // sorting two lists to ascending order
         endA.sort( );
         endB.sort( );
         // find endpair candidates
         for( itA = endA.begin( ); itA != endA.end( ); ++itA ){
            flag = 0;
            for( itB = endB.begin( ); itB != endB.end( ); ++itB ){
               if( *itA - *itB == gap ) {
                  tmppair.A = *itA;
                  tmppair.B = *itB;
                  endpairs.push_back( tmppair );
                  flag = 1;
                  break;
               }            
            }
         }
         //find the matching duplicates block
         for( itEP = endpairs.begin( ); itEP != endpairs.end( ); ++itEP ) {
            int NA = dS;
            int NB = (*itEP).B;
            
            duplicatesA2.push_back( NA );
            duplicatesB2.push_back( NB );
            
            while( NB < dT  ) {
               flag = 0;
               NA --;
               NB ++;
               if( !S.genes[NA].isBlastHit(NB) || ((int)S.genes[NA].BlastHits.size( ) <= 1 && (int)T.genes[NB].BlastHits.size( ) <= 1)  || S.genes[NA].Sign == T.genes[NB].Sign){
                  duplicatesA2.clear( );
                  duplicatesB2.clear( );
                  break;
               }
               duplicatesA2.push_back( NA );
               duplicatesB2.push_back( NB );
            }
            if( (int)duplicatesA2.size( ) > 0 )
               break;
         }

         for( int ii = 0; ii < (int)duplicatesA2.size( ); ii ++){
            
            int m = duplicatesA2[ii];
            int n = duplicatesB2[ii];
            if( matchesS[m] == n )
               continue;
            if( matchesS[m]!=-1 ){
               cout<<"ERROR@RELATEDGENOMES@SUBOPTIMALRULE2"<<endl;
               exit(1);
            } 
            cout<<m<<"-----"<<n<<endl;
            matchesS[m] = n;
            matchesT[n] = m;
            S.genes[m].Match = n;
            T.genes[n].Match = m;
            selectBlastHit( m,n );
            modified.push_back( m );
         }
         endA.clear( );
         endB.clear( );
         endpairs.clear( );       
      }     
    
      else{
         int dS = leftS + 1;
         int dT = leftT - 1;
         if( dS >= (int)S.genes.size( ) || dT < 0 )
            continue;
         int gap = dT + dS;
         if( S.genes[dS].ChromId != ChromidS || T.genes[dT].ChromId != ChromidT)
            continue;
         if(( (int)S.genes[dS].BlastHits.size( ) > 1 || (int)T.genes[dT].BlastHits.size( ) > 1) && S.genes[dS].ChromId == ChromidS && T.genes[dT].ChromId == ChromidT){
            for( int ii = 0; ii < (int)S.genes[dS].BlastHits.size( ); ii ++ ){
               int j = S.genes[dS].BlastHits[ii].ObjectIndex;
               if( S.genes[dS].Sign == T.genes[j].Sign && T.genes[j].ChromId == ChromidT )
                  endB.push_back( j );
            }
            for( int ii = 0; ii < (int)T.genes[dT].BlastHits.size( ); ii ++ ){
               int j = T.genes[dT].BlastHits[ii].ObjectIndex;
               if( T.genes[dT].Sign == S.genes[j].Sign && S.genes[j].ChromId == ChromidS )
                  endA.push_back( j );
            }            
         }

         // sorting two lists to ascending order
         endA.sort( );
         endB.sort( );
         // find endpair candidates
         list<int>::iterator itA;
         list<int>::iterator itB;
         for( itA = endA.begin( ); itA != endA.end( ); ++itA ){
            flag = 0;
            for( itB = endB.begin( ); itB != endB.end( ); ++itB ){
               if( *itA + *itB == gap ) {
                  tmppair.A = *itA;
                  tmppair.B = *itB;
                  endpairs.push_front( tmppair );
                  flag = 1;
                  break;
               }            
            }
         }

         //find the matching duplicates block
         list<EndPair>::iterator itEP;
         for( itEP = endpairs.begin( ); itEP != endpairs.end( ); ++itEP ) {
            int NA = dS;
            int NB = (*itEP).B;
            
            duplicatesA1.push_back( NA );
            duplicatesB1.push_back( NB );
            flag = 1;
            while( NB < dT  ) {
               flag = 0;
               NA ++;
               NB ++;
               if( !S.genes[NA].isBlastHit(NB) || ((int)S.genes[NA].BlastHits.size( ) <= 1 && (int)T.genes[NB].BlastHits.size( ) <= 1 ) || S.genes[NA].Sign != T.genes[NB].Sign){
                  duplicatesA1.clear( );
                  duplicatesB1.clear( );
                  break;
               }
               duplicatesA1.push_back( NA );
               duplicatesB1.push_back( NB );
            }
            if( (int)duplicatesA1.size( ) > 0 )
               break;
         }

         for( int ii = 0; ii < (int)duplicatesA1.size( ); ii ++){
            
            int m = duplicatesA1[ii];
            int n = duplicatesB1[ii];
            if( matchesS[m] == n )
               continue;
            if( matchesS[m]!=-1 ){
               cout<<"ERROR@RELATEDGENOMES@SUBOPTIMALRULE2"<<endl;
               exit(1);
            } 
             cout<<m<<"------"<<n<<endl;
            matchesS[m] = n;
            matchesT[n] = m;
            S.genes[m].Match = n;
            T.genes[n].Match = m;
            selectBlastHit( m,n );
            modified.push_back( m );
         }

         endA.clear( );
         endB.clear( );
         endpairs.clear( );
       
         //check backword direction
         dS = leftS - 1;
         dT = leftT + 1;
         if( dT >= (int)T.genes.size( ) || dS < 0 )
            continue;
         gap = dT + dS;
         if( S.genes[dS].ChromId != ChromidS || T.genes[dT].ChromId != ChromidT)
            continue;
         if( ((int)S.genes[dS].BlastHits.size( ) > 1 || (int)T.genes[dT].BlastHits.size( ) > 1) && S.genes[dS].ChromId == ChromidS && T.genes[dT].ChromId == ChromidT){
            for( int ii = 0; ii < (int)S.genes[dS].BlastHits.size( ); ii ++ ){
               int j = S.genes[dS].BlastHits[ii].ObjectIndex;
               if( S.genes[dS].Sign == T.genes[j].Sign && T.genes[j].ChromId == ChromidT )
                  endB.push_back( j );
            }
            for( int ii = 0; ii < (int)T.genes[dT].BlastHits.size( ); ii ++ ){
               int j = T.genes[dT].BlastHits[ii].ObjectIndex;
               if( T.genes[dT].Sign == S.genes[j].Sign && S.genes[j].ChromId == ChromidS )
                  endA.push_back( j );
            }            
         }
         
         // sorting two lists to ascending order
         endA.sort( );
         endB.sort( );
         // find endpair candidates
         for( itA = endA.begin( ); itA != endA.end( ); ++itA ){
            flag = 0;
            for( itB = endB.begin( ); itB != endB.end( ); ++itB ){
               if( *itA + *itB == gap ) {
                  tmppair.A = *itA;
                  tmppair.B = *itB;
                  endpairs.push_back( tmppair );
                  flag = 1;
                  break;
               }            
            }
         }
         
         //find the matching duplicates block
         for( itEP = endpairs.begin( ); itEP != endpairs.end( ); ++itEP ) {
            int NA = dS;
            int NB = (*itEP).B;
            
            duplicatesA2.push_back( NA );
            duplicatesB2.push_back( NB );
            flag = 1;
            while( NB > dT  ) {
               flag = 0;
               NA --;
               NB --;
               if( !S.genes[NA].isBlastHit(NB) || ((int)S.genes[NA].BlastHits.size( ) <= 1 && (int)T.genes[NB].BlastHits.size( ) <= 1)  || S.genes[NA].Sign != T.genes[NB].Sign){
                  duplicatesA2.clear( );
                  duplicatesB2.clear( );
                  break;
               }
               duplicatesA2.push_back( NA );
               duplicatesB2.push_back( NB );
            }
            if( (int)duplicatesA2.size( ) > 0 )
               break;
         }

         for( int ii = 0; ii < (int)duplicatesA2.size( ); ii ++){
      
            int m = duplicatesA2[ii];
            int n = duplicatesB2[ii];
            if( matchesS[m] == n )
               continue;
            if( matchesS[m]!=-1 ){
               cout<<"ERROR@RELATEDGENOMES@SUBOPTIMALRULE2"<<endl;
               exit(1);
            } 
            cout<<m<<"------"<<n<<endl;
            matchesS[m] = n;
            matchesT[n] = m;
            S.genes[m].Match = n;
            T.genes[n].Match = m;
            selectBlastHit( m,n );
            modified.push_back( m );
         }         
         endA.clear( );
         endB.clear( );
         endpairs.clear( );         
      }   
  }
}

/*Third sub-optimal rule*/
void RelatedGenomes::SuboptimalRule3( ){
 for( int i = 0; i < (int)S.genes.size( ); i++ ){
    //Find the left singleton pair
    if( (int)S.genes[i].BlastHits.size( ) != 1 )
      continue;
   
    int leftS = i;
    int leftT = S.genes[i].BlastHits[0].ObjectIndex;
    int rightS = 0;
    if( (int)T.genes[leftT].BlastHits.size( ) != 1 )
      continue;
    int ChromidS = S.genes[leftS].ChromId;
    int ChromidT = T.genes[leftT].ChromId;

    if( S.genes[leftS].Sign == T.genes[leftT].Sign ){ 
      int dS = leftS+1;
      int dT = leftT+1;
            
      while(  dS < (int)S.genes.size( ) && dT < (int)T.genes.size( ) && ((int)S.genes[dS].BlastHits.size( ) > 1 || (int)T.genes[dT].BlastHits.size( ) > 1 )&&S.genes[dS].isBlastHit(dT) && S.genes[dS].Sign == T.genes[dT].Sign && S.genes[dS].ChromId == ChromidS && T.genes[dT].ChromId == ChromidT){
	dS++;
	dT++;	  
      }
     
      if( dS == leftS+1 )
	continue;
     
      //update matches and modified
      rightS = dS - 1;
      //int rightT = dT;
      
      for(dS = leftS, dT = leftT; dS<= rightS; dS++, dT++ ){
	if( matchesS[dS] == dT )
	  continue;
	if( matchesS[dS]!=-1 ){
	  cout<<"ERROR@RELATEDGENOMES@SUBOPTIMALRULE3"<<endl;
	  exit(1);
	}
        cout <<dS<<"---"<<dT<<endl;
	matchesS[dS] = dT;
	matchesT[dT] = dS;
        S.genes[dS].Match = dT;
        T.genes[dT].Match = dS;
	modified.push_back( dS );
	selectBlastHit( dS, dT );
      }

      dS = leftS-1;
      dT = leftT-1;
      while(  dS >= 0 && dT >= 0 &&( (int)S.genes[dS].BlastHits.size( ) > 1 && (int)T.genes[dT].BlastHits.size( ) > 1) &&S.genes[dS].isBlastHit(dT) && S.genes[dS].Sign == T.genes[dT].Sign && S.genes[dS].ChromId == ChromidS && T.genes[dT].ChromId == ChromidT){
         dS--;
         dT--;	  
      }
     
      if( dS == leftS-1 )
         continue;
     
      //update matches and modified
      rightS = dS + 1;
      //int rightT = dT;
      
      for(dS = leftS, dT = leftT; dS>= rightS; dS--, dT-- ){
         if( matchesS[dS] == dT )
	  continue;
	if( matchesS[dS]!=-1 ){
	  cout<<"ERROR@RELATEDGENOMES@SUBOPTIMALRULE3"<<endl;
	  exit(1);
	}
          cout <<dS<<"---"<<dT<<endl;
	matchesS[dS] = dT;
	matchesT[dT] = dS;
        S.genes[dS].Match = dT;
        T.genes[dT].Match = dS;
	modified.push_back( dS );
	selectBlastHit( dS, dT );
      }
    }      
    
    else{
      int dS = leftS+1;
      int dT = leftT-1;
            
      while(  dS < (int)S.genes.size( ) && dT >= 0 &&( (int)S.genes[dS].BlastHits.size( ) > 1 || (int)T.genes[dT].BlastHits.size( ) > 1) && S.genes[dS].isBlastHit(dT) && S.genes[dS].Sign != T.genes[dT].Sign && S.genes[dS].ChromId == ChromidS && T.genes[dT].ChromId == ChromidT){
	dS++;
	dT--;	  
      }
     
      if( dS == leftS+1 )
	continue;
     
      //update matches and modified
      rightS = dS - 1;
      //int rightT = dT;
      
      for(dS = leftS, dT = leftT; dS<= rightS; dS++, dT-- ){
	if( matchesS[dS] == dT )
	  continue;
	if( matchesS[dS]!=-1 ){
	  cout<<"ERROR@RELATEDGENOMES@SUBOPTIMALRULE3"<<endl;
	  exit(1);
	}
//          cout <<dS<<"---"<<dT<<endl;
	matchesS[dS] = dT;
	matchesT[dT] = dS;
        S.genes[dS].Match = dT;
        T.genes[dT].Match = dS;
	modified.push_back( dS );
	selectBlastHit( dS, dT );
      }

      dS = leftS-1;
      dT = leftT+1;
            
      while(  dS >= 0 && dT < (int)T.genes.size( ) &&( (int)S.genes[dS].BlastHits.size( ) > 1 || (int)T.genes[dT].BlastHits.size( ) > 1) && S.genes[dS].isBlastHit(dT) && S.genes[dS].Sign != T.genes[dT].Sign && S.genes[dS].ChromId == ChromidS && T.genes[dT].ChromId == ChromidT){
	dS--;
	dT++;	  
      }
     
      if( dS == leftS-1 )
	continue;
     
      //update matches and modified
      rightS = dS + 1;
      //int rightT = dT;
      
      for(dS = leftS, dT = leftT; dS>= rightS; dS--, dT++ ){
	if( matchesS[dS] == dT )
	  continue;
	if( matchesS[dS]!=-1 ){
	  cout<<"ERROR@RELATEDGENOMES@SUBOPTIMALRULE3"<<endl;
	  exit(1);
	}
          cout <<dS<<"---"<<dT<<endl;
	matchesS[dS] = dT;
	matchesT[dT] = dS;
        S.genes[dS].Match = dT;
        T.genes[dT].Match = dS;
	modified.push_back( dS );
	selectBlastHit( dS, dT );
      }      
    }
    i = rightS - 1;
  } 
}


void RelatedGenomes::print( ){
   cout<<"------Genome S"<<endl;
   S.print( );
   cout<<"------Genome T"<<endl;
   T.print( );
   cout<<"------CompleteGraph"<<endl;
  G.print( );
  //cout<<"------PairmatchGraph"<<endl;
  //P.Print( );
   cout<<"------MatchesS"<<endl;
   for( int i = 0; i<(int)matchesS.size( ); i++ )
      cout << i << ": "<< matchesS[i]<<endl;
   cout<<"------MatchesT"<<endl;
   for( int i = 0; i<(int)matchesT.size( ); i++ )
      cout << i << ": "<< matchesT[i]<<endl;
  //cout<<"------Modified"<<endl;
  //for( int i = 0; i<(int)modified.size( ); i++ )
  // cout << modified[i]<<endl;

}

/*Print main orthologs assigned by MSOAR to file "outfile"*/
void RelatedGenomes::printresult( string outfile ){
   ofstream fout( outfile.c_str( ) );
   for( int i = 0; i < (int)matchesS.size( ); i++ ){
      int j = matchesS[i];
      if( j == -1 )
         continue;

      fout << S.genes[i].Accession_No<<"\t"<<S.genes[i].GeneSymbol<<"\t"<<T.genes[j].Accession_No<<"\t"<<T.genes[j].GeneSymbol<<endl;
   }
}

/*Check if two strings are the same discarding the case*/
bool RelatedGenomes::isSame( const char* s1, const char* s2 ){
   if( strlen( s1 ) != strlen( s2 )  )
      return false;
   int i;
   for( i = 0; i < (int)strlen( s1 ); i++ ){
      if( s1[i] != s2[i] && s1[i] != (s2[i] - 'a' + 'A') && s2[i] != (s1[i] - 'a' + 'A') )
         break;
   }
   if( i ==  (int)strlen( s1 ))
      return true;
   return false;
}

/*update related genomes based on MatchesS and MatchesT*/
void RelatedGenomes::Update( ){

   //Build two hashes for matched genes
   vector<int>hashS;
   vector<int>hashT;
   vector<int>tmp1;
   vector<int>tmp2;
   int index = 0;
   for( int i = 0; i < (int)matchesS.size( ); i++ ){
      int j = matchesS[i];
      if( j == -1 )
         hashS.push_back( -1 );
      else{
         hashS.push_back( index++ );
      }
   }
   index = 0;
   for( int i = 0; i < (int)matchesT.size( ); i++ ){
      int j = matchesT[i];
      if( j == -1 )
         hashT.push_back( -1 );
      else
         hashT.push_back( index++ );         
   }

   vector<Gene>A;
   A.clear( );
   for( int i = 0; i < (int)matchesS.size( ); i++ ){
      int j = matchesS[i];
      if( j == -1 )
         continue;
      for( int k = 0; k < (int)S.genes[i].BlastHits.size( ); k++ ){
         int obj = S.genes[i].BlastHits[k].ObjectIndex;
         int tmp =  hashT[obj];
         if( tmp  == -1 ){
            S.genes[i].deleteBlastHit( S.genes[i].BlastHits[k].ObjectIndex );
            k --;
         }
         else
            S.genes[i].BlastHits[k].ObjectIndex = tmp;
      }
      S.genes[i].Match = hashT[j];
      A.push_back( S.genes[i] );
      tmp1.push_back( S.genes[i].Match );
   }
 
   vector<Gene>B;
   B.clear( );
   for( int i = 0; i < (int)matchesT.size( ); i++ ){
      int j = matchesT[i];
      if( j == -1 )
         continue;
      for( int k = 0; k < (int)T.genes[i].BlastHits.size( ); k++ ){
         int tmp =  hashS[T.genes[i].BlastHits[k].ObjectIndex];
         if( tmp  == -1 ){
            T.genes[i].deleteBlastHit( T.genes[i].BlastHits[k].ObjectIndex );
            k --;
         }
         else
            T.genes[i].BlastHits[k].ObjectIndex = tmp;
      }
      T.genes[i].Match = hashS[j];
      B.push_back( T.genes[i] );
      tmp2.push_back( T.genes[i].Match );
   }
  
   matchesS.clear( );
   matchesT.clear( );
   for( int i = 0; i< (int)tmp1.size( ); i++ )
      matchesS.push_back( tmp1[i] );
   for( int i = 0; i< (int)tmp2.size( ); i++ )
      matchesT.push_back( tmp2[i] );


   //Construct Gene families
   GeneFamiliesConstruction( A, B );
  
   for( int i = 0; i < (int)A.size( ); i++ )
      A[i].selectBlastHit( A[i].Match );
   
   for( int i = 0; i < (int)B.size( ); i++ )
      B[i].selectBlastHit( B[i].Match );
   
   S.genes.clear( );
   S.chromsizes.clear( );
   for( int i = 0; i < (int)A.size( ); i++ ) 
      S.genes.push_back( A[i] );
   S.getChromSizes( );
   T.genes.clear( );
   T.chromsizes.clear( );
   for( int i = 0; i < (int)B.size( ); i++ ) 
      T.genes.push_back( B[i] );
   T.getChromSizes( );
   cout<<"Print out new Genome S"<<endl;
   S.print( );
   cout<<"Print out new Genome T"<<endl;
   T.print( ); 
   cout<<"Construct complete graph G"<<endl;
   G.ClearGraph( );
   G.ConstructGraph( S,T );
   cout<<"Print out the new completegraph"<<endl;                                                    
   G.print( );
   cout<<"Print out the gene families"<<endl;
   for( int i = 0; i < (int)GFs.size( ); i++ ){
      cout<<"GeneFamily"<<i<<":"<<endl;
      for( int j = 0; j < (int)GFs[i].size( ); j++){
         cout<<"\t"<<GFs[i][j].A<<"---"<<GFs[i][j].B<<endl;
      }
   }
   
}

/*Construct gene families after assigning main orthologs*/
void RelatedGenomes::GeneFamiliesConstruction( vector<Gene>A, vector<Gene>B ){
   GFs.clear( );
   
   if( A.size( ) != B.size( ) ){
      cout<<"Error@RELATEDGENOMES@GENEFAMILIESCONSTRUCTION"<<endl;
      exit( 1 );
   }
   
   vector<int>flag( (int)A.size( ), 0 );
   for( int i = 0; i < (int)A.size( ); i++ ){
      if( flag[i] != 0 )
         continue; //already visited

      queue<int>Q;//id Queue
      GeneFamily gf;
      EndPair ep;
      Q.push( i );
      flag[i] = 1;
      while( !Q.empty( ) ){
         int t = Q.front( );
         Q.pop( );
         ep.A = t;
         ep.B = A[t].Match;
         gf.push_back( ep );
         for( int j = 0; j < (int)A[t].BlastHits.size( ); j++ ){
            int index = A[t].BlastHits[j].ObjectIndex;
            if( index == A[t].Match || flag[B[index].Match] != 0)
               continue;
            Q.push( B[index].Match );
            flag[B[index].Match] = 1;
         }
         for( int j = 0; j < (int)B[ep.B].BlastHits.size( ); j++ ){
            int index = B[ep.B].BlastHits[j].ObjectIndex;
            if( index == t || flag[ index ] != 0 )
               continue;
            Q.push( index );
            flag[index] = 1;
         }
         
      }

      GFs.push_back( gf );
   }


   
}

/*Detect paralogs after main ortholog detection.            */
/*The core idea is to minimize reversal/duplication distance*/
void RelatedGenomes::ParalogsDetection( ){
   Paralogs.clear( );
   int old_dis = G.ApproximateDistance( );
   int new_dis = 0;
   for( int i = 0; i < (int)GFs.size( ); i++ ){
      if( (int)GFs[i].size( ) == 1 )
         continue;
      for( int j = 0; j < (int)GFs[i].size( ); j++ ){
         CompleteGraph G_copy = G;
         new_dis = G_copy.ApproximateDistance( GFs[i][j].A, GFs[i][j].B );
         cout <<"old_dis="<<old_dis<<"\t new_dis="<<new_dis<<"\t";
         cout<<"A:"<< GFs[i][j].A<<"\t B:"<< GFs[i][j].B<<endl;

         if( S.genes[GFs[i][j].A].Sign == T.genes[GFs[i][j].B].Sign ){
            if( old_dis - new_dis > WEIGHT_PARALOG ){
               EndPair ep;
               ep.A = GFs[i][j].A;
               ep.B = GFs[i][j].B;
               Paralogs.push_back( ep );
               matchesS[GFs[i][j].A] = -1;
               matchesT[GFs[i][j].B] = -1;
               S.genes[GFs[i][j].A].Match = -1;
               T.genes[GFs[i][j].B].Match = -1;
               G.deletePair( GFs[i][j].A, GFs[i][j].B );
               GFs[i].erase( GFs[i].begin( ) + j );
               j -- ;
               old_dis = new_dis;
            }
         }
         else{
            if( old_dis - new_dis >= WEIGHT_PARALOG ){
               EndPair ep;
               ep.A = GFs[i][j].A;
               ep.B = GFs[i][j].B;
               Paralogs.push_back( ep );
               matchesS[GFs[i][j].A] = -1;
               matchesT[GFs[i][j].B] = -1;
               S.genes[GFs[i][j].A].Match = -1;
               T.genes[GFs[i][j].B].Match = -1;
               G.deletePair( GFs[i][j].A, GFs[i][j].B );
               GFs[i].erase( GFs[i].begin( ) + j );
               j -- ;
               old_dis = new_dis;
            } 
          
         }
      }
   }
}

/*Print out the paralogs detected by MSOAR*/
void RelatedGenomes::ParalogsPrint( ){
   int TP = 0;
   int total = 0;
   cout <<"Paralog size = "<<Paralogs.size( )<<endl;
   for( int t = 0; t < (int)Paralogs.size( ); t++ ){
      int i = Paralogs[t].A;
      int j = Paralogs[t].B;
     
      total ++;
      cout << i <<"\t"<< j<<endl;
      cout << S.genes[i].GeneSymbol<<"\t"<< T.genes[j].GeneSymbol;
      if( isSame( S.genes[i].GeneSymbol.c_str(  ),T.genes[j].GeneSymbol.c_str() ) ){
         cout<<"\t  @";
         TP ++;
      }
      cout<<endl;
      cout << "------------------------"<<endl; 
   }
   cout <<"TP = " <<TP<<endl;
   cout <<"total= "<<total<<endl;
}

void RelatedGenomes::OrthologCandidate( ){
   for( int i = 0; i < (int)S.genes.size( ); i ++ ){
	for( int j = 0; j< (int)S.genes[i].BlastHits.size( ); j ++ ){
            cout<< S.genes[i].Accession_No <<"\t";
	    int k = S.genes[i].BlastHits[j].ObjectIndex;
	    cout << T.genes[k].Accession_No << "\n";
	}
   }   
}
