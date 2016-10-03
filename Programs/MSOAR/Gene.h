/*********************************************************************/
/*Object: gene 							     */
/*Each gene has a vector of BlastHits, as well as official gene      */
/*symbole, Accession_noChromosome id,     */
/*position in that chromosome, strand and match id. The match id is  */
/*the id of its orthologs.*/
/*********************************************************************/

#ifndef _GENE_H_
#define _GENE_H_

#include <string>
#include <vector>
#include <fstream>
#include "BlastHit.h"
using namespace std;

class Gene
{
 private: 
   
 public:
   
   string GeneSymbol;              // gene official symbol
   string Accession_No;            // gene accession no.
   int ChromId;                   // which chromosome this gene locates in Genome
   int PosonChrom;                // position in the chromosome of Genome
   int Sign;                      // the oritaton of the gene
   int Match;                     // the id of the match
   vector<BlastHit>BlastHits;     // blastp hits
  
   Gene( ){
      GeneSymbol="";
      Accession_No = "";
      ChromId = 0;
      PosonChrom = 0;
      Sign = 0;
      BlastHits.empty( );
      Match = -1;
   }
   
   Gene( const Gene& g ){
      GeneSymbol = g.GeneSymbol;
      Accession_No = g.Accession_No;
      ChromId = g.ChromId;
      PosonChrom = g.PosonChrom;
      Sign = g.Sign;
      BlastHits = g.BlastHits;
      Match = g.Match;  
   }

   Gene( string s, string a, int c, int p, int sign, int m, vector<BlastHit>h ){
      GeneSymbol = a;
      Accession_No = s;
      ChromId = c;
      PosonChrom = p;
      Sign = sign;
      BlastHits = h;
      Match = m;
   }
  
   /*Check if the gene "id" in another genome is one of the Blast Hits of this gene*/
   bool isBlastHit( int id ){
      for( int i = 0; i < (int)BlastHits.size( ); i ++){
         if( BlastHits[i].ObjectIndex == id )
            return true;
      }
      return false;
   }
   
   /*Delete one of the blast Hits*/
   void deleteBlastHit( int id ){
      for( int i = 0; i < (int)BlastHits.size( ); i ++){
         if( BlastHits[i].ObjectIndex == id ){
            BlastHits.erase( BlastHits.begin( ) + i );
            return;
         }
      }
      cout<<"ERROR@GENE@DELETEBLASTBLASTHIT id="<<id<<endl;
      print( );
      exit(1);
   }
   
   /*Selete one of the Blast Hits to be match and delete all the others*/
   void selectBlastHit( int id ){
      BlastHit tmp;
      for( int i = 0; i < (int)BlastHits.size( ); i ++){
         if( BlastHits[i].ObjectIndex == id ){
            tmp = BlastHits[i];
            break;
         }
      }
      if( tmp.ObjectIndex == -1 ){
         cout<<"ERROR@GENE@SELECTBLASTHIT id="<<id<<endl;
         print( );
         exit(1);
      }
      BlastHits.clear( );
      BlastHits.push_back( tmp );
   }

   double getBestEvalue( ){
      double min = 1.0;
      for( int i = 0; i< (int)BlastHits.size( ); i ++){
         if( BlastHits[i].Evalue <= min )
            min = BlastHits[i].Evalue;                 
      }
      return min;
   }

   double getBestBitScore( ){
      double max = 0.0;
      for( int i = 0; i< (int)BlastHits.size( ); i ++){
         if( BlastHits[i].BitScore >= max )
            max = BlastHits[i].BitScore;                 
      }
      return max;
   }

   const Gene& operator = ( const Gene& g ){
      GeneSymbol = g.GeneSymbol;
      Accession_No = g.Accession_No;
      ChromId = g.ChromId;
      PosonChrom = g.PosonChrom;
      Sign = g.Sign;
      BlastHits = g.BlastHits;
      Match = g.Match;
      return *this;
   }
   
   void print( ) {
      cout << GeneSymbol<< "\t";
      cout << Accession_No<< "\t";
      cout << ChromId << "\t";
      cout << PosonChrom << "\t"; 
      cout << Sign << "\t" ;
      cout << Match << "\n";
      for( int i = 0; i < (int)BlastHits.size( ); i++ )
         BlastHits[i].print( );      
   }
};

#endif
