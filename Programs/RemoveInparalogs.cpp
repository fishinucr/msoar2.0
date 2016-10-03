#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <vector>
#include <iomanip>
#include <sstream>

using namespace std;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		cout << "Usage: removeInparalogs <TAGs> <G1G2.blastp>" << endl;
		exit(1);
	}
	
	ifstream infile1(argv[1]), infile2(argv[2]);

	map<string, string> tag_name;
	
	string line;
	while( getline(infile1, line) ) 
	{
		stringstream ss;
		ss << line;
		string first_gene;
		ss >> first_gene;
		tag_name[first_gene] = first_gene;
		string gene;
		while( ss >> gene ) tag_name[gene] = first_gene;
	}
	infile1.close();

	////////////////////////////////////////////
	
	vector<string> representative;
	map<pair<string, string>, int> rep_score;
	vector<string> blastResult;

	while( getline(infile2, line) )
	{
		stringstream ss;
		ss << line;
		blastResult.push_back(line);
		string gene1, gene2;
		int bitScore;
		ss >> gene1 >> gene2 >> bitScore;
		pair<string, string> gene_pair = make_pair(gene1, gene2);

		if ( tag_name.count(gene1)!=0 and tag_name.count(gene2)==0 )
		{
			rep_score[make_pair(tag_name[gene1], gene2)] = max( bitScore, rep_score[make_pair(tag_name[gene1], gene2)] );
		}
		else if ( tag_name.count(gene1)==0 and tag_name.count(gene2)!=0 )
		{
			rep_score[make_pair(gene1, tag_name[gene2])] = max( bitScore, rep_score[make_pair(gene1, tag_name[gene2])] );
		}
		else if( tag_name.count(gene1)!=0 and tag_name.count(gene2)!=0 ) 
		{
			rep_score[make_pair(tag_name[gene1], tag_name[gene2])] = max( bitScore, rep_score[make_pair(tag_name[gene1], tag_name[gene2])] );
		}
	}
	infile2.close();

	ofstream outfile(argv[2]);

	for(int i=0; i<blastResult.size(); i++)
	{
		stringstream ss;
		ss << blastResult[i];
		string gene1, gene2;
		int bitScore;
		ss >> gene1 >> gene2 >> bitScore;

		if( tag_name.count(gene1)==0 and tag_name.count(gene2)==0 )
		{
			outfile << blastResult[i] << endl;
		}
		else if( tag_name.count(gene1)!=0 and tag_name.count(gene2)==0 )
		{
			if( bitScore == rep_score[make_pair(tag_name[gene1], gene2)] )
				outfile << blastResult[i] << endl;
		}
		else if( tag_name.count(gene1)==0 and tag_name.count(gene2)!=0 )
		{
			if( bitScore == rep_score[make_pair(gene1, tag_name[gene2])] )
				outfile << blastResult[i] << endl;
		}
		else if( tag_name.count(gene1)!=0 and tag_name.count(gene2)!=0 )
		{
			if( bitScore == rep_score[make_pair(tag_name[gene1], tag_name[gene2])] )
				outfile << blastResult[i] << endl;
		}
	}
	outfile.close();


	return 0;
}
