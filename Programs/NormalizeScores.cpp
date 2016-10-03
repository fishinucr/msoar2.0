#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>

using namespace std;

int main(int argc, char** argv)
{
	if(argc!=2)
	{
		cout << "Usage: a.out <G1.pep-G2.pep>" << endl;
		exit(1);
	}
	ifstream infile(argv[1]);

	string line;
	string gene1, gene2, firstGene="";
	double score, maxScore;

	while( getline(infile, line) )
	{
		stringstream ss(line);
		ss >> gene1 >> gene2 >> score;
		if( gene1 != firstGene )
		{
			firstGene = gene1;
			maxScore = score;
			cout << gene1 << "\t" << gene2 << "\t" << fixed << setprecision(1) << 100.0 <<endl;
		}
		else
		{
			score = (100.0*score)/maxScore;
			cout << gene1 << "\t" << gene2 << "\t" << fixed << setprecision(1) << score <<endl;
		}
	}
	infile.close();



	return 0;
}
