#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

int main(int argc, char** argv)
{
	if(argc!=6)
	{
		cout << "Usage: familySeparator <HM.cluster> <G1.pep> <G2.pep> <G1.nuc> <G2.nuc>" << endl;
		exit(1);
	}
	ifstream infile1(argv[1]), infile2(argv[2]), infile3(argv[3]), infile4(argv[4]), infile5(argv[5]);
	
	map<string, int> G1pos;
	map<string, int> G2pos;
	map<string, string> pep;
	map<string, string> nuc;

	int index=1;
	string name, seq;
	while(infile2 >> name >> seq)
	{
		name = name.substr(1);
		G1pos[name] = index++;
		pep[name] = seq;
	}
	infile2.close();

	index=100001;
	while(infile3 >> name >> seq)
	{
		name = name.substr(1);
		G2pos[name] = index++;
		pep[name] = seq;
	}
	infile3.close();

	while(infile4 >> name >> seq)
	{
		name = name.substr(1);
		nuc[name] = seq;
	}
	infile4.close();

	while(infile5 >> name >> seq)
	{
		name = name.substr(1);
		nuc[name] = seq;
	}
	infile5.close();
	
	///////////////////////////////////////////////////////////////
	
	vector<vector<string> > GeneFamily;

	string line;
	while( getline(infile1, line) )
	{
		stringstream ss;
		ss << line;
		vector <pair<int, string> > G1, G2;
		string gene;
		while( ss >> gene )
		{
			if( G1pos.count(gene) > 0 ) G1.push_back( make_pair( G1pos[gene], gene ) );
			else if( G2pos.count(gene) > 0 ) G2.push_back( make_pair( G2pos[gene], gene ) );
		}
		sort(G1.begin(), G1.end());
		sort(G2.begin(), G2.end());

		int n1=G1.size(), n2=G2.size();
		
		int i=0, j=0;

		int NumOfEach = 100;
		int group = (max(n1,n2)+ NumOfEach - 1)/NumOfEach;

		for(int k=0; k<group; k++)
		{
			vector <string> vf;
			for(int p=i; p<min(i+NumOfEach, n1); p++) vf.push_back( G1[p].second );
			for(int p=j; p<min(j+NumOfEach, n2); p++) vf.push_back( G2[p].second );
			GeneFamily.push_back(vf);
			i+=NumOfEach;
			j+=NumOfEach;
		}
	}

	///////////////////////////////////////////////////////////////

	int familyIndex = 1;
	for(int i=0; i<GeneFamily.size(); i++)
	{
		vector<int> vpp;
		for(int j=0; j<GeneFamily[i].size(); j++)	
		{
			string st = GeneFamily[i][j];
			if( G1pos.count(st) ) vpp.push_back( G1pos[st] );
			else if( G2pos.count(st) ) vpp.push_back( G2pos[st] );
		}
		sort(vpp.begin(), vpp.end());
		bool adjOK=false;
		for(int j=0; j<vpp.size()-1; j++) if( abs(vpp[j+1]-vpp[j])==1 ) { adjOK=true; break; }
		if( adjOK==true )
		{
			string familyName;
			stringstream ss;
			ss << "gf" << familyIndex;
			familyIndex++;
			ss >> familyName;
			string gfPep=familyName+".pep";
			string gfNuc=familyName+".nuc";
			string gfPos=familyName+".pos";
			familyName = "Families/";

			ofstream outfile1( ( familyName + gfPep).c_str() );
			if( outfile1.fail() ) { cout << "Fail open Pep" << endl; exit(1); }

			for(int j=0; j<GeneFamily[i].size(); j++)
			{
				string name = GeneFamily[i][j];
				outfile1 << ">" << name << endl;
				outfile1 << pep[name] << endl;
			}
			outfile1.close();
			
			ofstream outfile2( ( familyName + gfNuc).c_str() );
			if( outfile2.fail() ) { cout << "Fail open Nuc" << endl; exit(1); }

			for(int j=0; j<GeneFamily[i].size(); j++)
			{
				string name = GeneFamily[i][j];
				outfile2 << ">" << name << endl;
				outfile2 << nuc[name] << endl;
			}
			outfile2.close();
			
			ofstream outfile3( ( familyName + gfPos).c_str() );
			if( outfile3.fail() ) { cout << "Fail open Pos" << endl; exit(1); }

			for(int j=0; j<GeneFamily[i].size(); j++)
			{
				string name = GeneFamily[i][j];
				if( G1pos.count(name) ) outfile3 << G1pos[name] << "\t";
				else outfile3 << G2pos[name] << "\t";

				outfile3 << name << endl;
			}
			outfile3.close();
		}
	}
	cout << familyIndex-1 << endl;

	return 0;
}
