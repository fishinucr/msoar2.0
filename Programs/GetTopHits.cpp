#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <set>
#include <algorithm>

using namespace std;

int main(int argc, char** argv)
{
	if(argc!=10)
	{
		cout << "Usage: getTopHits <G1G2.blastp> <bitScore> <Evalue> <AlignLength> <TopNum> <-o G1-G2.top> <-o G2-G1.top> <G1.info> <G2.info>" << endl;
		exit(1);
	}

	int b_cut=atoi(argv[2]);
	double e_cut=atof(argv[3]);
	double len_cut=atof(argv[4]);
	int top_num=atoi(argv[5]);

	ifstream infile(argv[1]);
	ofstream outfile1(argv[6]), outfile2(argv[7]);

	ifstream infile2(argv[8]), infile3(argv[9]);
	string line5;
	map<string, int> mpName;
	while(getline(infile2, line5))
	{
		stringstream ss;
		ss<<line5;
		string tmp;
		ss>>tmp;
		mpName[tmp]=1;
	}
	infile2.close();

	while(getline(infile3, line5))
	{
		stringstream ss;
		ss<<line5;
		string tmp;
		ss>>tmp;
		mpName[tmp]=2;
	}
	infile3.close();


	string line;
	string gene1="ENSG0", gene2;
	double identity;
	double len, mismatches, gap, qstart, qend, sstart, send;
	double evalue;
	int bitscore;
	
	map<string, int> bitScoreMap;
	map<string, double> eValueMap;
	map<string, double> lenMap;

	int cnt=0;

	string genome1="";

	while( getline(infile, line) )
	{
		if(line[0]=='#')
		{
			vector<pair<int,string> > v;
			for(map<string,int>::iterator it=bitScoreMap.begin(); it!=bitScoreMap.end(); it++)
			{
				v.push_back(make_pair(it->second, it->first));
			}
			sort(v.begin(), v.end());
			reverse(v.begin(), v.end());

			for(int k=0; k<v.size() and k<top_num and v[k].first>=b_cut and eValueMap[v[k].second]<=e_cut and lenMap[v[k].second]>=(qend-qstart+1)*len_cut; k++)
			{
				if( mpName[gene1]==1 )
				{
					outfile1 << gene1 << "\t" << v[k].second << "\t" << v[k].first << "\t" << eValueMap[v[k].second] << endl;
				}
				else
				{
					outfile2 << gene1 << "\t" << v[k].second << "\t" << v[k].first << "\t" << eValueMap[v[k].second] << endl;
				}
			}

			bitScoreMap.clear();
			eValueMap.clear();
			lenMap.clear();

			cnt=1;
		}

		while(line[0]=='#') {getline(infile, line);}


		stringstream ss;
		ss << line;
		ss >> gene1 >> gene2 >> identity >> len >> mismatches >> gap >> qstart >> qend >> sstart >> send >> evalue >> bitscore;

		if(mpName[gene1]==mpName[gene2]) continue;

		if( bitscore < b_cut ) continue;
		if( evalue > e_cut ) continue;
		if( len < (qend-qstart)*len_cut ) continue;

		if( bitScoreMap.count(gene2) == 0 or bitScoreMap[gene2] < bitscore )
		{
			bitScoreMap[gene2]=bitscore;
			eValueMap[gene2]=evalue;
			lenMap[gene2]=len;
		}

		cnt++;
	}

	infile.close();
	outfile1.close();
	outfile2.close();


	return 0;
}
