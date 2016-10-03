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
	if(argc!=7)
	{
		cout << "Usage: a.out <G1_filtered_unique> <G2_filtered_unique> <Result> <G1-G2.score> <G2-G1.score> <cutoff>" << endl;
		exit(1);
	}
	ifstream infile1(argv[1]), infile2(argv[2]), infile3(argv[3]), infile4(argv[4]), infile5(argv[5]);
	double cutoff=atof(argv[6])*100.0;
	
	string line;
	map<string, int> pos;
	string name;

	map<int, string> mp1, mp2;

	int index=0;
	while( getline(infile1, line))
	{
		index++;
		stringstream ss;
		ss << line;
		ss >> name;
		pos[name] = index;
		mp1[index] = name;
	}
	infile1.close();

	index=0;
	while( getline(infile2, line))
	{
		index++;
		stringstream ss;
		ss << line;
		ss >> name;
		pos[name] = index;
		mp2[index] = name;
	}
	infile2.close();
	
	vector<pair<int,int> > v;
	while(getline(infile3, line))
	{
		stringstream ss;
		ss << line;
		string gene1, gene2;
		ss >> gene1 >> gene2;
		v.push_back(make_pair(pos[gene1], pos[gene2]));
	}
	infile3.close();

	sort(v.begin(), v.end());

	map<pair<string,string>, double> mp;
	
	while(getline(infile4, line))
	{
		stringstream ss;
		ss << line;
		string gene1, gene2;
		double score;
		ss >> gene1 >> gene2 >> score;
		mp[make_pair(gene1, gene2)]=score;
	}
	infile4.close();


	while(getline(infile5, line))
	{
		stringstream ss;
		ss << line;
		string gene1, gene2;
		double score;
		ss >> gene1 >> gene2 >> score;
		mp[make_pair(gene1, gene2)]=score;
	}
	infile5.close();


	cout << mp1[v[0].first] << "\t" << mp2[v[0].second] << endl;

	for(int i=1; i<v.size(); i++)
	{
		int gene1=v[i].first, gene2=v[i].second; 
		if( gene1-v[i-1].first==2 and gene2-v[i-1].second==2)
		{
			if(mp[make_pair(mp1[gene1-1], mp2[gene2-1])]>=cutoff-1e-5 or mp[make_pair(mp2[gene2-1], mp1[gene1-1])]>=cutoff-1e-5)
					cout << mp1[gene1-1] << "\t" << mp2[gene2-1] << endl;
		}
		else if( gene1-v[i-1].first==2 and gene2-v[i-1].second==-2)
		{
			if(mp[make_pair(mp1[gene1-1], mp2[gene2+1])]>=cutoff-1e-5 or mp[make_pair(mp2[gene2+1], mp1[gene1-1])]>=cutoff-1e-5)
			cout << mp1[gene1-1] << "\t" << mp2[gene2+1] << endl;
		}
		cout << mp1[gene1] << "\t" << mp2[gene2] << endl;
	}

	return 0;
}
