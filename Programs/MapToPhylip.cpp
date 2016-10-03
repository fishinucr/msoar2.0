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
	if(argc!=4)
	{
		cout << "Usage: mapToPhylip <Families/gf.codon> <-o gf.codon> <-o gf.map>" << endl;
		exit(1);
	}
	ifstream infile(argv[1]);
	ofstream outfile1(argv[2]), outfile2(argv[3]);

	string line;
	int numOfGenes, numOfAA;
	int index=0;
	
 	getline(infile, line);
	stringstream ss;
	ss<<line;
	ss>>numOfGenes>>numOfAA;

	outfile2<<setw(10)<<setfill('0')<<index++<<"\t"<<"OUTGROUPMY"<<endl;
	outfile1<<" "<<numOfGenes<<"   "<<numOfAA<<endl;
	while(getline(infile, line))
	{
		outfile2<<setw(10)<<setfill('0')<<index<<"\t"<<line<<endl;
		outfile1<<setw(10)<<setfill('0')<<index<<endl;
		index++;
		for(int i=0; i<(numOfAA+59)/60; i++)
		{
			getline(infile, line);
			outfile1<<line<<endl;
		}
	}
	
	infile.close();
	outfile1.close();
	outfile2.close();

	return 0;
}
