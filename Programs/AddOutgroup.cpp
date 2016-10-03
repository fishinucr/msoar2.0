#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>


using namespace std;

int main(int argc, char** argv)
{
	if(argc!=2) { "Input Error: "; return 0;}

	ifstream infile(argv[1]);
	int n;
	infile >> n;
	n++;
	vector <vector<double> > v(n, n);
	vector <string> vn(n);
	vn[0] = "0000000000";
	v[0][0]=0.0;

	double maxv=0.0;
	for(int i=1; i<n; i++)
	{
		infile >> vn[i];
		for(int j=1; j<n; j++)
		{
			infile >> v[i][j];
			maxv=max(maxv, v[i][j]);
		}
	}

	infile.close();

	for(int i=1; i<n; i++) v[i][0] = v[0][i] = min(99.0, 4*maxv);

	/////////////////////////////
	
	cout << right << setw(10) << n << endl;
	for(int i=0; i<n; i++)
	{
		cout << left << setw(15) << vn[i];
		for(int j=0; j<n; j++)
			cout << right << setw(11) << fixed << setprecision(6) << (v[i][j]>-1e-10?v[i][j]:2*maxv);
		cout << endl;
	}
	return 0;
}
