#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <fstream>
#include <set>

using namespace std;

struct Node
{
	int id;
	int property; // (property|1):hg; (property|2):mm; (property|3):speciation
	Node* parent;
	vector <Node*> child;
	Node(){}
	Node(int a, Node* b) { id=a; parent=b; property=0; }
};


map<string, Node*> gene_parent;


void Gene_Dating( Node* cur )
{
	for(int i=0; i < cur->child.size(); i++)
	{
		Node* p = cur->child[i];
		Gene_Dating( p );
		cur->property |= p->property;
	}
}

void Print_Tree( Node* cur )
{
	string prop = "Speciation";
 	if( cur->property != 3 )	prop = "Duplication";

	cout << cur->id << ": " << prop << endl;

	for(int i=0; i < cur->child.size(); i++)
	{
		Print_Tree( cur->child[i] );
	}
}

string LCA( string gene1, string gene2 )
{
	Node* p1 = gene_parent[gene1];
	while( p1 != NULL )
	{
		Node* p2 = gene_parent[gene2];
		while( p2 != NULL )
		{
			if( p2 == p1 )
			{
				if( p1->property == 3 ) return "Speciation";
				else	return "Duplication";
			}
			p2 = p2->parent;
		}
		p1 = p1->parent;
	}
	return "Speciation";
}

int main(int argc, char** argv)
{
	if(argc!=5)
	{
		cout << "Usage: a.out <gfxxx.tree> <gfxxx.pos> <gfxxx.map> <-o gfxxx.tags>" << endl;
		exit(1);
	}

	int id=0;
	Node* root = NULL;

	// Read Tree File
	string treefile = "";
	ifstream infile1(argv[1]);
	if(infile1.fail()) return 0;
	ifstream infile2(argv[2]);
	string line;
	while( getline(infile1, line) ) treefile += line;
	infile1.close();

	// Read Gene Position File
	map<int, string> gene_pos;
	map<string, int> gene_pos2;
	gene_pos[99999] = "OUTGROUPMY";
	gene_pos2["OUTGROUPMY"] = 99999;

	while( getline(infile2, line) )
	{
		stringstream ss;
		ss << line;
		string name;
		int pos;
		ss >> pos >> name;
		gene_pos[pos] = name;
		gene_pos2[name] = pos;
	}
	infile2.close();

	// Read Map File
	ifstream infile3(argv[3]);
	string newname, oldname;
	map<string, string> GeneNameConverter;
	while( infile3 >> newname >> oldname ) GeneNameConverter[newname] = oldname;

	/////////////////////////////////////////////////////////////////
	// Construct Gene Tree
	/////////////////////////////////////////////////////////////////
	
	Node* cur = root;
	int pos=-1;
	int prev_pos = pos;
	while( treefile.find_first_of(")(,", pos+1) < treefile.length() )
	{
		pos = treefile.find_first_of(")(,", pos+1);
		if( treefile[pos] == '(' )
		{
			Node* c = new Node(id++, cur);
			if( root == NULL )
			{
				root = c;
				cur = root;
			}
			else
			{
				cur->child.push_back(c);
				cur = c;
			}
		}
		else if ( treefile[pos] == ')' )
		{
			if( treefile[prev_pos] != ')')
			{
				string gene;
				stringstream ss;
				ss << treefile.substr(prev_pos+1, pos-prev_pos-1);
				ss >> gene;
				//gene = gene.substr(0, gene.find(':'));
				gene = GeneNameConverter[gene];
				gene_parent[gene] = cur;
				if( gene_pos2[gene] >= 100000 ) cur->property |= 2;
				else if( gene_pos2[gene] < 99999) cur->property |= 1;
			}
			cur = cur->parent;
		}
		else if( treefile[prev_pos] != ')' )
		{
			string gene;
		 	stringstream ss;
			ss << treefile.substr(prev_pos+1, pos-prev_pos-1);
			ss >> gene;
			//gene = gene.substr(0, gene.find(':'));
			gene = GeneNameConverter[gene];
			gene_parent[gene] = cur;
			if( gene_pos2[gene] >= 100000 ) cur->property |= 2;
			else if( gene_pos2[gene] < 99999) cur->property |= 1;
		}
		prev_pos = pos;
	}
	/*
	for(map<string, Node*>::iterator i=gene_parent.begin(); i!=gene_parent.end(); i++)
	{
		cout << (*i).first <<"'s parent: " << ((*i).second)->id << endl;
	}
	*/

	///////////////////////////////////////////////////////////////////
	// Differentiate Duplication and Speciation Events
	///////////////////////////////////////////////////////////////////
	
	Gene_Dating( root );
	
	/////// Print results
	
	// Print_Tree( root );

	//////////////////////////////////////////////////////////////////
	// Given Two Genes in TAGs, identify their LCA
	// corresponding to duplication or speciation event
	//////////////////////////////////////////////////////////////////
	

	/*
	vector <vector<string> > TAGs(1);
	TAGs[0].push_back((gene_parent.begin())->first);
	for(map<string, Node*>::iterator i=gene_parent.begin(); i!=gene_parent.end(); i++)
	{
		for(map<string, Node*>::iterator j=
		if( j == gene_parent.end() ) break;

		// Define TAGs as adjacent genes on genome with no spacers in between
		int adj_distance = gene_pos[ j->first ] - gene_pos[ i->first ];


		if( abs(adj_distance) <= 1 and LCA( (*i).first, (*j).first ) == "Duplication" )
		{
			TAGs[TAGs.size()-1].push_back((*j).first);
		}
		else
		{
			vector <string> newTag;
			newTag.push_back((*j).first);
			TAGs.push_back(newTag);
		}
	}

	*/

	vector <vector<string> > TAGs(1);
	TAGs[0].push_back((gene_pos.begin())->second);
	for(map<int, string>::iterator i=gene_pos.begin(); i!=gene_pos.end(); i++)
	{
		map<int, string>::iterator j=i; j++;
		if( j==gene_pos.end() ) break;

		if( j->first == i->first + 1 and LCA( i->second, j->second ) == "Duplication" )
			TAGs[TAGs.size()-1].push_back(j->second);
		else 
		{
			vector <string> newTag;
			newTag.push_back( j->second );
			TAGs.push_back(newTag);
		}
	}
	ofstream outfile(argv[4]);

	for(int i=0; i<TAGs.size();  i++)
	{
		if( TAGs[i].size() <= 1 ) continue;
		for(int j=0; j<TAGs[i].size(); j++)
		{
			outfile << TAGs[i][j] << " ";
		}
		outfile << endl;
	}

	return 0;

}
