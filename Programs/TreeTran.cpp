// =====================================================================================
// 
//       Filename:  main.cc
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  10/22/2009 09:08:53 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Jianxing Feng (feeldead), feeldead@gmail.com
//        Company:  THU
// 
// =====================================================================================
//

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stack>

using namespace std;

class TreeNode
{
public:
	int mOutGroup;
	string mID;
	vector<TreeNode*> mpChilds;
};

bool IsOutGroup(TreeNode* p_root, string outgroup_id)
{
	p_root->mOutGroup = -1;

	if (p_root->mpChilds.size() == 0)
	{
		if (p_root->mID == outgroup_id)
		{
			p_root->mOutGroup = -2;
			return true;
		}
		else
			return false;
	}

	for (unsigned i = 0; i < p_root->mpChilds.size(); i++)
	{
		if (IsOutGroup(p_root->mpChilds[i], outgroup_id))
			p_root->mOutGroup = i;
	}
	if (p_root->mOutGroup >= 0) return true;
	return false;
}

void PrintTree(TreeNode* p_root)
{
	if (p_root->mpChilds.size() == 0)
	{
		//cout << p_root->mOutGroup;
		cout << p_root->mID;
		return;
	}
	cout << "(";
	for (unsigned i = 0; i < p_root->mpChilds.size(); i++)
	{
		PrintTree(p_root->mpChilds[i]);
		if (i < p_root->mpChilds.size() - 1)
			cout << ",";
	}
	cout << ")";
	return;
}

void DeleteTree(TreeNode* p_root)
{
	for (unsigned i = 0; i < p_root->mpChilds.size(); i++)
		DeleteTree(p_root->mpChilds[i]);
	delete p_root;
}

// Return the depth of the outgroup
int PrintTreeOutGroupSub(TreeNode* p_root, int depth)
{
	int outgroup_depth = 0;
	if (-2 == p_root->mOutGroup)
		return depth;
	else if (-1 == p_root->mOutGroup)
		PrintTree(p_root);
	else if (p_root->mOutGroup >= 0)
	{
		outgroup_depth = PrintTreeOutGroupSub(p_root->mpChilds[p_root->mOutGroup], depth+1);

		int last_non_outgroup = p_root->mpChilds.size() - 1;
		if (last_non_outgroup == p_root->mOutGroup)
			last_non_outgroup --;

		// If 0-depth root has only two child, do not output an extra "("
		if (0 != depth || p_root->mpChilds.size() > 2)
			cout << "(";
		for (unsigned i = 0; i < p_root->mpChilds.size(); i++)
		{
			if (i != p_root->mOutGroup)
			{
				PrintTree(p_root->mpChilds[i]);
				if (i != last_non_outgroup || depth != 0)
					cout << ",";
			}
		}

		if (0 == depth)
		{
			int match_cnt = outgroup_depth;
			// If an extra "(" is output before, output an extra ")" too
			if (p_root->mpChilds.size() > 2)
				match_cnt++;
			for (int i = 0; i < match_cnt - 1; i++)
				cout << ")";
		}
	}

	return outgroup_depth;
}

void PrintTreeOutGroup(TreeNode* p_root, string outgroup_id)
{
	cout << "(";
	int depth = PrintTreeOutGroupSub(p_root, 0);
	cout << "," << outgroup_id << ")" << endl;
}

bool ReadTree (TreeNode* p_root, fstream& input)
{
	char c;
	input >> c;
	if (c == '(')
	{
		bool b_succ = false;
		while (!b_succ) 
		{
			TreeNode* p_node = new TreeNode;
			b_succ = ReadTree(p_node, input);
			p_root->mpChilds.push_back(p_node);
		}
		input >> c;
	}
	else
	{
		bool b_score = false;
		string id;
		while (c != ')' && c != ',')
		{
			if (c == ':') b_score = true;
			if (!b_score) id += c;
			input >> c;
		}
		p_root->mID = id;
		//cout << id << endl;
	}

	if (c == ',') return false;

	// Skip scores
	c = input.peek();
	while (c != ')' && c != ',' && c != EOF)
	{
		input >> c;
		c = input.peek();
	}

	if (p_root->mpChilds.size() == 1)
	{
		cerr << "WARNING : There is a node contains only one child. " << endl;
	}

	return true;
}

int 
main(int argc, char* argv[])
{
	if (argc <= 2)
	{
		cout << 											endl;
		cout << "        Usage : treetran <outgroup ID> <file name>" << endl;
		cout << 											endl;
		return 0;
	}

	string outgroup_id = argv[1];
	string input_file_name = argv[2];

	fstream infile;
	infile.open(input_file_name.data(), ios::in);
	if (!infile.is_open())
	{
		cerr << "File " << input_file_name << " can not be opened" << endl;
		return 0;
	}

	TreeNode* p_root = new TreeNode;

	ReadTree(p_root, infile);

	if (!IsOutGroup(p_root, outgroup_id))
		cerr << "ERROR: The outgroup does not exist in the tree." << endl;
	else
	{
		// PrintTree(p_root);
		// cout << endl;
		PrintTreeOutGroup(p_root, outgroup_id);
	}

	DeleteTree(p_root);

	infile.close();
	return 0;
}

