
#include <iostream>
#include <bitset>
#include <vector>
#include <string>
#include <map>
#include <regex>
#include "gurobi_c++.h"
#include "log.h"
#include "BooleanPolynomial.h"
#include "framework.h"
#include "flag.h"
#include "cipher.h"



using namespace std;



int main()
{

	

	string cipher_name = "trivium";
	set<int> cube_index;
	vector<BooleanPolynomial> noncube_exps;
	map<int, BooleanPolynomial> bit_conditions;
	int rounds;
	int r0, r1;
	int first_expand_step;
	int N;
	int fbound, cbound0, cbound1;
	int min_gap;
	int single_threads = 2;
	cipher* p_target_cipher = NULL;
	mode solver_mode = mode::OUTPUT_FILE;
	ThreadPool threadpool;

	if (cipher_name == "trivium")
	{
		
		cube_index = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 21, 24, 26, 28, 30,
32, 34, 36, 39, 41, 43, 45, 47, 49, 51, 54, 56, 58, 60, 62, 64, 66, 69, 71, 73, 75,
77, 79 };
		
		int noncube_size = 80 - cube_index.size();
		noncube_exps.resize(noncube_size);
		for(int i = 0; i < noncube_size; i++)
			noncube_exps[i] = BooleanPolynomial(288, "0");

		bit_conditions[62] = BooleanPolynomial(288, "0");
		bit_conditions[42] = BooleanPolynomial(288, "s40s41+s15");

		rounds = 852;
		r0 = 5;
		r1 = 20;
		first_expand_step = 300;
		N = 50000;
		fbound = 350;
		cbound0 = 0;
		cbound1 = 600;
		single_threads = 2;
		min_gap = 50;
		p_target_cipher = new cipher_trivium();
	}
	else
	{
		cerr << "The cipher is not defined." << endl;
		exit(-1);
	}



	framework MITM_framework(rounds, r0, r1, first_expand_step, N, fbound, cbound0, cbound1, single_threads,
		solver_mode, threadpool, *p_target_cipher, cube_index, min_gap, noncube_exps, bit_conditions);




	MITM_framework.start(); 

	MITM_framework.stop();

	
	if (solver_mode == mode::OUTPUT_FILE)
	{
		// if you want the exact superpoly, then set this variable to true; otherwise set it to false. We recommend true for kreyvium, but false for trivium and grain, otherwise you may encounter an out of memory (OOM) issue. 
		bool isAccurate = false;

		if (isAccurate)
		{
			MITM_framework.read_sols_and_output(true);
		}
		else
		{
			MITM_framework.read_sols_and_output(false);
			// vector<string> TERM_paths = { R"~(./TERM_841)~" };
			MITM_framework.analyze_superpoly_asLists(); 
		}
	}
	else if(solver_mode == mode::OUTPUT_EXP)
	{
		MITM_framework.analyze_superpoly();
	}


	

	delete p_target_cipher;
	
}

