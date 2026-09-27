
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
		
		cube_index = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 21,
24, 26, 28, 30, 32, 34, 36, 39, 41, 43, 45, 47, 49, 51, 
54, 56, 58, 60, 62, 64, 66, 69, 71, 73, 75, 77, 79 };
	
		

		rounds = 849;
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
	else if (cipher_name == "kreyvium")
	{
		set<int> noncube_index = { 66, 73, 85, 87 };
		for (int i = 0; i < 128; i++)
			if (noncube_index.find(i) == noncube_index.end())
				cube_index.emplace(i);

		rounds = 899;
		r0 = 5;
		r1 = 20;
		first_expand_step = 300;
		N = 15000;
		fbound = 320;
		cbound0 = 0;
		cbound1 = 420;
		min_gap = 100;
		p_target_cipher = new cipher_kreyvium();
	}
	else if (cipher_name == "grain128AEAD")
	{
		set<int> noncube_index = { 42,43 };
		for (int i = 0; i < 96; i++)
			if (noncube_index.find(i) == noncube_index.end())
				cube_index.emplace(i);

		rounds = 192;
		r0 = 1;
		r1 = 1;
		first_expand_step = 60;
		N = 15000;
		fbound = 90;
		cbound0 = 0;
		cbound1 = 150;
		min_gap = 20;
		p_target_cipher = new cipher_grain128AEAD();
	}
	else if (cipher_name == "acorn")
	{
		// the first 256 rounds
		set<int> noncube_index = { 2, 28 };
		for(int i = 0; i < 128; i++)
			if (noncube_index.find(i) == noncube_index.end())
				cube_index.emplace(i);

		acorn_basic basic_funcs(293 + 256);

		p_target_cipher = new cipher_acorn256(basic_funcs);

		// expand through the first 256 rounds
		auto initialcube = p_target_cipher->generate_cube(cube_index);
		auto initialstate = p_target_cipher->set_initial_state(initialcube);

		framework first256_framework(256, 0, 0, 0, 0, 0, 0, 0, single_threads,
			solver_mode, threadpool, *p_target_cipher, cube_index, 0);
		map<dynamic_bitset<>, BooleanPolynomial> expand_state_coeff;
		map<dynamic_bitset<>, ListsOfPolynomialsAsFactors> state_coeff_lists;
		first256_framework.noncallback_forwardexpand(0, 256, initialstate, expand_state_coeff, state_coeff_lists, 20, 2);
		if (expand_state_coeff.size() != 1)
		{
			cout << "Expanding to round 256 results in more than one monomial." << endl;
			exit(-1);
		}
		auto &state_coef = (*expand_state_coeff.begin());

		// start to recover the superpoly starting from round 256
		vector<Flag> round256_flag = first256_framework.getFlag(256);
		vector<BooleanPolynomial> round256_exps = first256_framework.getExps(256);
		delete p_target_cipher;
		basic_funcs.reset(293 + 129);

		p_target_cipher = new cipher_acorn(basic_funcs, round256_flag, round256_exps);
		
		auto& round256cube = state_coef.first;
		auto& round256coef = state_coef.second;
		if (!round256coef.isone())
		{
			cout << "Expanding to round 256 results in a non-zero coeffcient." << endl;
			exit(-1);
		}

		cube_index.clear();
		for (int i = 0; i < 293; i++)
			if (round256cube[i])
				cube_index.emplace(i);

		
		rounds = 776-256;
		r0 = 2;
		r1 = 5;
		first_expand_step = 100;
		N = 15000;
		fbound = 450-256;
		cbound0 = 0;
		cbound1 = 540-256;
		min_gap = 50;

	}
	else
	{
		cerr << "The cipher is not defined." << endl;
		exit(-1);
	}



	framework MITM_framework(rounds, r0, r1, first_expand_step, N, fbound, cbound0, cbound1, single_threads,
		solver_mode, threadpool, *p_target_cipher, cube_index, min_gap);



	// MITM_framework.debug_test();

	// MITM_framework.generateTikzCodes(250, 0, -40);

	MITM_framework.start();

	// MITM_framework.continue_after_failed(15, 305, false);

	MITM_framework.stop();

	
	
	if (solver_mode == mode::OUTPUT_FILE)
	{
		// if you want the exact superpoly, then set this variable to true; otherwise set it to false. We recommend true for kreyvium, but false for trivium and grain, otherwise you may encounter an out of memory (OOM) issue. 
		bool isAccurate = true;

		if (isAccurate)
		{
			MITM_framework.read_sols_and_output(true);
		}
		else
		{
			MITM_framework.read_sols_and_output(false);
			MITM_framework.analyze_superpoly_asLists();
		}
	}
	else if(solver_mode == mode::OUTPUT_EXP)
	{
		MITM_framework.analyze_superpoly();
	}


	

	delete p_target_cipher;
	
}

