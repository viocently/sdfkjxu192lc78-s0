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

		rounds = 850;
		r0 = 5;
		r1 = 20;
		first_expand_step = 300;
		N = 50000;
		fbound = 320;
		cbound0 = 0;
		cbound1 = 600;
		single_threads = 2;
		min_gap = 50;
		p_target_cipher = new cipher_trivium();
	}
	else if (cipher_name == "kreyvium")
	{
		set<int> noncube_index = { 66,73,106,110 };
		for (int i = 0; i < 128; i++)
			if (noncube_index.find(i) == noncube_index.end())
				cube_index.emplace(i);

		rounds = 897;
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
	else
	{
		cerr << "The cipher is not defined." << endl;
		exit(-1);
	}



	framework MITM_framework(rounds, r0, r1, first_expand_step, N, fbound, cbound0, cbound1, single_threads,
		solver_mode, threadpool, *p_target_cipher, cube_index, min_gap);

	// MITM_framework.compare_cube(8);

	MITM_framework.start();

	MITM_framework.stop();

	// MITM_framework.read_sols_and_output();
	
	delete p_target_cipher;
	
}

