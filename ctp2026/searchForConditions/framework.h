#pragma once
#ifndef _FRAMEWORK_H_
#define _FRAMEWORK_H_

#include<sstream>
#include<map>
#include <algorithm>
#include <random>
#include"cipher.h"
#include "thread_pool.h"
#include "deg.h"
#include "log.h"
#include "file_reader.h"
#include "listsOfPolynomials.h"

using namespace thread_pool;
enum class mode { OUTPUT_FILE, OUTPUT_EXP, NO_OUTPUT };
mutex solver_mutex0;
mutex solver_mutex1;
mutex expander_mutex0;
mutex expander_mutex1;
mutex reader_mutex;
mutex analyze_mutex;
mutex tester_mutex;

#define REF(x) std::ref(x)






// the framework
class framework
{
private:

	// the data structure and parameters
	map<node_pair, ListsOfPolynomialsAsFactors> P_asLists; // another representation of P
	map<node_pair, BooleanPolynomial> P;
	int rs, re;
	int r0, r1;
	int first_expand_step;
	int N;
	int fbound; // the number of rounds that choose the forward expansion first
	int cbound1; // the number of rounds that use callback for backward expansion
	int cbound0; // the number of rounds that use callback for forward expansion 
	int d; // 0 for forward, 1 for backward
	int single_threads;
	int min_gap;
	mode solver_mode;
	ThreadPool& threadpool;


	// target cipher
	cipher& target_cipher;
	int target_rounds;

	// precomputation
	vector<vector<vector<BooleanPolynomial>>> all_rounds_exps;
	vector<BooleanPolynomial> normal_exp;
	vector<vector<vector<Flag>>>  all_rounds_flags;

	vector<Flag> initial_flags;
	vector<BooleanPolynomial> initial_exps;
	dynamic_bitset<> initial_state;
	dynamic_bitset<> cube;
	vector<BooleanPolynomial> noncube_exps;
	map<int, BooleanPolynomial> bit_conditions;
	vector<BooleanPolynomial> bit_conditions_submap;

	// the final superpoly
	map<BooleanMonomial, int> superpoly;




public:
	framework(int rounds, int r0, int r1, int first_expand_step, int N, int fbound, int cbound0, int cbound1, int single_threads,
		mode solver_mode, ThreadPool& threadpool, cipher& target_cipher, const set<int>& cube_index, int min_gap, const vector<BooleanPolynomial> &noncube_exps, const map<int, BooleanPolynomial> & bit_conditions)
		: rs(0), re(rounds), target_rounds(rounds), r0(r0), r1(r1), first_expand_step(first_expand_step), N(N), fbound(fbound), cbound0(cbound0), cbound1(cbound1),
		single_threads(single_threads), solver_mode(solver_mode), threadpool(threadpool), target_cipher(target_cipher), min_gap(min_gap), noncube_exps(noncube_exps), bit_conditions(bit_conditions)
	{

		cube = target_cipher.generate_cube(cube_index);


		logger("Simplify target cipher according to bit operations.");
		target_cipher.simplify_cipher_under_bit_conditions(cube, noncube_exps, bit_conditions);
		logger("Simplifying target cipher finished.");

		initial_flags = target_cipher.set_initial_flag(cube, noncube_exps);
		initial_exps = target_cipher.set_initial_exps(cube, noncube_exps);
		initial_state = target_cipher.set_initial_state(cube, noncube_exps);
		normal_exp = target_cipher.normal_exp();

		// impose bit conditions on initial exps
		bit_conditions_submap = normal_exp;
		for (auto& fromto : bit_conditions)
			bit_conditions_submap[fromto.first] = fromto.second;
		for (auto& initial_exp : initial_exps)
			initial_exp = initial_exp.subs(bit_conditions_submap);

		for (int i = 0; i < target_cipher.statesize; i++)
			if (initial_exps[i].iszero())
				initial_flags[i].setflag("zero_c");
		//

		all_rounds_exps = target_cipher.generate_exps(0, rounds, initial_flags, initial_exps);


		target_cipher.calculate_flags(0, rounds, initial_flags, all_rounds_flags, true);
		auto last_flag = all_rounds_flags[rounds][0];
		auto last_exps = all_rounds_exps[rounds][0];
		vector<vector<vector<Flag>>> ks_flags;
		target_cipher.calculate_flags(rounds, rounds + 1, last_flag, ks_flags, false, target_cipher.outputks);
		vector<vector<vector<BooleanPolynomial>>> ks_exps = target_cipher.generate_exps(rounds, rounds + 1, last_flag, last_exps, target_cipher.outputks);

		all_rounds_flags[rounds] = ks_flags[0];
		all_rounds_exps[rounds] = ks_exps[0];
		all_rounds_flags.emplace_back(ks_flags[1]);
		all_rounds_exps.emplace_back(ks_exps[1]);

		cout << "Initialize framework complete." << endl;

	}




	void generateTikzCodes(int rounds, double xBasis, double yBasis)
	{
		FlagDisplay flagDisplay(all_rounds_flags[rounds][0], 5, 1, yBasis, xBasis);
		flagDisplay.generateTikzCodes(cout, rounds);
	}

	vector<Flag> getFlag(int round)
	{
		return all_rounds_flags[round][0];
	}

	vector<BooleanPolynomial> getExps(int round)
	{
		return all_rounds_exps[round][0];
	}



	void filterDeg()
	{
		if (rs == 0 && target_cipher.ciphername == "kreyvium")
			kreyvium_filterDeg();
	}

	// this function is specific to kreyvium
	void kreyvium_filterDeg()
	{
		string cubestr;
		to_string(cube, cubestr);

		string noncube_val_str;
		dynamic_bitset<> noncube_val(128);
		int j = 0;
		for (int i = 0; i < 128; i++)
			if (cube[i] == 0 && !noncube_exps[j++].iszero())
				noncube_val[i] = 1;

		to_string(noncube_val, noncube_val_str);
		bitset<128> kreyvium_cube(cubestr);
		bitset<128> kreyvium_noncube_val(noncube_val_str);

		int size0 = P.size();
		map<node_pair, BooleanPolynomial> tmp_P;
		map<node_pair, ListsOfPolynomialsAsFactors> tmp_P_asLists;
		for (auto& it : P)
		{
			auto& nodepair = it.first;
			auto& end_state = nodepair.second;
			auto& coef_lists = P_asLists.at(nodepair);
			string statestr;
			to_string(end_state, statestr);
			bitset<544> state(statestr);

			auto d = computeDegree(kreyvium_cube, kreyvium_noncube_val, re-128, state);
			if (d >= kreyvium_cube.count())
			{
				tmp_P.emplace(it);
				tmp_P_asLists[nodepair] = coef_lists;
			}
		}

		P = tmp_P;
		P_asLists = tmp_P_asLists;
		int size1 = P.size();

		logger(__func__ + string(": ") + to_string(size0) + string("\t") +
			to_string(size1));
	}


	void P_filter()
	{
		int size0 = P.size();
		map<node_pair, BooleanPolynomial> tmp_P;
		map<node_pair, ListsOfPolynomialsAsFactors> tmp_P_asLists;

		for (auto& nodes_coef : P)
		{
			auto& nodes = nodes_coef.first;
			auto& coef = nodes_coef.second;
			auto& coef_lists = P_asLists.at(nodes);
			if (!coef.iszero())
			{
				tmp_P[nodes] = coef;
				tmp_P_asLists[nodes] = coef_lists;
			}
		}

		P = tmp_P;
		P_asLists = tmp_P_asLists;
		int size1 = P.size();
		logger("Filter P: " + to_string(size0) + " " + to_string(size1));
	}




	void superpoly_filter()
	{
		int size0 = superpoly.size();
		map<BooleanMonomial, int> tmp_superpoly;
		for (auto& mon_cnt : superpoly)
		{
			auto& cnt = mon_cnt.second;
			if (cnt % 2)
				tmp_superpoly.emplace(mon_cnt);
		}

		superpoly = tmp_superpoly;
		int size1 = superpoly.size();
		logger("Filter superpoly: " + to_string(size0) + " " + to_string(size1));
	}

	inline void P_output()
	{
		logger("Output P to file.");
		string path = string("STATE/") + to_string(rs) + "_" + to_string(re) + "_P.txt";
		ofstream os;
		os.open(path, ios::out);
		for (auto& nodes_coef : P)
		{
			auto& nodes = nodes_coef.first;
			auto& coef = nodes_coef.second;
			auto& coef_lists = P_asLists.at(nodes);

			auto& start_state = nodes.first;
			auto& end_state = nodes.second;
			os << start_state << endl;
			os << end_state << endl;
			os << coef << endl;
			os << coef_lists << endl;
			os << endl;
		}
		os.close();
		logger("Output P to file finished.");
	}

	virtual void expander()
	{
		int execnt = 0;

		// if re is less than fbound, we choose to forward expand first
		if (re < fbound)
			d = 0;
		else
			d = 1;


		while (P.size() > 0 && (P.size() < N || execnt == 0) && re - rs > min_gap)
		{
			bool callback_flag0;
			bool callback_flag1;
			if (rs < cbound0)
				callback_flag0 = false;
			else
				callback_flag0 = true;


			if (re > cbound1)
				callback_flag1 = false;
			else
				callback_flag1 = true;

			if (d == 1)
			{
				map<node_pair, BooleanPolynomial> new_P;
				map<node_pair, ListsOfPolynomialsAsFactors> new_P_asLists;
				auto expander_time = target_cipher.set_expander_time(rs, re);
				backexpand(rs, re, r1, P, P_asLists, new_P, new_P_asLists, expander_time, single_threads, threadpool, callback_flag1);

				re -= r1;
				P = new_P;
				P_asLists = new_P_asLists;
			}
			else
			{
				map<node_pair, BooleanPolynomial> new_P;
				map<node_pair, ListsOfPolynomialsAsFactors> new_P_asLists;
				auto expander_time = target_cipher.set_expander_time(rs, re);
				forwardexpand(rs, re, r0, P, P_asLists, new_P, new_P_asLists, expander_time, single_threads, threadpool, callback_flag0);

				rs += r0;
				P = new_P;
				P_asLists = new_P_asLists;
			}


#ifdef _WIN32
#else
			showProcessMemUsage();
			malloc_trim(0);
			showProcessMemUsage();
#endif

			filterDeg();
			P_filter();

			// output current state and coeffs to file  for backup
			P_output();

			


			logger("Current size of P : " + to_string(P.size()) );
			logger("-------------------------------------------------------------");

			execnt++;
			if (re < fbound)
				d = (d + 1) % 2;


			
		}

	}

	virtual void solver()
	{
		logger("Start to solve P.");
		auto solver_time = target_cipher.set_solver_time(rs, re);

		

		
		
		
		
		map<node_pair,BooleanPolynomial> new_P;
		map<node_pair, ListsOfPolynomialsAsFactors> new_P_asLists;
		solve_nodes(rs, re, P, P_asLists, solver_time, single_threads, threadpool, new_P, new_P_asLists, superpoly);
		P = new_P;
		P_asLists = new_P_asLists;
		superpoly_filter();
		logger("Current superpoly size: " + to_string(superpoly.size()));
		logger("Current unsolved nodes: " + to_string(P.size()) );
	}

	virtual void expand_first()
	{
		while (true)
		{
			expander();






			solver();


#ifdef _WIN32
#else
			showProcessMemUsage();
			malloc_trim(0);
			showProcessMemUsage();
#endif

			if (P.size() == 0)
				break;


		}

		cout << "Success!" << endl;
	}

	

	virtual void solve_first()
	{
		solver();


#ifdef _WIN32
#else
		showProcessMemUsage();
		malloc_trim(0);
		showProcessMemUsage();
#endif

		if (P.size() == 0)
		{
			cout << "Success!" << endl;
			return;
		}

		expand_first();
	}

	virtual void start(int num_of_constraints)
	{

		double expander_time = target_cipher.set_expander_time(rs, re);
		first_expand(re, first_expand_step, expander_time, single_threads);
		re -= first_expand_step;

		filterDeg();
	
		P_output();

		// solve_first();

		tester(num_of_constraints);
	}

	virtual void tester(int num_of_constraints)
	{
		logger("Start to search optimal P constraints.");
		cout << "Start to search for optimal conditions." << endl;
		double solver_time = 0;
		search_optimal_constraints(rs, re, P, solver_time, single_threads, threadpool, num_of_constraints);
	}

	virtual map<BooleanMonomial, int> & retrive_superpoly()
	{
		return superpoly;
	}

	virtual void stop()
	{
		if (solver_mode == mode::OUTPUT_EXP)
		{
			cout << "Output superpoly to file." << endl;
			string path = string("TERM/") + string("superpoly.txt");
			ofstream os;
			os.open(path, ios::out | ios::app);
			for (auto& mon_cnt : superpoly)
				os << mon_cnt.first << endl;
			os.close();
			cout << "Output superpoly to file finished." << endl;
		}
	}





	// first back expand the output bit
	virtual void first_expand(int rounds, int step, double time, int threads)
	{
		int& statesize = target_cipher.statesize;
		auto start_flags = all_rounds_flags[rounds - step][0];
		if ( start_flags.size() != statesize)
		{
			cerr << __func__ << ": The number of start flags is invalid.";
			exit(-1);
		}

		// set env

		GRBEnv env = GRBEnv();

		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, MAX);
		env.set(GRB_IntParam_Threads, threads);

		GRBModel model = GRBModel(env);

		// set initial variables
		vector<GRBVar> start_vars(statesize);

		for (int i = 0; i < statesize; i++)
			if (start_flags[i] == "delta")
			{
				start_vars[i] = model.addVar(0, 1, 0, GRB_BINARY);
			}

		// build model round by round
		vector< vector<vector<pair<BooleanPolynomial, GRBVar>> > > rounds_p_maps;
		vector<GRBVar> end_vars(statesize);
		target_cipher.build_model(model, rounds-step, rounds, all_rounds_flags, start_vars, end_vars, rounds_p_maps, true);

		
		vector<GRBVar> ks_vars(statesize);
		vector< vector<vector<pair<BooleanPolynomial, GRBVar>>> > ks_p_maps;

		// impose constraints on the output bit
		target_cipher.build_model(model, rounds, rounds + 1, all_rounds_flags, end_vars, ks_vars, ks_p_maps, target_cipher.outputks,false);

		auto ks_flags = all_rounds_flags[rounds + 1][0];
		if (ks_flags[0] == "delta")
			model.addConstr(ks_vars[0] == 1);
		else
		{
			cerr << __func__ << ": The output is a constant under the chosen cube.";
			exit(-1);
		}
		

		vector<BooleanPolynomial> p_exps;
		vector<pair<int, int>> p_rms;
		vector<GRBVar> p_vars;
		
		for (int r = 0; r < step; r++)
		{
			int list_index = target_cipher.cal_list_index(r + rounds - step);

			for (int i = 0; i < target_cipher.updatelists[list_index].size(); i++)
			{
				for (auto& e_v : rounds_p_maps[r][i])
				{
					p_exps.emplace_back(e_v.first);
					p_vars.emplace_back(e_v.second);
					p_rms.emplace_back(pair(r, i));
				}
			}
		}

		
		for (int r = 0; r < 1; r++)
		{
			for (int i = 0; i < target_cipher.outputks.size(); i++)
			{
				for (auto& e_v : ks_p_maps[r][i])
				{
					p_exps.emplace_back(e_v.first);
					p_vars.emplace_back(e_v.second);
					p_rms.emplace_back(pair(step+r, i));
				}
			}
		}
		

		int p_num = p_vars.size();

		if (time > 0)
			model.set(GRB_DoubleParam_TimeLimit, time);

		model.optimize();

		if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
		{
			// this should not happen
			logger(__func__ + string(" : ") + to_string(rounds- step) + "-" + to_string(rounds) +
				string(" | Failed "));
		}
		else
		{
			int solCount = model.get(GRB_IntAttr_SolCount);
			if (solCount >= MAX)
			{
				cerr << "solCount value  is too big !" << endl;
				exit(-1);
			}

			double time = model.get(GRB_DoubleAttr_Runtime);
			logger(__func__ + string(" : ") + to_string(rounds - step) + "-" + to_string(rounds) + string(" | "
			) + to_string(p_num) + string(" | ") + to_string(time) + string(" | "
			) + to_string(solCount));

			if (solCount > 0)
			{

				// first collect solutions
				map<dynamic_bitset<>, map<dynamic_bitset<>, int> >  state_p_sols_counter;
				dynamic_bitset<> p_sol(p_num);
				dynamic_bitset<> state(statesize);

				for (int i = 0; i < solCount; i++)
				{
					model.set(GRB_IntParam_SolutionNumber, i);
					// read the input state
					for (int j = 0; j < statesize; j++)
						if (start_flags[j] == "delta")
							if (round(start_vars[j].get(GRB_DoubleAttr_Xn)) == 1)
								state[j] = 1;
							else
								state[j] = 0;
						else
							state[j] = 0;

					// read the p_sol
					for (int j = 0; j < p_num; j++)
					{
						if (round(p_vars[j].get(GRB_DoubleAttr_Xn)) == 1)
							p_sol[j] = 1;
						else
							p_sol[j] = 0;
					}

					state_p_sols_counter[state][p_sol]++;
				}

				for (auto& sp : state_p_sols_counter)
				{
					auto &state = sp.first;
					auto &p_sols_counter  = sp.second;
					
					vector<dynamic_bitset<>> p_sols;
					for (auto& p_cnt : p_sols_counter)
						if (p_cnt.second % 2)
							p_sols.emplace_back(p_cnt.first);

					if (p_sols.size() == 0)
						continue;

					ListsOfPolynomialsAsFactors cof_lists = get_cof_lists(rounds - step, rounds + 1, p_sols, p_exps, p_rms, 0, all_rounds_exps);
					P_asLists[pair(initial_state, state)] = cof_lists;


					auto expand_res = expand_cof(cof_lists);
					P[pair(initial_state, state)] = expand_res;
				}

			}
			else
			{
				;
			}
		}
	}

	virtual status callback_forwardexpand(int start, int end, int step, const dynamic_bitset<>& start_state, const dynamic_bitset<>& end_state,
		map<dynamic_bitset<>, BooleanPolynomial>& expand_state_coeff, map<dynamic_bitset<>, ListsOfPolynomialsAsFactors>& state_coeff_lists, double time, int threads)
	{

		int& statesize = target_cipher.statesize;
		auto start_flags = all_rounds_flags[start][0];
		if (start_flags.size() != statesize)
		{
			cerr << __func__ << ": The number of start flags is invalid.";
			exit(-1);
		}

		vector<Flag> middle_start_flags(start_flags);
		for (int i = 0; i < target_cipher.statesize; i++)
			if (start_flags[i] == "delta" && start_state[i] == 0)
				middle_start_flags[i] = "zero_c";

		vector<vector<vector<Flag>>> middle_rounds_flags;
		target_cipher.calculate_flags(start, end, middle_start_flags, middle_rounds_flags, true);

		// set env

		GRBEnv env = GRBEnv();

		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_IntParam_Threads, threads);
		if (target_cipher.ciphername == "kreyvium")
			env.set(GRB_IntParam_MIPFocus, 3);

		GRBModel model = GRBModel(env);

		// set initial variables
		vector<GRBVar> start_vars(statesize);

		

		target_cipher.set_callback_start_cstr(start, middle_start_flags, start_state, model, start_vars);

		vector< vector<vector<pair<BooleanPolynomial, GRBVar>> >> rounds_p_maps;
		vector<GRBVar> helper_midvars;
		target_cipher.callback_build_model_fast(model, start, start+step, middle_rounds_flags, start_vars, helper_midvars);
		vector<Flag> helper_midflags = middle_rounds_flags[start + step][0];
		vector<Flag> midflags = all_rounds_flags[start + step][0];

		GRBLinExpr obj = 0;
		vector<GRBVar> midvars(statesize);
		for (int i = 0; i < statesize; i++)
			if (midflags[i] == "delta")
			{
				midvars[i] = model.addVar(0, 1, 0, GRB_BINARY);
				if (helper_midflags[i] == "delta")
				{
					model.addConstr(midvars[i] >= helper_midvars[i]);
					obj += (midvars[i] - helper_midvars[i]);
				}
				else if (helper_midflags[i] == "zero_c")
					model.addConstr(midvars[i] == 0);
				else
					obj += midvars[i];
			}
		vector<GRBVar> end_vars(statesize);
		target_cipher.build_model_fast(model, start + step, end, all_rounds_flags, midvars, end_vars, rounds_p_maps, false);

		for (auto& round_p_maps : rounds_p_maps)
			for (auto& ele_p_map : round_p_maps)
				for (auto& p_var : ele_p_map)
					obj += p_var.second;


		model.setObjective(obj, GRB_MAXIMIZE);

		// impose constraints on final variables
		auto end_flags = all_rounds_flags[end][0];
		for (int i = 0; i < statesize; i++)
			if (end_state[i] == 1 && end_flags[i] != "delta")
			{
				cerr << __func__ << ": the end state is invalid";
				exit(-1);
			}
			else if (end_flags[i] == "delta")
				model.addConstr(end_vars[i] == end_state[i]);

		if (time > 0)
			model.set(GRB_DoubleParam_TimeLimit, time);

		vector<dynamic_bitset<> > callback_sols;
		vector<GRBVar>& callback_vars = midvars;
		vector<Flag>& callback_flags = midflags;
		ExpandCallback cb = ExpandCallback(callback_vars, callback_flags, callback_sols);
		model.setCallback(&cb);

		model.optimize();

		if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
		{
			// this may happen
			logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(start + step) +
				string(" | Failed "));

			return status::UNSOLVED;
		}
		else
		{
			int solCount = model.get(GRB_IntAttr_SolCount);
			if (solCount >= MAX)
			{
				cerr << "solCount value  is too big !" << endl;
				exit(-1);
			}

			int sol_num = callback_sols.size();

			double time = model.get(GRB_DoubleAttr_Runtime);
			logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(start + step) + string(" | "
			) + to_string(time) + string(" | "
			) + to_string(sol_num));


			// start to compute the coefficient
			for (auto& mid_state : callback_sols)
			{

				vector<dynamic_bitset<>> p_sols;
				vector<BooleanPolynomial> p_exps;
				vector<pair<int, int>> p_rms;
				set<int> end_constants;
				auto mid_status = target_cipher.solve_model(start,start+step, all_rounds_flags, start_state, mid_state, p_sols, p_exps, p_rms, end_constants, 1, 120, false);

				if (mid_status != status::SOLVED || end_constants.size() > 0)
				{
					cerr << __func__ << ": mid status error." << endl;
					exit(-1);
				}

				ListsOfPolynomialsAsFactors cof_lists = get_cof_lists(start, start + step, p_sols, p_exps, p_rms, 0, all_rounds_exps);
				state_coeff_lists[mid_state] = cof_lists;

				auto expand_res = expand_cof(cof_lists);
				expand_state_coeff[mid_state] = expand_res;
			}

			if (sol_num == 0)
				return status::NOSOLUTION;
			else
				return status::SOLVED;
		}


	}

	virtual status noncallback_forwardexpand(int start, int step, const dynamic_bitset<>& start_state,
		map<dynamic_bitset<>, BooleanPolynomial>& expand_state_coeff, map<dynamic_bitset<>, ListsOfPolynomialsAsFactors>& state_coeff_lists,
		double time, int threads)
	{
		int& statesize = target_cipher.statesize;
		auto start_flags = all_rounds_flags[start][0];
		if (start_flags.size() != statesize)
		{
			cerr << __func__ << ": The number of start flags is invalid.";
			exit(-1);
		}

		// set env

		GRBEnv env = GRBEnv();

		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, MAX);
		env.set(GRB_IntParam_Threads, threads);

		GRBModel model = GRBModel(env);

		// set initial variables
		vector<GRBVar> start_vars(statesize);

		target_cipher.set_start_cstr(start, start_flags, start_state, model, start_vars);

		// build model round by round
		vector< vector<vector<pair<BooleanPolynomial, GRBVar>>> > rounds_p_maps;
		vector<GRBVar> end_vars(statesize);
		target_cipher.build_model_fast(model, start, start + step, all_rounds_flags, start_vars, end_vars, rounds_p_maps, true);


		auto end_flags = all_rounds_flags[start+step][0];
		




		vector<BooleanPolynomial> p_exps;
		vector<pair<int, int>> p_rms;
		vector<GRBVar> p_vars;

		for (int r = 0; r < step; r++)
		{
			int list_index = target_cipher.cal_list_index(r + start);
			for (int i = 0; i < target_cipher.updatelists[list_index].size(); i++)
			{
				for (auto& e_v : rounds_p_maps[r][i])
				{
					p_exps.emplace_back(e_v.first);
					p_vars.emplace_back(e_v.second);
					p_rms.emplace_back(pair(r, i));
				}
			}
		}


		int p_num = p_vars.size();

		if (time > 0)
			model.set(GRB_DoubleParam_TimeLimit, time);

		model.optimize();


		if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
		{
			// this should not happen
			logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(start+ step) +
				string(" | Failed "));

			return status::UNSOLVED;
		}
		else
		{
			int solCount = model.get(GRB_IntAttr_SolCount);
			if (solCount >= MAX)
			{
				cerr << "solCount value  is too big !" << endl;
				exit(-1);
			}

			double time = model.get(GRB_DoubleAttr_Runtime);
			logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(start + step) + string(" | "
			) + to_string(time) + string(" | "
			) + to_string(solCount));

			if (solCount > 0)
			{

				// first collect solutions
				map<dynamic_bitset<>, map<dynamic_bitset<>, int> >  state_p_sols_counter;
				dynamic_bitset<> p_sol(p_num);
				dynamic_bitset<> state(statesize);

				for (int i = 0; i < solCount; i++)
				{
					model.set(GRB_IntParam_SolutionNumber, i);
					// read the input state
					for (int j = 0; j < statesize; j++)
						if (end_flags[j] == "delta")
							if (round(end_vars[j].get(GRB_DoubleAttr_Xn)) == 1)
								state[j] = 1;
							else
								state[j] = 0;
						else
							state[j] = 0;

					// read the p_sol
					for (int j = 0; j < p_num; j++)
					{
						if (round(p_vars[j].get(GRB_DoubleAttr_Xn)) == 1)
							p_sol[j] = 1;
						else
							p_sol[j] = 0;
					}

					state_p_sols_counter[state][p_sol]++;
				}

				for (auto& sp : state_p_sols_counter)
				{
					auto& state = sp.first;
					auto& p_sols_counter = sp.second;

					vector<dynamic_bitset<>> p_sols;
					for (auto& p_cnt : p_sols_counter)
						if (p_cnt.second % 2)
							p_sols.emplace_back(p_cnt.first);

					if (p_sols.size() == 0)
						continue;

					ListsOfPolynomialsAsFactors cof_lists = get_cof_lists(start, start + step, p_sols, p_exps, p_rms, 0, all_rounds_exps);
					state_coeff_lists[state] = cof_lists;


					auto expand_res = expand_cof(cof_lists);
					expand_state_coeff[state] = expand_res;
				}

				return status::SOLVED;
			}
			else
			{
				return status::NOSOLUTION;
			}
		}
	}

	virtual status callback_backexpand(int start, int end, int step, const dynamic_bitset<>& start_state, const dynamic_bitset<>& end_state,
		map<dynamic_bitset<>, BooleanPolynomial>& expand_state_coeff, map<dynamic_bitset<>, ListsOfPolynomialsAsFactors> & state_coeff_lists, double time, int threads)
	{
		int& statesize = target_cipher.statesize;
		auto start_flags = all_rounds_flags[start][0];
		if (start_flags.size() != statesize)
		{
			cerr << __func__ << ": The number of start flags is invalid.";
			exit(-1);
		}

		vector<Flag> middle_start_flags(start_flags);
		for (int i = 0; i < target_cipher.statesize; i++)
			if (start_flags[i] == "delta" && start_state[i] == 0)
				middle_start_flags[i] = "zero_c";

		vector<vector<vector<Flag>>> middle_rounds_flags;
		target_cipher.calculate_flags(start, end, middle_start_flags, middle_rounds_flags, true);

		// set env

		GRBEnv env = GRBEnv();

		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_IntParam_Threads, threads);
		if (target_cipher.ciphername == "kreyvium")
			env.set(GRB_IntParam_MIPFocus, 3);
		


		GRBModel model = GRBModel(env);

		// set initial variables
		vector<GRBVar> start_vars(statesize);

		target_cipher.set_callback_start_cstr(start, middle_start_flags, start_state, model, start_vars);


		vector< vector<vector<pair<BooleanPolynomial,GRBVar>> > > rounds_p_maps;

		vector<GRBVar> helper_midvars;
		int midround = end - step;
		if (target_cipher.ciphername == "acorn" && step < 30)
			midround = end - 30;

		target_cipher.callback_build_model_fast(model, start, midround, middle_rounds_flags, start_vars, helper_midvars);
		vector<Flag> helper_midflags = middle_rounds_flags[midround][0];
		vector<Flag> midflags = all_rounds_flags[midround][0];

		GRBLinExpr obj = 0;
		vector<GRBVar> midvars(statesize);
		for (int i = 0; i < statesize; i++)
			if (midflags[i] == "delta")
			{
				midvars[i] = model.addVar(0, 1, 0, GRB_BINARY);
				if (helper_midflags[i] == "delta")
				{
					model.addConstr(midvars[i] >= helper_midvars[i]);
					obj += (midvars[i] - helper_midvars[i]);
				}
				else if (helper_midflags[i] == "zero_c")
					model.addConstr(midvars[i] == 0);
				else
					obj += midvars[i];
			}

		rounds_p_maps.clear();
		vector<GRBVar> end_vars(statesize);
		vector<vector<GRBVar>> rounds_vars;
		if(target_cipher.ciphername == "acorn")
			target_cipher.build_model_fast_ret_round_variables(model, midround, end, all_rounds_flags, midvars, end_vars, rounds_p_maps, rounds_vars, false);
		else
			target_cipher.build_model_fast(model, midround, end, all_rounds_flags, midvars, end_vars, rounds_p_maps, false);

		for (auto& round_p_maps : rounds_p_maps)
			for (auto& ele_p_map : round_p_maps)
				for (auto& p_var : ele_p_map)
					obj += p_var.second;
		
		model.setObjective(obj, GRB_MAXIMIZE);

		// impose constraints on final variables
		auto end_flags = all_rounds_flags[end][0];
		for (int i = 0; i < statesize; i++)
			if (end_state[i] == 1 && end_flags[i] != "delta")
			{
				cerr << __func__ << ": the end state is invalid";
				exit(-1);
			}
			else if (end_flags[i] == "delta")
				model.addConstr(end_vars[i] == end_state[i]);

		if (time > 0)
			model.set(GRB_DoubleParam_TimeLimit, time);

		vector<dynamic_bitset<> > callback_sols;
		vector<GRBVar> callback_vars = midvars;
		vector<Flag> callback_flags = midflags;
		if (target_cipher.ciphername == "acorn")
		{
			callback_vars = rounds_vars[end - midround - step];
			callback_flags = all_rounds_flags[end - step][0];
		}

		ExpandCallback cb = ExpandCallback(callback_vars, callback_flags, callback_sols);
		model.setCallback(&cb);

		model.optimize();

		if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
		{
			// this may happen
			logger(__func__ + string(" : ") + to_string(end-step) + "-" + to_string(end) +
				string(" | Failed "));

			return status::UNSOLVED;
		}
		else
		{
			int solCount = model.get(GRB_IntAttr_SolCount);
			if (solCount >= MAX)
			{
				cerr << "solCount value  is too big !" << endl;
				exit(-1);
			}

			int sol_num = callback_sols.size();

			double time = model.get(GRB_DoubleAttr_Runtime);
			logger(__func__ + string(" : ") + to_string(end - step) + "-" + to_string(end) + string(" | "
			) + to_string(time) + string(" | "
			) + to_string(sol_num));


			// start to compute the coefficient
			for (auto& mid_state : callback_sols)
			{

				vector<dynamic_bitset<>> p_sols;
				vector<BooleanPolynomial> p_exps;
				vector<pair<int, int>> p_rms;
				set<int> end_constants;
				auto mid_status = target_cipher.solve_model(end-step, end, all_rounds_flags, mid_state,  end_state, p_sols, p_exps, p_rms, end_constants, 1, 120, false);

				if (mid_status != status::SOLVED || end_constants.size() > 0)
				{
					cerr << __func__<<": mid status error." << endl;
					exit(-1);
				}

				ListsOfPolynomialsAsFactors cof_lists = get_cof_lists(end - step, end, p_sols, p_exps, p_rms, 0, all_rounds_exps);
				state_coeff_lists[mid_state] = cof_lists;
				
				auto expand_res = expand_cof(cof_lists);
				expand_state_coeff[mid_state] = expand_res;
			}

			if (sol_num == 0)
				return status::NOSOLUTION;
			else
				return status::SOLVED;
		}


	}



	virtual status noncallback_backexpand(int rounds, int step, const dynamic_bitset<> & end_state, 
		map<dynamic_bitset<>, BooleanPolynomial> & expand_state_coeff, map<dynamic_bitset<>, ListsOfPolynomialsAsFactors> & state_coeff_lists,
		double time, int threads)
	{
		int& statesize = target_cipher.statesize;
		auto start_flags = all_rounds_flags[rounds - step][0];
		if (start_flags.size() != statesize)
		{
			cerr << __func__ << ": The number of start flags is invalid.";
			exit(-1);
		}

		// set env

		GRBEnv env = GRBEnv();

		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, MAX);
		env.set(GRB_IntParam_Threads, threads);

		GRBModel model = GRBModel(env);

		// set initial variables
		vector<GRBVar> start_vars(statesize);

		for (int i = 0; i < statesize; i++)
			if (start_flags[i] == "delta")
			{
				start_vars[i] = model.addVar(0, 1, 0, GRB_BINARY);
			}

		// build model round by round
		vector< vector<vector<pair<BooleanPolynomial, GRBVar>> > > rounds_p_maps;
		vector<GRBVar> end_vars(statesize);
		target_cipher.build_model_fast(model, rounds - step, rounds, all_rounds_flags, start_vars, end_vars, rounds_p_maps,true);


		auto end_flags = all_rounds_flags[rounds][0];
		for (int i = 0; i < statesize; i++)
			if (end_state[i] == 1 && end_flags[i] != "delta")
			{
				cerr << __func__ << ": the end state is invalid";
				exit(-1);
			}
			else if (end_flags[i] == "delta")
				model.addConstr(end_vars[i] == end_state[i]);




		vector<BooleanPolynomial> p_exps;
		vector<pair<int, int>> p_rms;
		vector<GRBVar> p_vars;

		for (int r = 0; r < step; r++)
		{
			int list_index = target_cipher.cal_list_index(r + rounds - step);
			for (int i = 0; i < target_cipher.updatelists[list_index].size(); i++)
			{
				for (auto& e_v : rounds_p_maps[r][i])
				{
					p_exps.emplace_back(e_v.first);
					p_vars.emplace_back(e_v.second);
					p_rms.emplace_back(pair(r, i));
				}
			}
		}


		int p_num = p_vars.size();

		if (time > 0)
			model.set(GRB_DoubleParam_TimeLimit, time);

		model.optimize();


		if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
		{
			// this should not happen
			logger(__func__ + string(" : ") + to_string(rounds - step) + "-" + to_string(rounds) +
				string(" | Failed "));

			return status::UNSOLVED;
		}
		else
		{
			int solCount = model.get(GRB_IntAttr_SolCount);
			if (solCount >= MAX)
			{
				cerr << "solCount value  is too big !" << endl;
				exit(-1);
			}

			double time = model.get(GRB_DoubleAttr_Runtime);
			logger(__func__ + string(" : ") + to_string(rounds - step) + "-" + to_string(rounds) + string(" | "
			) + to_string(time) + string(" | "
			) + to_string(solCount));

			if (solCount > 0)
			{

				// first collect solutions
				map<dynamic_bitset<>, map<dynamic_bitset<>, int> >  state_p_sols_counter;
				dynamic_bitset<> p_sol(p_num);
				dynamic_bitset<> state(statesize);

				for (int i = 0; i < solCount; i++)
				{
					model.set(GRB_IntParam_SolutionNumber, i);
					// read the input state
					for (int j = 0; j < statesize; j++)
						if (start_flags[j] == "delta")
							if (round(start_vars[j].get(GRB_DoubleAttr_Xn)) == 1)
								state[j] = 1;
							else
								state[j] = 0;
						else
							state[j] = 0;

					// read the p_sol
					for (int j = 0; j < p_num; j++)
					{
						if (round(p_vars[j].get(GRB_DoubleAttr_Xn)) == 1)
							p_sol[j] = 1;
						else
							p_sol[j] = 0;
					}

					state_p_sols_counter[state][p_sol]++;
				}

				for (auto& sp : state_p_sols_counter)
				{
					auto& state = sp.first;
					auto& p_sols_counter = sp.second;

					vector<dynamic_bitset<>> p_sols;
					for (auto& p_cnt : p_sols_counter)
						if (p_cnt.second % 2)
							p_sols.emplace_back(p_cnt.first);

					if (p_sols.size() == 0)
						continue;

					
					ListsOfPolynomialsAsFactors cof_lists = get_cof_lists(rounds - step, rounds + 1, p_sols, p_exps, p_rms, 0, all_rounds_exps);
					state_coeff_lists[state] = cof_lists;

					auto expand_res = expand_cof(cof_lists);		
					expand_state_coeff[state] = expand_res;
				}

				return status::SOLVED;
			}
			else
			{
				return status::NOSOLUTION;
			}
		}
	}

	virtual void forwardexpand_thread(int start, int end, int step, const dynamic_bitset<>& start_state, const dynamic_bitset<>& end_state, const BooleanPolynomial & cur_coef, const ListsOfPolynomialsAsFactors & cur_lists, 
		map<node_pair,BooleanPolynomial> & new_P, map<node_pair, ListsOfPolynomialsAsFactors> & new_P_asLists, 
		 double time, int threads, bool iscallback)
	{
		map<dynamic_bitset<>, BooleanPolynomial> expand_state_coeff;
		map<dynamic_bitset<>, ListsOfPolynomialsAsFactors> state_coeff_lists;

		if (iscallback)
		{
			auto expand_status = callback_forwardexpand(start, end, step, start_state, end_state, expand_state_coeff, state_coeff_lists, time, threads);
			if (expand_status == status::UNSOLVED)
			{
				auto expand_status2 = noncallback_forwardexpand(start, step, start_state, expand_state_coeff, state_coeff_lists, time, threads);
				if (expand_status2 == status::UNSOLVED)
				{
					cerr << __func__ << "Forwardexpand failed." << endl;
					exit(-1);
				}
			}
		}
		else
		{
			auto expand_status = noncallback_forwardexpand(start, step, start_state, expand_state_coeff, state_coeff_lists, time, threads);
			if (expand_status == status::UNSOLVED)
			{
				cerr << __func__ << "Forwardexpand failed." << endl;
				exit(-1);
			}
		}

		map<node_pair, BooleanPolynomial> part_P;
		map<node_pair, ListsOfPolynomialsAsFactors> part_P_asLists;
		for (auto& state_coef : expand_state_coeff)
		{
			auto& state = state_coef.first;
			auto& coef = state_coef.second;
			part_P[pair(state, end_state)] = coef * cur_coef;
		}

		for (auto& state_coef : state_coeff_lists)
		{
			auto& state = state_coef.first;
			auto& coef_lists = state_coef.second;
			part_P_asLists[pair(state, end_state)] = coef_lists * cur_lists;
		}



		// update the new_P
		{
			lock_guard<mutex> guard(expander_mutex0);
			for (auto& nodes_coef : part_P)
			{
				auto& nodes = nodes_coef.first;
				auto& coef = nodes_coef.second;
				auto it = new_P.find(nodes);
				if (it != new_P.end())
					it->second += coef;
				else
					new_P[nodes] = coef;
			}
		}

		// update the new_P_asLists
		{
			lock_guard<mutex> guard(expander_mutex1);
			for (auto& nodes_coef : part_P_asLists)
			{
				auto& nodes = nodes_coef.first;
				auto& coef_lists = nodes_coef.second;
				auto it = new_P_asLists.find(nodes);
				if (it != new_P_asLists.end())
					it->second += coef_lists;
				else
					new_P_asLists[nodes] = coef_lists;
			}
		}



	}



	virtual void forwardexpand(int start, int end, int step, const map<node_pair, BooleanPolynomial> & cur_P, const map<node_pair, ListsOfPolynomialsAsFactors>& cur_P_asLists, map<node_pair, BooleanPolynomial> & new_P, map<node_pair, ListsOfPolynomialsAsFactors> & new_P_asLists, double time, int threads, ThreadPool& thread_pool, bool iscallback)
	{
		vector<future<void>> futures;


		for (auto& nodes_coef : cur_P)
		{
			auto& nodepair = nodes_coef.first;
			auto& start_state = nodepair.first;
			auto& end_state = nodepair.second;
			auto& cur_coef = nodes_coef.second;
			auto& cur_coef_lists = cur_P_asLists.at(nodepair);
			futures.emplace_back(thread_pool.Submit(&framework::forwardexpand_thread, REF(*this), start, end, step, REF(start_state), REF(end_state),
				REF(cur_coef), REF(cur_coef_lists), REF(new_P), REF(new_P_asLists), time, threads, iscallback));
		}

		for (auto& it : futures)
			it.get();
	}

	virtual void backexpand_thread(int start, int end, int step, const dynamic_bitset<>& start_state, const dynamic_bitset<>& end_state,
		const BooleanPolynomial& cur_coef, const ListsOfPolynomialsAsFactors & cur_lists, map<node_pair, BooleanPolynomial>& new_P,
		map<node_pair, ListsOfPolynomialsAsFactors> & new_P_asLists, double time, int threads, bool iscallback)
	{
		map<dynamic_bitset<>, BooleanPolynomial> expand_state_coeff;
		map<dynamic_bitset<>, ListsOfPolynomialsAsFactors> state_coeff_lists;

		
		if (iscallback)
		{
			auto expand_status = callback_backexpand(start, end, step, start_state, end_state, expand_state_coeff, state_coeff_lists, time, threads);
			if (expand_status == status::UNSOLVED)
			{
				auto expand_status2 = noncallback_backexpand(end, step, end_state, expand_state_coeff, state_coeff_lists, time, threads);
				if (expand_status2 == status::UNSOLVED)
				{
					cerr << __func__ << "Backexpand failed." << endl;
					exit(-1);
				}
			}
		}
		else
		{
			auto expand_status = noncallback_backexpand(end, step, end_state, expand_state_coeff, state_coeff_lists, time, threads);
			if (expand_status == status::UNSOLVED)
			{
				cerr << __func__ << "Backexpand failed." << endl;
				exit(-1);
			}
		}

		map<node_pair, BooleanPolynomial> part_P;
		for (auto& state_coef : expand_state_coeff)
		{
			auto& state = state_coef.first;
			auto& coef = state_coef.second;
			part_P[pair(start_state, state)] = coef * cur_coef;
		}

		map<node_pair, ListsOfPolynomialsAsFactors> part_P_asLists;
		for (auto& state_coef : state_coeff_lists)
		{
			auto& state = state_coef.first;
			auto& coef_lists = state_coef.second;
			part_P_asLists[pair(start_state, state)] = coef_lists * cur_lists;
		}


		// update the new_P
		{
			lock_guard<mutex> guard(expander_mutex0);
			for (auto& nodes_coef : part_P)
			{
				auto& nodes = nodes_coef.first;
				auto& coef = nodes_coef.second;
				auto it = new_P.find(nodes);
				if (it != new_P.end())
					it->second += coef;
				else
					new_P[nodes] = coef;
			}
		}

		// update the new_P_asLists
		{
			lock_guard<mutex> guard(expander_mutex1);
			for (auto& nodes_coef : part_P_asLists)
			{
				auto& nodes = nodes_coef.first;
				auto& coef_lists = nodes_coef.second;
				auto it = new_P_asLists.find(nodes);
				if (it != new_P_asLists.end())
					it->second += coef_lists;
				else
					new_P_asLists[nodes] = coef_lists;
			}
		}





	}



	virtual void backexpand(int start, int end, int step, const map<node_pair, BooleanPolynomial>& cur_P, const map<node_pair, ListsOfPolynomialsAsFactors> & cur_P_asLists, map<node_pair, BooleanPolynomial>& new_P, 
		map<node_pair, ListsOfPolynomialsAsFactors> & new_P_asLists, double time, int threads, ThreadPool& thread_pool, bool iscallback)
	{
		vector<future<void>> futures;

	
		for (auto& nodes_coef : cur_P)
		{
			auto& nodepair = nodes_coef.first;
			auto& start_state = nodepair.first;
			auto& end_state = nodepair.second;
			auto& cur_coef = nodes_coef.second;
			auto& cur_coef_lists = cur_P_asLists.at(nodepair);
			futures.emplace_back(thread_pool.Submit(&framework::backexpand_thread, REF(*this), start, end, step, REF(start_state), REF(end_state),
				REF(cur_coef), REF(cur_coef_lists), REF(new_P), REF(new_P_asLists), time, threads, iscallback));
		}

		for (auto& it : futures)
			it.get();

	}



	

	// the function for output to file
	void print_sol(int start, int end, const dynamic_bitset<>& start_state, const dynamic_bitset<>& end_state,
		vector<dynamic_bitset<>> p_sols, vector<BooleanPolynomial>& p_exps,
		vector<pair<int, int>>& p_rms, const set<int> & end_constants) 
	{
		if (p_sols.size() == 0)
			return;

		int p_num = p_sols[0].size();

		stringstream ss;
		ss << this_thread::get_id();
		string thread_id;
		ss >> thread_id;

		string path = string("TERM/") + to_string(start) + "_" + to_string(end) + "_" + thread_id + string(".txt");
		ofstream os;
		os.open(path, ios::out | ios::app);
		os << "state from:" << start_state << endl;
		os << "state to:" << end_state << endl;

		dynamic_bitset<> p_sol_mask(p_num);
		for (auto& p_sol : p_sols)
		{
			p_sol_mask |= p_sol;
		}

		for(int i = 0;i < p_num;i++)
			if (p_sol_mask[i] == 1)
			{
				os << "p" << i<< "=" << p_rms[i].first << "-" << p_rms[i].second << "-" << p_exps[i] << " ";
			}

		os << endl;

		for (auto& p_sol : p_sols)
		{
			os << p_sol << endl;
		}

		os << "end constants : ";
		for (auto& i : end_constants)
			os << i << " ";
		os << endl;

		os << endl;
		os.close();
	}

	void print_sol(int start, int end, const ListsOfPolynomialsAsFactors& coeff_lists) 
	{
		

		stringstream ss;
		ss << this_thread::get_id();
		string thread_id;
		ss >> thread_id;

		string path = string("TERM/") + to_string(start) + "_" + to_string(end) + "_" + thread_id + string(".txt");
		ofstream os;
		os.open(path, ios::out | ios::app);
		coeff_lists.display(os);
		os.close();
	}

	void print_sol_debug(int start, int end, const dynamic_bitset<>& start_state,
		const dynamic_bitset<>& end_state, const ListsOfPolynomialsAsFactors& cur_lists, const ListsOfPolynomialsAsFactors & second_coeff_lists, const ListsOfPolynomialsAsFactors& coeff_lists)
	{
		stringstream ss;
		ss << this_thread::get_id();
		string thread_id;
		ss >> thread_id;

		string path = string("DEBUG/") + to_string(start) + "_" + to_string(end) + "_" + thread_id + string(".txt");
		ofstream os;
		os.open(path, ios::out | ios::app);
		os << start_state << endl;
		os << end_state << endl;
		cur_lists.display(os);
		second_coeff_lists.display(os);
		coeff_lists.display(os);
		os.close();
	}


	ListsOfPolynomialsAsFactors get_cof_lists(int start, int end,
		const vector<dynamic_bitset<>>& p_sols, const vector<BooleanPolynomial>& p_exps, const vector<pair<int, int>>& p_rms,
		int expandto, vector<vector<vector<BooleanPolynomial> >>& expandto_exps)
	{
		ListsOfPolynomialsAsFactors sols_lists(target_cipher.statesize);

		if (p_sols.size() == 0)
			return sols_lists;

		int p_num = p_sols[0].size();
		dynamic_bitset<> p_sol_mask(p_num);
		for (auto& p_sol : p_sols)
		{
			p_sol_mask |= p_sol;
		}

		vector<BooleanPolynomial> p_expands(p_num);

		for (int i = 0; i < p_num; i++)
		{
			if (p_sol_mask[i] == 1)
			{
				auto p_expand = p_exps[i].subs(expandto_exps[p_rms[i].first + start - expandto][p_rms[i].second]);
				p_expands[i] = p_expand;
			}

			
		}

		for (auto& p_sol : p_sols)
		{
			ListOfPolynomialsAsFactors sol_list(target_cipher.statesize);
			for (int i = 0; i < p_sol.size(); i++)
				if (p_sol[i] == 1)
				{
					sol_list.add(p_expands[i]);
				}

			sols_lists.add(sol_list);
		}

		return sols_lists;




	}


	// the function for expanding the coefficient
	BooleanPolynomial expand_cof(int start, int end, 
		const vector<dynamic_bitset<>> & p_sols, const vector<BooleanPolynomial>& p_exps, const vector<pair<int, int>>& p_rms,
		int expandto, vector<vector<vector<BooleanPolynomial> >> & expandto_exps )
	{
		ListsOfPolynomialsAsFactors lists = get_cof_lists(start, end, p_sols, p_exps, p_rms, expandto, expandto_exps);
		return lists.getSum();
	}

	BooleanPolynomial expand_cof(const ListsOfPolynomialsAsFactors& lists)
	{
		return lists.getSum();
	}

	// calculate the degree of non-cube vars for a nodepair
	virtual void cal_noncube_deg_nodepair(int start, int end, const vector<vector<vector<Flag>>>& rounds_flags, const dynamic_bitset<>& start_state,
		const dynamic_bitset<>& end_state, int threads, double time,
		const set<BooleanPolynomial>& excluded_expanded_p_exps, const set<BooleanPolynomial>& included_expanded_p_exps, 
		const set<int> & state_noncube_indices, int& ret_deg, string tag)
	{
		// this function assume that at round end the state bits are all delta bits

		if (start != 0)
		{
			throw Exception(__func__ + string(" : This function only works for start being 0."));
		}

		auto& start_flags = rounds_flags[start][0];
		auto& end_flags = rounds_flags[end][0];

		GRBEnv env = GRBEnv();

		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_Threads, threads);
		if (target_cipher.ciphername == "kreyvium")
			env.set(GRB_IntParam_MIPFocus, 3);

		GRBModel model = GRBModel(env);

		// set initial variables
		vector<GRBVar> start_vars(target_cipher.statesize);

		target_cipher.set_start_cstr(start, start_flags, start_state, model, start_vars);

		// from start to end
		vector< vector<vector<pair<BooleanPolynomial, GRBVar>> >> rounds_p_maps;

		vector<GRBVar> end_vars(target_cipher.statesize);
		target_cipher.build_model_fast(model, start, end, rounds_flags, start_vars, end_vars, rounds_p_maps, false, -1);

		// impose constraints on final variables
		for (int i = 0; i < target_cipher.statesize; i++)
			if (end_state[i] == 1 && end_flags[i] == "zero_c")
			{
				throw Exception(__func__ + string(": non delta bits in end state detected."));
			}

		for (int i = 0; i < target_cipher.statesize; i++)
			if (end_flags[i] == "delta")
			{
				model.addConstr(end_vars[i] == end_state[i]);
			}

		// collect p exps and p vars
		vector<GRBVar> p_vars;
		vector<BooleanPolynomial> p_exps;
		vector<pair<int, int>> p_rms;
		for (int r = 0; r < end - start; r++)
		{
			int list_index = target_cipher.cal_list_index(r + start);
			for (int i = 0; i < target_cipher.updatelists[list_index].size(); i++)
			{
				for (auto& e_v : rounds_p_maps[r][i])
				{
					p_exps.emplace_back(e_v.first);
					p_vars.emplace_back(e_v.second);
					p_rms.emplace_back(pair(r, i));
				}
			}

		}

		int p_num = p_vars.size();

		// calculate expanded p exps
		vector<BooleanPolynomial> expanded_p_exps;
		for (int i = 0; i < p_num; i++)
		{
			auto expanded_p_exp = p_exps[i].subs(all_rounds_exps[p_rms[i].first][p_rms[i].second]);
			expanded_p_exps.emplace_back(expanded_p_exp);

			if (expanded_p_exp.iszero())
				model.addConstr(p_vars[i] == 0);
		}

		// add constraints according to excluded and included expanded p exps
		vector<int> free_p_indices;
		vector<GRBVar> included_p_vars;
		for (int i = 0; i < p_num; i++)
		{
			if (excluded_expanded_p_exps.find(expanded_p_exps[i]) != excluded_expanded_p_exps.end())
			{
				model.addConstr(p_vars[i] == 0);
				// cout << "Exclude " << i << " " << expanded_p_exps[i] << endl;
			}
			else if (included_expanded_p_exps.find(expanded_p_exps[i]) != included_expanded_p_exps.end())
				included_p_vars.emplace_back(p_vars[i]);
			else
				free_p_indices.emplace_back(i);
		}

		if (included_p_vars.size() != 0)
		{
			GRBLinExpr sum = 0;
			for (auto& p_var : included_p_vars)
				sum += p_var;
			model.addConstr(sum >= 1);
		}


		// simplify expanded exps according to excluded and included p exps
		set<BooleanPolynomial> constrained_expanded_p_exps(excluded_expanded_p_exps);
		constrained_expanded_p_exps.insert(included_expanded_p_exps.begin(), included_expanded_p_exps.end());

		vector<BooleanPolynomial> cur_constrains_map = normal_exp;

		for (auto& expanded_p_exp : constrained_expanded_p_exps)
		{
			auto expanded_p_exp_after_sub = expanded_p_exp.subs(cur_constrains_map);

			map<int, vector<BooleanMonomial>> var_mons_distribution;
			for (auto& mon : expanded_p_exp_after_sub)
			{
				for (auto& v : mon.index())
					var_mons_distribution[v].emplace_back(mon);
			}

			bool is_noncube_involved = false;
			bool is_balanced = false;
			int balanced_v = 0;
			BooleanMonomial balanced_mon;
			for(auto& v_mons : var_mons_distribution)
				if (state_noncube_indices.find(v_mons.first) != state_noncube_indices.end())
				{
					is_noncube_involved = true;
					if (v_mons.second.size() == 1 && v_mons.second[0].count() == 1)
					{
						is_balanced = true;
						balanced_v = v_mons.first;
						balanced_mon = v_mons.second[0];
						break;
					}
				}

			if (is_noncube_involved)
				if (is_balanced)
				{
					cur_constrains_map[balanced_v] = expanded_p_exp_after_sub + balanced_mon;
				}
				else
					throw Exception(__func__ + string("Unbalanced constraints based on p detected."));
		}


		map<int, vector<GRBVar>> copy_to_per_noncube_var;
		for (auto& free_p_index : free_p_indices)
		{
			auto free_expanded_p_exp_after_sub = expanded_p_exps[free_p_index].subs(cur_constrains_map);
			// cout << expanded_p_exps[free_p_index] << endl;
			// cout << free_expanded_p_exp_after_sub << endl;
			// cout << endl;

			GRBLinExpr sum_GRBvars_all_mons = 0;
			for (auto& mon : free_expanded_p_exp_after_sub)
			{
				GRBVar GRBvar_this_mon = model.addVar(0, 1, 0, GRB_BINARY);
				sum_GRBvars_all_mons += GRBvar_this_mon;
				for(auto& v : mon.index())
					if (state_noncube_indices.find(v) != state_noncube_indices.end())
					{
						copy_to_per_noncube_var[v].emplace_back(GRBvar_this_mon);
					}
			}

			model.addConstr(sum_GRBvars_all_mons == p_vars[free_p_index]);
		}

		GRBLinExpr obj_exp = 0;
		for (auto& copy_to_this_noncube_var : copy_to_per_noncube_var)
		{
			GRBVar GRBvar_this_noncube_var = model.addVar(0, 1, 0, GRB_BINARY);
			obj_exp += GRBvar_this_noncube_var;
			if (copy_to_this_noncube_var.second.size() > 0)
				model.addGenConstrOr(GRBvar_this_noncube_var, &(copy_to_this_noncube_var.second[0]), copy_to_this_noncube_var.second.size());
			else
				model.addConstr(GRBvar_this_noncube_var == 0);
		}

		{
			lock_guard<mutex> guard(solver_mutex0);
			model.addConstr(obj_exp >= ret_deg + 1);
		}
		model.setObjective(obj_exp, GRB_MAXIMIZE);
		// 

		if (time > 0)
			model.set(GRB_DoubleParam_TimeLimit, time);


		model.optimize();

		if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
		{
			string startstr, endstr;
			to_string(start_state, startstr);
			to_string(end_state, endstr);

			double timecost = model.get(GRB_DoubleAttr_Runtime);
			logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(end) + " | "
				+ to_string(timecost) + string(" | Calculate noncube degree failed | ") + startstr + string(" | ") + endstr + tag);

			cal_noncube_deg_nodepair(start, end, rounds_flags, start_state, end_state, threads, time, excluded_expanded_p_exps, included_expanded_p_exps, state_noncube_indices, ret_deg, tag);
			return;
		}
		else if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
		{
			string startstr, endstr;
			to_string(start_state, startstr);
			to_string(end_state, endstr);
			double timecost = model.get(GRB_DoubleAttr_Runtime);
			logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(end) + " | "
				+ to_string(timecost) + string(" | No solution ") + tag);
			return;
		}
		else if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
		{
			string startstr, endstr;
			to_string(start_state, startstr);
			to_string(end_state, endstr);
			auto obj = model.getObjective();
			double timecost = model.get(GRB_DoubleAttr_Runtime);
			logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(end) + " | "
				+ to_string(timecost) + string(" | Solved | ") + to_string(round(obj.getValue())) + tag );

			lock_guard<mutex> guard(solver_mutex0);
			if (int(round(obj.getValue())) > ret_deg)
			{
				cout << "Update deg from " << ret_deg << " to " << int(round(obj.getValue())) << tag << endl;
				ret_deg = int(round(obj.getValue()));
			}


			return;
		}
		else
		{
			throw Exception(__func__ + string(": unknown model status detected."));
		}


	}

	// calculate the degree of p vars for a nodepair
	virtual void cal_p_deg_nodepair(int start, int end, const vector<vector<vector<Flag>>>& rounds_flags, const dynamic_bitset<>& start_state,
		const dynamic_bitset<>& end_state, int threads, double time, vector<BooleanPolynomial> & all_possible_expanded_p_exps,  const set<BooleanPolynomial>& excluded_expanded_p_exps, const set<BooleanPolynomial>& included_expanded_p_exps, int& ret_deg, bool set_all_possible_expanded_p_exps, string tag)
	{
		// this function assume that at round end the state bits are all delta bits

		if (start != 0)
		{
			throw Exception(__func__ + string(" : This function only works for start being 0."));
		}

		auto& start_flags = rounds_flags[start][0];
		auto& end_flags = rounds_flags[end][0];

		GRBEnv env = GRBEnv();

	    env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_Threads, threads);
		if (target_cipher.ciphername == "kreyvium")
			env.set(GRB_IntParam_MIPFocus, 3);

		GRBModel model = GRBModel(env);

		// set initial variables
		vector<GRBVar> start_vars(target_cipher.statesize);

		target_cipher.set_start_cstr(start, start_flags, start_state, model, start_vars);

		// from start to end
		vector< vector<vector<pair<BooleanPolynomial, GRBVar>> >> rounds_p_maps;

		vector<GRBVar> end_vars(target_cipher.statesize);
		target_cipher.build_model_fast(model, start, end, rounds_flags, start_vars, end_vars, rounds_p_maps, false, -1);

		// impose constraints on final variables
		for (int i = 0; i < target_cipher.statesize; i++)
			if (end_state[i] == 1 && end_flags[i] == "zero_c")
			{
				throw Exception(__func__ + string(": non delta bits in end state detected."));
			}

		for (int i = 0; i < target_cipher.statesize; i++)
			if (end_flags[i] == "delta")
			{
				model.addConstr(end_vars[i] == end_state[i]);
			}

		// collect p exps and p vars
		vector<GRBVar> p_vars;
		vector<BooleanPolynomial> p_exps;
		vector<pair<int, int>> p_rms;
		for (int r = 0; r < end - start; r++)
		{
			int list_index = target_cipher.cal_list_index(r + start);
			for (int i = 0; i < target_cipher.updatelists[list_index].size(); i++)
			{
				for (auto& e_v : rounds_p_maps[r][i])
				{
					p_exps.emplace_back(e_v.first);
					p_vars.emplace_back(e_v.second);
					p_rms.emplace_back(pair(r, i));
				}
			}

		}

		int p_num = p_vars.size();

		// calculate expanded p exps
		vector<BooleanPolynomial> expanded_p_exps;
		for (int i = 0; i < p_num; i++)
		{
			auto expanded_p_exp = p_exps[i].subs(all_rounds_exps[p_rms[i].first][p_rms[i].second]);
			expanded_p_exps.emplace_back(expanded_p_exp);
			if (expanded_p_exp.iszero())
				model.addConstr(p_vars[i] == 0);
		}

		if (set_all_possible_expanded_p_exps)
			all_possible_expanded_p_exps = expanded_p_exps;

		// add constraints according to excluded and included expanded p exps
		GRBLinExpr obj_p_vars = 0;
		vector<GRBVar> included_p_vars;
		for (int i = 0; i < p_num; i++)
		{
			if (excluded_expanded_p_exps.find(expanded_p_exps[i]) != excluded_expanded_p_exps.end())
			{
				model.addConstr(p_vars[i] == 0);
			}
			else if (included_expanded_p_exps.find(expanded_p_exps[i]) != included_expanded_p_exps.end())
				included_p_vars.emplace_back(p_vars[i]);
			else
				obj_p_vars += p_vars[i];
		}

		if (included_p_vars.size() != 0)
		{
			GRBLinExpr sum = 0;
			for (auto& p_var : included_p_vars)
				sum += p_var;
			model.addConstr(sum >= 1);
		}
		
		{
			lock_guard<mutex> guard(solver_mutex0);
			model.addConstr(obj_p_vars >= ret_deg + 1);
		}

		model.setObjective(obj_p_vars, GRB_MAXIMIZE);

		if (time > 0)
			model.set(GRB_DoubleParam_TimeLimit, time);


		model.optimize();

		if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
		{
			string startstr, endstr;
			to_string(start_state, startstr);
			to_string(end_state, endstr);

			double timecost = model.get(GRB_DoubleAttr_Runtime);
			logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(end) + " | " 
			+ to_string(timecost) + string(" | Calculate p degree failed | ") + startstr + string(" | ") + endstr + tag);

			return;
		}
		else if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
		{
			string startstr, endstr;
			to_string(start_state, startstr);
			to_string(end_state, endstr);
			double timecost = model.get(GRB_DoubleAttr_Runtime);
			logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(end) + " | "
				+ to_string(timecost) + string(" | No solution ") + tag);
			return;
		}
		else if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
		{
			string startstr, endstr;
			to_string(start_state, startstr);
			to_string(end_state, endstr);
			auto obj = model.getObjective();
			double timecost = model.get(GRB_DoubleAttr_Runtime);
			logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(end) + " | "
				+ to_string(timecost) + string(" | Solved | ") + to_string(round(obj.getValue())) + tag);

			lock_guard<mutex> guard(solver_mutex0);
			if (true)
			{
				if (int(round(obj.getValue())) > ret_deg)
				{
					cout << "Update deg from " << ret_deg << " to " << int(round(obj.getValue())) << tag << endl;
					ret_deg = int(round(obj.getValue()));
				}
			}

			
			return;
		}
		else
		{
			throw Exception(__func__ + string(": unknown model status detected."));
		}




	}

	virtual void get_all_possible_expanded_p_exps(vector<BooleanPolynomial> & all_possible_expanded_p_exps)
	{
		dynamic_bitset<> start_state(target_cipher.statesize);
		dynamic_bitset<> end_state(target_cipher.statesize);
		int ret_deg = 0;
		set<BooleanPolynomial> excluded_expanded_p_exps, included_expanded_p_exps;
		cal_p_deg_nodepair(0, target_rounds, all_rounds_flags, start_state, end_state, 2, 1, all_possible_expanded_p_exps, excluded_expanded_p_exps, included_expanded_p_exps, ret_deg, true, R"~(AllPossible)~");
	}


	virtual void search_optimal_constraints(int start, int end, const map<node_pair, BooleanPolynomial>& cur_P, 
		double time, int threads, ThreadPool& thread_pool, int num_of_constraints)
	{
		vector<BooleanPolynomial> all_possible_expanded_p_exps;
		get_all_possible_expanded_p_exps(all_possible_expanded_p_exps);
		map<BooleanPolynomial, int> expanded_p_exp_to_index;
		for (int i = 0; i < all_possible_expanded_p_exps.size(); i++)
			expanded_p_exp_to_index[all_possible_expanded_p_exps[i]] = i;


		logger("There are in total " + to_string(all_possible_expanded_p_exps.size()) + string(" expanded p exps."));

		for (auto& p_exp : all_possible_expanded_p_exps)
		{
			if (p_exp.isone() || p_exp.iszero())
			{
				cout << "Detect one or zero polynomial." << endl;
			}
		}

		set<int> state_noncube_indices;
		for (auto& noncube_exp : noncube_exps)
		{
			for (auto& mon : noncube_exp)
			{
				auto tmp = mon.index();
				state_noncube_indices.insert(tmp.begin(), tmp.end());
			}
		}

		int depth = num_of_constraints;
		int loop = 0;
		set<BooleanPolynomial> base_excluded_expanded_p_exps = {
		};

		while (loop < depth)
		{
			
			set<BooleanPolynomial> constrained_expanded_p_exps(base_excluded_expanded_p_exps);

			vector<BooleanPolynomial> cur_constrains_map = normal_exp;

			for (auto& expanded_p_exp : constrained_expanded_p_exps)
			{
				auto expanded_p_exp_after_sub = expanded_p_exp.subs(cur_constrains_map);

				map<int, vector<BooleanMonomial>> var_mons_distribution;
				for (auto& mon : expanded_p_exp_after_sub)
				{
					for (auto& v : mon.index())
						var_mons_distribution[v].emplace_back(mon);
				}

				bool is_noncube_involved = false;
				bool is_balanced = false;
				int balanced_v = 0;
				BooleanMonomial balanced_mon;
				for (auto& v_mons : var_mons_distribution)
					if (state_noncube_indices.find(v_mons.first) != state_noncube_indices.end())
					{
						is_noncube_involved = true;
						if (v_mons.second.size() == 1 && v_mons.second[0].count() == 1)
						{
							is_balanced = true;
							balanced_v = v_mons.first;
							balanced_mon = v_mons.second[0];
							break;
						}
					}

				if (is_noncube_involved)
					if (is_balanced)
					{
						cur_constrains_map[balanced_v] = expanded_p_exp_after_sub + balanced_mon;
					}
					else
						throw Exception(__func__ + string("Unbalanced constraints based on p detected."));
			}

			for (auto& expanded_p_exp : all_possible_expanded_p_exps)
			{
				if (expanded_p_exp.iszero() || expanded_p_exp.isone())
					continue;
				if (base_excluded_expanded_p_exps.find(expanded_p_exp) != base_excluded_expanded_p_exps.end())
					continue;

				auto expanded_p_exp_after_sub = expanded_p_exp.subs(cur_constrains_map);
				if (expanded_p_exp_after_sub.iszero())
					base_excluded_expanded_p_exps.emplace(expanded_p_exp);
			}
			

			logger("Start loop " + to_string(loop));
			logger("Current excluded p exps : ");
			for (auto& expanded_p_exp : base_excluded_expanded_p_exps)
			{
				stringstream ss;
				ss << expanded_p_exp;
				logger(ss.str());
			}

			set<BooleanPolynomial> have_tested_expanded_p_exps;
			int cur_optimal_deg3 = -1;
			int cur_optimal_deg4 = -1;
			BooleanPolynomial next_excluded_expanded_p_exp;
			string next_excluded_expanded_p_exp_str;


			// search for next excluded expanded p exp
			for (auto& expanded_p_exp : all_possible_expanded_p_exps)
			{


				if (expanded_p_exp.iszero() || expanded_p_exp.isone())
					continue;
				if (base_excluded_expanded_p_exps.find(expanded_p_exp) != base_excluded_expanded_p_exps.end())
					continue;

				bool is_noncube_involved = false;
				for (auto& mon : expanded_p_exp)
				{
					for (auto& v : mon.index())
						if (state_noncube_indices.find(v) != state_noncube_indices.end())
						{
							is_noncube_involved = true;
							break;
						}

					if (is_noncube_involved)
						break;
				}

				//if(!is_noncube_involved)
					// continue;

				int base_deg3 = cur_optimal_deg3;
				int base_deg4 = cur_optimal_deg4;

				if (have_tested_expanded_p_exps.find(expanded_p_exp) == have_tested_expanded_p_exps.end())
				{
					stringstream ss;
					ss << expanded_p_exp;

					logger("Test expanded p exp : " + ss.str());
					vector<future<void>> futures;


					set<BooleanPolynomial> excluded_expanded_p_exps(base_excluded_expanded_p_exps);
					set<BooleanPolynomial> included_expanded_p_exps = { expanded_p_exp };
					set<BooleanPolynomial> excluded_expanded_p_exps2(base_excluded_expanded_p_exps);
					excluded_expanded_p_exps2.emplace(expanded_p_exp);
					set<BooleanPolynomial> included_expanded_p_exps2;

					for (auto& nodes_coef : cur_P)
					{
						auto& nodepair = nodes_coef.first;
						auto& start_state = nodepair.first;
						auto& end_state = nodepair.second;

						futures.emplace_back(thread_pool.Submit(&framework::cal_p_deg_nodepair, REF(*this), start, end, REF(all_rounds_flags), REF(start_state), REF(end_state), threads, time, REF(all_possible_expanded_p_exps), REF(excluded_expanded_p_exps), REF(included_expanded_p_exps), REF(base_deg3), false,  R"~([1])~"));
						futures.emplace_back(thread_pool.Submit(&framework::cal_p_deg_nodepair, REF(*this), start, end, REF(all_rounds_flags), REF(start_state), REF(end_state), threads, time, REF(all_possible_expanded_p_exps), REF(excluded_expanded_p_exps2), REF(included_expanded_p_exps2), REF(base_deg4), false, R"~([0])~"));
					}

					for (auto& it : futures)
						it.get();

					if (cur_optimal_deg3 < base_deg3 && cur_optimal_deg4 == base_deg4)
					{
						logger("Better bound on base deg4 detected.");
						futures.clear();
						base_deg4 = 0;
						for (auto& nodes_coef : cur_P)
						{
							auto& nodepair = nodes_coef.first;
							auto& start_state = nodepair.first;
							auto& end_state = nodepair.second;
							futures.emplace_back(thread_pool.Submit(&framework::cal_p_deg_nodepair, REF(*this), start, end, REF(all_rounds_flags), REF(start_state), REF(end_state), threads, time, REF(all_possible_expanded_p_exps), REF(excluded_expanded_p_exps2), REF(included_expanded_p_exps2), REF(base_deg4), false, R"~([0])~"));
						}
						for (auto& it : futures)
							it.get();

						logger("Update the p or noncube deg to " + to_string(base_deg3) + " and " + to_string(base_deg4));
						cur_optimal_deg3 = base_deg3;
						cur_optimal_deg4 = base_deg4;
						next_excluded_expanded_p_exp = expanded_p_exp;
						next_excluded_expanded_p_exp_str = ss.str();
					}
					else if (cur_optimal_deg3 == base_deg3 && cur_optimal_deg4 == base_deg4)
					{
						logger("Better bound on base deg4 detected.");
						futures.clear();
						base_deg4 = -1;
						base_deg3 = -1;
						for (auto& nodes_coef : cur_P)
						{
							auto& nodepair = nodes_coef.first;
							auto& start_state = nodepair.first;
							auto& end_state = nodepair.second;
							futures.emplace_back(thread_pool.Submit(&framework::cal_p_deg_nodepair, REF(*this), start, end, REF(all_rounds_flags), REF(start_state), REF(end_state), threads, time, REF(all_possible_expanded_p_exps), REF(excluded_expanded_p_exps), REF(included_expanded_p_exps), REF(base_deg3), false, R"~([1])~"));
							futures.emplace_back(thread_pool.Submit(&framework::cal_p_deg_nodepair, REF(*this), start, end, REF(all_rounds_flags), REF(start_state), REF(end_state), threads, time, REF(all_possible_expanded_p_exps), REF(excluded_expanded_p_exps2), REF(included_expanded_p_exps2), REF(base_deg4), false, R"~([0])~"));
						}
						for (auto& it : futures)
							it.get();

						if (base_deg3 == cur_optimal_deg3 && base_deg4 < cur_optimal_deg4)
						{
							logger("Update the p or noncube deg to " + to_string(base_deg3) + " and " + to_string(base_deg4));
							cur_optimal_deg3 = base_deg3;
							cur_optimal_deg4 = base_deg4;
							next_excluded_expanded_p_exp = expanded_p_exp;
							next_excluded_expanded_p_exp_str = ss.str();
						}
					}
					else if (cur_optimal_deg3 == -1 && cur_optimal_deg4 == -1)
					{
						logger("deg4 is not initialized yet.");
						logger("Update the p or noncube deg to " + to_string(base_deg3) + " and " + to_string(base_deg4));
						cur_optimal_deg3 = base_deg3;
						cur_optimal_deg4 = base_deg4;
						next_excluded_expanded_p_exp = expanded_p_exp;
						next_excluded_expanded_p_exp_str = ss.str();
					}
					else if (cur_optimal_deg3 == -1 || cur_optimal_deg4 == -1)
						throw Exception("deg3 or deg4 is not initialized. Error detected.");
					else
						logger("The p or noncube degree is not updated.");


					have_tested_expanded_p_exps.emplace(expanded_p_exp);
				}

			}

			base_excluded_expanded_p_exps.emplace(next_excluded_expanded_p_exp);


			logger("loop " + to_string(loop) + " exclude " + next_excluded_expanded_p_exp_str + " with quotient deg " + to_string(cur_optimal_deg3) + " and remainder deg " + to_string(cur_optimal_deg4));
			logger("The index of excluded p exp is " + to_string(expanded_p_exp_to_index[next_excluded_expanded_p_exp]) );
			
			loop++;
		}

		cout << "The identified optimal constraints (set to 0) are : " << endl;
		for (auto& expanded_p_exp : base_excluded_expanded_p_exps)
		{
			cout << expanded_p_exp << endl;
		}
	}

	// the thread for solver
	virtual void solve_nodes_thread(int start, int end, const vector<Flag>& start_flags, const dynamic_bitset<>& start_state,
		const dynamic_bitset<>& end_state, const BooleanPolynomial & cur_coef,  const ListsOfPolynomialsAsFactors & cur_lists, double time, int threads, map<node_pair, BooleanPolynomial> & new_P, map<node_pair, ListsOfPolynomialsAsFactors> & new_P_asLists, map<BooleanMonomial, int>& sup_counter)
	{
		
		// generate local rounds flags for solver
		vector<Flag> middle_start_flags(start_flags);
		for (int i = 0; i < target_cipher.statesize; i++)
			if (start_flags[i] == "delta" && start_state[i] == 0)
				middle_start_flags[i] = "zero_c";

		vector<vector<vector<Flag>>> middle_rounds_flags;
		target_cipher.calculate_flags(start, end, middle_start_flags, middle_rounds_flags, true);

		// determine the midround
		int midround = start;
		while (midround < end)
		{
			bool is_full_delta = true;
			auto& this_round_flags = middle_rounds_flags[midround][0];
			for(auto & update_bit : target_cipher.update_bits)
				if (this_round_flags[update_bit] != "delta")
				{
					is_full_delta = false;
					break;
				}

			if (is_full_delta)
				break;
			else
				midround++;
		}

			

		vector<dynamic_bitset<>> p_sols;
		vector<BooleanPolynomial> p_exps;
		vector<pair<int, int>> p_rms;
		set<int> end_constants;
		if (midround >= end)
			// midround = start + (end - start) * 3 / 4;
			midround = end - 1;

		
		

		auto solver_status = target_cipher.two_stage_solve_model(start, end, middle_rounds_flags, start_state, end_state, p_sols, p_exps, p_rms, end_constants, threads, time, true, midround);

	


		if (solver_status == status::SOLVED)
		{
			
			if (solver_mode == mode::NO_OUTPUT)
			{
				; // do nothing
			}
			else if (solver_mode == mode::OUTPUT_FILE || solver_mode == mode::OUTPUT_EXP)
			{
				// print_sol(start, end, start_state, end_state, p_sols, p_exps, p_rms, end_constants);
				bool has_sols = (p_sols.size() > 0);
				if (!has_sols)
					return;



				bool has_constant = (end_constants.size() > 0);
				int p_num = p_sols[0].size();
				dynamic_bitset<> p_sol_mask(p_num);
				for (auto& p_sol : p_sols)
					p_sol_mask |= p_sol;

				// find the max number of rounds to generate middle exps for p_sols 
				int maxr = 0;
				if (has_constant)
					maxr = end - start;
				else
				{
					for (int i = p_num - 1; i >= 0; i--)
					{
						if (p_sol_mask[i] == 1)
						{
							maxr = p_rms[i].first;
							break;
						}
					}
				}

				auto middle_rounds_exps = target_cipher.generate_exps(start, start + maxr + 1, middle_start_flags, normal_exp);
				auto first_coeff_lists = get_cof_lists(start, end, p_sols, p_exps, p_rms, start, middle_rounds_exps);

				if (has_constant)
				{
					ListOfPolynomialsAsFactors constant_list(target_cipher.default_constant_1);
					for (auto& i : end_constants)
						constant_list.add(middle_rounds_exps[maxr][0][i]);


					first_coeff_lists *= constant_list;
				}

				ListsOfPolynomialsAsFactors second_coeff_lists(target_cipher.statesize);
				for (auto& coeff_list : first_coeff_lists)
				{
					ListOfPolynomialsAsFactors new_coeff_list(target_cipher.statesize);
					for (auto& poly : coeff_list)
					{
						auto newPoly = poly.subs(all_rounds_exps[start][0]);
						new_coeff_list.add(newPoly);
					}
					second_coeff_lists.add(new_coeff_list);
				}

				if (solver_mode == mode::OUTPUT_EXP)
				{
					/*
					auto contr2 = cur_lists * second_coeff_lists;
					contr2.filterLists();
					print_sol(contr2);
					*/

					auto contr = cur_coef * second_coeff_lists.getSum();
					lock_guard<mutex> guard(solver_mutex0);
					for (auto& mon : contr)
						sup_counter[mon] ++;
				}
				else
				{
					// print the solution
					auto contr = cur_lists * second_coeff_lists;
					contr.filterLists();
					print_sol(start, end, contr);
					print_sol_debug(start, end, start_state, end_state, cur_lists, second_coeff_lists, contr);
				}

			}
			
		}
		else if (solver_status == status::UNSOLVED)
		{


			lock_guard<mutex> guard(solver_mutex1);
			new_P[pair(start_state, end_state)] = cur_coef;
			new_P_asLists[pair(start_state, end_state)] = cur_lists;
		}
		
	}

	virtual void solve_nodes(int start, int end, const map<node_pair, BooleanPolynomial>& cur_P, 
		const map<node_pair, ListsOfPolynomialsAsFactors> & cur_P_asLists, double time, int threads, ThreadPool& thread_pool, map<node_pair, BooleanPolynomial>& new_P, map<node_pair, ListsOfPolynomialsAsFactors> & new_P_asLists, map<BooleanMonomial, int>& sup_counter)
	{
		
		auto start_flags = all_rounds_flags[start][0];

		vector<future<void>> futures;



		for (auto& nodes_coef : cur_P)
		{
			auto& nodepair = nodes_coef.first;
			auto& start_state = nodepair.first;
			auto& end_state = nodepair.second;
			auto& cur_coef = nodes_coef.second;
			auto& cur_lists = cur_P_asLists.at(nodepair);
			futures.emplace_back(thread_pool.Submit(&framework::solve_nodes_thread, REF(*this), start, end, REF(start_flags), REF(start_state),
				REF(end_state), REF(cur_coef), REF(cur_lists), time, threads, REF(new_P), REF(new_P_asLists), REF(sup_counter) ));
		}

		for (auto& it : futures)
			it.get();
		

	}

	void read_P(const string filepath, const FileReader& reader, map<node_pair, BooleanPolynomial> & new_P,
		map<node_pair,ListsOfPolynomialsAsFactors> & new_P_asLists)
	{
		fstream fs;
		fs.open(filepath);
		while (1)
		{
			node_pair state_pair;
			BooleanPolynomial coef;
			ListsOfPolynomialsAsFactors coef_lists;
			auto not_EOF = reader.read_one_pair(fs, state_pair, coef, coef_lists);
			if (not_EOF)
			{
				new_P[state_pair] = coef;
				new_P_asLists[state_pair] = coef_lists;
			}
			else
			{
				break;
			}
		}
		fs.close();
	}



	/**
	 * @brief This is the function for a single thread to read solutions from the file.
	 * @param filepath The path of the file
	 * @param reader The class used to read solutions from the file
	 * @param sup_counter The data structure used to count the occuring times of each monomial
	 * @param line_counter The data structure used to count the occuring times of each line
	 * @param isAccurate The parameter used to indicate whether we want to recover the concrete expression of the superpoly. For a massive superpoly, we tend to set isAccurate = False. 
	*/
	void read_sols_thread(const string filepath,  const FileReader & reader, map<BooleanMonomial, int> & sup_counter, map<string, int > & line_counter, bool isAccurate)
	{
		cout << "Start to read " << filepath << endl;


		fstream fs;
		fs.open(filepath);

		map<BooleanMonomial, int> thread_sup_counter;
		map<string, int> thread_line_counter;


		if (isAccurate)
		{
			while (1)
			{
				ListsOfPolynomialsAsFactors coef_lists;

				auto not_EOF = reader.read_lists_once(fs, coef_lists, target_cipher.statesize);

				if (!not_EOF)
					break;

				// coef_lists.filterLists();
				BooleanPolynomial coef = coef_lists.getSum();
				
					


				for (auto& mon : coef)
				{
					thread_sup_counter[mon]++;

				}


			}
		
		}
		else
		{
			string oneline;
			while (getline(fs, oneline))
			{
				if(!oneline.empty())
					thread_line_counter[oneline]++;
			}
		}






		fs.close();

		cout << "Read sol file complete : " <<filepath<<endl;

		lock_guard<mutex> guard(reader_mutex);
		if (isAccurate)
		{
			for (auto& mon_cnt : thread_sup_counter)
				if (mon_cnt.second % 2)
					sup_counter[mon_cnt.first]++;
		}
		else
		{
			for (auto& line_cnt : thread_line_counter)
				if (line_cnt.second % 2)
					line_counter[line_cnt.first]++;
					
		}


	}



	



	void read_sols(string TERM_path,  const FileReader& reader, map<BooleanMonomial, int> & sup_counter, map<string,int> & line_counter, bool isAccurate)
	{
		// first sort the files 
		using round_pair = pair<int, int>;

		vector<filesystem::path> sol_files;
		reader.getJustCurrentFilePaths(TERM_path, sol_files);

		map<round_pair, vector<string>> sort_sol_files;

		for (auto& sol_file : sol_files)
		{
			auto sol_filename = sol_file.filename().string();
			cout << "Detect sol_file: " << sol_filename << endl;
			smatch sol_filename_sm;
			auto match_status = regex_match(sol_filename, sol_filename_sm, sol_filename_regex);
			if (match_status)
			{
				int start = stoi(sol_filename_sm[1]);
				int end = stoi(sol_filename_sm[2]);
				sort_sol_files[pair(start, end)].emplace_back(sol_file.string() );
			}
		}

		vector<future<void>> futures;
		for (auto& rr_filepath : sort_sol_files)
		{
			auto& rr = rr_filepath.first;
			auto& sol_filepaths = rr_filepath.second;


			for (auto& sol_file : sol_filepaths)
			{
				futures.emplace_back(threadpool.Submit(&framework::read_sols_thread, this, sol_file, REF(reader),REF(sup_counter), REF(line_counter), isAccurate));
			}


		}


		for (auto& it : futures)
			it.get();

		#ifdef _WIN32
		#else
				showProcessMemUsage();
				malloc_trim(0);
				showProcessMemUsage();
		#endif
	}

	void analyze_superpoly(string TERM_path = R"~(./TERM)~")
	{

		// we only use single thread
		if (solver_mode == mode::OUTPUT_FILE || solver_mode == mode::OUTPUT_EXP)
		{
			map<BooleanMonomial, int> sup_counter;
			fstream fs;
			string path = TERM_path + string(R"(/superpoly.txt)");
			fs.open(path);
			string oneline;
			while (getline(fs, oneline))
			{
				if (!oneline.empty())
				{
					BooleanPolynomial poly(target_cipher.statesize, oneline);
					for (auto& mon : poly)
						sup_counter[mon]++;
				}
			}

			fs.close();
			

			if (true)
			{
				long long term_count = 0;
				int d = 0;
				map<int, int> key_bits_counter;
				vector<BooleanMonomial> mons_in_superpoly;

				if (target_cipher.ciphername != "acorn")
				{
					for (auto& mon_cnt : sup_counter)
						if (mon_cnt.second % 2)
						{
							term_count++;
							int mon_d = mon_cnt.first.count();
							if (mon_d > d)
								d = mon_d;

							for (auto& i : mon_cnt.first.index())
								key_bits_counter[i] ++;

							mons_in_superpoly.emplace_back(mon_cnt.first);
						};

					vector<int> involved_key_bits;
					for (auto& bit_cnt : key_bits_counter)
						involved_key_bits.emplace_back(bit_cnt.first);

					double balancedness = calculate_balancedness(involved_key_bits, mons_in_superpoly);




					cout << "The number of monomials appearing in the superpoly: " << term_count << endl;
					cout << "The algebraic degree of the superpoly: " << d << endl;
					cout << "The superpoly involves " << involved_key_bits.size() << " key bits: " << endl;
					for (auto& bit : involved_key_bits)
						cout << bit << " ";
					cout << endl;
					cout << "The balancedness is estimated to be " << balancedness << endl;
				}
				else
				{
					// first output superpoly of round 256
					for (auto& mon_cnt : sup_counter)
						if (mon_cnt.second % 2)
						{
							term_count++;
							int mon_d = mon_cnt.first.count();
							if (mon_d > d)
								d = mon_d;
						};


					cout << "The number of monomials appearing in the superpoly of round 256: " << term_count << endl;
					cout << "The algebraic degree of the superpoly of round 256: " << d << endl;

					// next transform the superpoly of round 256 into the final superpoly
					cipher_acorn& acorn = dynamic_cast<cipher_acorn&>(target_cipher);
					auto& round256_exps = acorn.round256exps;
					int round256_varsnum = round256_exps.size();

					map<BooleanMonomial, int> final_sup_counter;
					for (auto& mon_cnt : sup_counter)
						if (mon_cnt.second % 2)
						{
							set<BooleanPolynomial> expand_exps;
							for (auto& i : mon_cnt.first.index())
							{
								expand_exps.emplace(round256_exps[i]);
							}

							BooleanPolynomial expand_res = fastmul(round256_varsnum, expand_exps);
							for (auto& mon : expand_res)
								final_sup_counter[mon]++;
						}

					// output final superpoly to file
					term_count = 0;
					d = 0;
					for (auto& mon_cnt : final_sup_counter)
						if (mon_cnt.second % 2)
						{
							term_count++;
							int mon_d = mon_cnt.first.count();
							if (mon_d > d)
								d = mon_d;

							for (auto& i : mon_cnt.first.index())
								key_bits_counter[i] ++;

							mons_in_superpoly.emplace_back(mon_cnt.first);
						};


					vector<int> involved_key_bits;
					for (auto& bit_cnt : key_bits_counter)
						involved_key_bits.emplace_back(bit_cnt.first);

					double balancedness = calculate_balancedness(involved_key_bits, mons_in_superpoly);




					cout << "The number of monomials appearing in the superpoly: " << term_count << endl;
					cout << "The algebraic degree of the superpoly: " << d << endl;
					cout << "The superpoly involves " << involved_key_bits.size() << " key bits: " << endl;
					for (auto& bit : involved_key_bits)
						cout << bit << " ";
					cout << endl;
					cout << "The balancedness is estimated to be " << balancedness << endl;
				}

			}


		}
	}

	void read_sols_and_output(bool isAccurate = false, string TERM_path = R"~(./TERM)~")
	{
		if (solver_mode == mode::OUTPUT_FILE)
		{
			map<BooleanMonomial, int> sup_counter;
			map<string, int> line_counter;
			FileReader reader;

			read_sols(TERM_path, reader, sup_counter, line_counter, isAccurate);

			if (isAccurate)
			{
				long long term_count = 0;
				int d = 0;
				map<int, int> key_bits_counter;
				vector<BooleanMonomial> mons_in_superpoly;

				if (target_cipher.ciphername != "acorn")
				{
					cout << "Output superpoly to file." << endl;
					string path = TERM_path + string(R"(/superpoly.txt)");
					ofstream os;
					os.open(path, ios::out);
					for (auto& mon_cnt : sup_counter)
						if (mon_cnt.second % 2)
						{
							os << mon_cnt.first << endl;
							term_count++;
							int mon_d = mon_cnt.first.count();
							if (mon_d > d)
								d = mon_d;

							for (auto& i : mon_cnt.first.index())
								key_bits_counter[i] ++;

							mons_in_superpoly.emplace_back(mon_cnt.first);
						};
					os.close();
					cout << "Output superpoly to file finished." << endl;

					vector<int> involved_key_bits;
					for (auto& bit_cnt : key_bits_counter)
						involved_key_bits.emplace_back(bit_cnt.first);

					double balancedness = calculate_balancedness(involved_key_bits, mons_in_superpoly);




					cout << "The number of monomials appearing in the superpoly: " << term_count << endl;
					cout << "The algebraic degree of the superpoly: " << d << endl;
					cout << "The superpoly involves " << involved_key_bits.size() << " key bits: " << endl;
					for (auto& bit : involved_key_bits)
						cout << bit << " ";
					cout << endl;
					cout << "The balancedness is estimated to be " << balancedness << endl;
				}
				else
				{
					// first output superpoly of round 256
					cout << "Output superpoly of round 256 to file." << endl;
					string path = TERM_path + string(R"(/superpoly256.txt)");
					ofstream os;
					os.open(path, ios::out);
					for (auto& mon_cnt : sup_counter)
						if (mon_cnt.second % 2)
						{
							os << mon_cnt.first << endl;
							term_count++;
							int mon_d = mon_cnt.first.count();
							if (mon_d > d)
								d = mon_d;
						};
					os.close();
					cout << "Output superpoly of round 256 to file finished." << endl;


					cout << "The number of monomials appearing in the superpoly of round 256: " << term_count << endl;
					cout << "The algebraic degree of the superpoly of round 256: " << d << endl;

					// next transform the superpoly of round 256 into the final superpoly
					cipher_acorn& acorn = dynamic_cast<cipher_acorn&>(target_cipher);
					auto& round256_exps = acorn.round256exps;
					int round256_varsnum = round256_exps.size();

					map<BooleanMonomial, int> final_sup_counter;
					for (auto& mon_cnt : sup_counter)
						if (mon_cnt.second % 2)
						{
							set<BooleanPolynomial> expand_exps;
							for (auto& i : mon_cnt.first.index())
							{
								expand_exps.emplace(round256_exps[i]);
							}

							BooleanPolynomial expand_res = fastmul(round256_varsnum, expand_exps);
							for (auto& mon : expand_res)
								final_sup_counter[mon]++;
						}

					// output final superpoly to file
					term_count = 0;
					d = 0;
					cout << "Output superpoly to file." << endl;
					path = TERM_path + string(R"(/superpoly.txt)");
					os.open(path, ios::out);
					for (auto& mon_cnt : final_sup_counter)
						if (mon_cnt.second % 2)
						{
							os << mon_cnt.first << endl;
							term_count++;
							int mon_d = mon_cnt.first.count();
							if (mon_d > d)
								d = mon_d;

							for (auto& i : mon_cnt.first.index())
								key_bits_counter[i] ++;

							mons_in_superpoly.emplace_back(mon_cnt.first);
						};
					os.close();
					cout << "Output superpoly to file finished." << endl;


					vector<int> involved_key_bits;
					for (auto& bit_cnt : key_bits_counter)
						involved_key_bits.emplace_back(bit_cnt.first);

					double balancedness = calculate_balancedness(involved_key_bits, mons_in_superpoly);




					cout << "The number of monomials appearing in the superpoly: " << term_count << endl;
					cout << "The algebraic degree of the superpoly: " << d << endl;
					cout << "The superpoly involves " << involved_key_bits.size() << " key bits: " << endl;
					for (auto& bit : involved_key_bits)
						cout << bit << " ";
					cout << endl;
					cout << "The balancedness is estimated to be " << balancedness << endl;
				}

			}
			else
			{
				cout << "Output superpoly to file." << endl;
				string path = TERM_path + string(R"(/superpoly.txt)");
				ofstream os;
				os.open(path, ios::out);
				for(auto & line_cnt : line_counter)
					if (line_cnt.second % 2)
					{
						os << line_cnt.first << endl;
					}

				os.close();
				cout << "Output superpoly to file finished." << endl;
			}


		}
	}

	
	long long analyze_coef_list(const string& coef_list_str, map<string, set<int> > & vars_index_per_poly, set<int>& vars_index_this_list, int & d)
	{
		// if there is a poly like (0+0+0+0), this function can not identify it as 0.


		int i = 0;
		long long term_cnt = 1;

		string polystr;
		set<int> vars_index_this_poly;
		long long term_cnt_this_poly = 0;

		string indexstr;
		int poly_d = 0;
		int mon_d = 0;


		
		while (i < coef_list_str.size())
		{

			if (coef_list_str[i] == '(')
			{
				polystr = "";
				vars_index_this_poly.clear();
				poly_d = 0;
				term_cnt_this_poly = 0;
			}
			else if (coef_list_str[i] == ')')
			{
				term_cnt_this_poly++;
				term_cnt *= term_cnt_this_poly;

				if (!indexstr.empty())
				{
					vars_index_this_list.emplace(stoi(indexstr));
					vars_index_this_poly.emplace(stoi(indexstr));
					indexstr = "";
				}

				if (poly_d < mon_d)
				{
					poly_d = mon_d;
				}

				mon_d = 0;

				if (!polystr.empty())
				{
					vars_index_per_poly[polystr] = vars_index_this_poly;
					d += poly_d;
				}
			}
			else if (coef_list_str[i] == 's')
			{
				mon_d++;

				if (!indexstr.empty())
				{
					vars_index_this_list.emplace(stoi(indexstr));
					vars_index_this_poly.emplace(stoi(indexstr));
					indexstr = "";
				}

				polystr += coef_list_str[i];
			}
			else if (coef_list_str[i] == '+')
			{
				term_cnt_this_poly++;

				polystr += coef_list_str[i];
				if (poly_d < mon_d)
				{
					poly_d = mon_d;
				}

				mon_d = 0;

				if (!indexstr.empty())
				{
					vars_index_this_list.emplace(stoi(indexstr));
					vars_index_this_poly.emplace(stoi(indexstr));
					indexstr = "";
				}
			}
			else if (coef_list_str[i] == '1' && (coef_list_str[i - 1] == '(' || coef_list_str[i-1] == '+'))
			{
				polystr += coef_list_str[i];
			}
			else if (coef_list_str[i] == '0' )
			{
				if (coef_list_str[i - 1] == '+')
				{
					if (!polystr.empty() && polystr.back() == '+')
					{
						polystr.pop_back();
					}
					else if (polystr.empty() && coef_list_str[i+1] == '+')
						i++;
						
				}
				else if (coef_list_str[i - 1] == '(')
				{
					if (coef_list_str[i + 1] == ')')
					{
						d = 0;
						vars_index_per_poly.clear();
						vars_index_this_list.clear();
						return 0;
					}
					i++;
				}
				else
				{
					indexstr += coef_list_str[i];
					polystr += coef_list_str[i];
				}
			}
			else
			{
				indexstr += coef_list_str[i];
				polystr += coef_list_str[i];
			}

			i++;
		}

		return term_cnt;
		
	}


	void analyze_superpoly_asLists_thread(const vector<string>& lines, int from, int to, map<int, int>& lists_per_var, map<string, set<int> >& vars_index_per_poly, int& sup_deg, long long& term_cnt, vector<vector<string>> & polys_per_list)
	{

		map<int, int> thread_lists_per_var;
		map<string,set<int>> thread_vars_index_per_poly;
		vector<vector<string>> thread_polys_per_list;

		int thread_deg = 0;
		long long thread_term_cnt = 0;

		int i = from;
		while(i < to)
		{

			// 
			set<int> vars_index_this_list;
			map<string, set<int>> vars_index_per_poly_this_list;
			int deg_this_list = 0;
			thread_term_cnt += analyze_coef_list(lines[i], vars_index_per_poly_this_list, vars_index_this_list, deg_this_list);
			if (thread_deg < deg_this_list)
				thread_deg = deg_this_list;
			
			for (auto& var_index : vars_index_this_list)
				thread_lists_per_var[var_index]++;

			vector<string> polys_this_list;
			for (auto& poly_index : vars_index_per_poly_this_list)
			{
				thread_vars_index_per_poly.emplace(poly_index);
				polys_this_list.emplace_back(poly_index.first);
			}

			if(polys_this_list.size() > 0)
				thread_polys_per_list.emplace_back(polys_this_list);


			i++;
		}

		lock_guard<mutex> guard(reader_mutex);
		if (sup_deg < thread_deg)
		{
			sup_deg = thread_deg;
			cout << "Update the upper bound of the degree to " << sup_deg << endl;
		}

		for (auto& poly_index : thread_vars_index_per_poly)
			if (vars_index_per_poly.find(poly_index.first) == vars_index_per_poly.end())
				vars_index_per_poly.emplace(poly_index);
		cout << "Current collected polys : " << vars_index_per_poly.size() << endl;

		polys_per_list.insert(polys_per_list.end(), thread_polys_per_list.begin(), thread_polys_per_list.end());

		for (auto& var_cnt : thread_lists_per_var)
			lists_per_var[var_cnt.first] += var_cnt.second;

		term_cnt += thread_term_cnt;

		cout << "----------------------------------- thread over ----------------------------------------" << endl;
	}

	void calculate_superpoly_balancedness_thread(int thread_nr_tests, const vector<int>& involved_key_bits, const vector<BooleanMonomial>& mons_in_superpoly, int& res1_cnt)
	{

		if (mons_in_superpoly.size() == 0)
			return;

		int varsnum = mons_in_superpoly[0].size();

		random_device rd;
		mt19937 re(rd());
		bernoulli_distribution d(0.5);
		
		int thread_res1_cnt = 0;

		for (int i = 0; i < thread_nr_tests; i++)
		{
			// generate random values for the secret variables
			dynamic_bitset<> vals(varsnum);
			for (auto& bit : involved_key_bits)
				vals[bit] = d(re);

			int res = 0;
			for (auto& mon : mons_in_superpoly)
				res ^= mon.eval(vals);

			if (res == 1)
				thread_res1_cnt++;
		}

		cout << "calculate thread finished - " << thread_res1_cnt << endl;

		lock_guard<mutex> guard(reader_mutex);
		res1_cnt += thread_res1_cnt;
	}


	double calculate_balancedness(const vector<int>& involved_key_bits, const vector<BooleanMonomial>& mons_in_superpoly)
	{
		int res1_cnt = 0;
		int total_nr_tests = 1 << 15;
		int test_cnt = 0;
		int nr_tests_each_thread = 100;

		vector<future<void>> futures;
		while (test_cnt < total_nr_tests)
		{
			int thread_nr_tests = (test_cnt + nr_tests_each_thread < total_nr_tests) ? nr_tests_each_thread : total_nr_tests - test_cnt;
			futures.emplace_back(threadpool.Submit(&framework::calculate_superpoly_balancedness_thread, REF(*this), thread_nr_tests, REF(involved_key_bits), REF(mons_in_superpoly), REF(res1_cnt)));

			test_cnt += nr_tests_each_thread;
		}

		for (auto& it : futures)
			it.get();

		return (double)res1_cnt / (double)total_nr_tests;
	}

	void calculate_lists_balancedness_thread(int thread_nr_tests, const vector<int>& involved_key_bits, 
		const map<string, BooleanPolynomial>& polystrs_to_polys, 
		const vector<vector<string>>& polystrs_per_list, int& res1_cnt)
	{
		random_device rd;
		mt19937 re(rd());
		bernoulli_distribution d(0.5);

		int thread_res1_cnt = 0;

		for (int i = 0; i < thread_nr_tests; i++)
		{
			// generate random values for the secret variables
			dynamic_bitset<> vals(target_cipher.statesize);
			for (auto& bit : involved_key_bits)
				vals[bit] = d(re);

			/*
			while (polystrs_to_polys.at("s44+s42s43+s17").eval(vals) != 0 )
			{
				for (auto& bit : involved_key_bits)
					vals[bit] = d(re);
			}
			*/
			

			// generate the values for each poly under the values of secret variables
			map<string, int> polystr_vals;
			for (auto& polystr_poly : polystrs_to_polys)
			{
				polystr_vals[polystr_poly.first] = polystr_poly.second.eval(vals);
			}

			int res = 0;
			for (auto& polystrs_this_list : polystrs_per_list)
			{
				int res_this_list = 1;
				for (auto& polystr : polystrs_this_list)
					if (polystr_vals[polystr] == 0)
					{
						res_this_list = 0;
						break;
					}

				res ^= res_this_list;

			}

			if (res == 1)
				thread_res1_cnt++;
		}

		cout << "calculate thread finished - " << thread_res1_cnt << endl;

		lock_guard<mutex> guard(reader_mutex);
		res1_cnt += thread_res1_cnt;
	}

	void calculate_reducecost_thread(const vector<vector<string>>& polystrs_per_list, int from, int to, int guessed_num, 
		const map<string, vector<set<int>>> &  polystrs_to_coeff_vars, const map<string, vector<unsigned long long>> & polystrs_to_coeff_var_vals_num, const map<string, vector<int>>& polystrs_to_coeff_var_vals_num_exp, map<string, vector<dynamic_bitset<>>>& polystrs_to_unguessed_monreps, map<int, unsigned long long> & exp_counter, 
		int& max_deg_all_lists)
	{
		map<int, unsigned long long> thread_exp_counter;
		int thread_max_deg_all_lists = 0;
		
		for (int i = from; i < to; i++)
		{
			auto& polystrs_this_list = polystrs_per_list[i];
			int polystrs_num = polystrs_this_list.size();
			if (polystrs_num == 0)
				continue;

			vector<vector<set<int>>::const_iterator> coeff_vars_iters(polystrs_num);
			vector<vector<unsigned long long>::const_iterator> coeff_var_vals_num_iters(polystrs_num);
			vector<vector<int>::const_iterator> coeff_var_vals_num_exp_iters(polystrs_num);
			vector<vector<dynamic_bitset<>>::const_iterator> unguessed_iters(polystrs_num);

			for (int j = 0; j < polystrs_num; j++)
			{
				auto& polystr = polystrs_this_list[j];
				coeff_vars_iters[j] = polystrs_to_coeff_vars.at(polystr).cbegin();
				coeff_var_vals_num_iters[j] = polystrs_to_coeff_var_vals_num.at(polystr).cbegin();
				coeff_var_vals_num_exp_iters[j] = polystrs_to_coeff_var_vals_num_exp.at(polystr).cbegin();
				unguessed_iters[j] = polystrs_to_unguessed_monreps.at(polystr).cbegin();
			}


			while (1)
			{
				unsigned long long prod_of_coeff_var_vals_num = 1;
				int prod_exp_of_coeff_var_vals_num = 0;
				set<int> union_of_coeff_vars;
				int cur_coeff_vars_num = 0;
				dynamic_bitset<> prod_of_monreps(target_cipher.statesize);

				for (int j = 0; j < polystrs_num; j++)
				{
					auto& cur_coeff_vars = *(coeff_vars_iters[j]);
					auto& cur_coeff_var_vals_num = *(coeff_var_vals_num_iters[j]);
					auto& cur_coeff_var_vals_num_exp = *(coeff_var_vals_num_exp_iters[j]);

					prod_of_monreps |= *(unguessed_iters[j]);

					union_of_coeff_vars.insert(cur_coeff_vars.cbegin(), cur_coeff_vars.cend());


					if (cur_coeff_var_vals_num_exp != -1)
						prod_exp_of_coeff_var_vals_num += cur_coeff_var_vals_num_exp;
					else
						prod_of_coeff_var_vals_num *= cur_coeff_var_vals_num;
					


					cur_coeff_vars_num = union_of_coeff_vars.size();
					unsigned long long max_limit = (unsigned long long)1 << (cur_coeff_vars_num - prod_exp_of_coeff_var_vals_num);
					if (prod_of_coeff_var_vals_num > max_limit)
						prod_of_coeff_var_vals_num = max_limit;

					if (prod_of_coeff_var_vals_num == 0)
					{
						cout << "This should never happen." << endl;
						break;
					}
				}

				int deg_this_prod = prod_of_monreps.count();
				if (thread_max_deg_all_lists < deg_this_prod)
					thread_max_deg_all_lists = deg_this_prod;



				thread_exp_counter[guessed_num - cur_coeff_vars_num + prod_exp_of_coeff_var_vals_num] += (prod_of_coeff_var_vals_num * ((unsigned long long)polystrs_num-1));


				// increment the iterators
				int j = 0;
				for (j = 0; j < polystrs_num; j++)
				{
					auto& polystr = polystrs_this_list[j];
					++coeff_vars_iters[j];
					if (coeff_vars_iters[j] == polystrs_to_coeff_vars.at(polystr).cend())
					{
						coeff_vars_iters[j] = polystrs_to_coeff_vars.at(polystr).cbegin();
						coeff_var_vals_num_iters[j] = polystrs_to_coeff_var_vals_num.at(polystr).cbegin();
						coeff_var_vals_num_exp_iters[j] = polystrs_to_coeff_var_vals_num_exp.at(polystr).cbegin();
						unguessed_iters[j] = polystrs_to_unguessed_monreps.at(polystr).cbegin();
					}
					else
					{
						++coeff_var_vals_num_iters[j];
						++coeff_var_vals_num_exp_iters[j];
						++unguessed_iters[j];
						break;
					}
				}

				if (j == polystrs_num)
					break;
			}
		}

		lock_guard<mutex> guard(analyze_mutex);
		for (auto& exp_cnt : thread_exp_counter)
		{
			// cout << "Update " << exp_cnt.first << " from " << hex << showbase << setw(12)<< setfill('0')<< exp_counter[exp_cnt.first];
			// cout << " plus " << exp_cnt.second;
			exp_counter[exp_cnt.first] += exp_cnt.second;
			// cout << " to " << hex << showbase << setw(12) << setfill('0') << exp_counter[exp_cnt.first] << endl;
		}

		if (max_deg_all_lists < thread_max_deg_all_lists)
		{
			max_deg_all_lists = thread_max_deg_all_lists;
			// cout << "Update maximum degree after reduction to " << max_deg_all_lists << endl;
		}
		

		cout << "Thread over " << endl;
	}

	vector<BooleanPolynomial> generate_submap_from_excluded_polys(const set<BooleanPolynomial>& excluded_polys) 
	{
		set<BooleanPolynomial> constrained_expanded_p_exps(excluded_polys);

		vector<BooleanPolynomial> cur_constrains_map = normal_exp;

		for (auto& expanded_p_exp : constrained_expanded_p_exps)
		{
			auto expanded_p_exp_after_sub = expanded_p_exp.subs(cur_constrains_map);

			map<int, vector<BooleanMonomial>> var_mons_distribution;
			for (auto& mon : expanded_p_exp_after_sub)
			{
				for (auto& v : mon.index())
					var_mons_distribution[v].emplace_back(mon);
			}

			bool is_balanced = false;
			int balanced_v = 0;
			BooleanMonomial balanced_mon;
			for (auto& v_mons : var_mons_distribution)
			{
					if (v_mons.second.size() == 1 && v_mons.second[0].count() == 1)
					{
						is_balanced = true;
						balanced_v = v_mons.first;
						balanced_mon = v_mons.second[0];
						break;
					}
			}

			if (is_balanced)
			{
				cur_constrains_map[balanced_v] = expanded_p_exp_after_sub + balanced_mon;
			}
		}

		return cur_constrains_map;
	}


	unsigned long long calculate_reducecost(const vector<int>& guessed_key_bits, const vector<vector<string>>& polystrs_per_list)
	{
		dynamic_bitset<> guessed_mask(target_cipher.statesize);
		for (auto& bit : guessed_key_bits)
			guessed_mask[bit] = 1;
		dynamic_bitset<> unguessed_mask = ~guessed_mask;

		// we regard the guessed part of a monomial as the coefficient of the unguessed part
		map<string, vector<BooleanPolynomial>> polystrs_to_coeffs;
		map<string, vector<dynamic_bitset<>>> polystrs_to_unguessed_monreps;

		for (auto& polystrs_this_list : polystrs_per_list)
		{
			for (auto& polystr : polystrs_this_list)
			{
				if (polystrs_to_coeffs.find(polystr) == polystrs_to_coeffs.end())
				{
					map<BooleanMonomial, BooleanPolynomial> coeff_of_unguessed;
					BooleanPolynomial poly(target_cipher.statesize, polystr);
					// cout << "Current poly " << poly << endl;
					for (auto& mon : poly)
					{
						// cout << "Current monomial " << mon << endl;
						BooleanMonomial guessed_part = mon.intersect(guessed_mask);
						// cout << "guessed part " << guessed_part << endl;
						BooleanMonomial unguessed_part = mon.intersect(unguessed_mask);
						// cout << "unguessed part " << unguessed_part << endl;
						if (coeff_of_unguessed.find(unguessed_part) != coeff_of_unguessed.end())
							coeff_of_unguessed[unguessed_part] += guessed_part;
						else
						{
							coeff_of_unguessed[unguessed_part] = guessed_part;
						}
					}

					vector<BooleanPolynomial> coeffs_this_poly;
					vector<dynamic_bitset<>> unguessed_monreps_this_poly;
					for (auto& unguessed_coeff : coeff_of_unguessed)
					{
						coeffs_this_poly.emplace_back(unguessed_coeff.second);
						unguessed_monreps_this_poly.emplace_back(unguessed_coeff.first.rep());
					}

					polystrs_to_coeffs[polystr] = coeffs_this_poly;
					polystrs_to_unguessed_monreps[polystr] = unguessed_monreps_this_poly;
				}

				
			}
		}

		for (auto& poly_unguessed_monreps : polystrs_to_unguessed_monreps)
		{
			auto& poly = poly_unguessed_monreps.first;
			auto& monreps = poly_unguessed_monreps.second;
			// cout << poly << ": ";
			// for (auto& monrep : monreps)
				// cout << monrep.count() << " ";

			// cout << endl;
		}

		// cout << "After reduction, the max degree is " << max_deg_all_lists << endl;

		// record the set of variables involved in each coefficient and check if coefficient is balanced
		map<string, vector<set<int>>> polystrs_to_coeff_vars;
		map<string, vector<bool>> polystrs_to_coeff_is_balanced;
		for (auto& polystr_coeffs : polystrs_to_coeffs)
		{
			auto& polystr = polystr_coeffs.first;
			auto& coeffs = polystr_coeffs.second;
			vector<set<int>> vars_per_coeff;
			vector<bool> is_balanced_per_coeff; // If a variable appears as a single monomial in the coeffcient, then it is balanced; otherwise it is not
			for (auto& coeff : coeffs)
			{
				set<int> vars_this_coeff;
				map<int, vector<BooleanMonomial>> monomials_per_var; // monomials related to each variable
				for (auto& mon : coeff)
				{
					set<int> vars_this_mon = mon.index();
					for (auto& v : vars_this_mon)
					{
						vars_this_coeff.emplace(v);
						monomials_per_var[v].emplace_back(mon);
					}
				}

				vars_per_coeff.emplace_back(vars_this_coeff);

				bool is_balanced = false;
				for (auto& var_mons : monomials_per_var)
				{
					auto& mons_this_var = var_mons.second;
					if (mons_this_var.size() == 1 && mons_this_var[0].count() == 1)
					{
						is_balanced = true;
						break;
					}
					
				}

				is_balanced_per_coeff.emplace_back(is_balanced);
			}

			polystrs_to_coeff_vars[polystr] = vars_per_coeff;
			polystrs_to_coeff_is_balanced[polystr] = is_balanced_per_coeff;
		}

		// count the number of possible values on which the coefficient evaluates to 1 for each coefficient
		map<string, vector<unsigned long long>> polystrs_to_coeff_var_vals_num;
		for (auto& polystr_coeff_vars : polystrs_to_coeff_vars)
		{
			auto& polystr = polystr_coeff_vars.first;
			auto& coeffs = polystrs_to_coeffs[polystr];
			auto& vars_per_coeff = polystr_coeff_vars.second;
			auto& is_balanced_per_coeff = polystrs_to_coeff_is_balanced[polystr];

			int coeff_num = coeffs.size();

			// cout << "Current poly : " << polystr << endl;

			vector<unsigned long long> var_vals_num_per_coeff;
			for (int i = 0; i < coeff_num; i++)
			{
				auto& vars_this_coeff = vars_per_coeff[i];
				auto& cur_coeff = coeffs[i];
				// cout << "Current coeff : " << cur_coeff << endl;
				int vars_num_this_coeff = vars_this_coeff.size();
				// cout << "Current number of variables in this coeff : " << vars_num_this_coeff << endl;
				bool is_balanced_this_coeff = is_balanced_per_coeff[i];
				// cout << "Is balanced : " << ((is_balanced_this_coeff) ? ("Yes") : ("No")) << endl;
				unsigned long long var_vals_num_this_coeff = 0;
				if (is_balanced_this_coeff)
				{
					var_vals_num_this_coeff = (unsigned long long)1 << (vars_num_this_coeff - 1);
				}
				else if (vars_num_this_coeff <= 25)
				{
					
					for (unsigned long long vars_val = 0; vars_val < (unsigned long long)1 << vars_num_this_coeff; vars_val++)
					{
						dynamic_bitset<> vars_val_mask(target_cipher.statesize);
						auto iter = vars_this_coeff.begin();
						int cur_bit = 0;
						while (iter != vars_this_coeff.end())
						{
							vars_val_mask[*iter] = (vars_val >> cur_bit) & (unsigned long long)0x1;
							cur_bit++;
							iter++;
						}

						if (cur_coeff.eval(vars_val_mask) == 1)
							var_vals_num_this_coeff++;
					}
				}
				else
					throw Exception("The coefficient is not balanced and has too many variables.");

				// cout << "The number of values of variables that make this coefficient equal to 1 : " << var_vals_num_this_coeff << endl;
				var_vals_num_per_coeff.emplace_back(var_vals_num_this_coeff);

			}

			polystrs_to_coeff_var_vals_num[polystr] = var_vals_num_per_coeff;
		}

		// record the exponent of the number of possible values on which the coefficient evaluates to 1 for each coefficient
		// if the exponent is not an integer, we record it as -1
		map<string, vector<int>> polystrs_to_coeff_var_vals_num_exp;
		for (auto& polystr_coeff_var_vals_num : polystrs_to_coeff_var_vals_num)
		{
			auto& polystr = polystr_coeff_var_vals_num.first;
			auto& coeff_var_vals_num = polystr_coeff_var_vals_num.second;

			vector<int> coeff_var_vals_num_exp;
			for (auto& n : coeff_var_vals_num)
			{
				vector<int> indexes;
				for (int i = 0; i < 64; i++)
					if (((n >> i) & (unsigned long long)0x1) == 1)
						indexes.emplace_back(i);

				if (indexes.size() == 1)
					coeff_var_vals_num_exp.emplace_back(indexes[0]);
				else
					coeff_var_vals_num_exp.emplace_back(-1);
			}

			polystrs_to_coeff_var_vals_num_exp[polystr] = coeff_var_vals_num_exp;
		}

		vector<future<void>> futures;

		int start = 0;
		int prods_per_thread = 100000;
		int prods_num = polystrs_per_list.size();
		map<int, unsigned long long> exp_counter;
		int guessed_num = guessed_key_bits.size();
		int max_deg_all_lists = 0;
		while (start < prods_num)
		{
			int end = (start + prods_per_thread) < prods_num ? start + prods_per_thread : prods_num;
			futures.emplace_back(threadpool.Submit(&framework::calculate_reducecost_thread, REF(*this), REF(polystrs_per_list), start, end, guessed_num, REF(polystrs_to_coeff_vars), REF(polystrs_to_coeff_var_vals_num), REF(polystrs_to_coeff_var_vals_num_exp), REF(polystrs_to_unguessed_monreps), REF(exp_counter), REF(max_deg_all_lists)));


			start += prods_per_thread;
		}

		for (auto& it : futures)
			it.get();

		cout << "The maximum degree after reduction is " << max_deg_all_lists << endl;


		unsigned long long cost = 0;
		for (auto& exp_cnt : exp_counter)
		{
			auto& exp = exp_cnt.first;
			auto& cnt = exp_cnt.second;
			cost += (cnt >> (guessed_num - exp) )+(unsigned long long)(1);
		}


		return cost;

		
		

		
	}
	



	double calculate_balancedness(const vector<int>& involved_key_bits, const vector<vector<string>>& polystrs_per_list)
	{
		map<string, BooleanPolynomial> polystrs_to_polys;
		for (auto& polystrs_this_list : polystrs_per_list)
			for (auto& polystr : polystrs_this_list)
				if (polystrs_to_polys.find(polystr) == polystrs_to_polys.end())
				{
					polystrs_to_polys[polystr] = BooleanPolynomial(target_cipher.statesize, polystr);
				}



		int res1_cnt = 0;
		int total_nr_tests = 1 << 15;
		int test_cnt = 0;
		int nr_tests_each_thread = 100;

		vector<future<void>> futures;
		while (test_cnt < total_nr_tests)
		{
			int thread_nr_tests = (test_cnt + nr_tests_each_thread < total_nr_tests) ? nr_tests_each_thread : total_nr_tests - test_cnt;
			futures.emplace_back(threadpool.Submit(&framework::calculate_lists_balancedness_thread, REF(*this), thread_nr_tests, REF(involved_key_bits), REF(polystrs_to_polys), REF(polystrs_per_list), REF(res1_cnt)));

			test_cnt += nr_tests_each_thread;
		}

		for (auto& it : futures)
			it.get();

		return (double)res1_cnt / (double)total_nr_tests;
	}

	
	
	

	


	void analyze_superpoly_asLists(vector<string> TERM_paths = { R"~(./TERM)~" })
	{
		if (solver_mode == mode::OUTPUT_FILE)
		{
			fstream fs;
			string oneline;
			vector<string> lines;
			for (auto& TERM_path : TERM_paths)
			{
				string path = TERM_path + string(R"(/superpoly.txt)");
				fs.open(path, ios::in);

				while (getline(fs, oneline))
				{	

					if (!oneline.empty())
					{
						lines.emplace_back(oneline);
					}

				}

				fs.close();
				cout << "After reading " << TERM_path <<", Number of products: " << lines.size()<<endl;
			}

			cout << "Read in total " << lines.size() << " products (including constant 0 product)." << endl;

			vector<future<void>> futures;
			map<int, int> lists_per_var;
			map<string, set<int> > vars_index_per_poly;
			vector<vector<string>> polys_per_list;
			int sup_deg = 0;
			long long term_cnt = 0;


			

			int start = 0;
			int lines_per_thread = 100000;
			while (start < lines.size())
			{
				int end = (start + lines_per_thread) < lines.size() ? start + lines_per_thread : lines.size();
				futures.emplace_back(threadpool.Submit(&framework::analyze_superpoly_asLists_thread, REF(*this), REF(lines), start, end, REF(lists_per_var), REF(vars_index_per_poly), REF(sup_deg), REF(term_cnt), REF(polys_per_list)));

				start += lines_per_thread;
			}

			for (auto& it : futures)
				it.get();

			cout << "After removing constant 0 products, there are " << polys_per_list.size() << " products." << endl;

			// This is to remove s58 from the superpoly
			map<BooleanPolynomial, BooleanPolynomial> polys_to_sub;
			vector<BooleanPolynomial> submap = normal_exp;
			submap[58] = BooleanPolynomial(target_cipher.statesize, "0");

			for (auto& it : vars_index_per_poly)
			{
				BooleanPolynomial poly(target_cipher.statesize, it.first);
				if (poly.iszero())
					continue;
				if (polys_to_sub.find(poly) == polys_to_sub.end())
				{
					polys_to_sub[poly] = poly.subs(submap);
				}

			}

			

			string path = R"~(./TERM)~" + string(R"(/superpoly2.txt)");
			fs.open(path, ios::out);

			cout << "Start to output " << endl;
			ListsOfPolynomialsAsFactors lists_after_sub;

			int list_num = 0;
			for (auto& polys_this_list : polys_per_list)
			{
				bool isThisListZero = false;
				ListOfPolynomialsAsFactors list_after_sub;
				for (auto& polystr : polys_this_list)
				{
					BooleanPolynomial poly(target_cipher.statesize, polystr);
					BooleanPolynomial poly_after_sub = poly.subs(submap);
					if (poly_after_sub.iszero())
					{
						isThisListZero = true;
						break;
					}

					list_after_sub.add(poly_after_sub);
				}

				if (isThisListZero)
					continue;

				lists_after_sub.add(list_after_sub);
				cout << list_num++ << endl;
			}
			lists_after_sub.filterLists();
			lists_after_sub.display(fs);
			fs.close();

			exit(-1);
			// end



			long long num_of_factors = 0;
			long long max_num_of_factors_per_list = 0;
			for (auto& polys_this_list : polys_per_list)
			{
				num_of_factors += polys_this_list.size();
				if (polys_this_list.size() > max_num_of_factors_per_list)
					max_num_of_factors_per_list = polys_this_list.size();
			}

			cout << "There are in total " << num_of_factors << " factors in all products." << endl;
			cout << "There are in total " << vars_index_per_poly.size() << " polys that appear as factors of products. " << endl;

			// test if the polys appearing as factors are truly included by the set of all possible factors
			vector<BooleanPolynomial> all_possible_expanded_p_exps;
			// get_all_possible_expanded_p_exps(all_possible_expanded_p_exps);

			cout << "The number of all possible factors is " << all_possible_expanded_p_exps.size() << endl;

			for (auto& poly : all_possible_expanded_p_exps)
			{
				if(poly.isone())
					cout << "The constant one is included in the set of all possible factors." << endl;
			}

			for (auto& poly_vars : vars_index_per_poly)
			{
				BooleanPolynomial poly(target_cipher.statesize, poly_vars.first);
				if (find(all_possible_expanded_p_exps.begin(), all_possible_expanded_p_exps.end(), poly) == all_possible_expanded_p_exps.end())
					cout << "The polynomial " << poly << " is not in the set of all possible factors." << endl;
			}

			//

			map<int, int> polys_per_var;
			for (auto& poly_indices : vars_index_per_poly)
			{
				// cout << poly_indices.first << ": " << endl;
				for (auto& index : poly_indices.second)
				{
					// cout << index << " ";
					polys_per_var[index]++;
				}
				// cout << endl << endl;
			}
			// cout << endl;
			//cout << endl;





			cout << "The number of products related to each variable. " << endl;
			for (auto& index_cnt : lists_per_var)
				cout << index_cnt.first << " " << index_cnt.second << endl;

			cout << endl;
			cout << endl;

			vector<int> involved_key_bits;
			for (auto& index_cnt : lists_per_var)
				involved_key_bits.emplace_back(index_cnt.first);

			cout << "The superpoly involves " << involved_key_bits.size() << " key bits: " << endl;
			for (auto& bit : involved_key_bits)
				cout << bit << " ";

			cout << endl;

			cout << endl;

			vector<pair<int, int>> sorted_lists_per_var(lists_per_var.begin(), lists_per_var.end());
			sort(sorted_lists_per_var.begin(), sorted_lists_per_var.end(), [&](pair<int, int> a, pair<int, int> b)
				{return a.second > b.second; }
			);




			vector<int> guessed_key_bits;
			// int guessnum = modelGuessNum(vars_index_per_poly, target_cipher.keysize - 40, guessed_key_bits);
			for (int i = 0; i < int(involved_key_bits.size()) - 35; i++)
				guessed_key_bits.emplace_back(sorted_lists_per_var[i].first);
			cout << "Reduce the superpoly by guessing " << guessed_key_bits.size()<<" bits ";
			set<int> sorted_guessed_key_bits(guessed_key_bits.begin(), guessed_key_bits.end());
			for (auto& i : sorted_guessed_key_bits)
				cout << i << ",";
			cout << endl;

			unsigned long long cost = calculate_reducecost(guessed_key_bits, polys_per_list);
			
			
			cout << "After reduction, there are " << cost << " monomials remaining in the polynomial without considering cancellation." << endl;

			cout << endl;
			cout << endl;

			// calculate the number of lists each poly appear
			map<string, int> lists_num_per_polystr;
			set<string> excluded_polystrs = {  };
			set<BooleanPolynomial> excluded_polys;

			vector<vector<string>> polystrs_per_list_after_excluded = polys_per_list;
			int sub_depth = 0;

			while (1)
			{
				map<string, BooleanPolynomial> polystrs_to_polys_after_sub;
				vector<BooleanPolynomial> submap_from_excluded_polys = generate_submap_from_excluded_polys(excluded_polys);
				for (auto& poly_vars : vars_index_per_poly)
				{
					BooleanPolynomial poly(target_cipher.statesize, poly_vars.first);
					for(int i = 0; i < sub_depth; i++)
						poly = poly.subs(submap_from_excluded_polys);
					polystrs_to_polys_after_sub[poly_vars.first] = poly;
				}

				lists_num_per_polystr.clear();
				vector<vector<string>> next_polystrs_per_list_after_excluded;

				for (auto& polystrs_this_list : polystrs_per_list_after_excluded)
				{

					bool isThisProdExcluded = false;

					for (auto& polystr : polystrs_this_list)
					{
						if (excluded_polystrs.find(polystr) != excluded_polystrs.end())
						{
							isThisProdExcluded = true;
							break;
						}

						if (polystrs_to_polys_after_sub[polystr].iszero())
						{
							isThisProdExcluded = true;
							cout << polystr << " is zero after substitution." << endl;
							cout << "all polystrs in this list: ";
							for (auto& polystr : polystrs_this_list)
								cout << polystr << endl;
							break;
						}

					}


					if (!isThisProdExcluded)
					{
						next_polystrs_per_list_after_excluded.emplace_back(polystrs_this_list);
						for (auto& polystr : polystrs_this_list)
							lists_num_per_polystr[polystr] ++;
					}

				}

				cout << "Remaining lists num : " << next_polystrs_per_list_after_excluded.size() << endl;

				if (excluded_polystrs.size() == 15)
					break;

				if (next_polystrs_per_list_after_excluded.size() == 0)
					break;

				vector<pair<string, int>> vec_lists_num_per_poly(lists_num_per_polystr.begin(), lists_num_per_polystr.end());
				sort(vec_lists_num_per_poly.begin(), vec_lists_num_per_poly.end(), [&](pair<string, int> a, pair<string, int> b)
					{return a.second > b.second; }
				);


				cout << "Exclude poly " << vec_lists_num_per_poly[0].first << ", which occurs in " << vec_lists_num_per_poly[0].second << " lists." << endl;
				// check if the excluded poly string is in the set of all possible factors
				BooleanPolynomial excluded_poly(target_cipher.statesize, vec_lists_num_per_poly[0].first);
				if (find(all_possible_expanded_p_exps.begin(), all_possible_expanded_p_exps.end(), excluded_poly) == all_possible_expanded_p_exps.end())
					cout << "The excluded poly is not in the set of all possible factors." << endl;
				//

				excluded_polystrs.emplace(vec_lists_num_per_poly[0].first);
				excluded_polys.emplace(BooleanPolynomial(target_cipher.statesize, vec_lists_num_per_poly[0].first));



				
				polystrs_per_list_after_excluded  = next_polystrs_per_list_after_excluded;
			}



			double balancedness = calculate_balancedness(involved_key_bits, polys_per_list);
			cout << "The balancedness is estimated to be " << balancedness << endl;

			cout << endl;
			cout << endl;


			
			
			cout << "Analyzing the superpoly finished." << endl;

		}
	}

	void continue_after_failed(int continue_start, int continue_end, bool continue_solve)
	{
		if (solver_mode == mode::OUTPUT_FILE)
		{
			string TERM_path = R"~(./TERM)~";
			string STATE_path = R"~(./STATE)~";
			FileReader reader;


			// first sort the files 
			using round_pair = pair<int, int>;

			vector<filesystem::path> P_files;
			reader.getJustCurrentFilePaths(STATE_path, P_files);
			vector<filesystem::path> sol_files;
			reader.getJustCurrentFilePaths(TERM_path, sol_files);

			for (auto& P_file : P_files)
			{
				auto P_filename = P_file.filename().string();
				smatch P_filename_sm;
				auto match_status = regex_match(P_filename, P_filename_sm, P_filename_regex);
				if (match_status)
				{
					int start = stoi(P_filename_sm[1]);
					int end = stoi(P_filename_sm[2]);
					
					if (start == continue_start && end == continue_end)
					{
						map<node_pair, BooleanPolynomial> new_P;
						map<node_pair, ListsOfPolynomialsAsFactors> new_P_asLists;
						read_P(P_file.string(), reader, new_P, new_P_asLists);
						P = new_P;
						P_asLists = new_P_asLists;
						cout << "Read P file complete: " << P_filename << endl;
					}
					else if (start >= continue_start && end <= continue_end)
					{
						remove(P_file);
						cout << "Delete P file: " << P_filename << endl;
					}
				}
			}

			for (auto& sol_file : sol_files)
			{
				auto sol_filename = sol_file.filename().string();
				smatch sol_filename_sm;
				auto match_status = regex_match(sol_filename, sol_filename_sm, sol_filename_regex);
				if (match_status)
				{
					int start = stoi(sol_filename_sm[1]);
					int end = stoi(sol_filename_sm[2]);
					if (start >= continue_start && end <= continue_end)
					{
						remove(sol_file);
						cout << "Delete sol_file: " << sol_filename << endl;
					}
				}

			}

			rs = continue_start;
			re = continue_end;

			logger("Continue after failed :" + to_string(continue_start) + "-" + to_string(continue_end));
			logger("Current P size :" + to_string(P.size()));

			if (continue_solve)
				solve_first();
			else
				expand_first();

			
		}
	}
};

#endif

