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
	int d;
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

	vector<Flag> initial_flag;
	vector<BooleanPolynomial> initial_exps;
	dynamic_bitset<> initial_state;
	dynamic_bitset<> cube;

	// the final superpoly
	map<BooleanMonomial, int> superpoly;

	int debug_mode = 0;


public:
	framework(int rounds, int r0, int r1, int first_expand_step, int N, int fbound, int cbound0, int cbound1, int single_threads,
		mode solver_mode, ThreadPool& threadpool, cipher& target_cipher, const set<int>& cube_index, int min_gap)
		: rs(0), re(rounds), target_rounds(rounds), r0(r0), r1(r1), first_expand_step(first_expand_step), N(N), fbound(fbound), cbound0(cbound0), cbound1(cbound1),
		single_threads(single_threads), solver_mode(solver_mode), threadpool(threadpool), target_cipher(target_cipher), min_gap(min_gap)
	{
		cube = target_cipher.generate_cube(cube_index);
		initial_flag = target_cipher.set_initial_flag(cube);
		initial_exps = target_cipher.set_initial_exps(cube);
		initial_state = target_cipher.set_initial_state(cube);
		normal_exp = target_cipher.normal_exp();

		all_rounds_exps = target_cipher.generate_exps(0, rounds, initial_flag, initial_exps, target_cipher.updatelist);


		target_cipher.calculate_flags(0, rounds, initial_flag, all_rounds_flags, true, target_cipher.updatelist);
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

	void debug_test()
	{
		// we write code here to debug 
		
		debug_mode = 1;
		
		/*
		string filepath = R"(./STATE_DEBUG/8_235_P.txt)";
		map<node_pair, BooleanPolynomial> new_P;
		map<node_pair, ListsOfPolynomialsAsFactors> new_P_asLists;
		FileReader reader;
		read_P(filepath, reader, new_P, new_P_asLists);
		for (auto& ss_coef : new_P)
		{
			stringstream ss;
			ss << ss_coef.second;
			if (ss.str() == "s111s167+s66s167")
			{
				P.emplace(ss_coef);
				P_asLists[ss_coef.first] = new_P_asLists[ss_coef.first];
			}
		}

		rs = 8;
		re = 235;

		filterDeg();

		P_output();

		solve_first();
		*/

		
		string from_str = "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011111111111111111111111111111111111111111111111111111111111111111101111110111111111110101111111111111111111111111111111111111111000000000000000000000000000000000000000000000000000000000000000000011111111111111111111111111111111111111110101111111111101111110111111111111111111111111111111111111111111111111111111111111111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
		string to_str = "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110000000000110000000011000000000000011011000000110000001100000000111100000000000111001100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110000000010000000000";



		dynamic_bitset<> from(from_str), to(to_str);

		ListOfPolynomialsAsFactors constant_1_list(target_cipher.default_constant_1);

		P[pair(from, to)] = target_cipher.default_constant_1;
		P_asLists[pair(from, to)] = ListsOfPolynomialsAsFactors(constant_1_list);

		rs = 0;
		re = 319;
		r0 = 1;
		r1 = 20;

		fbound = 300;



		filterDeg();

		P_output();

		expand_first();
		
		
		
		
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

	/*
	void reset_cube(const dynamic_bitset<>& new_cube)
	{
		cube = new_cube;
		initial_flag = target_cipher.set_initial_flag(cube);
		initial_exps = target_cipher.set_initial_exps(cube);
		initial_state = target_cipher.set_initial_state(cube);

		all_rounds_exps = target_cipher.generate_exps(0, target_rounds, initial_flag, initial_exps, target_cipher.updatelist);

		all_rounds_flags.clear();
		all_rounds_flags.shrink_to_fit();
		target_cipher.calculate_flags(0, target_rounds, initial_flag, all_rounds_flags, true, target_cipher.updatelist);
		auto last_flag = all_rounds_flags[target_rounds][0];
		auto last_exps = all_rounds_exps[target_rounds][0];
		vector<vector<vector<Flag>>> ks_flags;
		target_cipher.calculate_flags(target_rounds, target_rounds + 1, last_flag, ks_flags, false, target_cipher.outputks);

		vector<vector<vector<BooleanPolynomial>>> ks_exps = target_cipher.generate_exps(target_rounds, target_rounds + 1, last_flag, last_exps, target_cipher.outputks);

		all_rounds_flags[target_rounds] = ks_flags[0];
		all_rounds_exps[target_rounds] = ks_exps[0];
		all_rounds_flags.emplace_back(ks_flags[1]);
		all_rounds_exps.emplace_back(ks_exps[1]);

		cout << "target_roundsset complete." << endl;
	}

	void reset_cube(const set<int> & new_cube_index)
	{
		cube = target_cipher.generate_cube(new_cube_index);
		initial_flag = target_cipher.set_initial_flag(cube);
		initial_exps = target_cipher.set_initial_exps(cube);
		initial_state = target_cipher.set_initial_state(cube);

		all_rounds_exps = target_cipher.generate_exps(0, target_rounds, initial_flag, initial_exps, target_cipher.updatelist);

		all_rounds_flags.clear();
		all_rounds_flags.shrink_to_fit();
		target_cipher.calculate_flags(0, target_rounds, initial_flag, all_rounds_flags, true, target_cipher.updatelist);
		auto last_flag = all_rounds_flags[target_rounds][0];
		auto last_exps = all_rounds_exps[target_rounds][0];
		vector<vector<vector<Flag>>> ks_flags;
		target_cipher.calculate_flags(target_rounds, target_rounds + 1, last_flag, ks_flags, false, target_cipher.outputks);

		vector<vector<vector<BooleanPolynomial>>> ks_exps = target_cipher.generate_exps(target_rounds, target_rounds + 1, last_flag, last_exps, target_cipher.outputks);

		all_rounds_flags[target_rounds] = ks_flags[0];
		all_rounds_exps[target_rounds] = ks_exps[0];
		all_rounds_flags.emplace_back(ks_flags[1]);
		all_rounds_exps.emplace_back(ks_exps[1]);

		cout << "target_roundsset complete." << endl;
	}

	void compare_cube(int forward_rounds)
	{
		cout << "The number of forward rounds is set to " << forward_rounds << endl;

		auto state0 = target_cipher.set_initial_state(cube);
		dynamic_bitset<> tmp_end_state;
		BooleanPolynomial cur_coef(target_cipher.default_constant_1);

		map<node_pair, BooleanPolynomial> new_P0;
		double expander_time = 2000;
		forwardexpand_thread(0, 0, forward_rounds, state0, tmp_end_state, cur_coef, new_P0, expander_time, 2, false);
		set<dynamic_bitset<>> forward_states0;
		for (auto& nodes_coef : new_P0)
			forward_states0.emplace(nodes_coef.first.first);

		auto cur_flags0 = all_rounds_flags[forward_rounds][0];

		int min_weight0 = 9999;
		for (auto& state : forward_states0)
			if (state.count() < min_weight0)
				min_weight0 = state.count();


		for (int i = 0; i < cube.size(); i++)
		{
			if (true)
			{
				dynamic_bitset<> cube1(cube);


				framework f1(*this);
				f1.reset_cube(cube1);
				auto state1 = target_cipher.set_initial_state(cube1);
				map<node_pair, BooleanPolynomial> new_P1;
				f1.forwardexpand_thread(0, 0, forward_rounds, state1, tmp_end_state, cur_coef, new_P1, expander_time, 2, false);
				set<dynamic_bitset<>> forward_states1;
				for (auto& nodes_coef : new_P1)
					forward_states1.emplace(nodes_coef.first.first);

				auto cur_flags1 = f1.all_rounds_flags[forward_rounds][0];
				bool is_flags_equal = true;
				for(int i = 0;i < target_cipher.statesize; i++)
					if (cur_flags0[i] != cur_flags1[i])
					{
						is_flags_equal = false;
						break;
					}

				int c0 = 0, c1 = 0, c01 = 0;
				for (auto& state : forward_states0)
					if (forward_states1.find(state) == forward_states1.end())
						c0++;
					else
						c01++;

				c1 = forward_states1.size() - c01;

				int min_weight1 = 9999;
				for (auto& state : forward_states1)
					if (state.count() < min_weight1)
						min_weight1 = state.count();

				cout << "Flip cube bit " << i << endl;
				cout << "The number of states generated by original cube: " << forward_states0.size() << " with minimum hamming weight " << min_weight0<<endl;
				cout << "The number of states generated after flipping: " << forward_states1.size() << " with minimum hamming weight " << min_weight1 << endl;
				cout << "The number of states in the intersection: " << c01 << endl;
				cout << "Whether the flags are equal : " << is_flags_equal << endl;
				cout << "----------------------------------------------------------" << endl;
			}
		}
	}
	*/


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
		bitset<128> kreyvium_cube(cubestr);

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

			auto d = computeDegree(kreyvium_cube, re, state);
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

	virtual void start()
	{

		double expander_time = target_cipher.set_expander_time(rs, re);
		first_expand(re, first_expand_step, expander_time, single_threads);
		re -= first_expand_step;

		filterDeg();
	
		P_output();

		solve_first();


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
		target_cipher.build_model(model, rounds-step, rounds, all_rounds_flags, start_vars, end_vars, rounds_p_maps, target_cipher.updatelist, true);

		
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
			for (int i = 0; i < target_cipher.updatelist.size(); i++)
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
		target_cipher.calculate_flags(start, end, middle_start_flags, middle_rounds_flags, true, target_cipher.updatelist);

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
		expandcallback cb = expandcallback(callback_vars, callback_flags, callback_sols);
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
			for (int i = 0; i < target_cipher.updatelist.size(); i++)
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
		target_cipher.calculate_flags(start, end, middle_start_flags, middle_rounds_flags, true, target_cipher.updatelist);

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

		expandcallback cb = expandcallback(callback_vars, callback_flags, callback_sols);
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
			for (int i = 0; i < target_cipher.updatelist.size(); i++)
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
		target_cipher.calculate_flags(start, end, middle_start_flags, middle_rounds_flags, true, target_cipher.updatelist);

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
			midround = start + (end - start) * 3 / 4;

		

		

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

				auto middle_rounds_exps = target_cipher.generate_exps(start, start + maxr + 1, middle_start_flags, normal_exp, target_cipher.updatelist);
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
		if (solver_mode == mode::OUTPUT_FILE)
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

	/*
	int modelGuessNum(const map<string, set<int> >& vars_index_per_poly, int guessNum)
	{
		GRBEnv env = GRBEnv();

		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_Threads, 2);

		GRBModel model = GRBModel(env);

		set<int> vars_indices;
		for (auto& poly_vars : vars_index_per_poly)
			for (auto& index : poly_vars.second)
				vars_indices.emplace(index);

		vector<GRBVar> guessvars(vars_indices.size());
		GRBLinExpr guessnum = 0;
		for (auto& v : guessvars)
		{
			v = model.addVar(0, 1, 0, GRB_BINARY);
			guessnum += v;
		}
		model.addConstr(guessnum == guessNum);

		map<int, int> vars_order;
		int order = 0;
		for (auto& index : vars_indices)
		{
			vars_order[index] = order;
			order++;
		}

		GRBLinExpr guesspolys = 0;
		for (auto& poly_vars : vars_index_per_poly)
		{
			GRBVar polyvar = model.addVar(0, 1, 0, GRB_BINARY);
			vector<GRBVar> reqvars;
			for (auto& index : poly_vars.second)
				reqvars.emplace_back(guessvars[vars_order[index]]);

			model.addGenConstrAnd(polyvar, &(reqvars[0]), reqvars.size());
			guesspolys += polyvar;
		}

		model.setObjective(guesspolys, GRB_MAXIMIZE);

		model.optimize();

		return round(model.get(GRB_DoubleAttr_ObjVal));
	}
	*/
	

	


	void analyze_superpoly_asLists(string TERM_path = R"~(./TERM)~")
	{
		if (solver_mode == mode::OUTPUT_FILE)
		{
			fstream fs;
			string path = TERM_path + string(R"(/superpoly.txt)");
			fs.open(path);


			vector<future<void>> futures;
			map<int, int> lists_per_var;
			map<string, set<int> > vars_index_per_poly;
			vector<vector<string>> polys_per_list;
			int sup_deg = 0;
			long long term_cnt = 0;



			string oneline;
			vector<string> lines;
			while (getline(fs, oneline))
			{
				if (!oneline.empty())
				{
					lines.emplace_back(oneline);
				}
			}
			cout << "Read in total " << lines.size() << " lists." << endl;

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

			cout << "The upperbound for the number of monomials appears in the superpoly : " << term_cnt << endl;

			map<int, int> polys_per_var;
			cout << "The distribution of variable indices in " << vars_index_per_poly.size() << " polys. " << endl;
			for (auto& poly_indices : vars_index_per_poly)
			{
				cout << poly_indices.first << ": " << endl;
				for (auto& index : poly_indices.second)
				{
					cout << index << " ";
					polys_per_var[index]++;
				}
				cout << endl << endl;
			}
			cout << endl;
			cout << endl;


			cout << "The number of polys related to each variable. " << endl;
			for (auto& index_cnt : polys_per_var)
			{
				cout << index_cnt.first << " " << index_cnt.second << endl;
			}

			cout << endl;
			cout << endl;
			

			cout << "The number of lists related to each variable. " << endl;
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

			double balancedness = calculate_balancedness(involved_key_bits, polys_per_list);
			cout << "The balancedness is estimated to be " << balancedness << endl;

			cout << endl;
			cout << endl;

			
			
			cout << "Analyzing the superpoly finished." << endl;

			fs.close();
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

