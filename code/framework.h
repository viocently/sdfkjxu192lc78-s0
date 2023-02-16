#ifndef __FRAMEWORK_H_
#define __FRAMEWORK_H_
#include<sstream>
#include<map>
#include"cipher.h"
#include "thread_pool.h"
#include "deg.h"
#include "log.h"

using namespace thread_pool;
enum class mode { OUTPUT_FILE, OUTPUT_EXP, NO_OUTPUT };
mutex solver_mutex0;
mutex solver_mutex1;
mutex expander_mutex0;
mutex expander_mutex1;

#define REF(x) std::ref(x)

// the callback class
class expandcallback : public GRBCallback
{
public:
	vector<dynamic_bitset<> >* psolutions = NULL;
	vector<GRBVar> solvars;
	vector<Flag> *psolflags;


	expandcallback(vector<GRBVar>& solvars, vector<Flag> & solflags, vector<dynamic_bitset<> >& solutions)
	{
		this->solvars = solvars;
		this->psolutions = &solutions;
		this->psolflags = &solflags;
	}

protected:
	void callback()
	{
		try {
			if (where == GRB_CB_MIPSOL)
			{
				// read solution
				int solvar_num = solvars.size();
				dynamic_bitset<> sol(solvar_num);
				for (int i = 0; i < solvar_num; i++)
					if ((*psolflags)[i] == "delta" && round(getSolution(solvars[i])) == 1)
						sol[i] = 1;

				// add lazy constr
				GRBLinExpr excludeCon = 0;
				for (int i = 0; i < solvar_num; i++)
					if((*psolflags)[i] == "delta")
						if (sol[i] == 1)
							excludeCon += (1 - solvars[i]);
						else
							excludeCon += solvars[i];
				addLazy(excludeCon >= 1);
	
				(*psolutions).emplace_back(sol);

				// cout << "find one solution by callback: " << sol << endl;
			}
		}
		catch (GRBException e) {
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Error during callback" << endl;
		}
	}
};



// the framework
class framework
{
private:

	// the data structure and parameters
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
	cipher & target_cipher;
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
		for (auto& it : P)
		{
			auto& node_pair = it.first;
			auto& end_state = node_pair.second;
			string statestr;
			to_string(end_state, statestr);
			bitset<544> state(statestr);

			auto d = computeDegree(kreyvium_cube, re, state);
			if (d >= kreyvium_cube.count())
				tmp_P.emplace(it);
		}

		P = tmp_P;
		int size1 = P.size();

		logger(__func__ + string(": ") + to_string(size0) + string("\t") +
			to_string(size1));
	}


	inline void P_filter()
	{
		int size0 = P.size();
		map<node_pair, BooleanPolynomial> tmp_P;

		for (auto& nodes_coef : P)
		{
			auto& nodes = nodes_coef.first;
			auto& coef = nodes_coef.second;
			if (!coef.iszero())
				tmp_P[nodes] = coef;
		}

		P = tmp_P;
		int size1 = P.size();
		logger("Filter P: " + to_string(size0) + " " + to_string(size1));
	}

	inline void P_filter(map<node_pair, BooleanPolynomial> & P)
	{
		int size0 = P.size();
		map<node_pair, BooleanPolynomial> tmp_P;

		for (auto& nodes_coef : P)
		{
			auto& nodes = nodes_coef.first;
			auto& coef = nodes_coef.second;
			if (!coef.iszero())
				tmp_P[nodes] = coef;
		}

		P = tmp_P;
		int size1 = P.size();
		logger("Filter P: " + to_string(size0) + " " + to_string(size1));
	}

	inline void superpoly_filter()
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

			auto& start_state = nodes.first;
			auto& end_state = nodes.second;
			os << start_state << endl;
			os << end_state << endl;
			os << coef << endl;
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
				auto expander_time = target_cipher.set_expander_time(rs, re);
				backexpand(rs, re, r1, P, new_P, expander_time, single_threads, threadpool, callback_flag1);

				re -= r1;
				P = new_P;
			}
			else
			{
				map<node_pair, BooleanPolynomial> new_P;
				auto expander_time = target_cipher.set_expander_time(rs, re);
				forwardexpand(rs, re, r0, P, new_P, expander_time, single_threads, threadpool, callback_flag0);

				rs += r0;
				P = new_P;
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
		solve_nodes(rs, re, P, solver_time, single_threads, threadpool, new_P, superpoly);
		P = new_P;
		superpoly_filter();
		logger("Current superpoly size: " + to_string(superpoly.size()));
		logger("Current unsolved nodes: " + to_string(P.size()) );
	}

	virtual void start()
	{

		double expander_time = target_cipher.set_expander_time(rs, re);
		first_expand(re, first_expand_step, expander_time, single_threads);
		re -= first_expand_step;

		filterDeg();
		
		
		
		
		

		P_output();

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

					cout << e_v.first << endl;
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

					auto expand_res = expand_cof(rounds - step, rounds + 1, p_sols, p_exps, p_rms, 0, all_rounds_exps);
					//cout << "State:" << state << endl;
					//cout << expand_res << endl;

					// insert into P
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
		map<dynamic_bitset<>, BooleanPolynomial>& expand_state_coeff, double time, int threads)
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
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_IntParam_Threads, threads);
		// env.set(GRB_IntParam_Presolve, 2);

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

				auto expand_res = expand_cof(start, start + step, p_sols, p_exps, p_rms, 0, all_rounds_exps);
				expand_state_coeff[mid_state] = expand_res;
			}

			if (sol_num == 0)
				return status::NOSOLUTION;
			else
				return status::SOLVED;
		}


	}

	virtual status noncallback_forwardexpand(int start, int step, const dynamic_bitset<>& start_state,
		map<dynamic_bitset<>, BooleanPolynomial>& expand_state_coeff, double time, int threads)
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



					auto expand_res = expand_cof(start, start+step, p_sols, p_exps, p_rms, 0, all_rounds_exps);

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
		map<dynamic_bitset<>, BooleanPolynomial>& expand_state_coeff, double time, int threads)
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
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_IntParam_Threads, threads);
		// env.set(GRB_IntParam_Presolve, 2);

		GRBModel model = GRBModel(env);

		// set initial variables
		vector<GRBVar> start_vars(statesize);

		target_cipher.set_callback_start_cstr(start, middle_start_flags, start_state, model, start_vars);


		vector< vector<vector<pair<BooleanPolynomial,GRBVar>> > > rounds_p_maps;

		vector<GRBVar> helper_midvars;
		target_cipher.callback_build_model_fast(model, start, end-step, middle_rounds_flags, start_vars, helper_midvars);
		vector<Flag> helper_midflags = middle_rounds_flags[end-step][0];
		vector<Flag> midflags = all_rounds_flags[end-step][0];

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
		target_cipher.build_model_fast(model, end-step, end, all_rounds_flags, midvars, end_vars, rounds_p_maps, false);

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

				auto expand_res = expand_cof(end - step, end, p_sols, p_exps, p_rms, 0, all_rounds_exps);
				expand_state_coeff[mid_state] = expand_res;
			}

			if (sol_num == 0)
				return status::NOSOLUTION;
			else
				return status::SOLVED;
		}


	}



	virtual status noncallback_backexpand(int rounds, int step, const dynamic_bitset<> & end_state, 
		map<dynamic_bitset<>, BooleanPolynomial> & expand_state_coeff, double time, int threads)
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

					

					auto expand_res = expand_cof(rounds - step, rounds + 1, p_sols, p_exps, p_rms, 0, all_rounds_exps);
					
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

	virtual void forwardexpand_thread(int start, int end, int step, const dynamic_bitset<>& start_state, const dynamic_bitset<>& end_state,
		const BooleanPolynomial & cur_coef, map<node_pair,BooleanPolynomial> & new_P,
		 double time, int threads, bool iscallback)
	{
		map<dynamic_bitset<>, BooleanPolynomial> expand_state_coeff;

		if (iscallback)
		{
			auto expand_status = callback_forwardexpand(start, end, step, start_state, end_state, expand_state_coeff, time, threads);
			if (expand_status == status::UNSOLVED)
			{
				auto expand_status2 = noncallback_forwardexpand(start, step, start_state, expand_state_coeff, time, threads);
				if (expand_status2 == status::UNSOLVED)
				{
					cerr << __func__ << "Forwardexpand failed." << endl;
					exit(-1);
				}
			}
		}
		else
		{
			auto expand_status = noncallback_forwardexpand(start, step, start_state, expand_state_coeff, time, threads);
			if (expand_status == status::UNSOLVED)
			{
				cerr << __func__ << "Forwardexpand failed." << endl;
				exit(-1);
			}
		}

		map<node_pair, BooleanPolynomial> part_P;
		for (auto& state_coef : expand_state_coeff)
		{
			auto& state = state_coef.first;
			auto& coef = state_coef.second;
			part_P[pair(state, end_state)] = coef * cur_coef;
		}

		// update the new_P
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



	virtual void forwardexpand(int start, int end, int step, const map<node_pair, BooleanPolynomial> & cur_P,
		map<node_pair, BooleanPolynomial> & new_P,
		double time, int threads, ThreadPool& thread_pool, bool iscallback)
	{
		vector<future<void>> futures;


		for (auto& nodes_coef : cur_P)
		{
			auto& start_state = nodes_coef.first.first;
			auto& end_state = nodes_coef.first.second;
			auto& cur_coef = nodes_coef.second;
			futures.emplace_back(thread_pool.Submit(&framework::forwardexpand_thread, REF(*this), start, end, step, REF(start_state), REF(end_state),
				REF(cur_coef), REF(new_P), time, threads, iscallback));
		}

		for (auto& it : futures)
			it.get();
	}

	virtual void backexpand_thread(int start, int end, int step, const dynamic_bitset<>& start_state, const dynamic_bitset<>& end_state,
		const BooleanPolynomial& cur_coef, map<node_pair, BooleanPolynomial>& new_P,
		double time, int threads, bool iscallback)
	{
		map<dynamic_bitset<>, BooleanPolynomial> expand_state_coeff;

		
		if (iscallback)
		{
			auto expand_status = callback_backexpand(start, end, step, start_state, end_state, expand_state_coeff, time, threads);
			if (expand_status == status::UNSOLVED)
			{
				auto expand_status2 = noncallback_backexpand(end, step, end_state, expand_state_coeff, time, threads);
				if (expand_status2 == status::UNSOLVED)
				{
					cerr << __func__ << "Backexpand failed." << endl;
					exit(-1);
				}
			}
		}
		else
		{
			auto expand_status = noncallback_backexpand(end, step, end_state, expand_state_coeff, time, threads);
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

		// update the new_P
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



	virtual void backexpand(int start, int end, int step, const map<node_pair, BooleanPolynomial>& cur_P,
		map<node_pair, BooleanPolynomial>& new_P,
		double time, int threads, ThreadPool& thread_pool, bool iscallback)
	{
		vector<future<void>> futures;

	
		for (auto& nodes_coef : cur_P)
		{
			auto& start_state = nodes_coef.first.first;
			auto& end_state = nodes_coef.first.second;
			auto& cur_coef = nodes_coef.second;
			futures.emplace_back(thread_pool.Submit(&framework::backexpand_thread, REF(*this), start, end, step, REF(start_state), REF(end_state),
				REF(cur_coef), REF(new_P), time, threads, iscallback));
		}

		for (auto& it : futures)
			it.get();

	}



	

	// the function for output to file
	inline void print_sol(int start, int end, const dynamic_bitset<>& start_state, const dynamic_bitset<>& end_state,
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

	// the function for expanding the coefficient
	inline BooleanPolynomial expand_cof(int start, int end, 
		const vector<dynamic_bitset<>> & p_sols, const vector<BooleanPolynomial>& p_exps, const vector<pair<int, int>>& p_rms,
		int expandto, vector<vector<vector<BooleanPolynomial> >> & expandto_exps )
	{
		if (p_sols.size() == 0)
			return target_cipher.default_constant_0;

		
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
				// cout << "At " << p_rms[i].first << "-" << p_rms[i].second << endl;
				// cout << p_exps[i] << " sub to " << p_expands[i] << endl;
			}
		}


		vector<BooleanPolynomial> expand_sols;


		for (auto& p_sol : p_sols)
		{
			// cout << p_sol.count() << endl;

			vector<BooleanPolynomial> expand_comps;
			for (int i = 0; i < p_sol.size(); i++)
				if (p_sol[i] == 1)
				{
					expand_comps.emplace_back(p_expands[i]);
					//cout << p_expands[i] << " ";
				}

			//cout << endl;

			BooleanPolynomial expand_sol = fastmul(target_cipher.statesize, expand_comps);
			//cout << expand_sol << endl;
			expand_sols.emplace_back(expand_sol);
			
		}

		auto expand_res = fastsum(target_cipher.statesize, expand_sols);
		return expand_res;
	}

	// the thread for solver
	virtual void solve_nodes_thread(int start, int end, const vector<Flag>& start_flags, const dynamic_bitset<>& start_state,
		const dynamic_bitset<>& end_state, const BooleanPolynomial & cur_coef, double time, int threads, map<node_pair, BooleanPolynomial> & new_P,
		map<BooleanMonomial, int>& sup_counter)
	{
		// generate local rounds flags for solver
		vector<Flag> middle_start_flags(start_flags);
		for (int i = 0; i < target_cipher.statesize; i++)
			if (start_flags[i] == "delta" && start_state[i] == 0)
				middle_start_flags[i] = "zero_c";

		vector<vector<vector<Flag>>> middle_rounds_flags;
		target_cipher.calculate_flags(start, end, middle_start_flags, middle_rounds_flags, true, target_cipher.updatelist);

			

		vector<dynamic_bitset<>> p_sols;
		vector<BooleanPolynomial> p_exps;
		vector<pair<int, int>> p_rms;
		set<int> end_constants;
		auto solver_status = target_cipher.solve_model(start, end, middle_rounds_flags, start_state, end_state, p_sols, p_exps, p_rms, end_constants, threads, time, true);

		if (solver_status == status::SOLVED)
		{
			if (solver_mode == mode::NO_OUTPUT)
			{
				; // do nothing
			}
			else if (solver_mode == mode::OUTPUT_FILE)
			{
				print_sol(start, end, start_state, end_state, p_sols, p_exps, p_rms, end_constants);
			}
			else if (solver_mode == mode::OUTPUT_EXP)
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
				
				BooleanPolynomial second_expand_res;
				// vector<BooleanPolynomial> subexps;
				// BooleanPolynomial first_expand_res;
				{
					vector<BooleanPolynomial> subexps;
					{
						BooleanPolynomial first_expand_res;
						{
							auto middle_rounds_exps = target_cipher.generate_exps(start, start + maxr + 1, middle_start_flags, normal_exp, target_cipher.updatelist);
							first_expand_res = expand_cof(start, end, p_sols, p_exps, p_rms, start, middle_rounds_exps);
							if (has_constant)
							{
								auto end_constants_expand = target_cipher.default_constant_1;
								for (auto& i : end_constants)
								{
									end_constants_expand *= middle_rounds_exps[maxr][0][i];
									// cout << i << endl;
								}
								first_expand_res *= end_constants_expand;
							}
						}

						//cout << "From " << first_expand_res << endl;
						for (auto& mon : first_expand_res)
						{
							auto subexp = mon.subs(all_rounds_exps[start][0]);
							// cout << mon << " Sub To " << subexp << endl;
							if (!subexp.iszero())
								subexps.emplace_back(subexp);
						}
					}

					second_expand_res = fastsum(target_cipher.statesize, subexps);
				}
	

				cout << __func__ << " : To multiply " << second_expand_res.moncnt() << " with " << cur_coef.moncnt() << endl;


				auto contr = second_expand_res * cur_coef;

				

				//cout << "To " << second_expand_res << endl;
				lock_guard<mutex> guard(solver_mutex0);
				for (auto& mon : contr)
					sup_counter[mon] ++;

			}
		}
		else if (solver_status == status::UNSOLVED)
		{
			lock_guard<mutex> guard(solver_mutex1);
			new_P[pair(start_state, end_state)] = cur_coef;
		}
	}

	virtual void solve_nodes(int start, int end, const map<node_pair, BooleanPolynomial>& cur_P,
		double time, int threads, ThreadPool& thread_pool, map<node_pair, BooleanPolynomial>& new_P, map<BooleanMonomial, int>& sup_counter)
	{
		auto start_flags = all_rounds_flags[start][0];

		vector<future<void>> futures;



		for (auto& nodes_coef : cur_P)
		{
			auto& start_state = nodes_coef.first.first;
			auto& end_state = nodes_coef.first.second;
			auto& cur_coef = nodes_coef.second;
			futures.emplace_back(thread_pool.Submit(&framework::solve_nodes_thread, REF(*this), start, end, REF(start_flags), REF(start_state),
				REF(end_state), REF(cur_coef), time, threads, REF(new_P), REF(sup_counter) ));
		}

		for (auto& it : futures)
			it.get();

	}
};

#endif
