#ifndef __CIPHER_H_
#define __CIPHER_H_

#include<iostream>
#include<vector>
#include<map>
#include "gurobi_c++.h"
#include"BooleanPolynomial.h"
#include"flag.h"
#include "log.h"
#include "comparator.h"

using namespace std;

#define MAX 200000000

enum class status {UNSOLVED, SOLVED, NOSOLUTION};





class cipher
{
public:
	string ciphername;
	int statesize;
	int keysize;
	int ivsize;


	using update_element = map<int, BooleanPolynomial>;
	/* tips: the type update_element doesn't allow for example s12->1, s12->2. You may split this into two steps*/

	using new_monrep = pair<BooleanMonomial, BooleanMonomial>;
	using new_polyrep = map<BooleanMonomial, BooleanPolynomial>;

	vector<update_element> updatelist;
	vector<update_element> outputks;

	BooleanPolynomial default_constant_1;
	BooleanPolynomial default_constant_0;

	// for resolving the updatelist
	vector<vector<vector<int>>> round_COPY_v_i;
	vector<vector<vector<int>>> round_AND_v_i;
	vector<vector<vector<int>>> round_XOR_v_i;
	vector<map<int, int>> round_rev_COPY_v_i;
	vector<vector<int>> round_update_in;
	vector<vector<int>> round_update_out;
	vector<int> round_COPY_v_num;

	vector<map<int, int> > round_moves;
	vector<map<int, dynamic_bitset<>>  > round_AND_rep;
	



	cipher()
	{
		ciphername = "undefined";
		statesize = -1;
		keysize = -1;
		ivsize = -1;
		// define_updatelist();
	}

	// these member functions should be redefined for each cipher
	virtual void define_updatelist() = 0;

	virtual void define_outputks() = 0;

	virtual dynamic_bitset<> generate_cube(const set<int>& cube_index) = 0;

	virtual dynamic_bitset<> set_initial_state(const dynamic_bitset<> & cube) = 0;

	virtual vector<Flag> set_initial_flag(const dynamic_bitset<>& cube) = 0;

	virtual vector<BooleanPolynomial>  set_initial_exps(const dynamic_bitset<>& cube) = 0;

	virtual double set_solver_time(int rs, int re) = 0;

	virtual double set_expander_time(int rs, int re) = 0;



public:
	virtual void set_start_cstr(int start, const vector<Flag>& start_flags, const dynamic_bitset<>& start_state, GRBModel& model, vector<GRBVar>& start_vars)
	{
		start_vars.resize(statesize);

		for (int i = 0; i < statesize; i++)
			if (start_flags[i] == "delta")
			{
				start_vars[i] = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(start_vars[i] == start_state[i]);
			}
	}

	virtual void set_callback_start_cstr(int start, const vector<Flag>& start_flags, const dynamic_bitset<>& start_state, GRBModel& model, vector<GRBVar>& start_vars)
	{
		start_vars.resize(statesize);

		for (int i = 0; i < statesize; i++)
			if (start_flags[i] == "delta")
			{
				start_vars[i] = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(start_vars[i] == start_state[i]);
			}
	}

	virtual void resolve_updatelist()
	{
		for (auto& ele_map : updatelist)
		{
			vector<vector<int>> COPY_v_i(statesize);
			map<int, int> rev_COPY_v_i;
			vector<vector<int>> AND_v_i;
			vector<vector<int>> XOR_v_i(statesize);

			int mon_i = 0;
			int copy_v_i = 0;


			set<int> have_determined;
			map<int, dynamic_bitset<> > AND_rep;
			for (auto& fromto : ele_map)
			{
				auto& from = fromto.second;
				auto& to = fromto.first;

				have_determined.emplace(to);

				for (auto& mon : from)
				{
					XOR_v_i[to].emplace_back(mon_i);

					vector<int> this_mon_v_i;
					for (auto& i : mon.index())
					{
						this_mon_v_i.emplace_back(copy_v_i);
						COPY_v_i[i].emplace_back(copy_v_i);
						rev_COPY_v_i[copy_v_i] = i;
						copy_v_i++;
					}

					AND_rep[mon_i] = mon.rep();

					AND_v_i.emplace_back(this_mon_v_i);
					mon_i++;
				}
			}

			// the others are considered unchanged
			for (int i = 0; i < statesize; i++)
			{
				auto it = have_determined.find(i);
				if (it == have_determined.end())
				{
					XOR_v_i[i].emplace_back(mon_i);
					COPY_v_i[i].emplace_back(copy_v_i);
					rev_COPY_v_i[copy_v_i] = i;
					vector<int> this_mon_v_i = { copy_v_i };
					AND_v_i.emplace_back(this_mon_v_i);

					dynamic_bitset<> mon_rep(statesize);
					mon_rep[i] = 1;
					AND_rep[mon_i] = mon_rep;

					mon_i++;
					copy_v_i++;
				}
			}

			round_COPY_v_num.emplace_back(copy_v_i);



			// analyze
			map<int, int> back_COPY_equals;
			for (int i = 0; i < statesize; i++)
			{
				if (COPY_v_i[i].size() == 1)
				{
					back_COPY_equals[COPY_v_i[i][0]] = i;
				}
			}

			map<int, int> back_AND_equals;
			for (int i = 0; i < AND_v_i.size(); i++)
			{
				if (AND_v_i[i].size() == 1)
				{
					back_AND_equals[i]= AND_v_i[i][0];
				}
			}

			map<int, int> back_XOR_equals;
			for (int i = 0; i < statesize; i++)
			{
				if (XOR_v_i[i].size() == 1)
				{
					back_XOR_equals[i] = XOR_v_i[i][0];
				}
			}

			map<int, int> move;
			for (auto& i_j : back_XOR_equals)
			{
				auto& i = i_j.first;
				auto& j = i_j.second;
				if (back_AND_equals.find(j) != back_AND_equals.end())
				{
					auto& k = back_AND_equals[j];
					if (back_COPY_equals.find(k) != back_COPY_equals.end())
					{
						auto& l = back_COPY_equals[k];
						move[l] = i;
					}
				}
			}

			round_moves.emplace_back(move);

			// update COPY_v_i, ..., XOR_v_i and rev_COPY_v_i
			for (auto& i_j : move)
			{
				auto& i = i_j.first;
				auto& j = i_j.second;
				rev_COPY_v_i.erase(COPY_v_i[i][0]);

				COPY_v_i[i].clear();
				for (auto& k : XOR_v_i[j])
				{
					AND_v_i[k].clear();
					AND_rep.erase(k);
				}
				XOR_v_i[j].clear();
			}

			round_COPY_v_i.emplace_back(COPY_v_i);
			round_AND_v_i.emplace_back(AND_v_i);
			round_XOR_v_i.emplace_back(XOR_v_i);
			round_rev_COPY_v_i.emplace_back(rev_COPY_v_i);
			round_AND_rep.emplace_back(AND_rep);

			vector<int> update_in;
			dynamic_bitset<> update_out_mask(statesize);
			update_out_mask.flip();
			for (int i = 0; i < statesize; i++)
				if (move.find(i) == move.end())
					update_in.emplace_back(i);
				else
					update_out_mask[move[i]] = 0;
			round_update_in.emplace_back(update_in);

			vector<int> update_out;
			for (int i = 0; i < statesize; i++)
				if (update_out_mask[i])
					update_out.emplace_back(i);

			round_update_out.emplace_back(update_out);



			// 
			
		}
	}

	// this is the common part for computing the flags
	virtual void calculate_onemap_flag(const vector<Flag>& prev_flag, vector<Flag>& next_flag, const update_element & ele_map) const
	{
		if (prev_flag.size() != statesize)
			return;

		next_flag.resize(statesize);
		
		for (int i = 0; i < statesize; i++)
		{
			next_flag[i] = prev_flag[i];
		}

		for (auto& fromto : ele_map)
		{
			auto from = fromto.second;
			auto to = fromto.first;
			Flag fromflag("zero_c");
			for (auto& mon : from)
			{
				Flag monflag("one_c_1");
				auto indexes = mon.index();


				for (auto& index : indexes)
				{
					// cout << index << "|";
					// cout << prev_flag[index] << "    ";
					monflag *= prev_flag[index];
				}
				
				fromflag += monflag;

				// cout << endl;
			
			}

			// cout << "t flag :" << fromflag << endl;
			next_flag[to] = fromflag;
		}
	}

	virtual void calculate_oneround_flag(const vector<Flag>& prev_flag, vector<Flag> & next_flag, vector<vector<Flag>> & round_flags,
		const vector<update_element> & updatefuncs) const
	{
		if (prev_flag.size() != statesize)
			return;

		auto tmp_prev_flag = prev_flag;
		for (auto& ele_map : updatefuncs)
		{
			round_flags.emplace_back(tmp_prev_flag);
			
			vector<Flag> tmp_next_flag;
			calculate_onemap_flag(tmp_prev_flag, tmp_next_flag, ele_map);
			tmp_prev_flag = tmp_next_flag;
		}

		next_flag = tmp_prev_flag;
	}

	// only for round functions
	virtual void build_round_model_fast(GRBModel& model, vector<GRBVar>& in_vars, vector<GRBVar>& out_vars,
		const vector<vector<Flag>>& round_flags, vector<vector<pair<BooleanPolynomial, GRBVar>>>& p_maps)
	{
		int map_num = updatelist.size();
		vector<GRBVar> map_in_vars(in_vars);
		for (int map_i = 0; map_i < map_num; map_i++)
		{
			vector<GRBVar> map_out_vars(statesize);

			vector<pair<BooleanPolynomial, GRBVar> > p_map;

			auto& in_flags = round_flags[map_i];
			auto& COPY_v_i = round_COPY_v_i[map_i];
			auto& rev_COPY_v_i = round_rev_COPY_v_i[map_i];
			auto& AND_v_i = round_AND_v_i[map_i];
			auto& XOR_v_i = round_XOR_v_i[map_i];
			auto& move = round_moves[map_i];
			auto& update_out = round_update_out[map_i];
			auto& update_in = round_update_in[map_i];
			auto& AND_rep = round_AND_rep[map_i];
			auto& COPY_v_num = round_COPY_v_num[map_i];

			for (auto & fromto : move)
			{
				auto& from = fromto.first;
				auto& to = fromto.second;
				map_out_vars[to] = map_in_vars[from];
			}

			dynamic_bitset<> one_c_mask(statesize);
			dynamic_bitset<> delta_mask(statesize);
			dynamic_bitset<> zero_c_mask(statesize);

			for (int j = 0;j< statesize;j++)
				if (in_flags[j] == "delta")
					delta_mask[j] = 1;
				else if (in_flags[j] == "zero_c")
					zero_c_mask[j] = 1;
				else
					one_c_mask[j] = 1;

			map<int, GRBVar> AND_vars;
			for (auto& j : update_out)
			{
				map<dynamic_bitset<>, int> mon_rep_i;
				map<dynamic_bitset<>, vector<dynamic_bitset<>>> mon_cof;
				for (auto& mon_i : XOR_v_i[j])
				{
					auto& mon_rep = AND_rep[mon_i];
					auto zero_c_part = mon_rep & zero_c_mask;
					if (zero_c_part.count() == 0)
					{
						auto delta_part = mon_rep & delta_mask;
						auto one_c_part = mon_rep & one_c_mask;
						auto iter = mon_cof.find(delta_part);
						if (iter == mon_cof.end())
						{
							mon_rep_i[delta_part] = mon_i;
							mon_cof[delta_part] = { one_c_part };
						}
						else
							iter->second.emplace_back(one_c_part);
					}
				}

				// constant 0
				if (mon_cof.size() == 0)
					continue;

				// non-zero constant
				if (mon_cof.size() == 1 && (*mon_cof.begin()).first.count() == 0)
					continue;

				// XOR
				map_out_vars[j] = model.addVar(0, 1, 0, GRB_BINARY);

				GRBLinExpr xor_exp = 0;
				for (auto& rep_i : mon_rep_i)
				{
					auto& rep = rep_i.first;
					auto& i = rep_i.second;
					AND_vars[i] = model.addVar(0, 1, 0, GRB_BINARY);
					xor_exp += AND_vars[i];
					if (mon_cof[rep].size() == 1 && mon_cof[rep][0].count() == 0)
						continue;
					else
						p_map.emplace_back(pair(BooleanPolynomial(statesize, mon_cof[rep]), AND_vars[i]));
				}

				model.addConstr(xor_exp == map_out_vars[j]);

			}


			dynamic_bitset<> COPY_v_delta_mask(COPY_v_num);
			for (auto& j : update_in)
				if(in_flags[j] == "delta")
					for (auto& k : COPY_v_i[j])
						COPY_v_delta_mask[k] = 1;

			// AND
			map<int, GRBVar> COPY_vars;
			for (auto& i_v : AND_vars)
			{
				auto& i = i_v.first;
				auto& v = i_v.second;
				for (auto& j : AND_v_i[i])
					if (COPY_v_delta_mask[j])
						COPY_vars[j] = v;
			}

			// COPY
			map<int, vector<GRBVar>> COPY_from_to;
			for (auto& i_v : COPY_vars)
			{
				auto& i = i_v.first;
				auto& v = i_v.second;
				auto& from = rev_COPY_v_i[i];
				COPY_from_to[from].emplace_back(v);
			}

			for(auto& j : update_in)
				if (in_flags[j] == "delta")
				{
					auto &from_v = map_in_vars[j];
					auto &to_vs = COPY_from_to[j];
					if (to_vs.size() == 0)
						model.addConstr(from_v == 0);
					else
						model.addGenConstrOr(from_v, &(to_vs[0]), to_vs.size());
				}


			p_maps.emplace_back(p_map);
			map_in_vars = map_out_vars;






		}

		out_vars = map_in_vars;
	}

	// CBDP based propagation
	virtual void callback_build_round_model_fast(GRBModel& model, vector<GRBVar>& in_vars, vector<GRBVar>& out_vars,
		const vector<vector<Flag>>& round_flags)
	{
		int map_num = updatelist.size();
		vector<GRBVar> map_in_vars(in_vars);
		for (int map_i = 0; map_i < map_num; map_i++)
		{
			vector<GRBVar> map_out_vars(statesize);


			auto& in_flags = round_flags[map_i];
			auto& COPY_v_i = round_COPY_v_i[map_i];
			auto& rev_COPY_v_i = round_rev_COPY_v_i[map_i];
			auto& AND_v_i = round_AND_v_i[map_i];
			auto& XOR_v_i = round_XOR_v_i[map_i];
			auto& move = round_moves[map_i];
			auto& update_out = round_update_out[map_i];
			auto& update_in = round_update_in[map_i];
			auto& AND_rep = round_AND_rep[map_i];
			auto& COPY_v_num = round_COPY_v_num[map_i];

			for (auto& fromto : move)
			{
				auto& from = fromto.first;
				auto& to = fromto.second;
				map_out_vars[to] = map_in_vars[from];
			}

			dynamic_bitset<> one_c_mask(statesize);
			dynamic_bitset<> delta_mask(statesize);
			dynamic_bitset<> zero_c_mask(statesize);

			for (int j = 0; j < statesize; j++)
				if (in_flags[j] == "delta")
					delta_mask[j] = 1;
				else if (in_flags[j] == "zero_c")
					zero_c_mask[j] = 1;
				else
					one_c_mask[j] = 1;

			map<int, GRBVar> AND_vars;
			for (auto& j : update_out)
			{
				map<dynamic_bitset<>, int> mon_rep_i;
				map<dynamic_bitset<>, vector<dynamic_bitset<>>> mon_cof;
				for (auto& mon_i : XOR_v_i[j])
				{
					auto& mon_rep = AND_rep[mon_i];
					auto zero_c_part = mon_rep & zero_c_mask;
					if (zero_c_part.count() == 0)
					{
						auto delta_part = mon_rep & delta_mask;
						auto one_c_part = mon_rep & one_c_mask;
						auto iter = mon_cof.find(delta_part);
						if (iter == mon_cof.end())
						{
							mon_rep_i[delta_part] = mon_i;
							mon_cof[delta_part] = { one_c_part };
						}
						else
							iter->second.emplace_back(one_c_part);
					}
				}

				// constant 0
				if (mon_cof.size() == 0)
					continue;

				// non-zero constant
				if (mon_cof.size() == 1 && (*mon_cof.begin()).first.count() == 0)
					continue;

				// XOR
				map_out_vars[j] = model.addVar(0, 1, 0, GRB_BINARY);

				GRBLinExpr xor_exp = 0;
				for (auto& rep_i : mon_rep_i)
				{
					auto& rep = rep_i.first;
					auto& i = rep_i.second;
					AND_vars[i] = model.addVar(0, 1, 0, GRB_BINARY);
					xor_exp += AND_vars[i];
					if (mon_cof[rep].size() == 1 && mon_cof[rep][0].count() == 0)
						continue;
				}

				model.addConstr(xor_exp == map_out_vars[j]);

			}


			dynamic_bitset<> COPY_v_delta_mask(COPY_v_num);
			for (auto& j : update_in)
				if (in_flags[j] == "delta")
					for (auto& k : COPY_v_i[j])
						COPY_v_delta_mask[k] = 1;

			// AND
			map<int, GRBVar> COPY_vars;
			for (auto& i_v : AND_vars)
			{
				auto& i = i_v.first;
				auto& v = i_v.second;
				vector<GRBVar> from_vars;
				for (auto& j : AND_v_i[i])
					if (COPY_v_delta_mask[j])
					{
						auto from_v = model.addVar(0, 1, 0, GRB_BINARY);
						from_vars.emplace_back(from_v);
						COPY_vars[j] = from_v;
					}

				if (from_vars.size() == 0)
					model.addConstr(v == 0);
				else
					model.addGenConstrOr(v, &(from_vars[0]), from_vars.size());
			}

			// COPY
			map<int, vector<GRBVar>> COPY_from_to;
			for (auto& i_v : COPY_vars)
			{
				auto& i = i_v.first;
				auto& v = i_v.second;
				auto& from = rev_COPY_v_i[i];
				COPY_from_to[from].emplace_back(v);
			}

			for (auto& j : update_in)
				if (in_flags[j] == "delta")
				{
					auto& from_v = map_in_vars[j];
					auto& to_vs = COPY_from_to[j];
					if (to_vs.size() == 0)
						model.addConstr(from_v == 0);
					else
					{
						GRBLinExpr xor_exp = 0;
						for (auto& v : to_vs)
							xor_exp += v;
						model.addConstr(from_v == xor_exp);
					}
				}


			map_in_vars = map_out_vars;






		}

		out_vars = map_in_vars;
	}

	virtual void calculate_flags(int start, int end, const vector<Flag>& start_flag, vector<vector<vector<Flag>>> & rounds_flags, bool padding,
		const vector<update_element>& updatefuncs) const
	{
		if (start_flag.size() != statesize)
			return;

		if (padding)
		{
			for (int r = 0; r < start; r++)
			{
				vector<vector<Flag>>  padding_round_flags;
				rounds_flags.emplace_back(padding_round_flags);
			}
		}

		vector<Flag> cur_flag = start_flag;
		for (int r = start; r < end; r++)
		{


			vector<vector<Flag> > round_flags;
			vector<Flag> next_flag;
			calculate_oneround_flag(cur_flag, next_flag, round_flags, updatefuncs);
			rounds_flags.emplace_back(round_flags);


			cur_flag = next_flag;
		}

		// the flags of last round
		vector<vector<Flag>> last_flag = { cur_flag };
		rounds_flags.emplace_back(last_flag);
	}

	// this is the common part for constructing MILP models
		// subfunction
	inline new_monrep divide_mon(const BooleanMonomial& mon, const vector<Flag> & flags)
	{
		new_monrep new_mon;

		int flagnum = flags.size();


		map<int, int> submap;
		dynamic_bitset<> one_c_k_mask(flagnum);
		dynamic_bitset<> delta_mask(flagnum);
		for (int i = 0; i < flagnum; i++)
		{
			if (flags[i] == "zero_c")
				submap[i] = 0;
			else if (flags[i] == "one_c_1")
				submap[i] = 1;
			else if (flags[i] == "one_c_k")
				one_c_k_mask[i] = 1;
			else if (flags[i] == "delta")
				delta_mask[i] = 1;
			else
			{
				cerr << __func__ << ": Undefined flag." << endl;
				exit(-1);
			}
		}

		
		BooleanMonomial tmpmon(mon);
		tmpmon.subs(submap);
		new_mon.first = tmpmon.intersect(one_c_k_mask);
		new_mon.second = tmpmon.intersect(delta_mask);

		return new_mon;
		
	}



	// this modeling method is too general
	virtual void build_onemap_model(GRBModel& model, vector<GRBVar>& in_vars, vector<GRBVar>& out_vars, 
		const update_element& ele_map, const vector<Flag> & in_flags, 
		vector<pair<BooleanPolynomial, GRBVar>> & p_map)
	{

		if (in_vars.size() != statesize || out_vars.size() != statesize)
		{
			cerr << __func__ << ": The number of input or output vars is invalid.";
			exit(-1);
		}

		if (in_flags.size() != statesize )
		{
			cerr << __func__ << ": The number of input flags is invalid.";
			exit(-1);
		}

		
		set<int> have_determined;


		vector<vector<int>> XOR_v_i(statesize);

		vector<vector<int>> COPY_v_i(statesize);

		int v_i = 0;

		vector<int> c_v_i;

		vector<pair<BooleanPolynomial, int> > p_map_helper;



		for (auto& fromto : ele_map)
		{
			// cout << "From " << fromto.second << " To " << fromto.first << endl;
			new_polyrep new_poly;
			for (auto& mon : fromto.second)
			{
				auto new_mon = divide_mon(mon, in_flags);
				// cout << new_mon.first << "*" << new_mon.second << endl;

				if (new_mon.second.iszero())
					continue;

				if (new_poly.find(new_mon.second) == new_poly.end())
					new_poly[new_mon.second] = BooleanPolynomial(new_mon.first);
				else
					new_poly[new_mon.second] += new_mon.first;
			}

			new_polyrep tmp_new_poly;
			for (auto& new_term : new_poly)
				if (!new_term.second.iszero())
					tmp_new_poly.emplace(new_term);
			new_poly = tmp_new_poly;

			// cout << "New poly size : " << new_poly.size() << endl;

			// check if new_poly is constant
			if (new_poly.size() == 1 && (*new_poly.begin()).first.isone())
			{
				// non-zero constant
				;
			}
			else if(new_poly.size() == 0)
			{
				// constant 0
				;
			}
			else
			{

				for (auto& new_term : new_poly)
				{

					XOR_v_i[fromto.first].emplace_back(v_i);
					if (!new_term.second.isone())
					{
						p_map_helper.emplace_back(pair(new_term.second, v_i));
						cout << "p constant: "<<new_term.second << endl;
					}

					if (new_term.first.isone())
						c_v_i.emplace_back(v_i);
					else
					{
						for (auto& index : new_term.first.index())
						{
							COPY_v_i[index].emplace_back(v_i);
							//cout << "Tag 1:" << v_i << " in " << index << endl;
						}
					}
					v_i++;
				}
			}

			// 
			have_determined.emplace(fromto.first);
		}
		

		// We assume the other vars remain unchanged
		for(int i = 0;i <statesize;i++)
			if (have_determined.find(i) == have_determined.end())
			{
				if (in_flags[i] == "delta")
				{
					XOR_v_i[i].emplace_back(v_i);
					COPY_v_i[i].emplace_back(v_i);
					//cout << "Tag 2:"<<v_i << endl;
					v_i++;
				}
			}

		// Turn v_i into MILP variables
		// Add constraints for COPY and XOR
		map<int, GRBVar> v_map;
		for (auto& each_v_i : c_v_i)
		{
			v_map[each_v_i] = model.addVar(0, 1, 0, GRB_BINARY);
		}

		for (int i = 0; i < statesize; i++)
			if (COPY_v_i[i].size() == 1)
			{
				int this_v_i = COPY_v_i[i][0];
				auto it = v_map.find(this_v_i);
				if (it == v_map.end())
				{
					//cout << "This_v_i: " << this_v_i << endl;
					v_map[this_v_i] = in_vars[i];
				}
				else
					model.addConstr(it->second == in_vars[i]);

			}
			else if (COPY_v_i[i].size() > 1)
			{
				vector<GRBVar> COPY_vars;
				for (auto& each_v_i : COPY_v_i[i])
				{
					auto it = v_map.find(each_v_i);
					if (it == v_map.end())
					{
						GRBVar a = model.addVar(0, 1, 0, GRB_BINARY);
						v_map[each_v_i] = a;
						//cout << "Each_v_i: " << each_v_i << endl;
						COPY_vars.emplace_back(a);
					}
					else
						COPY_vars.emplace_back(it->second);
				}

				model.addGenConstrOr(in_vars[i], &(COPY_vars[0]), COPY_vars.size());
			}
			else if (COPY_v_i[i].size() == 0 && in_flags[i] == "delta")
				model.addConstr(in_vars[i] == 0);


		for (int i = 0; i < statesize; i++)
		{
			if (XOR_v_i[i].size() == 1)
			{
				int this_v_i = XOR_v_i[i][0];
				out_vars[i] = v_map[this_v_i];
			}
			else if(XOR_v_i[i].size() > 1)
			{
				vector<GRBVar> XOR_vars;
				for (auto& each_v_i : XOR_v_i[i])
					XOR_vars.emplace_back(v_map[each_v_i]);

				GRBLinExpr xor_exp = 0;
				for (auto& v : XOR_vars)
					xor_exp += v;
				
				out_vars[i] = model.addVar(0, 1, 0, GRB_BINARY);
				model.addConstr(out_vars[i] == xor_exp);
			}
		}

		// set p_map
		for (auto& p_vi : p_map_helper)
			p_map.emplace_back(pair(p_vi.first, v_map[p_vi.second]));
	}

	virtual void build_oneround_model(GRBModel& model, vector<GRBVar>& in_vars, vector<GRBVar>& out_vars,
		const vector<vector<Flag>> & round_flags, vector<vector<pair<BooleanPolynomial, GRBVar>> >& p_maps,
		const vector<update_element> & updatefuncs)
	{
		if (in_vars.size() != statesize)
		{
			cerr << __func__ << ": The number of input vars is invalid.";
			exit(-1);
		}

		vector<GRBVar> ele_in_vars(in_vars);
		int ele_in_flags_i = 0;
		


		for (auto& ele_map : updatefuncs)
		{
			vector<Flag> ele_in_flags = round_flags[ele_in_flags_i++];

			vector<GRBVar> ele_out_vars(statesize);

			vector<pair<BooleanPolynomial, GRBVar> > p_map;

			build_onemap_model(model, ele_in_vars, ele_out_vars, ele_map, ele_in_flags, p_map);

			p_maps.emplace_back(p_map);

			ele_in_vars = ele_out_vars;
		}

		out_vars = ele_in_vars;
	}

	virtual void build_model(GRBModel & model, int start, int end, const vector<vector<vector<Flag>>> & rounds_flags, const vector<GRBVar>& start_vars,
		vector<GRBVar> & end_vars, vector<vector<vector<pair<BooleanPolynomial, GRBVar>> > >& rounds_p_maps,
		const vector<update_element> & updatefuncs, bool add_obj)
	{
		if (start_vars.size() != statesize)
		{
			cerr << __func__ << ": The number of input vars is invalid.";
			exit(-1);
		}

		vector<GRBVar> round_in_vars(start_vars);

		for (int r = start; r < end; r++)
		{
			vector<vector<Flag>>   round_flags = rounds_flags[r];

			vector<GRBVar> round_out_vars(statesize);

			vector<vector<pair<BooleanPolynomial, GRBVar>> > p_maps;

			if (add_obj)
			{
				if (r == (start + end) / 2)
				{
					GRBLinExpr obj = 0;
					for (int i = 0; i < statesize; i++)
						if (round_flags[0][i] == "delta")
							obj += round_in_vars[i];
					model.setObjective(obj, GRB_MAXIMIZE);
				}
			}

			build_oneround_model(model, round_in_vars, round_out_vars, round_flags,  p_maps, updatefuncs);



			rounds_p_maps.emplace_back(p_maps);

			round_in_vars = round_out_vars;
		}

		end_vars = round_in_vars;
	}

	virtual void build_model_fast(GRBModel& model, int start, int end, const vector<vector<vector<Flag>>>& rounds_flags, const vector<GRBVar>& start_vars,
		vector<GRBVar>& end_vars, vector<vector<vector<pair<BooleanPolynomial, GRBVar>> > >& rounds_p_maps, bool add_obj)
	{
		if (start_vars.size() != statesize)
		{
			cerr << __func__ << ": The number of input vars is invalid.";
			exit(-1);
		}

		vector<GRBVar> round_in_vars(start_vars);

		for (int r = start; r < end; r++)
		{
			vector<vector<Flag>>   round_flags = rounds_flags[r];

			vector<GRBVar> round_out_vars(statesize);

			vector<vector<pair<BooleanPolynomial, GRBVar>> > p_maps;

			if (add_obj)
			{
				if (r == (start + end) / 2)
				{
					GRBLinExpr obj = 0;
					for (int i = 0; i < statesize; i++)
						if (round_flags[0][i] == "delta")
							obj += round_in_vars[i];
					model.setObjective(obj, GRB_MAXIMIZE);
				}
			}

			build_round_model_fast(model, round_in_vars, round_out_vars, round_flags, p_maps);



			rounds_p_maps.emplace_back(p_maps);

			round_in_vars = round_out_vars;
		}

		end_vars = round_in_vars;
	}

	/*
	virtual void callback_build_model(GRBModel& model, int start, int end, const vector<vector<vector<Flag>>>& rounds_flags, const vector<GRBVar>& start_vars,
		vector<GRBVar>& end_vars, vector<vector<map<BooleanPolynomial, GRBVar>> >& rounds_p_maps,
		const vector<update_element>& updatefuncs, bool add_obj, int midround, vector<GRBVar> & midvars,vector<Flag> & midflags)
	{
		if (start_vars.size() != statesize)
		{
			cerr << __func__ << ": The number of input vars is invalid.";
			exit(-1);
		}

		vector<GRBVar> round_in_vars(start_vars);

		for (int r = start; r < end; r++)
		{
			vector<vector<Flag>>   round_flags = rounds_flags[r];

			vector<GRBVar> round_out_vars(statesize);

			vector<map<BooleanPolynomial, GRBVar>> p_maps;

			if (add_obj)
			{
				if (r == midround)
				{
					GRBLinExpr obj = 0;
					for (int i = 0; i < statesize; i++)
						if (round_flags[0][i] == "delta")
							obj += round_in_vars[i];
					model.setObjective(obj, GRB_MAXIMIZE);
				}
			}

			if (r == midround)
			{
				midvars = round_in_vars;
				midflags = round_flags[0];
			}

			build_oneround_model(model, round_in_vars, round_out_vars, round_flags, p_maps, updatefuncs);



			rounds_p_maps.emplace_back(p_maps);

			round_in_vars = round_out_vars;
		}

		end_vars = round_in_vars;
	}
	*/



	virtual void callback_build_model_fast(GRBModel& model, int start, int end, const vector<vector<vector<Flag>>>& rounds_flags, const vector<GRBVar>& start_vars,
		vector<GRBVar>& end_vars)
	{
		if (start_vars.size() != statesize)
		{
			cerr << __func__ << ": The number of input vars is invalid.";
			exit(-1);
		}

		vector<GRBVar> round_in_vars(start_vars);

		for (int r = start; r < end; r++)
		{
			vector<vector<Flag>>   round_flags = rounds_flags[r];

			vector<GRBVar> round_out_vars(statesize);



			callback_build_round_model_fast(model, round_in_vars, round_out_vars, round_flags);




			round_in_vars = round_out_vars;
		}

		end_vars = round_in_vars;
	}
	

	virtual status solve_model(int start, int end, const vector<vector<vector<Flag>>>& rounds_flags, const dynamic_bitset<>& start_state,
		const dynamic_bitset<>& end_state, vector<dynamic_bitset<>>& p_sols,  vector<BooleanPolynomial> &  p_exps, 
		vector<pair<int,int>> & p_rms, set<int> & end_constants, int threads, double time, bool outputflag)
	{
		auto start_flags = rounds_flags[start][0];
		auto end_flags = rounds_flags[end][0];

		if (start_state.size() != statesize || end_state.size() != statesize || start_flags.size() != statesize || end_flags.size() != statesize)
		{
			cerr << __func__ << ": The number of input state, output state, input flags or output flags is invalid.";
			exit(-1);
		}

		for (int i = 0; i < statesize; i++)
		{
			if (start_state[i] == 1 && start_flags[i] != "delta")
			{
				cerr << __func__ << ": The input state doesn't match the input flags";
				exit(-1);
			}
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

		set_start_cstr(start, start_flags, start_state, model, start_vars);

		// build model round by round
		vector< vector<vector<pair<BooleanPolynomial, GRBVar>> >> rounds_p_maps;
		vector<GRBVar> end_vars(statesize);
		build_model_fast(model, start, end, rounds_flags, start_vars, end_vars, rounds_p_maps, true);
		// build_model(model, start, end, rounds_flags, start_vars, end_vars, rounds_p_maps, updatelist,  true);

		// impose constraints on final variables
		for (int i = 0; i < statesize; i++)
			if (end_state[i] == 1 && end_flags[i] == "zero_c")
			{
				return status::NOSOLUTION;
			}
			else if (end_state[i] == 1 && end_flags[i] == "one_c_k")
			{
				end_constants.emplace(i);
			}
		

		for (int i = 0; i < statesize; i++)
			if(end_flags[i] == "delta")
				model.addConstr(end_vars[i] == end_state[i]);


		
		vector<GRBVar> p_vars; 
		for (int r = 0; r < end-start; r++)
		{
			for (int i = 0; i < updatelist.size(); i++)
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

		// cout << "P variables: " << p_cnt << endl;

		if (time > 0)
			model.set(GRB_DoubleParam_TimeLimit, time);

		model.optimize();

		if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
		{

			logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(end) +
				string(" | EXPAND "));

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
			
			if (outputflag)
			{
				double time = model.get(GRB_DoubleAttr_Runtime);
				logger(__func__ + string(" : ") + to_string(start) + "-" + to_string(end) + string(" | "
				) + to_string(p_num) + string(" | ") + to_string(time) + string(" | "
				) + to_string(solCount));
			}

			if (solCount > 0)
			{
				map<dynamic_bitset<>, int> p_sols_counter;
				dynamic_bitset<> p_sol(p_num);

				for (int i = 0; i < solCount; i++)
				{
					model.set(GRB_IntParam_SolutionNumber, i);
					for (int j = 0; j < p_num; j++)
					{
						if (round(p_vars[j].get(GRB_DoubleAttr_Xn)) == 1)
							p_sol[j] = 1;
						else
							p_sol[j] = 0;
					}

					p_sols_counter[p_sol]++;
				}

				for (auto& p_sol_count : p_sols_counter)
					if (p_sol_count.second % 2)
						p_sols.emplace_back(p_sol_count.first);

				return status::SOLVED;
			}
			else
			{
				return status::NOSOLUTION;
			}
		}
	}



	

	void calculate_onemap_exps(const vector<BooleanPolynomial>& prev_exps, 
		vector<BooleanPolynomial>& next_exps, const vector<Flag> & next_flag,
		const update_element& ele_map) const
	{
		BooleanPolynomial default_exp(default_constant_1);

		next_exps = prev_exps;

		for (auto& fromto : ele_map)
		{
			if (next_flag[fromto.first] == "delta" || next_flag[fromto.first] == "zero_c")
				continue;

			auto from = fromto.second;
			BooleanPolynomial new_exp(default_constant_0);
			for (auto& mon : from)
			{
				BooleanPolynomial new_mon(default_exp);
				for (auto& i : mon.index())
					new_mon *= prev_exps[i];

				new_exp += new_mon;
			}

			next_exps[fromto.first] = new_exp;

			// cout << fromto.first << " : " << new_exp << " to " << fromto.second << endl;
		}
	}


	vector<BooleanPolynomial> normal_exp()
	{
		vector<BooleanPolynomial> normal_exp(statesize);
		for (int i = 0; i < statesize; i++)
			normal_exp[i] = BooleanPolynomial(statesize, "s" + to_string(i));
		return normal_exp;
	}

	vector<vector<vector<BooleanPolynomial>>> generate_exps(int start, int end, const vector<Flag>& start_flag, 
		const vector<BooleanPolynomial>& start_exps, const vector<update_element> & updatefuncs)
	{
		vector<vector<vector<BooleanPolynomial>>> rounds_exps;
		vector<BooleanPolynomial> cur_exps(start_exps);
		vector<Flag> cur_flag(start_flag);


		for (int r = start; r < end; r++)
		{
			vector<vector<BooleanPolynomial>> round_exps;
			for (auto& ele_map : updatefuncs)
			{
				for (int i = 0; i < statesize; i++)
					if (cur_flag[i] == "delta" || cur_flag[i] == "zero_c")
						cur_exps[i] = default_constant_0;

				vector<Flag> next_flag;
				calculate_onemap_flag(cur_flag, next_flag, ele_map);

				round_exps.emplace_back(cur_exps);
				vector<BooleanPolynomial> next_exps;
				calculate_onemap_exps(cur_exps, next_exps, next_flag, ele_map);
				cur_exps = next_exps;
				cur_flag = next_flag;
			}

			rounds_exps.emplace_back(round_exps);
		}

		for (int i = 0; i < statesize; i++)
			if (cur_flag[i] == "delta" || cur_flag[i] == "zero_c")
				cur_exps[i] = default_constant_0;
		vector<vector<BooleanPolynomial>> last_round_exps = { cur_exps };
		rounds_exps.emplace_back(last_round_exps);

		return rounds_exps;
	}

	




	
};


class cipher_grain128AEAD : public cipher
{
public:
	cipher_grain128AEAD()
	{
		ciphername = "grain128AEAD";
		statesize = 256; // b + s
		keysize = 128;
		ivsize = 96;

		default_constant_0 = BooleanPolynomial(statesize, "0");
		default_constant_1 = BooleanPolynomial(statesize, "1");

		define_updatelist();
		define_outputks();
		resolve_updatelist();
	}

private:
	BooleanPolynomial g()
	{
		string g_str = "s0+s26+s56+s91+s96+s3s67+s11s13+s17s18+s27s59+s40s48+s61s65+s68s84+s88s92s93s95+s22s24s25+s70s78s82";
		BooleanPolynomial g(statesize, g_str);
		return g;
	}

	BooleanPolynomial f()
	{
		string f_str = "s135+s166+s198+s209+s224";
		BooleanPolynomial f(statesize, f_str);
		return f;
	}

	BooleanPolynomial z()
	{
		string h_str = "s12s136+s141s148+s95s170+s188s207+s12s95s222";
		string z_str = h_str + "+" + "s221+s2+s15+s36+s45+s64+s73+s89";
		BooleanPolynomial z(statesize, z_str);
		return z;
	}

	void define_updatelist()
	{
		
		update_element ele0_map;
		ele0_map[0] = g();
		BooleanPolynomial s0(statesize, "s128");
		ele0_map[128] = z() + s0;

		updatelist.emplace_back(ele0_map); // after which b0 = g, s0 = s0 + z


		update_element ele1_map;
		BooleanPolynomial new_b0(statesize, "s128+s0");
		BooleanPolynomial new_s0(statesize, "s128");
		new_s0 += f();
		ele1_map[0] = new_b0;
		ele1_map[128] = new_s0;

		updatelist.emplace_back(ele1_map);

		update_element ele2_map;
		for (int i = 0; i < 128; i++)
		{
			int from = (i + 1) % 128;
			ele2_map[i] = BooleanPolynomial(statesize, "s" + to_string(from));
			from += 128;
			ele2_map[i+128] = BooleanPolynomial(statesize, "s" + to_string(from));
		}

		updatelist.emplace_back(ele2_map);
		
		
		
	};

	void define_outputks()
	{
		update_element ele0_map;
		ele0_map[0] = z();

		for (int i = 1; i < statesize; i++)
		{
			ele0_map[i] = default_constant_0;
		}
		outputks.emplace_back(ele0_map);

		// cout << "outputKS" << endl;
	}

public:
	dynamic_bitset<> generate_cube(const set<int>& cube_index)
	{
		int max_i = (*cube_index.rbegin());
		int min_i = (*cube_index.begin());
		if (max_i >= ivsize || min_i <= -1)
		{
			cout << "Max index : " << max_i << " Min index : " << min_i << " ivsize : " << ivsize << endl;
			cerr << __func__ << ": The cube indexes are out of range";
			exit(-1);
		}

		dynamic_bitset<> cube(ivsize);

		for (auto& i : cube_index)
			cube[i] = 1;

		return cube;

	}

	vector<BooleanPolynomial>  set_initial_exps(const dynamic_bitset<>& cube)
	{
		vector<BooleanPolynomial> initial_exps(statesize);
		for (auto& exp : initial_exps)
			exp = BooleanPolynomial(statesize, "0");

		for (int i = 0; i < 128; i++)
		{
			initial_exps[i] = BooleanPolynomial(statesize, "s" + to_string(i));
		}

		for(int i = 0;i < 96;i++)
			if (cube[i] == 1)
				initial_exps[i + 128] = BooleanPolynomial(statesize, "s" + to_string(i + 128));

		for (int i = 96; i < 127; i++)
			initial_exps[i+128] = BooleanPolynomial(statesize, "1");

		return initial_exps;
	}

	vector<Flag> set_initial_flag(const dynamic_bitset<>& cube)
	{
		if (cube.size() != ivsize)
		{
			cerr << __func__ << ": The cube size is not correct";
			exit(-1);
		}

		vector<Flag> initial_flag(statesize, string("zero_c"));
		for (int i = 0; i < keysize; i++)
			initial_flag[i].setflag("one_c_k");
		for (int i = 0; i < ivsize; i++)
			if (cube[i] == 1)
				initial_flag[128 + i].setflag("delta");

		for (int i = 96; i < 127; i++)
			initial_flag[i+128].setflag("one_c_1");

		return initial_flag;

	}

	dynamic_bitset<> set_initial_state(const dynamic_bitset<>& cube)
	{
		if (cube.size() != ivsize)
		{
			cerr << __func__ << ": The cube size is not correct";
			exit(-1);
		}

		dynamic_bitset<> initial_state(statesize);
		for (int i = 0; i < ivsize; i++)
		{
			initial_state[128 + i] = cube[i];
		}

		return initial_state;

	}

	double set_solver_time(int rs, int re)
	{
		double time;
		int gap = re - rs;
		if (gap > 120)
			time = 60;
		else if (gap > 110)
			time = 120;
		else if (gap > 100)
			time = 180;
		else if (gap > 33)
			time = 360;
		else if (gap > 10)
			time = 360;
		else
			time = 0;
		return time;
	}

	double set_expander_time(int rs, int re)
	{
		return 3000;
	}



};

class cipher_kreyvium : public cipher
{
public:
	cipher_kreyvium()
	{
		ciphername = "kreyvium";
		statesize = 288+256; // state + iv + key
		keysize = 128;
		ivsize = 128;

		default_constant_0 = BooleanPolynomial(statesize, "0");
		default_constant_1 = BooleanPolynomial(statesize, "1");

		define_updatelist();
		define_outputks();
		resolve_updatelist();
	}

private:
	void define_updatelist()
	{
		// non-linear part
		update_element ele0_map;
		BooleanPolynomial t1(statesize, "s65 + s90s91 + s92 + s170 + s288");
		BooleanPolynomial t2(statesize, "s161 + s174s175 + s176 + s263");
		BooleanPolynomial t3(statesize, "s242 + s285s286 + s287 + s68 + s416");

		ele0_map[93] = t1;
		ele0_map[177] = t2;
		ele0_map[0] = t3;

		// linera part
		// state
		for (int i = 0; i < 288; i++)
		{
			if (i == 92 || i == 176 || i == 287)
				continue;
			BooleanPolynomial expi(statesize, "s" + to_string(i));

			int nexti = (i + 1) % 288;
			ele0_map[nexti] = expi;
		}

		 // iv
		for (int i = 0; i < 128; i++)
		{
			int from_iv = 288 + ((i + 1) % 128);
			BooleanPolynomial expi(statesize, "s" + to_string(from_iv));
			ele0_map[288 + i] = expi;
		}

		// key
		for (int i = 0; i < 128; i++)
		{
			int from_key = 288 + 128 + ((i + 1) % 128);
			BooleanPolynomial expi(statesize, "s" + to_string(from_key));
			ele0_map[288 + 128 + i] = expi;
		}

		updatelist.emplace_back(ele0_map);
	}

	void define_outputks()
	{
		update_element ele0_map;
		BooleanPolynomial ks(statesize, "s65+s92+s161+s176+s242+s287");
		ele0_map[0] = ks;

		for (int i = 1; i < statesize; i++)
		{
			ele0_map[i] = default_constant_0;
		}
		outputks.emplace_back(ele0_map);

		// cout << "outputKS" << endl;
	}

public:

	dynamic_bitset<> generate_cube(const set<int>& cube_index)
	{
		int max_i = (*cube_index.rbegin());
		int min_i = (*cube_index.begin());
		if (max_i >= ivsize || min_i <= -1)
		{
			cout << "Max index : " << max_i << " Min index : " << min_i << " ivsize : " << ivsize << endl;
			cerr << __func__ << ": The cube indexes are out of range";
			exit(-1);
		}

		dynamic_bitset<> cube(ivsize);

		for (auto& i : cube_index)
			cube[i] = 1;

		return cube;

	}

	dynamic_bitset<> set_initial_state(const dynamic_bitset<>& cube)
	{
		if (cube.size() != ivsize)
		{
			cerr << __func__ << ": The cube size is not correct";
			exit(-1);
		}

		dynamic_bitset<> initial_state(statesize);
		for (int i = 0; i < ivsize; i++)
		{
			initial_state[93 + i] = cube[i];
			initial_state[288 + 127 - i] = cube[i];
		}

		return initial_state;

	}

	vector<Flag> set_initial_flag(const dynamic_bitset<>& cube)
	{
		if (cube.size() != ivsize)
		{
			cerr << __func__ << ": The cube size is not correct";
			exit(-1);
		}

		vector<Flag> initial_flag(statesize, string("zero_c"));
		for (int i = 0; i < 93; i++)
			initial_flag[i].setflag("one_c_k");
		for (int i = 0; i < keysize; i++)
			initial_flag[288 + 128 + i].setflag("one_c_k");

		for (int i = 0; i < ivsize; i++)
			if (cube[i] == 1)
			{
				initial_flag[93 + i].setflag("delta");
				initial_flag[288 + 127 - i].setflag("delta");
			}

		for (int i = 93+128; i < 287; i++)
			initial_flag[i].setflag("one_c_1");

		return initial_flag;

	}

	vector<BooleanPolynomial>  set_initial_exps(const dynamic_bitset<>& cube)
	{
		vector<BooleanPolynomial> initial_exps(statesize);
		for (auto& exp : initial_exps)
			exp = BooleanPolynomial(statesize, "0");

		// state
		for (int i = 0; i < 93; i++)
		{
			initial_exps[i] = BooleanPolynomial(statesize, "s" + to_string(i));
		}

		for(int i = 0;i < 128; i++)
			if(cube[i] == 1)
				initial_exps[i+93] = BooleanPolynomial(statesize, "s" + to_string(i+128));


		for (int i = 93+128; i < 287; i++)
			initial_exps[i] = BooleanPolynomial(statesize, "1");

		// iv
		for (int i = 0; i < 128; i++)
		{
			if(cube[i] == 1)
				initial_exps[288+127 - i] = BooleanPolynomial(statesize, "s" + to_string(128+i));
		}

		// key
		for (int i = 0; i < 128; i++)
		{
			initial_exps[288 + 128+127 - i] = BooleanPolynomial(statesize, "s" + to_string(i));
		}

		return initial_exps;
	}

	double set_solver_time(int rs, int re)
	{
		double time;
		int gap = re - rs;
		if (gap > 600)
			time = 40;
		else if (gap > 500)
			time = 128;
		else if (gap > 400)
			time = 160;
		else if (gap > 300)
			time = 320;
		else if (gap > 250)
			time = 640;
		else if (gap > 100)
			time = 1200;
		else if (gap > 20)
			time = 3600;
		else
			time = 0;
		return time;
	}

	double set_expander_time(int rs, int re)
	{
		return 14400;
	}

	void set_start_cstr(int start, const vector<Flag>& start_flags, const dynamic_bitset<>& start_state, GRBModel& model, vector<GRBVar>& start_vars)
	{
		if (start > 0)
		{
			start_vars.resize(statesize);

			for (int i = 0; i < statesize; i++)
				if (start_flags[i] == "delta")
				{
					start_vars[i] = model.addVar(0, 1, 0, GRB_BINARY);
					model.addConstr(start_vars[i] == start_state[i]);
				}
		}
		else
		{
			// in this case start_state = initial_state
			start_vars.resize(statesize);
			for(int i = 0;i < 128; i++)
				if (start_state[i + 93])
				{
					start_vars[i+93] = model.addVar(0, 1, 0, GRB_BINARY);
					start_vars[288+127 - i] = model.addVar(0, 1, 0, GRB_BINARY);
					model.addConstr(start_vars[i + 93] + start_vars[288 + 127 - i] >= 1);
				}
				
		}
	}

	void set_callback_start_cstr(int start, const vector<Flag>& start_flags, const dynamic_bitset<>& start_state, GRBModel& model, vector<GRBVar>& start_vars)
	{
		if (start > 0)
		{
			start_vars.resize(statesize);

			for (int i = 0; i < statesize; i++)
				if (start_flags[i] == "delta")
				{
					start_vars[i] = model.addVar(0, 1, 0, GRB_BINARY);
					model.addConstr(start_vars[i] == start_state[i]);
				}
		}
		else
		{
			// in this case start_state = initial_state
			start_vars.resize(statesize);
			for (int i = 0; i < 128; i++)
				if (start_state[i + 93])
				{
					start_vars[i + 93] = model.addVar(0, 1, 0, GRB_BINARY);
					start_vars[288 + 127 - i] = model.addVar(0, 1, 0, GRB_BINARY);
					model.addConstr(start_vars[i + 93] + start_vars[288 + 127 - i] == 1);
				}

		}
	}
};











// this class should have only one instance (singleton pattern)
class cipher_trivium : public cipher
{
public:
	cipher_trivium()
	{
		ciphername = "trivium";
		statesize = 288;
		keysize = 80;
		ivsize = 80;

		default_constant_0 = BooleanPolynomial(statesize, "0");
		default_constant_1 = BooleanPolynomial(statesize, "1");

		define_updatelist();
		define_outputks();
		resolve_updatelist();
	}

private:
	void define_updatelist()
	{
		
		// non-linear part
		update_element ele0_map;
		BooleanPolynomial t1(statesize, "s65 + s90s91 + s92 + s170");
		BooleanPolynomial t2(statesize, "s161 + s174s175 + s176 + s263");
		BooleanPolynomial t3(statesize, "s242 + s285s286 + s287 + s68");

		ele0_map[93] = t1;
		ele0_map[177] = t2;
		ele0_map[0] = t3;



		// linera part
		for (int i = 0; i < 288; i++)
		{
			if (i == 92 || i == 176 || i == 287)
				continue;
			BooleanPolynomial expi(statesize, "s" + to_string(i));

			int nexti = (i + 1) % 288;
			ele0_map[nexti] = expi;
		}
		updatelist.emplace_back(ele0_map);
		
	}

	void define_outputks()
	{
		update_element ele0_map;
		BooleanPolynomial ks(statesize, "s65+s92+s161+s176+s242+s287");
		ele0_map[0] = ks;

		for (int i = 1; i < statesize; i++)
		{
			ele0_map[i] = default_constant_0;
		}
		outputks.emplace_back(ele0_map);

		// cout << "outputKS" << endl;
	}

public:

	dynamic_bitset<> generate_cube(const set<int>& cube_index)
	{
		int max_i = (*cube_index.rbegin());
		int min_i = (*cube_index.begin());
		if (max_i >= ivsize || min_i <= -1)
		{
			cout << "Max index : " << max_i << " Min index : " << min_i << " ivsize : " << ivsize << endl;
			cerr << __func__ << ": The cube indexes are out of range";
			exit(-1);
		}

		dynamic_bitset<> cube(ivsize);

		for (auto& i : cube_index)
			cube[i] = 1;

		return cube;
	
	}

	dynamic_bitset<> set_initial_state(const dynamic_bitset<>& cube)
	{
		if (cube.size() != ivsize)
		{
			cerr << __func__ << ": The cube size is not correct";
			exit(-1);
		}

		dynamic_bitset<> initial_state(statesize);
		for (int i = 0; i < ivsize; i++)
		{
			initial_state[93 + i] = cube[i];
		}

		return initial_state;

	}

	vector<Flag> set_initial_flag(const dynamic_bitset<>& cube)
	{
		if (cube.size() != ivsize)
		{
			cerr << __func__ << ": The cube size is not correct";
			exit(-1);
		}

		vector<Flag> initial_flag(statesize, string("zero_c"));
		for (int i = 0; i < keysize; i++)
			initial_flag[i].setflag("one_c_k");
		for (int i = 0; i < ivsize; i++)
			if (cube[i] == 1)
				initial_flag[93 + i].setflag("delta");

		for (int i = 285; i < 288; i++)
			initial_flag[i].setflag("one_c_1");

		return initial_flag;

	}

	vector<BooleanPolynomial>  set_initial_exps(const dynamic_bitset<>& cube)
	{
		vector<BooleanPolynomial> initial_exps(statesize);
		for (auto& exp : initial_exps)
			exp = BooleanPolynomial(statesize, "0");

		for (int i = 0; i < 80; i++)
		{
			initial_exps[i] = BooleanPolynomial(statesize, "s" + to_string(i));
			if(cube[i] == 1)
				initial_exps[i + 93] = BooleanPolynomial(statesize, "s" + to_string(i + 93));
		}

		for (int i = 285; i < 288; i++)
			initial_exps[i] = BooleanPolynomial(statesize, "1");

		return initial_exps;
	}

	double set_solver_time(int rs, int re)
	{
		double time;
		int gap = re - rs;
		if (gap > 600)
			time = 40;
		else if (gap > 500)
			time = 80;
		else if (gap > 400)
			time = 160;
		else if (gap > 300)
			time = 320;
		else if (gap > 250)
			time = 640;
		else if (gap > 100)
			time = 1200;
		else if (gap > 20)
			time = 3600;
		else
			time = 0;
		return time;
	}

	double set_expander_time(int rs, int re)
	{
		return 500;
	}

	




};

#endif
