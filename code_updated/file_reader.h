#pragma once
#ifndef __FILE_READER_H
#define __FILE_READER_H


#include<map>
#include<vector>
#include<filesystem>
#include<string>
#include"BooleanPolynomial.h"
#include"dynamic_bitset.hpp"
#include"thread_pool.h"

using namespace std;
using namespace thread_pool;

const regex sol_filename_regex(R"~((\d+)_(\d+)_\d+.txt)~");
const regex P_filename_regex(R"~((\d+)_(\d+)_P.txt)~");
const regex sol_start_state_regex(R"~(state from:([01]+))~");
const regex sol_end_state_regex(R"~(state to:([01]+))~");
const regex sol_p_regex(R"~(p(\d+)=(\d+)-(\d+)-([a-zA-Z\d\+]+))~");
const regex integer_regex(R"~((\d+))~");

void repeat_match(string& instr, vector<vector<string>>& smatches, const regex pattern)
{
	string& subject = instr;
	smatch sm;
	while (regex_search(subject, sm, pattern))
	{
		vector<string> ssmatches;
		for (auto& ssmatch : sm)
			ssmatches.emplace_back(ssmatch.str());
		smatches.emplace_back(ssmatches);
		subject = sm.suffix().str();
	}
}

class FileReader
{
public:

	void getJustCurrentFilePaths(string path, vector<filesystem::path>& files) const
	{
		for (const auto& entry : filesystem::directory_iterator(path))
			files.emplace_back(entry.path());
	}

	bool read_one_pair(fstream& fs, node_pair& state_pair, BooleanPolynomial& coef) const
	{
		string oneline;
		if (getline(fs, oneline))
		{
			state_pair.first = dynamic_bitset<>(oneline);
			getline(fs, oneline);
			state_pair.second = dynamic_bitset<>(oneline);
			getline(fs, oneline);
			coef = BooleanPolynomial(state_pair.first.size(), oneline);
			getline(fs, oneline);
			return true;
		}
		else
			return false;
	}


	bool read_one_sol(fstream& fs, dynamic_bitset<>& start_state, dynamic_bitset<>& end_state,
		vector<dynamic_bitset<>>& p_sols, vector<BooleanPolynomial>& p_exps,
		vector<pair<int, int>>& p_rms, set<int>& end_constants) const
	{

		string oneline;

		if (getline(fs, oneline))
		{
			smatch start_state_sm;
			regex_match(oneline, start_state_sm, sol_start_state_regex);
			auto start_state_str = start_state_sm[1].str();
			start_state = dynamic_bitset<>(start_state_str);

			getline(fs, oneline);
			smatch end_state_sm;
			regex_match(oneline, end_state_sm, sol_end_state_regex);
			auto end_state_str = end_state_sm[1].str();
			end_state = dynamic_bitset<>(end_state_str);

			int statesize = end_state.size();

			getline(fs, oneline);
			smatch p_sol_sm;
			vector<vector<string>> p_sol_sms;
			repeat_match(oneline, p_sol_sms, sol_p_regex);

			while (getline(fs, oneline))
			{
				if (oneline[0] != 'e')
				{
					dynamic_bitset<> p_sol(oneline);
					p_sols.emplace_back(p_sol);
				}
				else
					break;
			}

			if (p_sols.size() == 0)
				return true;

			int p_num = p_sols[0].size();

			p_exps.resize(p_num);
			p_rms.resize(p_num);

			for (auto& sm : p_sol_sms)
			{
				int p_index = stoi(sm[1]);
				int p_round = stoi(sm[2]);
				int p_mapi = stoi(sm[3]);
				BooleanPolynomial p_exp(statesize, sm[4]);

				p_exps[p_index] = p_exp;
				p_rms[p_index] = pair(p_round, p_mapi);
			}

			vector<vector<string>> int_sms;
			repeat_match(oneline, int_sms, integer_regex);
			for (auto& sm : int_sms)
			{
				end_constants.emplace(stoi(sm[1]));
			}

			getline(fs, oneline);
			return true;
		}
		else
			return false;

	}





};

#endif
