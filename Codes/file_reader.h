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
#include"listsOfPolynomials.h"

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

	bool read_one_pair(istream& is, node_pair& state_pair, BooleanPolynomial& coef, ListsOfPolynomialsAsFactors & coef_lists) const
	{
		string oneline;
		while (getline(is, oneline))
		{
			if (!oneline.empty())
				break;
		}

		if (!is)
			return false;

		state_pair.first = dynamic_bitset<>(oneline);
		getline(is, oneline);
		state_pair.second = dynamic_bitset<>(oneline);
		getline(is, oneline);
		coef = BooleanPolynomial(state_pair.first.size(), oneline);
		coef_lists.load(state_pair.first.size(), is);

		return true;
	}


	bool read_lists_once(istream& is, ListsOfPolynomialsAsFactors& coef_lists, int varsNum) const
	{
		return coef_lists.load(varsNum, is);
	}






};

#endif