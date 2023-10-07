#pragma once
#ifndef __LISTSOFPOLYNOMIALS_H_
#define __LISTSOFPOLYNOMIALS_H_
#include <iostream>
#include <vector>
#include <set>
#include "BooleanPolynomial.h"
#include "log.h"

using namespace std;

class ListOfPolynomialsAsFactors
{
private:
	set<BooleanPolynomial> listOfPolynomials;
	int varsNum;

public:
	ListOfPolynomialsAsFactors(int aVarsNum = 0) : varsNum(aVarsNum) {};

	ListOfPolynomialsAsFactors(const BooleanPolynomial& aPoly)
	{
		if (aPoly.isnull())
		{
			logger(__func__ + string(": try to add a null polynomial."));
			exit(-1);
		}

		varsNum = aPoly.size();
		if (!aPoly.isone())
			listOfPolynomials.emplace(aPoly);
	}

	void add(const BooleanPolynomial& aPoly)
	{
		if (aPoly.isnull())
		{
			logger(__func__ + string(": try to add a null polynomial."));
			exit(-1);
		}


		if (listOfPolynomials.size() == 0 && varsNum == 0)
		{
			varsNum = aPoly.size();
		}
		else
		{
			if (varsNum != aPoly.size())
			{
				logger(__func__ + string(": mismatch in the number of variables."));
				exit(-1);
			}
		}

		if (!aPoly.isone())
			listOfPolynomials.emplace(aPoly);

	}

	BooleanPolynomial getProd() const
	{
		BooleanPolynomial prod = fastmul(varsNum, listOfPolynomials);
		return prod;
	}

	int eval(const map<int, int>& vals) const
	{
		int res = 1;
		for (auto& poly : listOfPolynomials)
		{
			res *= poly.eval(vals);
			if (res == 0)
				return 0;
		}

		return res;
	}


	

	int maxDeg() const
	{
		int max_deg = 0;
		for (auto& poly : listOfPolynomials)
			max_deg += poly.deg();
		
		return max_deg;
	}

	int maxNrTerms() const
	{
		int max_num_terms = 1;
		for (auto& poly : listOfPolynomials)
			max_num_terms *= poly.moncnt();
		return max_num_terms;
	}

	int getVarsNum() const { return varsNum; }

	void merge(ListOfPolynomialsAsFactors& otherList)
	{
		if (varsNum != otherList.varsNum)
		{
			logger(__func__ + string(": try to add a null polynomial."));
			exit(-1);
		}

		listOfPolynomials.merge(otherList.listOfPolynomials);
	}

	void merge(const ListOfPolynomialsAsFactors& otherList)
	{
		if (varsNum != otherList.varsNum)
		{
			logger(__func__ + string(": try to add a null polynomial."));
			exit(-1);
		}

		listOfPolynomials.insert(otherList.listOfPolynomials.begin(), otherList.listOfPolynomials.end());
	}

	friend ListOfPolynomialsAsFactors operator * (const ListOfPolynomialsAsFactors& l0, const ListOfPolynomialsAsFactors& l1)
	{
		if (l0.varsNum != l1.varsNum)
		{
			logger(__func__ + string(": mismatch in the number of variables."));
			exit(-1);
		}

		ListOfPolynomialsAsFactors newList(l0.varsNum);
		set_union(l0.listOfPolynomials.begin(), l0.listOfPolynomials.end(), l1.listOfPolynomials.begin(), l1.listOfPolynomials.end(), insert_iterator(newList.listOfPolynomials, newList.listOfPolynomials.end()));

		return newList;
	}

	friend bool operator< (const ListOfPolynomialsAsFactors& l0, const ListOfPolynomialsAsFactors& l1)
	{
		if (l0.varsNum < l1.varsNum)
			return true;
		else if (l0.varsNum > l1.varsNum)
			return false;

		auto it0 = l0.listOfPolynomials.begin();
		auto it1 = l1.listOfPolynomials.begin();

		while (it0 != l0.listOfPolynomials.end() && it1 != l1.listOfPolynomials.end())
		{
			if (*it0 < *it1)
				return true;
			else if (*it0 > * it1)
				return false;

			it0++;
			it1++;
		}

		if (it0 == l0.listOfPolynomials.end() && it1 != l1.listOfPolynomials.end())
			return true;

		return false;
	}

	void display(ostream& fout) const
	{
		if (listOfPolynomials.size() == 0)
		{
			fout << "(1)" << endl;
			return;
		}

		for (auto& poly : listOfPolynomials)
			fout << "(" << poly << ")";

		fout << endl;
	}

	void load(int aVarsNum, const string& str)
	{
		varsNum = aVarsNum;
		listOfPolynomials.clear();

		if (str == "(1)")
			listOfPolynomials.clear();
		else
		{
			string polystr;
			for (int i = 0; i < str.size(); i++)
			{
				if (str[i] == '(')
					polystr = "";
				else if (str[i] == ')')
				{
					BooleanPolynomial poly(aVarsNum, polystr);
					listOfPolynomials.emplace(poly);
				}
				else
				{
					polystr += str[i];
				}
			}
		}

	}




	auto begin() const
	{
		return listOfPolynomials.begin();
	}

	auto end() const
	{
		return listOfPolynomials.end();
	}
};

class ListsOfPolynomialsAsFactors
{
private:
	vector<ListOfPolynomialsAsFactors> listsOfPolynomials;
	int varsNum;
public:
	ListsOfPolynomialsAsFactors(int aVarsNum = 0) : varsNum(aVarsNum) {};

	ListsOfPolynomialsAsFactors(const ListOfPolynomialsAsFactors& aList)
	{
		varsNum = aList.getVarsNum();
		listsOfPolynomials.emplace_back(aList);
	}

	void add(const ListOfPolynomialsAsFactors& list)
	{
		int listVarsNum = list.getVarsNum();

		if (listsOfPolynomials.size() == 0 && varsNum == 0)
			varsNum = listVarsNum;
		else
		{
			if (listVarsNum != varsNum)
			{
				logger(__func__ + string(": mismatch in the number of variables."));
				exit(-1);
			}
		}

		listsOfPolynomials.emplace_back(list);
	}

	/**
	 * @brief This function is used to remove those lists of polynomials that occur an even-number of times
	*/
	void filterLists()
	{
		map<ListOfPolynomialsAsFactors, int> listCounter;
		for (auto& list : listsOfPolynomials)
			listCounter[list]++;

		listsOfPolynomials.clear();
		for (auto& list_cnt : listCounter)
			if (list_cnt.second % 2)
				listsOfPolynomials.emplace_back(list_cnt.first);
	}

	int eval(const map<int, int>& vals) const
	{
		int res = 0;
		for (auto& list : listsOfPolynomials)
		{
			res ^= list.eval(vals);
		}
		return res;
	}

	BooleanPolynomial getSum() const
	{
		vector<BooleanPolynomial> prodsOfLists;
		for (auto& list : listsOfPolynomials)
			prodsOfLists.emplace_back(list.getProd());

		return fastsum(varsNum, prodsOfLists);
	}

	int maxDeg() const
	{
		int max_deg = 0;
		for (auto& list : listsOfPolynomials)
		{
			int tmp_deg = list.maxDeg();
			if (tmp_deg > max_deg)
				max_deg = tmp_deg;
		}

		return max_deg;
	}

	int maxNrTerms() const
	{
		int max_num_terms = 0;
		for (auto& list : listsOfPolynomials)
			max_num_terms += list.maxNrTerms();
		return max_num_terms;
	}

	void display(ostream& fout) const
	{
		for (auto& list : listsOfPolynomials)
			list.display(fout);

		fout << endl;
	}

	bool load(int aVarsNum, istream& is)
	{
		varsNum = aVarsNum;
		listsOfPolynomials.clear();

		string oneline;

		// omit the first empty lines
		while (getline(is, oneline))
		{
			if (!oneline.empty())
				break;
		}

		if (!is)
		{
			return false;
		}

		ListOfPolynomialsAsFactors list;
		list.load(aVarsNum, oneline);
		listsOfPolynomials.emplace_back(list);

		while (getline(is, oneline))
		{
			if (oneline.empty())
				break;

			list.load(aVarsNum, oneline);
			listsOfPolynomials.emplace_back(list);
		}

		return true;
	}


	friend ostream& operator << (ostream& os, const ListsOfPolynomialsAsFactors& ls)
	{
		ls.display(os);
		return os;
	}

	friend ListsOfPolynomialsAsFactors operator * (const ListsOfPolynomialsAsFactors& ls0, const ListsOfPolynomialsAsFactors& ls1)
	{
		if (ls0.varsNum != ls1.varsNum)
		{
			logger(__func__ + string(": mismatch in the number of variables."));
			exit(-1);
		}


		ListsOfPolynomialsAsFactors newLists(ls0.varsNum);
		for (auto& l0 : ls0.listsOfPolynomials)
			for (auto& l1 : ls1.listsOfPolynomials)
			{
				ListOfPolynomialsAsFactors newList = l0 * l1;
				newLists.add(newList);
			}

		return newLists;
	}

	void operator *=(const ListOfPolynomialsAsFactors& otherList)
	{
		if (varsNum != otherList.getVarsNum())
		{
			logger(__func__ + string(": mismatch in the number of variables."));
			exit(-1);
		}

		for (auto& list : listsOfPolynomials)
		{
			list.merge(otherList);
		}
	}

	void operator += (const ListsOfPolynomialsAsFactors& otherLists)
	{
		if (varsNum != otherLists.varsNum)
		{
			logger(__func__ + string(": mismatch in the number of variables."));
			exit(-1);
		}

		listsOfPolynomials.insert(listsOfPolynomials.end(), otherLists.listsOfPolynomials.begin(), otherLists.listsOfPolynomials.end());
	}

	auto begin() const
	{
		return listsOfPolynomials.begin();
	}

	auto end() const
	{
		return listsOfPolynomials.end();
	}



};
#endif