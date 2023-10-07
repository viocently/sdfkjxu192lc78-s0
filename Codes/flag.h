#pragma once
#ifndef __FLAG_H_
#define __FLAG_H_

#include <iostream>
#include <string>
#include <map>



using namespace std;


enum class FLAG { ZERO_C, ONE_C_1, ONE_C_K, DELTA, UNDEFINED };
map<string, FLAG> str_FLAG_map = {
	{"zero_c" , FLAG::ZERO_C},
	{"one_c_k", FLAG::ONE_C_K},
	{"one_c_1", FLAG::ONE_C_1},
	{"delta", FLAG::DELTA}
};

map<FLAG, int> add_priority = {
	{FLAG::ZERO_C, 0},
	{FLAG::ONE_C_1, 1},
	{FLAG::ONE_C_K, 2},
	{FLAG::DELTA, 3},
	{FLAG::UNDEFINED, 4}
};

map<FLAG, int> mul_priority = {
	{FLAG::ONE_C_1, 0},
	{FLAG::ONE_C_K, 1},
	{FLAG::DELTA, 2},
	{FLAG::ZERO_C, 3},
	{FLAG::UNDEFINED, 4}
};

map<FLAG, string> FLAG_str_map = {
	{FLAG::ZERO_C, "zero_c"},
	{FLAG::ONE_C_K, "one_c_k"},
	{FLAG::ONE_C_1, "one_c_1"},
	{FLAG::DELTA, "delta"},
	{FLAG::UNDEFINED, "undefined"}
};



class Flag
{
private:
	FLAG flag;
public:
	Flag() : flag(FLAG::UNDEFINED) {};
	Flag(const string s)
	{
		if (str_FLAG_map.find(s) != str_FLAG_map.end())
			flag = str_FLAG_map[s];
		else
			flag = FLAG::UNDEFINED;
	}

	void operator= (const string s)
	{
		if (str_FLAG_map.find(s) != str_FLAG_map.end())
			flag = str_FLAG_map[s];
		else
			flag = FLAG::UNDEFINED;
	}

	void operator *= (const Flag& rhs)
	{
		auto lp = mul_priority[flag];
		auto rp = mul_priority[rhs.flag];
		if (rp > lp)
			flag = rhs.flag;
	}

	void operator += (const Flag& rhs)
	{
		
		if (flag == FLAG::ONE_C_1 &&
			rhs.flag == FLAG::ONE_C_1)
		{
			flag = FLAG::ZERO_C;
			return;
		}
		
		


		auto lp = add_priority[flag];
		auto rp = add_priority[rhs.flag];

		if (rp > lp)
			flag = rhs.flag;
	}

	void setflag(const string s)
	{
		if (str_FLAG_map.find(s) != str_FLAG_map.end())
			flag = str_FLAG_map[s];
		else
			flag = FLAG::UNDEFINED;
	}

	string str() const
	{
		return FLAG_str_map[flag];
	}

	friend ostream& operator << (ostream& os, const Flag& flag)
	{
		os << flag.str();
		return os;
	}

	friend bool operator == (const Flag & lhs, const Flag & rhs) 
	{
		return lhs.flag == rhs.flag;
	}

	friend bool operator != (const Flag& lhs, const Flag& rhs)
	{
		return lhs.flag != rhs.flag;
	}

	friend bool operator < (const Flag& lhs, const Flag& rhs)
	{
		return lhs.flag < rhs.flag;
	}

	friend bool operator > (const Flag& lhs, const Flag& rhs)
	{
		return lhs.flag > rhs.flag;
	}

	friend Flag operator + (const Flag & lhs, const Flag & rhs) 
	{
		
		if (lhs.flag == FLAG::ONE_C_1 &&
			rhs.flag == FLAG::ONE_C_1)
			return Flag("zero_c");
		
		

		auto lp = add_priority[lhs.flag];
		auto rp = add_priority[rhs.flag];

		return (lp < rp) ? rhs : lhs;
	}

	friend Flag operator * (const Flag& lhs, const Flag& rhs)
	{

		auto lp = mul_priority[lhs.flag];
		auto rp = mul_priority[rhs.flag];

		return (lp < rp) ? rhs : lhs;
	}

	friend bool operator ==(const Flag& lhs, const string rhs)
	{
		return (FLAG_str_map[lhs.flag]) == rhs;
	}

	friend bool operator !=(const Flag& lhs, const string rhs)
	{
		return (FLAG_str_map[lhs.flag]) != rhs;
	}
};

/**
 * @brief This class is used to generate the basic tikz code for identifying the delta, 1_c and 0_c bits
*/
class FlagDisplay
{
private:
	vector<Flag> stateFlags;
	int heightUnit;
	int widthUnit;
	int heightBasis;
	int widthBasis;

public:

	FlagDisplay(const vector<Flag>& aStateFlags, int aHeightUnit, int aWidthUnit, int aHeightBasis, int aWidhtBasis) :
		stateFlags(aStateFlags), heightUnit(aHeightUnit), widthUnit(aWidthUnit),
		heightBasis(aHeightBasis), widthBasis(aWidhtBasis)
	{}

	string coordinate(double x, double y)
	{
		return "(" + to_string(x) + "," + to_string(y) + ")";
	}

	string getColor(const Flag& flag)
	{
		if (flag == "zero_c")
			return "gray";
		else if (flag == "delta")
			return "green!40";
		else if (flag == "one_c_1" || flag == "one_c_k")
			return "red!40";
		else
			return "black!40";
	}

	void generateTikzCodes(ostream & os, int rounds)
	{
		// first draw the rectangle as the state
		int stateSize = stateFlags.size();

		int curX = widthBasis;
		int curFlagIndex = 0;
		string curColor = getColor(stateFlags[curFlagIndex]);
		
		while (curFlagIndex < stateSize)
		{
			int nextFlagIndex = curFlagIndex + 1;
			while (nextFlagIndex < stateSize && stateFlags[nextFlagIndex] == stateFlags[curFlagIndex])
				nextFlagIndex++;
			int nextX = nextFlagIndex * widthUnit + widthBasis;

			os << R"(\draw [dashed, fill=)" << curColor << "] " << coordinate(curX, heightBasis) << " rectangle " << coordinate(nextX, heightBasis + heightUnit)
				<< ";" << endl;

			if (nextFlagIndex >= stateSize)
				break;

			curX = nextX;
			curFlagIndex = nextFlagIndex;
			curColor = getColor(stateFlags[curFlagIndex]);
		}

		os << R"(\draw )" + coordinate(widthBasis, heightBasis) + " rectangle " + coordinate(widthUnit * stateSize + widthBasis, heightBasis + heightUnit) + ";" << endl;

		os << R"(\node[scale=15] at )" << coordinate(-15, heightBasis + double(heightUnit) / 2) << R"({Round $)" << rounds << R"($};)" << endl;

	}
};
#endif