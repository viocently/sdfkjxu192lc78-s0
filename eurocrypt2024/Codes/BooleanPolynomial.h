#pragma once
#ifndef __BOOLEANPOLYNOMIAL_H_
#define __BOOLEANPOLYNOMIAL_H_
#include<iostream>
#include<map>
#include<vector>
#include<set>
#include<string>
#include<regex>
#include<cctype>
#include "boost/dynamic_bitset.hpp"
#include "comparator.h"

using namespace std;
using namespace boost;

// A personal implementation of Boolean Monomial and Polynomial
// Not so professional  : )


static string remove_blank(const string instr)
{
    string subject = instr;
    subject.erase(remove_if(subject.begin(), subject.end(), ::isspace), subject.end());
    return subject;
}

static bool find_repeat_patterns(string & instr, vector<vector<string>>& smatches, const regex pattern)
{
    string & subject = instr;
    smatch sm;
    while (regex_search(subject, sm, pattern, regex_constants::match_continuous))
    {
        vector<string> ssmatches;
        for (auto& ssmatch : sm)
            ssmatches.emplace_back(ssmatch.str());
        smatches.emplace_back(ssmatches);
        subject = sm.suffix().str();
    }

    if (!subject.empty())
    {
        return false;
    }
    else
        return true;
}


namespace BooleanResolution
{
    enum class SYMBOL_TYPE {ALPHA, ZERO, ONE, UNDEFINED};
    using SYMBOL_VALUE = char;

    class symbol
    {
    public:
        SYMBOL_VALUE val;
        SYMBOL_TYPE type;

        symbol() : type(SYMBOL_TYPE::UNDEFINED), val('?') {};
        symbol(const char v)
        {
            switch (v)
            {
            case '0':
                type = SYMBOL_TYPE::ZERO;
                val = '0';
                break;
            case '1':
                type = SYMBOL_TYPE::ONE;
                val = '1';
                break;
            default:
                if (isalpha(v))
                {
                    type = SYMBOL_TYPE::ALPHA;
                    val = v;
                    break;
                }
                else
                {
                    type = SYMBOL_TYPE::UNDEFINED;
                    val = '?';
                    break;
                }

            }
        }
        symbol(const symbol& s) : type(s.type), val(s.val) {};

        void setval(const char v)
        {
            switch (v)
            {
            case '0':
                type = SYMBOL_TYPE::ZERO;
                val = '0';
                break;
            case '1':
                type = SYMBOL_TYPE::ONE;
                val = '1';
                break;
            default:
                if (isalpha(v))
                {
                    type = SYMBOL_TYPE::ALPHA;
                    val = v;
                    break;
                }
                else
                {
                    type = SYMBOL_TYPE::UNDEFINED;
                    val = '?';
                    break;
                }

            }
        }
    };

    const regex MonomialRegex(R"~(([a-zA-Z])(\d+))~");

    using MonomialInf = pair<symbol, set<int>>;

    MonomialInf resolve_mon_exp(const string monexp)
    {
        MonomialInf moninf;
        string exp = remove_blank(monexp);

        if (exp == "0")
        {
            moninf.first.setval('0');
            return moninf;
        }
        else if (exp == "1")
        {
            moninf.first.setval('1');
            return moninf;
        }
        else
        {
            vector<vector<string>> smatches;
            auto isvalid = find_repeat_patterns(exp, smatches, MonomialRegex);


            symbol s;
            if (isvalid)
            {
                set<int> vi;
                for (auto& sm : smatches)
                {
                    if (s.type == SYMBOL_TYPE::UNDEFINED)
                    {
                        vi.emplace(stoi(sm[2]));
                        s.setval(sm[1][0]);
                    }
                    else
                        if (s.val == sm[1][0])
                        {
                            vi.emplace(stoi(sm[2]));
                            continue;
                        }
                        else
                            return moninf;
                }

                moninf.first = s;
                moninf.second = vi;
            }

            return moninf;
        }
    }

    using PolyInf = vector<MonomialInf>;

    PolyInf resolve_poly_exp(const string polyexp)
    {
        string exp = remove_blank(polyexp);

        regex plus_regex(R"~(\+)~");
        sregex_token_iterator rend;

        sregex_token_iterator rbegin(exp.begin(), exp.end(), plus_regex, -1);

        PolyInf moninfs;

        for (auto i = rbegin; i != rend; i++)
        {
            string monexp = i->str();


            if (monexp.empty())
            {
                moninfs.clear();
                return moninfs;
            }

            
            auto moninf = resolve_mon_exp(monexp);
            auto& monsymbol = moninf.first;
            if (monsymbol.type == SYMBOL_TYPE::UNDEFINED)
            {
                moninfs.clear();
                return moninfs;
            }

            moninfs.emplace_back(moninf);
        }

        return moninfs;

    }



};



class BooleanPolynomial;

class BooleanMonomial
{
private:
    friend class BooleanPolynomial;
    
    int coefficient = 1;

    dynamic_bitset<> monrep;



public:
    BooleanMonomial() {};

    BooleanMonomial(const BooleanMonomial& _mon) : monrep(_mon.monrep), coefficient(_mon.coefficient) {};

    BooleanMonomial(int ini_size, int c = 1) { monrep.resize(ini_size);  coefficient = c; }

private:
    BooleanMonomial(const dynamic_bitset<> rep) : monrep(rep) {}

    BooleanMonomial(int ini_size, const set<int>& true_index)
    {
        monrep.resize(ini_size);
        for (auto& i : true_index)
        {
            if(i < ini_size)
                monrep[i] = 1;
            else
            {
                monrep.resize(0);
                break;
            }
        }
    }



    BooleanMonomial(int ini_size, const string monexp)
    {

        
        string exp = remove_blank(monexp);
        using BooleanResolution::MonomialInf;
        using BooleanResolution::resolve_mon_exp;
        using BooleanResolution::SYMBOL_TYPE;
        MonomialInf moninf = resolve_mon_exp(monexp);

        auto monsymbol = moninf.first;
        if (monsymbol.type == SYMBOL_TYPE::ALPHA)
        {
            monrep.resize(ini_size);
            for (auto& i : moninf.second)
            {
                if (i < ini_size)
                    monrep[i] = 1;
                else
                {
                    monrep.resize(0);
                    break;
                }
            }
        }
        else if (monsymbol.type == SYMBOL_TYPE::ONE)
        {
            monrep.resize(ini_size);
        }
        else if (monsymbol.type == SYMBOL_TYPE::ZERO)
        {
            monrep.resize(ini_size);
            coefficient = 0;
        }

        
        
    }

public:
    int size() const { return monrep.size(); } 

    bool isnull() const { return monrep.size() == 0; }

    bool iszero() const { return coefficient == 0 && monrep.size() != 0; }

    bool isone() const
    {
        return coefficient == 1 && monrep.none() == true;
    }

    void resize(int n)
    {
        monrep.resize(n);
    }

    dynamic_bitset<> rep() const
    {
        return monrep;
    }

    set<int> index() const
    {
        set<int> indexes;
        if (coefficient == 0)
            return indexes;

        for (int i = 0; i < monrep.size(); i++)
            if (monrep[i])
                indexes.emplace(i);

        return indexes;
    }

    int eval(const dynamic_bitset<>& evalmap) const
    {
        if (coefficient == 0)
            return 0;

        if (evalmap.size() != monrep.size())
        {
            cerr << __func__ << "The size of evalmap does not match the number of variables in this monomial." << endl;
            exit(-1);
        }


        dynamic_bitset<> evalval = evalmap & monrep;

        if (evalval == monrep)
            return 1;
        else
            return 0;
    }

    int eval(const map<int, int>& evalmap) const
    {
        if (coefficient == 0)
            return 0;

        for (int i = 0; i < monrep.size(); i++)
        {
            if (monrep[i] == 1)
            {
               const auto & it = evalmap.find(i);
               if (it == evalmap.end())
                   return 0;
               else if (it->second == 0)
                   return 0;
            }
        }

        return 1;
        
    }

    void subs(const map<int, int>& submap) 
    {

        for (auto& var_val : submap)
        {
            int val = var_val.second % 2;
            int var = var_val.first;
            try
            {
                if (var >= monrep.size())
                    throw string("The input substitute table is invalid.");
            }
            catch (string s)
            {
                cerr << __func__<<":"<< s << endl;
                exit(-1);
            }

            if(monrep[var] == 1)
                if (val == 0)
                {
                    monrep.reset();
                    coefficient = 0;
                }
                else
                {
                    monrep[var] = 0;
                }


        }
    }

    BooleanPolynomial subs(const vector<BooleanPolynomial>& submap) const;

   


    BooleanMonomial intersect(const dynamic_bitset<>& mask)
    {
        BooleanMonomial _res;

        if (mask.size() != monrep.size())
            return _res;

        _res.resize(monrep.size());

        if (coefficient == 0)
        {
            _res.coefficient = 0;
            return _res;
        }
        else
        {
            _res.monrep = monrep & mask;
            return _res;
        }

    }


    bool operator[] (int i) const
    {
        return monrep[i] & (coefficient != 0);
    }

    int count() const
    {
        if (coefficient == 0)
            return 0;

        return monrep.count();
    }

    void print(ostream& os, const char c = 's') const
    {
        // symbol is 's' by default
        int len = this->size();
        if (len == 0)
            os << "null";
        else if (coefficient == 0)
            os << "0";
        else
        {
            if (this->count() > 0)
            {
                for (int i = 0; i < len; i++)
                    if ((*this)[i])
                        os << c << i;
            }
            else
                os << "1";

        }
    }


    void operator *= (const BooleanMonomial& _mon)
    {

        int _size0 = this->size();
        int _size1 = _mon.size();
        if (_size0 == _size1)
        {
            if(!_mon.iszero())
                monrep |= _mon.monrep;
            else
            {
                coefficient = 0;
                monrep.reset();
            }
        }
    }

    void operator /= (const BooleanMonomial& _mon)
    {
        if (_mon.iszero())
            return;

        int _size0 = this->size();
        int _size1 = _mon.size();
        if (_size0 == _size1 &&  !_mon.iszero())
        {
            auto checkflag = ((monrep & _mon.monrep) == _mon.monrep);
            if (checkflag)
                monrep = (monrep & (~_mon.monrep));
        }
    }

    friend bool operator == (const BooleanMonomial& _mon0, const BooleanMonomial& _mon1)
    {

        return ((_mon0.iszero() && _mon1.iszero()) && (_mon0.size() == _mon1.size())) || 
            ( (_mon0.monrep == _mon1.monrep) && (!_mon0.iszero() && !_mon1.iszero()) );
    }


    friend bool operator != (const BooleanMonomial& _mon0, const BooleanMonomial& _mon1)
    {
        return !(_mon0 == _mon1);
    }

    friend BooleanMonomial operator * (const BooleanMonomial& _mon0, const BooleanMonomial& _mon1)
    {
        int _size0 = _mon0.size();
        int _size1 = _mon1.size();
        if (_size0 == _size1)
        {
            if (_mon1.iszero() || _mon0.iszero())
            {
                return BooleanMonomial(_size0, 0);
            }
            else
            {
                BooleanMonomial _res(_mon0);
                _res *= _mon1;
                return _res;
            }
        }

        return BooleanMonomial();
    }


    friend BooleanMonomial operator / (const BooleanMonomial& _mon0, const BooleanMonomial& _mon1)
    {
        int _size0 = _mon0.size();
        int _size1 = _mon1.size();
        if (_size0 == _size1 && !_mon1.iszero())
        {
            if (_mon0.iszero())
            {
                return BooleanMonomial(_size0, 0);
            }
            else
            {
                BooleanMonomial _res(_mon0);
                _res /= _mon1;
                return _res;
            }
        }

        return BooleanMonomial();
    }

    friend bool operator < (const BooleanMonomial& _mon0, const BooleanMonomial& _mon1)
    {
        auto _size0 = _mon0.size();
        auto _size1 = _mon1.size();
        if (_size0 < _size1)
            return true;
        else if (_size0 > _size1)
            return false;
        else
        {
            if (!_mon0.iszero() && !_mon1.iszero())
            {
                for (int i = 0; i < _size0; i++)
                    if (_mon0[i] < _mon1[i])
                        return true;
                    else if (_mon0[i] > _mon1[i])
                        return false;


                return false;
            }
            else if (_mon0.iszero() && !_mon1.iszero())
                return true;
            else
                return false;
        }
    }

    friend bool operator > (const BooleanMonomial& _mon0, const BooleanMonomial& _mon1)
    {
        return _mon1 < _mon0;
    }

    
    friend BooleanPolynomial operator +(const BooleanMonomial& _mon0, const BooleanMonomial& _mon1);


    friend ostream& operator << (ostream& os, const BooleanMonomial& mon)
    {
        mon.print(os);

        return os;
    }
    
    
};





class BooleanPolynomial
{
private:
    vector<BooleanMonomial> polyrep;

    char symbol = 's';

    int varsnum = 0;


public:
    BooleanPolynomial() {};
    BooleanPolynomial(const BooleanMonomial& mon, const char s = 's') {
        polyrep.emplace_back(mon);
        varsnum = mon.size();
        symbol = s;
    }


    BooleanPolynomial(const BooleanPolynomial& poly)
    {
        polyrep = poly.polyrep;
        varsnum = poly.varsnum;
        symbol = poly.symbol;
    }
    
    BooleanPolynomial(int ini_size, vector<dynamic_bitset<>>& mon_reps)
    {
        varsnum = ini_size;
        for (auto& mon_rep : mon_reps)
            polyrep.emplace_back(BooleanMonomial(mon_rep));
    }

    
    BooleanPolynomial(int ini_size, const string polyexp, const char s = 's')
    {
        if (ini_size <= 0)
            return;

        string exp = remove_blank(polyexp);

        using BooleanResolution::PolyInf;
        using BooleanResolution::resolve_poly_exp;
        using BooleanResolution::SYMBOL_TYPE;
        using BooleanResolution::SYMBOL_VALUE;


        PolyInf moninfs = resolve_poly_exp(exp);
        if (moninfs.size() == 0)
            return;

        
        SYMBOL_VALUE prev_symbol = '?';

        for (auto& moninf : moninfs)
        {
            auto& monsymbol = moninf.first;
            if (monsymbol.type == SYMBOL_TYPE::ONE)
            {
                BooleanMonomial mon(ini_size);
                polyrep.emplace_back(mon);
                symbol = s;
                continue;
            }
            else if (monsymbol.type == SYMBOL_TYPE::ZERO)
            {
                symbol = s;
                continue;
            }
            else
            {
                if (prev_symbol == '?')
                    prev_symbol = monsymbol.val;
                else
                    if (prev_symbol != monsymbol.val)
                    {
                        polyrep.clear();
                        polyrep.shrink_to_fit();
                        return;
                    }
            }

            BooleanMonomial mon(ini_size, moninf.second);
            if (mon.isnull())
            {
                polyrep.clear();
                polyrep.shrink_to_fit();
                return;
            }
            else
            {
                polyrep.emplace_back(mon);
            }
        }

        if (prev_symbol == '?')
            symbol = s;
        else
            symbol = char(prev_symbol);

        varsnum = ini_size;

        
        sort(polyrep.begin(), polyrep.end());
    }
    



    void resize(int n)
    {
        varsnum = n;
        for (auto& mon : polyrep)
            mon.resize(n);
    }

    int size() const
    {
        return varsnum;
    }

    int moncnt() const
    {
        return polyrep.size();
    }

    int deg() const
    {
        int deg = 0;
        for (auto& mon : polyrep)
            if (mon.count() > deg)
                deg = mon.count();

        return deg;
    }

    

    BooleanPolynomial & emplace_back(const BooleanMonomial& mon)
    {
        if (mon.size() == varsnum)
            polyrep.emplace_back(mon);

        return *this;
    }

    void reduce()
    {
       
        map<BooleanMonomial, int> monmap;
        for (auto& mon : polyrep)
            if(!mon.iszero())
                monmap[mon]++;

        polyrep.clear();
        polyrep.shrink_to_fit();
        for (auto& moncnt : monmap)
            if (moncnt.second % 2)
                polyrep.emplace_back(moncnt.first);
        
    }




    bool iszero() const { return (polyrep.size() == 0 && varsnum != 0) || (polyrep.size() == 1 && polyrep[0].iszero()); } 

    bool isone() const { return polyrep.size() == 1 && polyrep[0].isone(); }

    bool isnull() const { return varsnum == 0; }

    void operator = (const BooleanMonomial& _mon)
    {
        polyrep.emplace_back(_mon);
        varsnum = _mon.size();
        symbol = 's';
    }

    void operator += (const BooleanMonomial& _mon)
    {
        if (varsnum != _mon.size())
            return;

        this->emplace_back(_mon);
        this->reduce();
    }

    void operator += (const BooleanPolynomial& _poly)
    {
        if (_poly.varsnum != varsnum || _poly.symbol != symbol)
            return;

        map<BooleanMonomial, int> monmap;
        for (auto& mon : polyrep)
            monmap[mon]++;
        for (auto& mon : _poly.polyrep)
            monmap[mon]++;

        polyrep.clear();
        polyrep.shrink_to_fit();
        for (auto& moncnt : monmap)
            if (moncnt.second % 2)
                polyrep.emplace_back(moncnt.first);

    }

    void operator *= (const BooleanPolynomial& _poly)
    {
        if (_poly.varsnum != varsnum || _poly.symbol != symbol || _poly.isone())
            return;

        map<BooleanMonomial, int> monmap;
        for (auto& mon0 : polyrep)
            for (auto& mon1 : _poly.polyrep)
                monmap[mon0 * mon1] ++;
	
	
	    polyrep.clear();
	    polyrep.shrink_to_fit();
	    //cout << polyrep.capacity()<<endl;

        for (auto& moncnt : monmap)
            if (moncnt.second % 2)
                polyrep.emplace_back(moncnt.first);

    }

    void operator *= (const BooleanMonomial& _mon)
    {
        if (varsnum != _mon.size() || _mon.isone())
            return;

        map<BooleanMonomial, int> monmap;
        for (auto& mon0 : polyrep)
            monmap[mon0 * _mon] ++;

        polyrep.clear();
        polyrep.shrink_to_fit();
        for (auto& moncnt : monmap)
            if (moncnt.second % 2)
                polyrep.emplace_back(moncnt.first);

    }

    auto begin() const
    {
        return polyrep.begin();
    }

    auto end() const
    {
        return polyrep.end();
    }

    BooleanMonomial operator[] (int i) const
    {
        return polyrep[i];
    }

    int eval(const map<int, int>& evalmap) const
    {
        int res = 0;
        for (auto& mon : polyrep)
        {
            res ^= mon.eval(evalmap);
        }

        return res;
    }

    int eval(const dynamic_bitset<>& evalmap) const
    {
        int res = 0;
        for (auto& mon : polyrep)
        {
            res ^= mon.eval(evalmap);
        }

        return res;
    }

    void subs(const map<int, int>& submap)
    {
        for (auto& mon : polyrep)
        {
            mon.subs(submap);
        }

        this->reduce();
    }

    BooleanPolynomial subs(const vector<BooleanPolynomial>& submap) const
    {

        vector<BooleanPolynomial> sum_mons;
        for (auto& mon : polyrep)
        {
            auto part_res = mon.subs(submap);
            sum_mons.emplace_back(part_res);
        }
        return fastsum(varsnum, sum_mons);
    }


    
    friend BooleanPolynomial operator + (const BooleanPolynomial& _poly0, const BooleanPolynomial& _poly1)
    {
        BooleanPolynomial _res;
        if (_poly0.symbol != _poly1.symbol)
            return _res;
        int _size0 = _poly0.size();
        int _size1 = _poly1.size();
        if (_size0 == _size1)
        {
            _res = _poly0;
            _res += _poly1;
        }
        return _res;
    }

    friend BooleanPolynomial fastsum(int varsnum, const vector<BooleanPolynomial>& polys)
    {
        BooleanPolynomial _res;
        _res.varsnum = varsnum;
        if (polys.size() == 0)
            return _res;

        _res.symbol = polys[0].symbol;



        map<BooleanMonomial, int> monmap;
        for (auto& poly : polys)
            for (auto& mon : poly.polyrep)
                monmap[mon]++;
        
        for (auto& moncnt : monmap)
            if (moncnt.second % 2)
                _res.polyrep.emplace_back(moncnt.first);


        return _res;
    }

    friend BooleanPolynomial fastmul(int varsnum, const set<BooleanPolynomial>& polys)
    {
        BooleanMonomial mon_1(varsnum);
        BooleanPolynomial _res(mon_1);
        if (polys.size() == 0)
            return _res;

        auto cur_res = polys;

        while (cur_res.size() != 1)
        {
            auto& first_poly = (*cur_res.begin());
            if (first_poly.iszero())
            {
                _res.polyrep.clear();
                //_res.polyrep.shrink_to_fit();
                return _res;
            }

            set<BooleanPolynomial> next_res;
            auto it0 = cur_res.begin();
            auto it1 = cur_res.rbegin();
            while ((*it0) < (*it1))
            {
                auto new_poly = (*it0) * (*it1);
                next_res.emplace(new_poly);

                it0++;
                it1++;
            }

            if ((*it0) == (*it1))
                next_res.emplace(*it0);

            cur_res = next_res;
        }

        _res = (*cur_res.begin());
        return _res;
    }

    friend BooleanPolynomial fastmul(int varsnum, const vector<BooleanPolynomial>& polys)
    {
        set<BooleanPolynomial> polysSet(polys.begin(), polys.end());

        return fastmul(varsnum, polysSet);
    }


    friend BooleanPolynomial operator * (const BooleanPolynomial& _poly0, const BooleanPolynomial&  _poly1)
    {
        BooleanPolynomial _res;
        if (_poly0.symbol != _poly1.symbol)
            return _res;
        int _size0 = _poly0.size();
        int _size1 = _poly1.size();
        if (_size0 == _size1)
        {
            _res = _poly0;
            _res *= _poly1;
        }
        return _res;
    }

    friend bool operator < (const BooleanPolynomial& _poly0, const BooleanPolynomial& _poly1)
    {
        try
        {
            if (_poly0.symbol != _poly1.symbol || _poly0.size() != _poly1.size())
            {
                throw string("Try to compare two polynomials with different symbols or variables");
            }
        }
        catch(string s)
        {
            cerr << __func__ << ":" << s << endl;
            exit(1);
        }
        int _cnt0 = _poly0.moncnt();
        int _cnt1 = _poly1.moncnt();
        if (_cnt0 < _cnt1)
            return true;
        else if (_cnt0 > _cnt1)
            return false;
        else
        {
            for (int i = 0; i < _cnt0; i++)
            {
                if (_poly0[i] < _poly1[i])
                    return true;
                else if (_poly0[i] > _poly1[i])
                    return false;
                
            }

            return false;
        }
    }

    friend bool operator > (const BooleanPolynomial& _poly0, const BooleanPolynomial& _poly1)
    {
        return _poly1 < _poly0;
    }

    friend bool operator ==(const BooleanPolynomial& _poly0, const BooleanPolynomial& _poly1)
    {
        try
        {
            if (_poly0.symbol != _poly1.symbol || _poly0.size() != _poly1.size())
            {
                throw string("Try to compare two polynomials with different symbols or variables");
            }
        }
        catch (string s)
        {
            cerr << __func__ << ":" << s << endl;
            exit(1);
        }

        int _cnt0 = _poly0.moncnt();
        int _cnt1 = _poly1.moncnt();
        if (_cnt0 != _cnt1)
            return false;

        for (int i = 0; i < _cnt0; i++)
            if (_poly0[i] != _poly1[i])
                return false;
        return true;
    }

    friend bool operator !=(const BooleanPolynomial& _poly0, const BooleanPolynomial& _poly1)
    {
        return !(_poly0 == _poly1);
    }


    friend ostream& operator << (ostream& os, const BooleanPolynomial& poly)
    {
        // symbol is 's' by default
        if (poly.size() == 0)
            os << "null";
        else if (poly.iszero())
            os << "0";
        else
        {
            string sep = "";
            for (auto& mon : poly.polyrep)
            {
                os << sep;
                mon.print(os, poly.symbol);
                sep = "+";
            }

        }

        return os;
    }

};



// 
BooleanPolynomial operator +(const BooleanMonomial& _mon0, const BooleanMonomial& _mon1)
{
    BooleanPolynomial _res;
    auto _size0 = _mon0.size();
    auto _size1 = _mon1.size();
    if (_size0 == _size1)
    {
        _res.resize(_size0);
        if (_mon0 != _mon1)
            (_res.emplace_back(_mon0)).emplace_back(_mon1);
    }

    _res.reduce();

    return _res;

}

BooleanPolynomial BooleanMonomial::subs(const vector<BooleanPolynomial>& submap) const
{
    if (monrep.size() != submap.size())
    {
        cerr << __func__ << ":" << "The input substitute table is invalid." << endl;
        cerr << "Monrep size: " << monrep.size() << endl;
        cerr << "Submap size: " << submap.size() << endl;
        exit(-1);
    }

    BooleanPolynomial res(monrep.size(), "1");

    for (auto& i : this->index())
        res *= submap[i];


    return res;

}

#endif
