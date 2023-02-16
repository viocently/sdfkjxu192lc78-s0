#pragma once
#ifndef __NODE_H__
#define __NODE_H__
#include<vector>
#include<bitset>
#include<map>
#include"log.h"
#include"dynamic_bitset.hpp"

using namespace std;
using namespace boost;

enum STATUS { SOLUTION, NOSOLUTION, EXPAND };



template<size_t N>
bool CMP( const bitset<N> & a, const bitset<N> & b )
{
        for ( int i = 0; i < N; i++ )
            if ( a[i] < b[i] ) return true;
            else if ( a[i] > b[i] ) return false;
        return false; // equal
}

template<size_t N>
struct CMPS
{
    bool operator()( const bitset<N> & a, const bitset<N> & b ) const
    {
        for ( int i = 0; i < N; i++ )
            if ( a[i] < b[i] ) return true;
            else if ( a[i] > b[i] ) return false;
        return false; // equal
    }
};

template<size_t N>
bool cmps(const bitset<N>& a, const bitset<N>& b) 
{
    for (int i = 0; i < N; i++)
        if (a[i] < b[i]) return true;
        else if (a[i] > b[i]) return false;
    return false; // equal
}

template<size_t N>
struct VCMPS
{
    bool operator()(const vector<bitset<N>>& a, const vector<bitset<N>>& b) const
    {
        if (a.size() < b.size())
            return true;
        else if (a.size() > b.size())
            return false;

        for (int i = 0; i < a.size(); i++)
            if (cmps<N>(a[i],b[i]) ) return true;
            else if (cmps<N>(b[i],a[i])) return false;
        return false; // equal
    }
};


struct DCMPS
{
    bool operator()(const dynamic_bitset<> & a, const dynamic_bitset<> & b) const
    {
        if (a.size() < b.size())
            return true;
        else if (a.size() > b.size())
            return false;

        for (int i = 0; i < a.size(); i++)
            if (a[i] < b[i]) return true;
            else if (a[i] > b[i]) return false;
        return false; // equal
    }
};


struct STAIRS
{
    int nextStep;
	float time;
};

struct Node
{
    int _rnd;
    STAIRS stair;    
    bitset<288> _vector;
    STATUS _status;
    map< bitset<80>, int, CMPS<80> > _solutions;
    vector<Node> _child;
};

struct cmpNode
{
    bool operator() ( const Node & a, const Node & b ) const 
    {
        if ( a._rnd < b._rnd )
            return true;
        else if ( a._rnd > b._rnd )
            return false;
        return CMP( a._vector, b._vector );
    }
};

template<int N>
bool vecCmp(vector<bitset<N> > _Left, vector<bitset<N> > _Right)
{
    if (_Left.size() == _Right.size())
    {
        for (int i = 0; i < _Left.size(); i++)
        {
            if (_Left[i] == _Right[i])
                continue;
            else
                return false;
        }
        return true;
    }
    return false;
}

template<int n>
vector<bitset<n>> vecAdd(vector<bitset<n>>& _left, vector<bitset<n>>& _right)
{
    map<bitset<n>, int, CMPS<n>> resMap;
    vector<bitset<n>> res;
    for (auto& x : _left)
        resMap[x]++;
    for (auto& y : _right)
        resMap[y] ++;
    
    for (auto& x : resMap)
        if (x.second % 2)
            res.emplace_back(x.first);

    return res;
}

template<int n>
vector<bitset<n>> vecMul (vector<bitset<n>>& _left, vector<bitset<n>>& _right)
{
    map<bitset<n>, int, CMPS<n>> resMap;
    vector<bitset<n>> res;

    for (auto& x : _left)
        for (auto& y : _right)
            resMap[x | y] ++;

    for (auto& x : resMap)
        if (x.second % 2)
            res.emplace_back(x.first);

    return res;
}

template<int n>
void filterMap(map<bitset<n>, int, CMPS<n>>& mp)
{
    int size0 = mp.size();
    map< bitset<n>, int, CMPS<n>> tmp(mp);
    mp.clear();

    for (auto& it : tmp)
    {
        if (it.second % 2 == 1)
            mp[it.first] = it.second;
    }
    int size1 = mp.size();

    logger(__func__ + string(": ") + to_string(size0) + string("\t") +
        to_string(size1));
}


#endif
