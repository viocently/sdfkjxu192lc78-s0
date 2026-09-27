#pragma once
#include <iostream>
#include <bitset>
#include "boost/dynamic_bitset.hpp"

using namespace std;
using namespace boost;

template<int N>
class bitset_cmp
{
public:
    bool operator()(const bitset<N>& _l, const bitset<N>& _r) const
    {

        for (int i = 0; i < N; i++)
            if (_l[i] < _r[i])
                return true;
            else if (_l[i] > _r[i])
                return false;

        return false;
    }
};

class dynamic_bitset_cmp
{
public:
    bool operator() (const dynamic_bitset<>& _l, const dynamic_bitset<>& _r) const
    {
        int _lSize = _l.size();
        int _rSize = _r.size();
        if (_lSize < _rSize)
            return true;
        else if (_lSize > _rSize)
            return false;

        int _Size = _lSize;
        for (int i = 0; i < _Size; i++)
            if (_l[i] < _r[i])
                return true;
            else if (_l[i] > _r[i])
                return false;

        return false;
    }
};

template<int N>
bool operator< (const bitset<N>& _l, const bitset<N>& _r) 
{

    for (int i = 0; i < N; i++)
        if (_l[i] < _r[i])
            return true;
        else if (_l[i] > _r[i])
            return false;

    return false;
}

bool operator< (const dynamic_bitset<>& _l, const dynamic_bitset<>& _r) 
{
    int _lSize = _l.size();
    int _rSize = _r.size();
    if (_lSize < _rSize)
        return true;
    else if (_lSize > _rSize)
        return false;

    int _Size = _lSize;
    for (int i = 0; i < _Size; i++)
        if (_l[i] < _r[i])
            return true;
        else if (_l[i] > _r[i])
            return false;

    return false;
}

using node_pair = pair<dynamic_bitset<>, dynamic_bitset<>>;


inline bool partial_less(const dynamic_bitset<> _l, const dynamic_bitset<>& _r)
{
    return ((_l & _r) == _l);
}

