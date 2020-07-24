#ifndef _SIMPLE_COUNTER_
#define _SIMPLE_COUNTER_

#include <cassert>
#include <iostream>

#include <queue>
#include <algorithm>
#include <sstream>
#include <cstring>

#include <mpi.h>

#include "../common/stdlib_compatibility.hpp"

using namespace std;

template <class T, class H>
class SimpleCount
{
public:
    typedef U_MAP< T, int, H> Map;
    typedef typename Map::value_type Value;
    typedef typename Map::iterator Iterator;
    typedef std::vector<T> Vector;
    SimpleCount(size_t totrack):maxsize(totrack),inserted(0),reduced(0),erased(0),reduceTime(0.0),mergeTime(0.0),incrementTime(0.0),insertTime(0.0),findTime(0.0) {
//        mapmajorityset.reserve(totrack);
    }
    
    template<typename InputIterator>
    void PushAll(InputIterator first, InputIterator last)
    {
        while (first!=last)
        {
            Push(*first);
            ++first;
        }
    }
    
    void Push(T item)
    {
        double now = MPI_Wtime();
        auto fentry = mapmajorityset.find(item); // FIXME This take a huge amount of time when kmer_len is small (15, 13, 11, ....)
        double now2 = MPI_Wtime();
        if (now2>now) findTime += now2-now;
        if( fentry != mapmajorityset.end())
        {
            ++(fentry->second);
            now = MPI_Wtime();
            if (now>now2) incrementTime += now - now2;
        }
        else
        {
            mapmajorityset.insert(make_pair(item,1));
            inserted++;
            now = MPI_Wtime();
            if (now>now2) insertTime += now - now2;
        }
        if(mapmajorityset.size() >= maxsize)
        {
            Reduce(1);
        }
    }
    void Reduce(int decrement) {
            double now = MPI_Wtime();

            assert(decrement > 0);
            vector< Iterator > victims;
            
            for (auto mit=mapmajorityset.begin(); mit!=mapmajorityset.end(); ++mit)
            {
                mit->second -= decrement; // first decrement
                assert(mit->second >= 0);
                
                // map::erase = Iterators, pointers and references referring to elements removed by the function
                // are invalidated. All other iterators, pointers and references keep their validity.
                if(mit->second <= 0) { victims.push_back(mit); }
            }
            for(auto qit=victims.begin(); qit!=victims.end(); ++qit)
            {
                mapmajorityset.erase(*qit);
            }
            size_t numErased = victims.size();
            erased += numErased;
            reduced++;

            reduceTime += MPI_Wtime() - now;
    }
    std::ostream & PrintAll(std::ostream &os)
    {
        os << "We have " << mapmajorityset.size() << " entries " << endl;
        for (auto mit=mapmajorityset.begin(); mit!=mapmajorityset.end(); ++mit)
        {
            os << mit->first << "(" << mit->second << ") ";
        }
        os << endl;
        return os;
    }
    std::ostream & PrintTop(std::ostream &os, int currank)
    {
        os << currank << " has " << mapmajorityset.size() << " entries but only printing top 40" << endl;
        vector< pair<int,T> > majoritysorted;
        for (auto mit=mapmajorityset.begin(); mit!=mapmajorityset.end(); ++mit)
        {
            majoritysorted.push_back(make_pair(mit->second, mit->first));
        }
        sort(majoritysorted.begin(), majoritysorted.end(), std::greater<pair<int,T>>());
        
        int i=0;
        for (auto mit=majoritysorted.begin(); mit!=majoritysorted.end() && i < 40; ++mit, ++i)
        {
            os << mit->second << "(" << mit->first << ") ";
        }
        os << endl;
        return os;
    }
    
    void MergeableSummary(vector<T> & keys, vector<int> & counts)
    {
        double now = MPI_Wtime();
        for_each(mapmajorityset.begin(),mapmajorityset.end(),
                 [&counts, &keys](pair<const T,int> & p){ counts.push_back(p.second); keys.push_back(p.first); });
        mergeTime += MPI_Wtime() - now;
    }
    
    bool IsMember(T item)
    {
        return (mapmajorityset.find(item) != mapmajorityset.end());
    }
    
    void CreateIndex()
    {
        assert(mapmajorityset.size() <= maxsize);
        int idx = 0;
        internal_v.resize( mapmajorityset.size(), T());
        for(auto iter = mapmajorityset.begin(); iter != mapmajorityset.end(); iter++) {
            internal_v[idx] = iter->first;
            ++idx;
        }
        internal_m.clear();
        sort( internal_v.begin(), internal_v.end() );
        for(idx = 0; idx < internal_v.size(); idx++) {
            internal_m[ internal_v[idx] ] = idx;
        }
    }
    void Clear()
    {
        Map().swap(mapmajorityset);
        Map().swap(internal_m);
        Vector().swap(internal_v);
    }
    
    //! make sure CreateIndex() is called first
    size_t FindIndex(T item) // finds the location of the item in internal_ vector, given the value
    {
        int dummy = 0;
        assert( internal_v.size() > 0 || mapmajorityset.empty() );
        assert( internal_v.size() <= maxsize );
        assert( internal_v.size() == internal_m.size() );
        auto i = internal_m.find(item);
        if (i != internal_m.end()) {
          assert(i->second < maxsize);
          return i->second;
        } else {
          return maxsize;
        }
    }
    
    //! make sure CreateIndex() is called first
    T Get(size_t index) // finds the location of the item in internal_ vector, given the value
    {
        assert(index < internal_v.size());
        return internal_v[index];
    }
    
    size_t Size()
    {
        return mapmajorityset.size();
    }

    size_t GetMaxSize()
    {
        return maxsize;
    }

    size_t Grow()
    {
        maxsize *= 2;
        mapmajorityset.reserve(maxsize);
        return maxsize;
    }

    void Shrink()
    {
        maxsize /= 2;
    }

    size_t GetNumReduced() { return reduced; }
    size_t GetNumErased() { return erased; }
    size_t GetNumInserted() { return inserted; }
    double GetReduceTime() { return reduceTime; }
    double GetMergeTime() { return mergeTime; }
    double GetIncrementTime() { return incrementTime; }
    double GetFindTime() { return findTime; }
    double GetInsertTime() { return insertTime; }

    std::string getPerformanceStats() {
        std::stringstream ss;
        ss << "Reduced " << reduced << " times, Erased " << erased << " out of " << inserted << " inserted kmers: reduce=" << reduceTime << "  merge=" << mergeTime << " insert=" << insertTime << " find=" << findTime << " increment=" << incrementTime;
        std::string s = ss.str();
        return s;
    }
    
    void Assign(T * lhskeys, int * lhscounts, size_t lhssize)   // inverse of MergeableSummary
    {
        double now = MPI_Wtime();
        mapmajorityset.clear(); // delete existing map but retain the maxsize
        for(int i=0; i< lhssize; ++i)
        {
            mapmajorityset.insert(make_pair(lhskeys[i], lhscounts[i]));
        }
        assert(mapmajorityset.size() == lhssize);
        mergeTime += MPI_Wtime() - now;
    }

    
    void Merge(T * lhskeys, int * lhscounts, size_t lhssize, bool allowExpansion = false)
    {
        double now = MPI_Wtime();
        int minweight;
        for (int i = 0; i < lhssize; i++)
        {
            while(lhscounts[i] > 0)
            {
                auto it = mapmajorityset.find(lhskeys[i]);
                if (it != mapmajorityset.end()) // X[z] is the monitored element of a counter c
                {
                    it->second = it->second + lhscounts[i];
                    lhscounts[i] = 0;   // done with this element
                }
                else    // not monitored yet
                {
                    if(allowExpansion || mapmajorityset.size() < maxsize)    // there is space for him
                    {
                        mapmajorityset.insert(make_pair(lhskeys[i],lhscounts[i]));
                        lhscounts[i] = 0;   // done with this element
                    }
                    else
                    {
                        // we have to find minimum weight (count) not the minimum key, which is linear if implemented with minelement
                        // however, this is OK because the "else" clause potentially touches every element in the list anyway
                        auto minentry = min_element(mapmajorityset.begin(), mapmajorityset.end(),
                                        [] (Value & l, Value & r) { return l.second < r.second; });
                        minweight = minentry->second;
                        
                        int decrementval = std::min(minweight, lhscounts[i]);   // because counts can not go negative
                        assert(decrementval > 0);
                        Reduce(decrementval);
                        lhscounts[i] -= decrementval;
    
                    }
                }
            }
          }
          mergeTime += MPI_Wtime() - now;
	}

    Map mapmajorityset;
    size_t maxsize;
private:
    size_t inserted, reduced, erased;
    double reduceTime, mergeTime, incrementTime, insertTime, findTime;
    Map internal_m;
    Vector internal_v;
};




#endif
