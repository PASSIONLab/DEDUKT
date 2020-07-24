#ifndef _EFFICIENT_COUNTER_
#define _EFFICIENT_COUNTER_

#include <cassert>
#include <iostream>

#include "../common/stdlib_compatibility.hpp"

#include <list>
#include <set>
#include <algorithm>
#include <sstream>
#include <mpi.h>
#include <cstring>


using namespace std;


/* 
 ** To support decrementing all counters at once in constant time, we store the counters
 ** in sorted order using a differential encoding. That is, each counter actually
 ** only stores how much larger it is compared to the next smallest counter.
 ** Now incrementing and decrementing counters requires them to move significantly in the
 ** total order; to support these operations, we coalesce equal counters
 ** (differentials of zero) into common groups.
 */
template <class T>
class Group
{
public:
    Group(T newele, int mydiff)
    {
        elements.insert(newele);
        diff = mydiff;
    }
    std::set<T> elements;
    int diff;   // diff from the previous group
};


/*
 ** The overall structure is a doubly linked list of groups, ordered by counter value.
 ** Each group represents a collection of equal counters, consisting of two parts:
 **      (1) a doubly linked list of counters
 **          (in no particular order, because they all have the same value), and
 **      (2) the difference in value between these counters and the counters in the previous group, or,
 **         for the first group, the value itself (i.e., diff from a dummy of count 0).
 */
template <class T, class H>
class EfficientCount
{
    using GroupList = list< Group<T> >;
    using GroupListItr = typename GroupList::iterator;
    using Map2Chain = U_MAP< T, GroupListItr , H>;
    
private:
    Map2Chain mapmajorityset;
    GroupList groupchain; // list guarantee: inserting/deleting elements doesn't invalidate pointers, references, and iterators to other elements
    size_t maxsize;
    
public:
    EfficientCount(size_t totrack):maxsize(totrack),inserted(0),reduced(0),erased(0),reduceTime(0.0),mergeTime(0.0) {}
    
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
        auto fentry = mapmajorityset.find(item);
        if( fentry != mapmajorityset.end()) // the element already exists
        {
            GroupListItr groupptr = fentry->second;  // ListItr is bidirectional
            auto entryptr = (groupptr->elements).find(item);     // entryptr is of type set<T>::iterator
            if(groupptr->elements.size() > 1)   // the group was not a singleton
            {
                groupptr->elements.erase(entryptr);
                
                GroupListItr nextgroupptr = groupptr;
                nextgroupptr++; // move to the next group
                if(nextgroupptr != groupchain.end())
                {
                    if(nextgroupptr->diff == 1) // either the next group accepts this new element (diff = 1)
                    {
                        (nextgroupptr->elements).insert(item);
                        fentry->second  = nextgroupptr;
                    }
                    else //  or we create a new group with this new element
                    {
                        // list::insert -> the container is extended by inserting new elements *before* the element at the specified position
                        GroupListItr nEltItr= groupchain.insert(nextgroupptr, Group<T>(item, 1));
                        (nextgroupptr->diff)--; // decrement the diff of the next group
                        fentry->second = nEltItr;    // update the iterator in the map
                    }
                }
                else    // the next group doesn't exist (end of list); so we create it at the end of the list
                {
                    groupchain.push_back(Group<T>(item, 1));
                    fentry->second = std::prev(groupchain.end());  // update the iterator in the map
                }
            }
            else    // the group is a singleton
            {
                GroupListItr nextgroupptr = groupptr;
                nextgroupptr++; // move to the next group
                if(nextgroupptr != groupchain.end())
                {
                    if(nextgroupptr->diff == 1) // merge with the next group
                    {
                        (nextgroupptr->elements).insert(item);
                        nextgroupptr->diff += groupptr->diff;   // groupptr will be deleted in a moment so nextgroupptr's diff increases by that much
                        groupchain.erase(groupptr);
                        fentry->second = nextgroupptr; // update the iterator in the map
                    }
                    else    // just adjust diff scores
                    {
                        ++(groupptr->diff);         // no need to update the iterator
                        --(nextgroupptr->diff);
                    }
                }
                else    // the next group doesn't exist
                {
                    ++(groupptr->diff);           // no need to update the iterator
                }
            }
        }
        else    // item doesn't exist
        {
            
            GroupListItr groupptr = groupchain.begin();
            if(groupptr != groupchain.end())
            {
                if(groupptr->diff == 1) // we are assuming the first entry being diff=1 means its count is one (i.e., diff from a dummy of count 0)
                {
                    (groupptr->elements).insert(item);
                    mapmajorityset.insert(make_pair(item,groupptr));  // nothign to delete in the map
                }
                else    // there is no group with count==1
                {
                    --(groupptr->diff); // decrement that guy's diff as we will insert a new entry with "diff=1" before it
                    groupchain.push_front(Group<T>(item, 1));
                    mapmajorityset.insert(make_pair(item, groupchain.begin())); // nothign to delete in the map
                }
            }
            else    // the whole thing was empty
            {
                groupchain.push_front(Group<T>(item, 1));
                mapmajorityset.insert(make_pair(item, groupchain.begin())); // nothign to delete in the map
            }
        }
        // original rule: "decrement every counter only if the new item is not monitored and there is no space for it"
        // we are implementing a slightly different rule here as mapmajorityset can only overflow after insertions
        // so we are inserting even if there is no place for it originally, and then dealing with it here
        // ReduceByOne actually might not change the size of mapmajorityset if the diff of the first entry is larger than one
        if(mapmajorityset.size() >= maxsize)
        {
            ReduceByOne();
        }
    }
    //! @pre{never called on an empty set since it is only called by Push, hence no need to check for emptiness}
    void ReduceByOne()
    {
            double now = MPI_Wtime();
            GroupListItr groupptr = groupchain.begin();
            size_t numErased;
            if(groupptr->diff == 1)
            {
                numErased = groupptr->elements.size();
                for(auto toerase : groupptr->elements)
                {
                    mapmajorityset.erase(toerase);
                }
                groupchain.pop_front();
            }
            else
            {
                --(groupptr->diff);
            }
        
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
        
        for (auto gr : groupchain ) // gr is of type
        {
            for (auto entry : gr.elements)  // entry is of type T
            {
                majoritysorted.push_back(make_pair(gr.diff, entry));
            }
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
       /* double now = MPI_Wtime();
        for_each(mapmajorityset.begin(),mapmajorityset.end(),
                 [&counts, &keys](pair<const T,int> & p){ counts.push_back(p.second); keys.push_back(p.first); });
        mergeTime += MPI_Wtime() - now; */
    }
    
    bool IsMember(T item)
    {
        return (mapmajorityset.find(item) != mapmajorityset.end());
    }
    
    void CreateIndex()
    {
        /*
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
        }*/
    }
    void Clear()
    {
        /*
        Map().swap(mapmajorityset);
        Map().swap(internal_m);
        Vector().swap(internal_v);*/
    }
    
    //! make sure CreateIndex() is called first
    size_t FindIndex(T item) // finds the location of the item in internal_ vector, given the value
    {
        /*
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
        }*/
        return 1;
    }
    
    //! make sure CreateIndex() is called first
   /* T Get(size_t index) // finds the location of the item in internal_ vector, given the value
    {
        assert(index < internal_v.size());
        return internal_v[index];
    }*/

    size_t GetNumReduced() { return reduced; }
    size_t GetNumErased() { return erased; }
    double GetReduceTime() { return reduceTime; }
    double GetMergeTime() { return mergeTime; }

    std::string getPerformanceStats()
    {
        std::stringstream ss;
        ss << "Reduced " << reduced << " times, Erased " << erased;
        ss << " out of " << inserted << " inserted kmers: " << reduceTime << " s. Merge: " << mergeTime;
        std::string s = ss.str();
        return s;
    }
    
    void Assign(T * lhskeys, int * lhscounts, size_t lhssize)   // inverse of MergeableSummary
    {
        /*
        double now = MPI_Wtime();
        mapmajorityset.clear(); // delete existing map but retain the maxsize
        for(int i=0; i< lhssize; ++i)
        {
            mapmajorityset.insert(make_pair(lhskeys[i], lhscounts[i]));
        }
        assert(mapmajorityset.size() == lhssize);
        mergeTime += MPI_Wtime() - now;*/
    }

    
    void Merge(T * lhskeys, int * lhscounts, size_t lhssize, bool allowExpansion = false)
    {
        /*
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
          mergeTime += MPI_Wtime() - now; */
	}

   
private:
    size_t inserted, reduced, erased;
    double reduceTime, mergeTime;
    //Map internal_m;
    //Vector internal_v;
};




#endif
