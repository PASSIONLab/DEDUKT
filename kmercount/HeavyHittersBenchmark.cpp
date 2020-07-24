#include "EfficientCount.h"
#include "SimpleCount.h"
#include <iomanip>
#include <random>
#include <map>
#include <iostream>
#include <sstream>
#include <string>
#include "../common/hash_funcs.h"

using namespace std;

static int HLL_HASH_SEED=313;

struct DumbHash {
    uint64_t operator()(const int64_t &km) const {
        //uint64_t hash;
        return MurmurHash3_x64_64(&km, sizeof(km));//, HLL_HASH_SEED, (void*)&hash);
        //return hash;
    }
};

int main(int argc, char * argv[])
{
    int64_t count;
    if (argc > 1)
        count = atoi(argv[1]);
    else
        count = 1000000;
    
    double const exp_dist_mean   = 1.8;
    double const exp_dist_lambda = 1.0 / exp_dist_mean;
    
    random_device rd_exp;
    random_device rd_norm1;
    random_device rd_norm2;
    
    exponential_distribution<> exp_d (exp_dist_lambda);
    normal_distribution<> norm_d1 (5,1); // normal_distribution( mean = 0.0, stddev = 1.0 );
    normal_distribution<> norm_d2 (7.5,1); // normal_distribution( mean = 0.0, stddev = 1.0 );

    
    mt19937 exp_gen (rd_exp ());
    mt19937 norm_gen1(rd_norm1());
    mt19937 norm_gen2(rd_norm2());
   
    
    std::map<int64_t, int64_t> histogram;
    std::map<int64_t, int64_t> exacts;
    vector<int64_t> all;
    
    
    for (int64_t i =0; i < count; ++i)
    {
        int64_t value_exp = static_cast<int64_t>(exp_d (exp_gen) * 1000);
        int64_t value_norm1 = static_cast<int64_t>(norm_d1 (norm_gen1) * 1000);
        int64_t value_norm2 = static_cast<int64_t>(norm_d2 (norm_gen2) * 1000);
        
        if(value_exp > 0)
        {
            ++histogram[value_exp/100];
            ++exacts[value_exp];
            all.push_back(value_exp);
        }
        if(value_norm1 > 0)
        {
            ++histogram[value_norm1/100];
            ++exacts[value_norm1];
            all.push_back(value_norm1);

        }
        if(value_norm2 > 0)
        {
            ++histogram[value_norm2/100];
            ++exacts[value_norm2];
            all.push_back(value_norm2);
        }
    }
    
    cout << "Smallest number in the list is " << exacts.begin()->first << endl;
    cout << "Largest number in the list is " << exacts.rbegin()->first << endl;
    
    
    SimpleCount<int64_t, DumbHash> HH_ver1(count/1000);
    double now = MPI_Wtime();
    HH_ver1.PushAll(all.begin(), all.end());
    cout << "Simple Count took " << MPI_Wtime() - now << " seconds" << endl;
    string ss1 = HH_ver1.getPerformanceStats();
    cout << ss1 << endl;
    HH_ver1.PrintTop(cout, 0);
    
    
    EfficientCount<int64_t, DumbHash> HH_ver2(count/1000);
    double now2 = MPI_Wtime();
    HH_ver2.PushAll(all.begin(), all.end());
    cout << "Efficient Count took " << MPI_Wtime() - now2 << " seconds" << endl;
    string ss2 = HH_ver2.getPerformanceStats();
    cout << ss2 << endl;
    HH_ver2.PrintTop(cout, 0);

    std::priority_queue<pair<int64_t, int64_t>> pq1;
    for (auto& elt : exacts) {
        pq1.push( make_pair(elt.second, elt.first));
    }
    int success = 0;
    int failure = 0;
    std::priority_queue<pair<int64_t, int64_t>> pq2 = pq1; // deep copy
    for(int i=0; i< 100; ++i)
    {
        pair<int64_t, int64_t> revent = pq1.top();
        cout << revent.first << " -> "  << revent.second << endl;
        
        if(HH_ver1.IsMember(revent.second)) success++;
        else failure++;
        pq1.pop();
    }
    cout << "Out of 100 top entries, found " << success << " of them in SimpleCount" << endl;
    
    success = 0;
    failure = 0;
    for(int i=0; i< 100; ++i)
    {
        pair<int64_t, int64_t> revent = pq2.top();
        cout << revent.first << " -> "  << revent.second << endl;
        
        if(HH_ver2.IsMember(revent.second)) success++;
        else failure++;
        pq2.pop();
    }
    cout << "Out of 100 top entries, found " << success << " of them in EfficientCount" << endl;
    
    for (auto& h : histogram) {
        if (h.first < 0) {
            continue;
        }
        cout << setprecision (2) << std::fixed;
        
        cout << h.first <<  " -> ";
        
        cout << string (h.second/(count/1000), '.') << std::endl;
        // fill constructor: string (size_t n, char c);
        // Fills the string with n consecutive copies of character c.
        
        if (h.second/(count/1000) == 0)
            break;
    }
}
