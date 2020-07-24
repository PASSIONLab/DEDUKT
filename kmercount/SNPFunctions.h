#ifndef SNP_FUNCTIONS_H
#define SNP_FUNCTIONS_H

#include "Kmer.hpp"

struct SNPdata
{
	Kmer::MERARR karr;
	char extA;
	char extB;

	bool operator > (const SNPdata & rhs) const
	{ return (karr > rhs.karr); }
	bool operator < (const SNPdata & rhs) const
	{ return (karr < rhs.karr); }
	bool operator == (const SNPdata & rhs) const
	{ return (karr == rhs.karr); }
};


void OpenDebugFile(string prefix, ofstream & output) 
{
	stringstream ss;
	string rank;
	ss << myrank;
	ss >> rank;
	string ofilename = prefix;
	ofilename += rank;
	output.open(ofilename.c_str(), ios_base::app );
}

//! Merge is simple because the entries are guarenteed to be distinct
//! Arrs2Merge is a vector whose elements are sorted vectors within themselves
//! Arrs2Merge doesn't contain the lexicographically larger "twin" k-mer
//! So this function also creates the twins and merges them on the fly
void MergeAll(vector< vector<kmerpack> > & Arrs2Merge, vector< pair<Kmer,int> > & MergedKMers)
{
	int hsize =  Arrs2Merge.size();		
	size_t mergedsize;
	size_t repind = 0;
	size_t twnind;
	if(hsize == 1)
	{
		mergedsize = Arrs2Merge[0].size();
		twnind = mergedsize;
		MergedKMers.resize(mergedsize*2);
		for(auto itr=Arrs2Merge[0].begin(); itr != Arrs2Merge[0].end(); ++itr)
		{
			Kmer kmer(itr->arr);
			MergedKMers[repind++] = make_pair(kmer, itr->count);	
			MergedKMers[twnind++] = make_pair(kmer.twin(), itr->count);
		}
	}
	else if(hsize > 1)
	{
		pair<Kmer::MERARR, int> * heap = new pair<Kmer::MERARR, int> [hsize];	// int -> source id
		vector<size_t> curptr(hsize, 0);

		mergedsize = 0;
		for(int i=0; i< hsize; ++i)
		{
			mergedsize += Arrs2Merge[i].size();
			heap[i] = make_pair(Arrs2Merge[i][0].arr, i);
		}	
		make_heap(heap, heap+hsize, greater< pair<Kmer::MERARR, int > >());	// to ensure a min-heap
		MergedKMers.resize(mergedsize*2);	// twice as big for the twins
		twnind = mergedsize;

		while(hsize > 0)
		{
			pop_heap(heap, heap + hsize, greater< pair<Kmer::MERARR, int> >());         // result is stored in heap[hsize-1]
			int source = heap[hsize-1].second;

			kmerpack kpack = Arrs2Merge[source][curptr[source]++];
			Kmer kmer(kpack.arr);
			MergedKMers[repind++] = make_pair(kmer, kpack.count);
			MergedKMers[twnind++] = make_pair(kmer.twin(), kpack.count);

			if(curptr[source] != Arrs2Merge[source].size())	// That array has not been depleted
			{
				heap[hsize-1] = make_pair(Arrs2Merge[source][curptr[source]].arr, source);
				push_heap(heap, heap+hsize, greater< pair<Kmer::MERARR, int> >());
			}
			else
			{
				--hsize;
			}
		}
		delete [] heap;
	
		for(size_t i=0; i<hsize; ++i)
			vector<kmerpack>().swap(Arrs2Merge[i]);
		vector< vector<kmerpack> >().swap(Arrs2Merge);
	}
	if(mergedsize != repind)	cout << "Counting is wrong" << endl;

	if(myrank == 0)
	{
		for(int j=1; j< mergedsize; ++j)
		{
			if(MergedKMers[j] < MergedKMers[j-1])
			{
				cout << "Sorting issue: "<< endl;
				cout << "Loc[" << j-1 << "]:" << MergedKMers[j-1].first.toString(); 
				cout << ", Loc[" << j << "]:" << MergedKMers[j].first.toString() << endl;
			}
		}
	}
	if(!is_sorted(MergedKMers.begin(), MergedKMers.begin()+mergedsize)) cout << "Original section not sorted" << endl;	
	sort(MergedKMers.begin()+mergedsize, MergedKMers.end());	// sort the twin range
	inplace_merge (MergedKMers.begin(),MergedKMers.begin()+mergedsize,MergedKMers.end());
	if(!is_sorted(MergedKMers.begin(), MergedKMers.end())) cout << "Whole thing is not sorted" << endl;	
}

#endif
