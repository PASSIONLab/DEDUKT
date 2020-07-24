
#ifndef _FRIENDS_MPI_H_
#define _FRIENDS_MPI_H_
#include "SimpleCount.h"
#include "MPIType.h"

extern int nprocs;
extern int myrank;

struct UFX2ReduceObj
{
    UFX2ReduceObj(): count(0)
    {
        array<int,4> nulltuple = {0,0,0,0};
        ACGTleft = nulltuple;
        ACGTrigh = nulltuple;
    }
    array<int,4> ACGTleft;
    array<int,4> ACGTrigh;
    int count;
};

static bool UFX2ReduceObjHasCount(const struct UFX2ReduceObj &a) {
    return (a.count > 0 || a.ACGTleft[0] || a.ACGTleft[1] || a.ACGTleft[2] || a.ACGTleft[3] || a.ACGTrigh[0] || a.ACGTrigh[1] || a.ACGTrigh[2] || a.ACGTrigh[3]);
}

struct UFXReducer
{
    const UFX2ReduceObj operator()(const UFX2ReduceObj & lUFX, const UFX2ReduceObj & rUFX)
    {
        UFX2ReduceObj retobj;
        transform (lUFX.ACGTleft.begin(), lUFX.ACGTleft.end(), rUFX.ACGTleft.begin(), retobj.ACGTleft.begin(), plus<int>());
        transform (lUFX.ACGTrigh.begin(), lUFX.ACGTrigh.end(), rUFX.ACGTrigh.begin(), retobj.ACGTrigh.begin(), plus<int>());
        retobj.count = lUFX.count + rUFX.count;
        return retobj;
    }
};

static void MPI_UFXReduce(void * invec, void * inoutvec, int * len, MPI_Datatype *datatype)
{
    UFXReducer doReduce;
    UFX2ReduceObj * inveccast = (UFX2ReduceObj *) invec;
    UFX2ReduceObj * inoutveccast = (UFX2ReduceObj *) inoutvec;
    for (int i=0; i<*len; i++ )
        inoutveccast[i] = doReduce(inveccast[i], inoutveccast[i]); // call UFXReducer::operator()
}



vector< filedata >  GetFiles(char * filename)
{
        int64_t totalsize = 0;
        int numfiles = 0;
        vector< filedata > filesview;
        if(myrank == 0)
        {
            filedata fdata;
            ifstream allfiles(filename);
                if(!allfiles.is_open()) {
                        cerr << "Could not open " << filename << endl;
                        MPI_Abort(MPI_COMM_WORLD, 1);
                }
            allfiles.getline(fdata.filename,MAX_FILE_PATH);
            while(!allfiles.eof())
            {
                    struct stat st;
                    stat(fdata.filename, &st);
                    fdata.filesize = st.st_size;
            
                    filesview.push_back(fdata);
                    cout << filesview.back().filename << " : " << filesview.back().filesize / (1024*1024) << " MB" << endl;
                    allfiles.getline(fdata.filename,MAX_FILE_PATH);
                    totalsize += fdata.filesize;
                    numfiles++;
            }
        }
        MPI_Bcast(&numfiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&totalsize, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        if(myrank != 0)
            filesview.resize(numfiles);
    
        MPI_Datatype MPI_filedata;
        MPI_Type_contiguous(sizeof(filedata), MPI_CHAR, &MPI_filedata);
        MPI_Type_commit(&MPI_filedata);
        MPI_Bcast(filesview.data(), numfiles, MPI_filedata, 0, MPI_COMM_WORLD);
    return filesview;
}

template <typename T, typename H>
void ParallelAllReduce(SimpleCount<T, H> & counter)
{
    double t01 = MPI_Wtime();

    // merging, first reduce to power of two, then butterfly pattern
    // pow is greatest power of two smaller than size
    int pow = 1;
    int a = nprocs;
    while (a > 1){
        a >>= 1;
        pow <<= 1;
    }
#ifdef DEBUG
    if(myrank == 0) printf("pow = %d, myrank = %d\n", pow, myrank);
#endif
    if(myrank >= pow) {
        //fprintf(stderr, "Thread %d: Sending to %d (%d)\n", myrank, myrank - pow, pow);
        vector<T> keys;
        vector<int> counts;
        counter.MergeableSummary(keys, counts);
        size_t length = counts.size();
        MPI_Send(&length, 1, MPIType<size_t>(), myrank - pow, 0, MPI_COMM_WORLD);
        
        MPI_Send(keys.data(), keys.size(), MPIType<T>(), myrank - pow, 0, MPI_COMM_WORLD);
        MPI_Send(counts.data(), counts.size(), MPIType<int>(), myrank - pow, 0, MPI_COMM_WORLD);
    } else if(myrank <= nprocs - pow - 1) {
        //fprintf(stderr, "Thread %d: Receiving from %d (%d)\n", myrank, pow+myrank, pow);
        MPI_Status status;
        size_t recvsize;
        MPI_Recv(&recvsize, 1,  MPIType<size_t>(), pow + myrank, 0, MPI_COMM_WORLD, &status);
        vector<T> remotekeys(recvsize);
        vector<int> remotecounts(recvsize);
        
        MPI_Recv(remotekeys.data(), recvsize, MPIType<T>(), pow + myrank, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(remotecounts.data(), recvsize, MPIType<int>(), pow + myrank, 0, MPI_COMM_WORLD, &status);
        
        counter.Merge(remotekeys.data(), remotecounts.data(), recvsize);
        //if(myrank == 0) counter.PrintAll();
        
    }
    //fprintf(stderr, "Thread %d: Barrier (first)\n", myrank);
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    while(pow > 1) { //Butterfly time 
        if(myrank < pow / 2) {
            //fprintf(stderr, "Thread %d: Receiving from %d (%d)\n", myrank, pow/2+myrank, pow);
#ifdef DEBUG
            if(myrank  == 0) printf("myrank %d is receiving from myrank %d\n", myrank, pow / 2 + myrank);
#endif
            
            MPI_Status status;
            size_t recvsize;
            MPI_Recv(&recvsize, 1,  MPIType<size_t>(), pow / 2 + myrank, 0, MPI_COMM_WORLD, &status);
            vector<T> remotekeys(recvsize);
            vector<int> remotecounts(recvsize);
            
            MPI_Recv(remotekeys.data(), recvsize, MPIType<T>(), pow / 2 + myrank, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(remotecounts.data(), recvsize, MPIType<int>(), pow / 2 + myrank, 0, MPI_COMM_WORLD, &status);
            
            counter.Merge(remotekeys.data(), remotecounts.data(), remotekeys.size());
            //if(myrank == 0) counter.PrintAll();
        } else if(myrank < pow) {
            //fprintf(stderr, "Thread %d: Sending to %d (%d)\n", myrank, myrank - pow/2, pow);
            vector<T> keys;
            vector<int> counts;
            counter.MergeableSummary(keys, counts);
            size_t length = counts.size();
            MPI_Send(&length, 1, MPIType<size_t>(), myrank - pow / 2, 0, MPI_COMM_WORLD);
            
            MPI_Send(keys.data(), length, MPIType<T>(), myrank - pow / 2, 0, MPI_COMM_WORLD);
            MPI_Send(counts.data(), length, MPIType<int>(), myrank - pow / 2, 0, MPI_COMM_WORLD);
        }
        //fprintf(stderr, "Thread %d: Barrier (%d)\n", myrank, pow);
        pow >>= 1;  // divide pow by 2
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    if(myrank == 0) {
        vector<T> keys;
        vector<int> counts;
        counter.MergeableSummary(keys, counts);
        size_t length = counts.size();
        MPI_Bcast(&length, 1, MPIType<size_t>(), 0, MPI_COMM_WORLD);
        MPI_Bcast(keys.data(), length, MPIType<T>(), 0, MPI_COMM_WORLD);
        MPI_Bcast(counts.data(), length, MPIType<int>(), 0, MPI_COMM_WORLD);
    } else {
        size_t recvsize;
        MPI_Bcast(&recvsize, 1, MPIType<size_t>(), 0, MPI_COMM_WORLD);

        vector<T> remotekeys(recvsize);
        vector<int> remotecounts(recvsize);
        MPI_Bcast(remotekeys.data(), recvsize, MPIType<T>(), 0, MPI_COMM_WORLD);
        MPI_Bcast(remotecounts.data(), recvsize, MPIType<int>(), 0, MPI_COMM_WORLD);
        counter.Assign(remotekeys.data(), remotecounts.data(), recvsize);
    }
    //counter.PrintTop(myrank);

    double t02 = MPI_Wtime();
    if(myrank == 0)
    cout << "ParallelAllReduce done in " << t02-t01 << " seconds" << endl;

}
#endif

