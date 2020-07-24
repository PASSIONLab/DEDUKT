#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <iostream>
#include <string>
#include <sstream>
#include <mpi.h>

#include "Friends.h"
#include "FriendsMPI.h"
using namespace std;

int myrank;
int nprocs;

void unlink_ufxfiles(char * ufxname)
{
    	stringstream ss;
    	string rank;
    	ss << myrank;
    	ss >> rank;
    	string ofilename(ufxname);
    	unsigned found = ofilename.find_last_of("/");
    	ofilename = ofilename.substr(found+1);
    	string fullname = ofilename + rank;
    	string sizefile = fullname + string(".entries");    // filename<pid>.entries

    	if(shm_unlink(sizefile.c_str()) == -1) {
		printf("Error, shm_unlink failed %d %s when trying to unlink %s\n", errno, strerror(errno), sizefile.c_str());		
    	}
    	else {
		if(myrank == 0)	cout << "shm_unlink successfully unlinked " << sizefile.c_str() << endl;
	}
    
    	if(shm_unlink(fullname.c_str()) == -1) {
		printf("Error, shm_unlink failed %d %s when trying to unlink %s\n", errno, strerror(errno), fullname.c_str());
	}
	else {
		if(myrank == 0) cout << "shm_unlink successfully unlinked " << fullname.c_str() << endl;
	}
}


void unlink_fastqfiles(char * fastqname)
{
    	stringstream ss;
    	string rank;
   	ss << myrank;
    	ss >> rank;
    	string sfastqname(fastqname);
    	unsigned found = sfastqname.find_last_of("/");
    	sfastqname = sfastqname.substr(found+1);
    	string fullfastq = sfastqname + rank;
    	if(shm_unlink(fullfastq.c_str()) == -1) {
		printf("Error, shm_unlink failed %d %s when trying to unlink %s\n", errno, strerror(errno), fullfastq.c_str());
        }
        else {
                if(myrank == 0) cout << "shm_unlink successfully unlinked " << fullfastq.c_str() << endl;
        }
}


int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	if(argc < 3)
	{
		if(myrank  == 0)
		{
			cout << "Usage: ./unlinker fastq.files ufx" << endl;
			cout << "Example: ./unlinker humans.txt humans.txt.ufx.bin" << endl;
			cout << "humans.txt is a text file that contains the paths to all the fastq files used during the meraculous runs" << endl;
		}
		MPI_Finalize(); 
		return -1;
	}
	vector<filedata> allfiles = GetFiles(argv[1]);
	for(auto itr= allfiles.begin(); itr != allfiles.end(); itr++)
	{
		unlink_fastqfiles(itr->filename);	
	}
	unlink_ufxfiles(argv[2]);
	return 0;
}
    


