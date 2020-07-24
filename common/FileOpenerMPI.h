/*
 * FileOpenerMPI.h
 *
 *  Created on: Nov 21, 2017
 *      Author: mellis
 */

#ifndef COMMON_FILEOPENERMPI_H_
#define COMMON_FILEOPENERMPI_H_

#include <vector>

using namespace std;

struct filenamesize
{
    char filename[MAX_FILE_PATH];
    size_t filesize;
};

vector< filenamesize >  loadAllFiles(char * filename);



#endif /* COMMON_FILEOPENERMPI_H_ */
