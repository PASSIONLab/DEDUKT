/*
 * edlib_kernel.h
 *
 *  Created on: Dec 2, 2017
 *      Author: mellis
 */

#ifndef EDLIB_INTERFACE_H_
#define EDLIB_INTERFACE_H_

#include "../../edlib/include/edlib.h"

/*
 * code adapted from longreads/mtspgemm2017/multiply.h
 */
#define pCORR .7225
#define pWRON .2275
#define minLEN 17
/* Compute potential overlap */
void getOverlapLen(PosInRead & posInA, PosInRead & posInB, std::string & seqA, std::string seqB, int & len, int & left)
{
    /* Left margin passed as argument */
    int right;           /* Define the right margin */
    int overlaplength;   /* Expected overlap region length */
    int temp1, temp2;

    if (posInA <= posInB) { left = posInA; }
    else { left = posInB; }

    temp1 = seqA.length() - posInA;   /* 1st value referred to A's  */
    temp2 = seqB.length() - posInB;  /* 2nd value referred to B's  */

    if (temp1 <= temp2) { right = temp1; }
    else { right = temp2; }

    len = left+right; /* Estimated overlap */

}

/* EDLIB Local Alignment */
/*
 * EdlibAlignResult stores alignment result
 * edlibAlignConfig needs to be set-up
 * edlibAlign(char * read1, int lenRead1, char * read2, int lenRead2, const EdlibConfig config)
 */
bool edlibOp(PosInRead & posInA, PosInRead & posInB, std::string & seqA, std::string & seqB, int & len, int & ed)
{
    bool align;
    int left = 0;
    len=0;
    /* Compute the overlap potential length and obtain substring of
    overlaplength from the considered pair and pass them to alignop */
    getOverlapLen(posInA, posInB, seqA, seqB, len, left);
    /* Obtain row and col substr of overlaplength */
    std::string substrA = seqA.substr(posInA-left, len);
    std::string substrB = seqB.substr(posInB-left, len);

    if(len >= minLEN)
    {
    		// TODO: why bother? commented-out process produces the same ACTG substring right? unless reversal happens? */
        Kmer kmerfromA;
        Kmer kmerfromB;

        kmerfromA.set_kmer(substrA.c_str());
        kmerfromB.set_kmer(substrB.c_str());

        kmerfromA = kmerfromA.rep();
        kmerfromB = kmerfromB.rep();

        char *read1 = new char[len+1];
        char *read2 = new char[len+1];

        kmerfromA.toString(read1);
        kmerfromB.toString(read2);

        /* In len we expect len*pow(pCORR, 2) correct base-pairs, so here we compute the maximum
        edit distance as len - mean + 2 standard deviation = len - (len*pow(pCORR, 2) + 2 * qrt(len*pow(pCORR, 2)*(1-pow(pCORR, 2)) */
		int maxed = len-len*pCORR+2*sqrt(len*pCORR*pWRON);

        EdlibAlignResult result = edlibAlign(read1, len, read2, len, edlibNewAlignConfig(maxed, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
        delete [] read1;
        delete [] read2;

        if(result.editDistance != -1)
        {
            align = true;
            ed = result.editDistance;
        }
        else align = false;

        edlibFreeAlignResult(result);
    }
    else { align = false; }

    return align;
}



#endif /* EDLIB_INTERFACE_H_ */
