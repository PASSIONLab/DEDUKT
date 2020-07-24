/*
 * seqan_kernel.h
 *
 *  Created on: Dec 2, 2017
 *      Author: mellis
 */

#ifndef SEQAN_INTERFACE_H_
#define SEQAN_INTERFACE_H_
//#include "Aligner.h"

#include <string>

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/seeds.h>


// TODO parameterize scoring scheme
// TODO include scoring scheme and gapped (versus ungapped) info. in output
seqan::Score<int, seqan::Simple> scoringScheme(1, -1, -1); // per Giulia's specification
const int dropThreshold = 7;

int seqan_align(string read1, string read2) {
	seqan::Align<seqan::String<char> > algn_obj;
	seqan::resize(seqan::rows(algn_obj), 2);
	seqan::assignSource(seqan::row(algn_obj, 0), read1);
	seqan::assignSource(seqan::row(algn_obj, 1), read2);

	/*
	 * See http://seqan.readthedocs.io/en/master/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html#local-alignments
	 * and http://seqan.readthedocs.io/en/master/Tutorial/DataStructures/Alignment/ScoringSchemes.html#linear-gap-model
	 * and http://docs.seqan.de/seqan/2.1.0/group_AlignmentAlgorithmTags.html#AlignmentAlgorithmTags%23DynamicGaps
	 */
	seqan::Score<int> scoringScheme(3, -3, -2, -2);
	int score = seqan::localAlignment(algn_obj, scoringScheme, seqan::DynamicGaps());

	return score;
}

double seqan_seedExtend(const string read1, const string read2, const PosInRead pos1, const PosInRead pos2,
		int kmer_len, align::alignment_data& data)
{
	double seedextend_time;
	int longest_extension;

	seqan::Dna5String seqH = read1;
	seqan::Dna5String seqV = read2;

	seqan::Seed<seqan::Simple> seed(pos1, pos2, pos1 + kmer_len, pos2 + kmer_len);
	seqan::Dna5String seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
	seqan::Dna5String seedV = infix(seqV, beginPositionV(seed), endPositionV(seed));

	seqan::Dna5StringReverseComplement twin(seedH);

	if(twin == seedV) {
		seqan::Dna5StringReverseComplement twinRead(seqH);
		PosInRead pos1_reverse = read1.length()-pos1-kmer_len;
		seqan::Seed<seqan::Simple> seed2(pos1_reverse, pos2, pos1_reverse + kmer_len, pos2 + kmer_len);
		seedextend_time = MPI_Wtime();
		longest_extension = seqan::extendSeed(seed2, twinRead, seqV, seqan::EXTEND_BOTH, scoringScheme,
				dropThreshold, seqan::GappedXDrop());
		seedextend_time = MPI_Wtime() - seedextend_time;
		data.reset_positions(seqan::beginPositionH(seed2), seqan::endPositionH(seed2), seqan::beginPositionV(seed2), seqan::endPositionV(seed2));
	} else {
		seedextend_time = MPI_Wtime();
		longest_extension = seqan::extendSeed(seed, seqH, seqV, seqan::EXTEND_BOTH, scoringScheme,
				dropThreshold, seqan::GappedXDrop());
		seedextend_time = MPI_Wtime() - seedextend_time;
		data.reset_positions(seqan::beginPositionH(seed), seqan::endPositionH(seed), seqan::beginPositionV(seed), seqan::endPositionV(seed));
	}

	data.length_r1 = read1.length();
	data.length_r2 = read2.length();
	data.score = longest_extension;

	return seedextend_time;
}




#endif /* SEQAN_INTERFACE_H_ */
