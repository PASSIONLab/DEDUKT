/*
 * read_compression.cpp
 *
 *  Created on: Jun 21, 2019
 *      Author: mellis
 */
#include "read_compression.h"

Encoder3bit::Encoding Encoder3bit::encode(std::string longread) {
	ASSERT(!longread.empty(), "should not be trying to encode an empty string!");
	Encoding enc = Encoding(longread);
	int str_i = 0;
	for (int i = 0; i < enc.num_chunks; i++) {
		for (int j = 0; j < LTTRS_PER_CHUNK && str_i < enc.len_orig; j++, str_i++) { // stuff 21 letters into last 63 bits
			shift_and_encode(enc.chunks[i], longread[str_i]);
		}
	}
	return enc;
}

void Encoder3bit::decode(Encoder3bit::Encoding& enc_read, char* decoding) {
	int str_i = enc_read.len_orig-1;
	uint8_t switch_code = 0;
	for (int i = 0; i < enc_read.num_chunks; i++) {
		for (int j = 0; j < LTTRS_PER_CHUNK && str_i >= 0; j++, str_i--) {
			switch_code = shift_and_decode(enc_read.chunks[i]);
			switch(switch_code) {
				case 1: { decoding[str_i] = 'A'; break; } // A, 001
				case 3: { decoding[str_i] = 'C'; break; } // C, 011
				case 4: { decoding[str_i] = 'T'; break; } // T, 100
				case 6: { decoding[str_i] = 'N'; break; } // N, 110
				case 7: { decoding[str_i] = 'G'; break; } // G, 111
				default: // oh noes!
				{
					ASSERT(0, "We tried to decode something we weren't supposed to read or the data is corrupt!\n");
				}
			}
		}
	}
	decoding[enc_read.len_orig] = '\0';
}



