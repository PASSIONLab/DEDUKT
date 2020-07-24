/*
 * read_compression.h
 *
 *  Created on: Jun 21, 2019
 *      Author: mellis
 */

#ifndef ALIGN_READ_COMPRESSION_H_
#define ALIGN_READ_COMPRESSION_H_

#include <iostream> // cout ...
#include <cstring> // memset
#include <bitset>
#include "../common/defines.h"

class Encoder3bit {

public:
	static const int MASK = 7;
	static const int BITS_PER_LTTR = 3;
	static const int LTTRS_PER_CHUNK = 21;

	class Encoding {
		public:
			int len_orig;
			int num_chunks;
			uint64_t* chunks = NULL;
			Encoding(std::string longread) {
				this->len_orig = longread.length();
				this->num_chunks = num_uint64s(this->len_orig);
				if (this->num_chunks > 0) {
					this->chunks = new uint64_t[this->num_chunks];
					std::memset(this->chunks, 0, this->num_chunks*sizeof(uint64_t));
				}
			}
			~Encoding() {
				if (chunks != NULL) { delete chunks; }
			}
	};

	static inline int num_bits(const std::string read) { return read.length()*BITS_PER_LTTR; }
	static inline int num_uint64s(const int num_letters) { return (num_letters+LTTRS_PER_CHUNK-1)/LTTRS_PER_CHUNK; }
	static inline int num_uint64s(const std::string read) { return num_uint64s(read.length()); }
	static inline void shift_and_encode(uint64_t& encoding, const char letter) { encoding = (encoding << BITS_PER_LTTR) | (letter&MASK); }
	static inline uint8_t shift_and_decode(uint64_t& encoding) { uint8_t toreturn = encoding&MASK; (encoding >>= BITS_PER_LTTR); return toreturn; }
	Encoding encode (std::string read);
	void decode (Encoding &read, char* decoding);
};


#endif /* ALIGN_READ_COMPRESSION_H_ */
