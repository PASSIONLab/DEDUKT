/*
 * test_compression.cpp
 *
 *  Created on: Jun 21, 2019
 *      Author: mellis
 */
#include <cassert>
#include <bitset>
#include "read_compression.h"

using namespace std;

void test_num_bits() {
	if (Encoder3bit::num_bits("AAA") == 3*Encoder3bit::BITS_PER_LTTR) {
		cout << __FUNCTION__ << " passed!" << endl;
	}
	else { cout << __FUNCTION__ << " failed!" << endl; }
}

void test_num_uint64s() {
	if ( Encoder3bit::num_uint64s(4) == 1
		&& Encoder3bit::num_uint64s(22) == 2
		&& Encoder3bit::num_uint64s(35) == 2
		&& Encoder3bit::num_uint64s(42) == 2
		&& Encoder3bit::num_uint64s(0) == 0
	) { cout << __FUNCTION__ << " passed!" << endl; }
	else { cout << __FUNCTION__ << " failed!" << endl; }
}

void test_encode_decode_one_helper(char letter, int expected) {
	uint64_t code = 0;
	Encoder3bit::shift_and_encode(code, letter);
	cout << bitset<64>(code) << " is the coding of "<< letter << endl;
	uint8_t decode = Encoder3bit::shift_and_decode(code);
	cout << "  " << bitset<8>(decode) << " is the decoding of "<< letter << endl;
	if (decode != expected) {
		cout << "  Failed on " << letter << endl;
	}
	else {
		cout << "  Passed on " << letter << endl;
	}
}

void test_encode_decode_one() {
	cout << "Testing shift_and_encode and shift_and_decode" << endl;
	test_encode_decode_one_helper('A', 1);
	test_encode_decode_one_helper('C', 3);
	test_encode_decode_one_helper('T', 4);
	test_encode_decode_one_helper('N', 6);
	test_encode_decode_one_helper('G', 7);
}

void test_single_chunk_encode(Encoder3bit* encoder) {
	string test_read = "GGGGGGGGGGGGGGGGGGGGG";
	bitset<64> expected("0111111111111111111111111111111111111111111111111111111111111111");
	Encoder3bit::Encoding code = encoder->encode(test_read);
	bitset<64> result(code.chunks[0]);
	if (code.len_orig == 21 && code.num_chunks == 1 && result == expected) {
		cout << __FUNCTION__ << " passed!" << endl;
	}
	else {
		cout << __FUNCTION__ << " failed!" << endl;
	}
}

void test_two_chunk_encode(Encoder3bit* encoder) {
	string test_read = "GGGGGGGGGGGGGGGGGGGGGG";
	bitset<64> expected("0000000000000000000000000000000000000000000000000000000000000111");
	Encoder3bit::Encoding code = encoder->encode(test_read);
	bitset<64> result(code.chunks[1]);
	if (code.len_orig == 22 && code.num_chunks == 2 && result == expected) {
		cout << __FUNCTION__ << " passed!" << endl;
	}
	else {
		cout << __FUNCTION__ << " failed!" << endl;
	}
}

/*
 * something broken with this test...
 */
void test_single_chunk_decode(Encoder3bit* decoder) {
	cout << "Testing " << __FUNCTION__ << endl;
	string expected = "GGGGGGGGGGGGGGGGGGGGG";
	uint64_t test_code = 9223372036854775807;
	cout << "  test code is: " << bitset<64>(test_code) << endl;
	Encoder3bit::Encoding t(expected);
	t.chunks[0] = test_code;
	char result[t.len_orig+1];
	decoder->decode(t, result);
	if (result == expected) {
		cout << __FUNCTION__ << " passed!" << endl;
	}
	else {
		cout << __FUNCTION__ << " failed!" << endl;
		cout << result << " was result" << endl;
		cout << expected << " was expected" << endl;
	}
}

void test_two_chunk_decode(Encoder3bit* decoder) {
	cout << "Testing " << __FUNCTION__ << endl;
	string expected = "GGGGGGGGGGGGGGGGGGGGGG";

	uint64_t test_code_one = 9223372036854775807;
	uint64_t test_code_two = 7;
	cout << "  test code one is: " << bitset<64>(test_code_one) << endl;
	cout << "  test code two is: " << bitset<64>(test_code_two) << endl;

	Encoder3bit::Encoding t(expected);
	t.chunks[0] = test_code_one;
	t.chunks[1] = test_code_two;

	char result[t.len_orig+1];
	decoder->decode(t, result);
	if (result == expected) {
		cout << __FUNCTION__ << " passed!" << endl;
	}
	else {
		cout << __FUNCTION__ << " failed!" << endl;
		cout << result << " was result" << endl;
		cout << expected << " was expected" << endl;
	}
}

void test_encode_and_decode(Encoder3bit* encoder) {
	string test = "AAAAGTCAAGGCCNAAAT";
	Encoder3bit::Encoding t = encoder->encode(test);
	char decoded[t.len_orig+1];
	encoder->decode(t, decoded);
	if (test.compare(decoded) == 0) {
		cout << __FUNCTION__ << " passed!" << endl;
	}
	else {
		cout << __FUNCTION__ << " failed!" << endl;
		cout << test << " was test string " << endl;
		cout << decoded << " was decoded result" << endl;
		cout << "bits are: " << endl;
		for (int i = 0; i < t.num_chunks; i++) {
			cout << bitset<64>(t.chunks[i]) << "  ";
		}
		cout << endl;
	}
}

int main(void) {
	test_num_bits();
	test_num_uint64s();
	test_encode_decode_one();
	Encoder3bit* encoder = new Encoder3bit();
	test_single_chunk_encode(encoder);
	test_two_chunk_encode(encoder);
	test_single_chunk_decode(encoder);
	test_two_chunk_decode(encoder);
	test_encode_and_decode(encoder);
	delete encoder;
	return 0;
}


