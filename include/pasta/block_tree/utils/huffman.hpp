# pragma once

#include <vector>
#include <queue>
#include <assert.h>

// TODO possibly change code_length to missing_bits

namespace pasta {

template <typename input_type, typename size_type> class HuffmanCode {
public:

	uint8_t max_code_length;
	size_type bit_length;
	std::vector<uint64_t> bits;

	size_type sample_pos;
	std::vector<size_type> samples;
	
	class HuffmanEntry {
	public:
		input_type letter;
		uint8_t code_word_length;
	};

	std::vector<uint64_t> encode_table;
	std::vector<HuffmanEntry> decode_table;

	class HuffmanNode {
	public:
		input_type letter;
		size_type occs;
		HuffmanNode* left_child;
		HuffmanNode* right_child;

		~HuffmanNode() {
			delete left_child;
			delete right_child;
		}

		uint8_t get_depth() {
			uint8_t d1 = (left_child == nullptr) ? 0 : left_child->get_depth();
			uint8_t d2 = (right_child == nullptr) ? 0 : right_child->get_depth();
			return (d1 > d2) ? (d1 + 1) : (d2 + 1);
		}
	};

	class HuffmanNodeCompare {
	public:
    	bool operator() (HuffmanNode* a, HuffmanNode* b) {
			if (a->occs > b->occs) return true;
			return false;
    	}
	};

	HuffmanCode(std::vector<input_type> &text, size_type sample_pos) {
		if (text.empty()) return;
		this->sample_pos = sample_pos;
		create_from_text(text);
	}

	void create_from_text(std::vector<input_type> &text) {
		// use encode_table for storing frequencies
		encode_table.resize(256,0);
		std::fill(encode_table.begin(), encode_table.end(), 0);

		for (size_type i = 0; i < text.size(); i++) {
			assert(text[i] < 256); // only small alphabets
			encode_table[text[i]]++;
		}

		std::priority_queue<HuffmanNode*, std::vector<HuffmanNode*>, HuffmanNodeCompare> pq;

		for (size_t i = 0; i < encode_table.size(); i++) {
			if (encode_table[i] == 0) continue;
			HuffmanNode* n = new HuffmanNode();
			n->letter = i;
			n->occs = encode_table[i];
			n->left_child = nullptr;
			n->right_child = nullptr;
			pq.push(n);
		}

		// reset encode table
		std::fill(encode_table.begin(), encode_table.end(), 0);

		while (pq.size() > 1) {
			HuffmanNode* n1 = pq.top();
			pq.pop();

			HuffmanNode* n2 = pq.top();
			pq.pop();

			HuffmanNode* parent = new HuffmanNode();

			parent->occs = n1->occs + n2->occs;
			parent->left_child = n1;
			parent->right_child = n2;

			pq.push(parent);
		}

		HuffmanNode* root = pq.top();

		// find depth for longest bit length
		max_code_length = root->get_depth()-1;

		if (max_code_length > 64) throw std::invalid_argument("Huffman code would require more than 64 bits!");

		// calculate size of decode table
		decode_table.resize(1 << max_code_length);

		// build encode and decode tables
		traverse_tree(root,0,0);

		delete root; // recursively deletes the entire tree

		// write bits

		bits.push_back(0);
		size_type current_bit_offset = 0;

		for (size_type i = 0; i < text.size(); i++) {

			if (sample_pos > 0 && (i % sample_pos == 0)) {
				samples.push_back((bits.size()-1)*64 + current_bit_offset);
			}

			uint64_t code = encode_table[text[i]];
			size_type length = decode_table[code].code_word_length;

			uint64_t left_code = code << (64 - max_code_length);

			bits.back() |= left_code >> current_bit_offset;

			if (64-current_bit_offset < length) {
				bits.push_back(0);
				bits.back() |= left_code << (64 - current_bit_offset);
				current_bit_offset = length - (64-current_bit_offset);
			} else {
				current_bit_offset += length;
			}
		}

		// Done
	}

	void traverse_tree(HuffmanNode* root, uint64_t current_code, uint8_t current_bit_length) {

		// base case
		if (root->left_child == nullptr && root->right_child == nullptr) {

			for (size_t i = 0; i < (1ul << (max_code_length - current_bit_length)); i++) {
				decode_table[current_code + i].code_word_length = current_bit_length;
				decode_table[current_code + i].letter = root->letter;
			}

			encode_table[root->letter] = current_code;

		} else {
			assert(root->left_child != nullptr && root->right_child != nullptr);
			assert(current_bit_length < max_code_length);

			traverse_tree(root->left_child, current_code, current_bit_length+1);

			traverse_tree(root->right_child, current_code | (1ull << (max_code_length - current_bit_length - 1)), current_bit_length+1);
		}

	}

	std::vector<input_type>* decode(size_type start_index, size_type num_elements) {
		
		std::vector<input_type>* decoded = new std::vector<input_type>();

		size_type block_index = start_index / 64;
		size_type bit_index = start_index % 64;

		while (num_elements > 0) {
			uint64_t block = bits[block_index];
			uint64_t symbol = block << bit_index;

			uint64_t next_block = 0;
			if (block_index < bits.size()-1) next_block = bits[block_index+1];

			symbol |= next_block >> (64-bit_index);

			symbol = symbol >> (64 - max_code_length);

			decoded->push_back(decode_table[symbol].letter);

			bit_index += decode_table[symbol].code_word_length;
			if (bit_index >= 64) block_index++;
			bit_index = bit_index % 64;

			num_elements--;
		}

		return decoded;
	}

	std::vector<input_type>* access(size_type index, size_type num_elements) {
		size_type closest_sample = index / sample_pos;

		size_type extended_num_elements = num_elements + (index % sample_pos);

		std::vector<input_type>* output = decode(samples[closest_sample], extended_num_elements);

		if (index % sample_pos == 0) return output;

		return new std::vector<input_type>(output->end()-num_elements, output->end());
	}

	int64_t print_space_usage() {
		return -1;
	}

};
}
