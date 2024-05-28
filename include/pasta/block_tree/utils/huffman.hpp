# pragma once

#include <vector>
#include <queue>

// TODO possibly change code_length to missing_bits

namespace pasta {

template <typename input_type, typename size_type> class HuffmanCode {
public:

	uint8_t max_code_length;
	std::vector<uint64_t> bits;
	
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

		uint8_t get_depth() {
			uint8_t d1 = (left_child == nullptr) ? 0 : left_child->get_depth();
			uint8_t d2 = (right_child == nullptr) ? 0 : right_child->get_depth();
			return (d1 > d2) ? (d1 + 1) : (d2 + 1);
		}
	};

	class HuffmanNodeCompare {
	public:
    	int operator() (HuffmanNode* a, HuffmanNode* b) {
        	return a->occs - b->occs;
    	}
	};

	HuffmanCode(std::vector<input_type> &text) {
		if (text.empty()) return;
		create_from_text(text);
	}

	void create_from_text(std::vector<input_type> &text) {
		// use encode_table for storing frequencies
		encode_table.resize(256,0);
		std::fill(encode_table.begin(), encode_table.end(), 0);

		for (size_type i = 0; i < text.size(); i++) {
			encode_table[text[i]]++;
		}

		std::priority_queue<HuffmanNode*, std::vector<HuffmanNode*>, HuffmanNodeCompare> pq;

		for (size_t i = 0; i < encode_table.size(); i++) {
			HuffmanNode n = new HuffmanNode();
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
			pq.top();

			HuffmanNode* parent = new HuffmanNode();

			parent->occs = n1->occs + n2->occs;
			parent->left_child = n1;
			parent->right_child = n2;

			pq.push(parent);
		}

		HuffmanNode* root = pq.top();

		// find depth for longest bit length
		max_code_length = root->get_depth()-1;
		// calculate size of decode table
		decode_table.resize(1 << max_code_length);

		traverse_tree(root);
	}

	void traverse_tree(HuffmanNode* root, uint64_t current_code, uint8_t current_bit_length) {

		// base case
		if (root->left_child == nullptr && root->right_child == nullptr) {

			for (size_t i = 0; i < (1 << (max_code_length - current_bit_length)); i++) {
				decode_table[current_code + i].code_word_length = current_bit_length;
				decode_table[current_code + i].letter = root->letter;
			}

			encode_table[root->letter] = current_code;

		} else {
			assert(root->left_child != nullptr && root->right_child != nullptr);

			traverse_tree(root->left_child, current_code, current_bit_length+1);

			traverse_tree(root->right_child, current_code & (1 << (max_code_length - current_bit_length - 1)), current_bit_length+1);
		}

	}

	std::vector<input_type>* decode(size_type start_index, size_type num_elements) {
		return nullptr;
	}

	int64_t print_space_usage() {
		return 0;
	}

};
}
