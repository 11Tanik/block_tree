/*******************************************************************************
 * This file is for measuring the block trees
 ******************************************************************************/

#include <string>
#include <cstdint>
#include <iostream>
#include <random>
#include <fstream>
#include <iostream>
#include <bitset>
#include <assert.h>
#include <pasta/block_tree/utils/huffman.hpp>
#include <pasta/block_tree/construction/block_tree_lpf.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/int_vector.hpp>

bool check_correct = true;

double calculate_entropy(std::vector<uint8_t> text) {
	std::vector<int32_t> freqs;
	freqs.resize(256);
	for (size_t i = 0; i < 256; i++) {
		freqs[i] = 0;
	}

	for (size_t i = 0; i < text.size(); i++) {
		freqs[text[i]]++;
	}

	double entropy = 0;
	for (size_t i = 0; i < 256; i++) {
		if (freqs[i] == 0) continue;
		double freq = (double) freqs[i] / (double) text.size();
		entropy += freq * log2(freq);
	}

	return entropy * -1.0;
}

void measure_for_text(std::string filename, int32_t tau, int32_t max_leaf_length, bool s_equal_z) {
	// read text
	std::ifstream file(filename);
	std::ostringstream ss;
	ss << file.rdbuf();
	const std::string& s = ss.str();
	std::vector<uint8_t> text(s.begin(), s.end());

	// construct block tree
	auto* bt = pasta::make_block_tree_lpf<uint8_t, int32_t>(text, tau, max_leaf_length, s_equal_z);

	// check text for correctness
	if (check_correct) {
		for (size_t i = 0; i < text.size(); ++i) {
			assert(bt->access(i) == text[i]);
			if (bt->access(i) != text[i]) throw std::runtime_error("block tree failed check");
		}
	}

	// look at pointers
	/*
	std::vector<std::vector<size_t>> distrib(bt->block_tree_pointers_.size(), std::vector<size_t>(64, 0));
	size_t unnecessary_bits = 0;
	size_t num_pointers = 0;
	for (size_t l = 0; l < bt->block_tree_pointers_.size(); l++) {
		auto &lvl = *bt->block_tree_pointers_[l];
		size_t bytes = sdsl::size_in_bytes(lvl);

		size_t width = (8*bytes) / lvl.size();
		if (width > 64) std::cout << "illegal width " << width << "\n";
		for (size_t e = 0; e < bt->block_tree_pointers_[l]->size(); e++) {
			size_t ptr = lvl[e];
			size_t nb = 1;
			if (sizeof(ptr) >= 8u && (ptr >> 32u)) {nb += 32; ptr >>= 32u;}
			if (sizeof(ptr) >= 4u && (ptr >> 16u)) {nb += 16; ptr >>= 16u;}
			if (sizeof(ptr) >= 2u && (ptr >> 8u)) {nb += 8; ptr >>= 8u;}
			if (ptr >> 4u) {nb += 4; ptr >>= 4u;}
			if (ptr >> 2u) {nb += 2; ptr >>= 2u;}
			if (ptr >> 1u) {nb += 1; ptr >>= 1u;}

			if (width < nb || 0 >= nb) std::cout << "illegal nb " << nb << " for value " << lvl[e] << " and width " << width << "\n";

			distrib[l][nb-1]++;

			unnecessary_bits += (width - nb);

			num_pointers++;
		}
	}
	*/
	//double avg_per_ptr = (double) unnecessary_bits / (double) num_pointers;
	
	// extract leaves
	/*
	std::vector<uint8_t> leave_text = {};
  	for (size_t i = 0; i < bt->compressed_leaves_.size(); ++i) {
		uint8_t letter = bt->compressed_leaves_[i];
		leave_text.push_back(letter);
  	}
	*/

	// calculate size
	int64_t bt_base_space = bt->print_space_usage();
	//int64_t leaves_space = sdsl::size_in_bytes(bt->compressed_leaves_);
	//int64_t num_leave_chars = bt->compressed_leaves_.size();

	//double leave_entropy = calculate_entropy(leave_text);
	//double text_entropy = calculate_entropy(text);

	//int64_t optimal_entropy_encoding = (leave_entropy * num_leave_chars) / 8.0;

	//int64_t bt_entropy_space = bt_base_space - leaves_space + optimal_entropy_encoding;

	bt->huffman_compress_leaves();

	if (check_correct) {
		for (size_t i = 0; i < text.size(); ++i) {
			assert(bt->access(i) == text[i]);
			if (bt->access(i) != text[i]) throw std::runtime_error("wavelet block tree failed at access");
		}
	}

	int64_t bt_wavelet_space = bt->print_space_usage();

	std::cout << filename << ", " << tau << ", " << max_leaf_length << ", " << bt_base_space << ", " << bt_wavelet_space << "\n";

	/*
	for (size_t i = 0; i < distrib[0].size(); i++) {
		std::cout << filename << ", " << bt->block_tree_types_[0]->size() << ", " << tau << ", " << max_leaf_length << ", " << bt_base_space;
		for (size_t l = 0; l < bt->block_tree_pointers_.size(); l++) {
			std::cout << ", " << distrib[l][i];
		}
		std::cout << "\n";
	}
	*/

  	// Clean-up
  	delete bt;
}

int32_t main(int argc, char* argv[])
{
	if (argc != 2) throw std::runtime_error("No text given.");
	std::string filename = argv[1];
	
	std::cout << "text, tau, max_leaf_size, size, wt_size\n";
	for (int32_t maxLS = 4; maxLS <= 4; maxLS *= 2) {
		for (int32_t tau = 2; tau <= 2; tau *= 2) {
			measure_for_text(filename, tau, maxLS, true);
		}
	}
  
  return 0;
}
