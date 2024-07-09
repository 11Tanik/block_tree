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

bool check_correct = false;

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
	
	// extract leaves
	std::vector<uint8_t> leave_text = {};

  	for (size_t i = 0; i < bt->compressed_leaves_.size(); ++i) {
		uint8_t letter = bt->compressed_leaves_[i];
		leave_text.push_back(letter);
  	}

	// check text for correctness
	if (check_correct) {
		for (size_t i = 0; i < text.size(); ++i) {
			assert(bt->access(i) == text[i]);
			if (bt->access(i) != text[i]) throw std::runtime_error("block tree failed check");
		}
	}

	// calculate size
	int64_t bt_base_space = bt->print_space_usage();
	int64_t leaves_space = sdsl::size_in_bytes(bt->compressed_leaves_);
	int64_t num_leave_chars = bt->compressed_leaves_.size();

	double leave_entropy = calculate_entropy(leave_text);
	//double text_entropy = calculate_entropy(text);

	int64_t optimal_entropy_encoding = (leave_entropy * num_leave_chars) / 8.0;

	int64_t bt_entropy_space = bt_base_space - leaves_space + optimal_entropy_encoding;

	bt->huffman_compress_leaves();

	if (check_correct) {
		for (size_t i = 0; i < text.size(); ++i) {
			assert(bt->access(i) == text[i]);
			if (bt->access(i) != text[i]) throw std::runtime_error("wavelet block tree failed at access");
		}
	}

	int64_t bt_wavelet_space = bt->print_space_usage();

	std::cout << filename << ", " << tau << ", " << max_leaf_length << ", " << bt_base_space << ", " << leaves_space << ", " << optimal_entropy_encoding << ", " << bt_entropy_space << ", " << bt_wavelet_space << "\n";

  	// Clean-up
  	delete bt;
}

int32_t main()
{

	std::vector<std::string> files;
	files.push_back("./testtext.txt");
	files.push_back("./english.50MB");
	files.push_back("./dna.50MB");
	files.push_back("./dblp.xml.50MB");
	files.push_back("./proteins.50MB");
	files.push_back("./sources.50MB");
	
	std::cout << "text, tau, max_leaf_size, bt base, leaves, entropy leaves, bt entropy, bt wavelet\n";
	for (int32_t maxLS = 2; maxLS <= 64; maxLS *= 2) {
		for (int32_t tau = 2; tau <= 8; tau *= 2) {
			for (auto s : files) {
				measure_for_text(s, tau, maxLS, true);
			}
		}
	}
  
  return 0;
}

/******************************************************************************/
