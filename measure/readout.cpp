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

bool check_correct = false;
bool short_output = true;

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
	if (!short_output) std::cout << "Building block tree for " << filename << " ...\n";
	auto* bt = pasta::make_block_tree_lpf<uint8_t, int32_t>(text, tau, max_leaf_length, s_equal_z);
	if (!short_output) std::cout << "... Done\n";
	
	std::vector<uint8_t> leave_text = {};

	if (!short_output) std::cout << "Extracting leaves ...\n";
  	for (size_t i = 0; i < bt->compressed_leaves_.size(); ++i) {
		uint8_t letter = bt->compressed_leaves_[i];
		leave_text.push_back(letter);
  	}
	if (!short_output) std::cout << "... Done\n";

	// check text for correctness
	if (check_correct) {
		if (!short_output) std::cout << "Checking correctness ...\n";
		for (size_t i = 0; i < text.size(); ++i) {
			assert(bt->access(i) == text[i]);
		}
		if (!short_output) std::cout << "... Done\n";
	}

	// calculate size

	int64_t total_space = bt->print_space_usage();
	int64_t leave_space = sdsl::size_in_bytes(bt->compressed_leaves_);
	int64_t num_leave_chars = bt->compressed_leaves_.size();
	int64_t leave_space_bits = bt->compressed_leaves_.bit_size();
	int64_t char_width_in_leaf = leave_space_bits / num_leave_chars;
	int64_t leave_bit_size = (char_width_in_leaf * num_leave_chars) / 8.0;

	if(!short_output) {
		std::cout << "Total space of block tree: " << total_space << "\n";
		std::cout << "Space of leaves: " << leave_space << "\n";
		std::cout << "Space of optimally aligned leaves: " << leave_bit_size << "\n";
		std::cout << "Bits per character in leaves: " << char_width_in_leaf << "\n";
		std::cout << "Number of characters in leaves: " << num_leave_chars << "\n";
	}

	double leave_entropy = calculate_entropy(leave_text);
	double text_entropy = calculate_entropy(text);

	if (!short_output) {
		std::cout << "Entropy of leaves: " << leave_entropy << "\n";
		std::cout << "Entropy of text: " << text_entropy << "\n";
	}

	int64_t optimal_entropy_encoding = (leave_entropy * num_leave_chars) / 8.0;

	if (!short_output) {
		std::cout << "Space of leaves with optimal entropy encoding: " << optimal_entropy_encoding << "\n";
		std::cout << "------------------------\n";
	}

	int64_t max_entropy_savings = leave_space - optimal_entropy_encoding;
	double max_percent_savings = 100.0 * ((double) max_entropy_savings) / ((double) total_space);

	if (!short_output) {
		std::cout << "Maximum space reduction achievable with entropy encoding:\n";
		std::cout << "Absolute: " << max_entropy_savings << " bytes\n";
		std::cout << "Ratio: " << max_percent_savings << " %\n";
	}

	int64_t final_space = total_space - leave_space + optimal_entropy_encoding;

	if (!short_output) {
		std::cout << "Final space: " << final_space << "\n";

		std::cout << "========================\n";
	}

	bt->huffman_compress_leaves();

	if (check_correct) {
		for (size_t i = 0; i < text.size(); ++i) {
			assert(bt->access(i) == text[i]);
			if (bt->access(i) != text[i]) throw std::runtime_error("huffman block tree failed check");
		}
	}

	int64_t praxis_huffman_space = bt->print_space_usage();
	int64_t praxis_optimal_huffman = praxis_huffman_space - bt->huffman_compressed_leaves->print_space_usage() + bt->huffman_compressed_leaves->print_space_usage_without_sampels();

	std::cout << filename << ", " << tau << ", " << max_leaf_length << ", " << total_space << ", " << leave_space << ", " << optimal_entropy_encoding << ", " << final_space << ", " << praxis_huffman_space << ", " << praxis_optimal_huffman << ", " << bt->huffman_compressed_leaves->max_sampled_span << ", " << bt->huffman_compressed_leaves->decode_table.size() << "\n";

  	// Clean-up
  	delete bt;
}

int32_t main()
{

	std::vector<std::string> files;
	//files.push_back("./english.50MB");
	files.push_back("./dna.50MB");
	//files.push_back("./dblp.xml.50MB");
	//files.push_back("./proteins.50MB");
	//files.push_back("./sources.50MB");
	
	std::cout << "text, tau, max_leaf_size, total space, leaf space, entropy encoded leaf space, total entropy encoded space, current implementation space, achievable space, longest huffman encoded leaf, #elements in decode_table\n";
	for (int32_t maxLS = 2; maxLS <= 32; maxLS *= 2) {
		for (int32_t tau = 2; tau <= 8; tau *= 2) {
			for (auto s : files) {
				measure_for_text(s, tau, maxLS, true);
			}
		}
	}
	
	/*
	std::ifstream file("./dna.50MB");
	std::ostringstream ss;
	ss << file.rdbuf();
	const std::string& s = ss.str();
	std::vector<uint8_t> text(s.begin(), s.end());

	size_t sample_pos = 0;

	pasta::HuffmanCode<uint8_t,size_t>* code = new pasta::HuffmanCode<uint8_t,size_t>(text,sample_pos);

	for(size_t i = 0; i < code->decode_table.size(); i++) {
		std::cout << "Index " << i << ": " << code->decode_table[i].letter << " with bit length " << (int) code->decode_table[i].code_word_length << "\n";
	}

	std::cout << "Space of HuffmanCode without sampels: " << code->print_space_usage_without_sampels() << "\n";

	std::vector<uint8_t>* output = code->decode(0, text.size());

	for (size_t i = 0; i < text.size(); i++) {
		if (text[i] != output->at(i)) throw std::runtime_error("huffman failed");
	}
	*/
  
  return 0;
}

/******************************************************************************/
