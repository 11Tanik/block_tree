/*******************************************************************************
 * This file is for measuring the block trees
 ******************************************************************************/

#include <chrono>
#include <cstddef>
#include <iomanip>
#include <numeric>
#include <vector>
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
std::random_device dev;
std::mt19937 rng(dev());
const int64_t repetitions = 10000;

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

double calculate_entropy(sdsl::int_vector<> text) {
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

int64_t time_access(pasta::BlockTreeLPF<uint8_t,int32_t>* bt, size_t max) {
	rng.seed(4389762);
    std::uniform_int_distribution<std::mt19937::result_type>  access_dist(0, max-1);

	std::vector<size_t> pos(repetitions);

	for (size_t i = 0; i < max; i++) {
		pos.push_back(access_dist(rng));
	}

	const auto start_access_time = std::chrono::steady_clock::now();
	for (size_t i = 0; i < pos.size(); i++) {
		bt->access(pos[i]);
	}
	const auto end_access_time = std::chrono::steady_clock::now();
	const auto access_time = std::chrono::duration_cast<std::chrono::microseconds>(end_access_time - start_access_time);
	return access_time.count() / repetitions;
}

int64_t time_rank(pasta::BlockTreeLPF<uint8_t,int32_t>* bt, size_t max) {
	rng.seed(4389762);
    std::uniform_int_distribution<std::mt19937::result_type>  rank_dist(0, max-1);
	std::uniform_int_distribution<std::mt19937::result_type>  symbol_dist(0, 255);

	std::vector<size_t> pos(repetitions);
	std::vector<uint8_t> symbols(repetitions);

	for (size_t i = 0; i < max; i++) {
		pos.push_back(rank_dist(rng));
		symbols.push_back(symbol_dist(rng));
	}

	const auto start_rank_time = std::chrono::steady_clock::now();
	for (size_t i = 0; i < pos.size(); i++) {
		bt->rank(symbols[i],pos[i]);
	}
	const auto end_rank_time = std::chrono::steady_clock::now();
	const auto rank_time = std::chrono::duration_cast<std::chrono::microseconds>(end_rank_time - start_rank_time);
	return rank_time.count() / repetitions;

}

void measure_for_text(std::string filename, int32_t tau, int32_t max_leaf_length, bool s_equal_z) {
	// read text
	std::ifstream file(filename);
	std::ostringstream ss;
	ss << file.rdbuf();
	const std::string& s = ss.str();
	std::vector<uint8_t> text(s.begin(), s.end());

	// construct block tree
	const auto start_bt_construction_time = std::chrono::steady_clock::now();
	auto* bt = pasta::make_block_tree_lpf<uint8_t, int32_t>(text, tau, max_leaf_length, s_equal_z);
	const auto end_bt_construction_time = std::chrono::steady_clock::now();

	const auto bt_construction_time = std::chrono::duration_cast<std::chrono::microseconds>(end_bt_construction_time - start_bt_construction_time);

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

	double leave_entropy = calculate_entropy(bt->compressed_leaves_);
	double text_entropy = calculate_entropy(text);

	int64_t bt_access_time = time_access(bt, text.size());
	int64_t bt_rank_time = time_rank(bt, text.size());

	// compress leaves with wavelet tree
	auto start_bt_wt_construction_time = std::chrono::steady_clock::now();
	bt->huffman_compress_leaves();
	auto end_bt_wt_construction_time = std::chrono::steady_clock::now();

	const auto bt_wt_construction_time = std::chrono::duration_cast<std::chrono::microseconds>(end_bt_wt_construction_time - start_bt_wt_construction_time);

	if (check_correct) {
		for (size_t i = 0; i < text.size(); ++i) {
			assert(bt->access(i) == text[i]);
			if (bt->access(i) != text[i]) throw std::runtime_error("wavelet block tree failed at access");
		}
	}

	int64_t bt_wavelet_space = bt->print_space_usage();
	int64_t wt_space = sdsl::size_in_bytes(bt->wavelet_leaves);

	int64_t bt_wt_access_time = time_access(bt, text.size());
	int64_t bt_wt_rank_time = time_rank(bt, text.size());

	std::cout << filename
			<< ", " << text_entropy
			<< ", " << tau
			<< ", " << max_leaf_length
			<< ", " << bt_construction_time.count()
			<< ", " << bt_access_time
			<< ", " << bt_rank_time
			<< ", " << bt_base_space
			<< ", " << leaves_space
			<< ", " << bt_wt_construction_time.count()
			<< ", " << bt_wt_access_time
			<< ", " << bt_wt_rank_time
			<< ", " << bt_wavelet_space
			<< ", " << wt_space
			<< ", " << num_leave_chars
			<< ", " << leave_entropy
			<< "\n";

  	// Clean-up
  	delete bt;
}

int32_t main(int argc, char* argv[])
{
	if (argc != 2) throw std::runtime_error("No text given.");
	std::string filename = argv[1];
	
	std::cout << "text, entropy, tau, max_leaf_size, bt_construction, bt_access, bt_rank, bt_size, leaves_size, bt_wt_construction, bt_wt_access, bt_wt_rank, bt_wt_size, wt_size, num_leaves_chars, leaves_entropy\n";
	for (int32_t maxLS = 2; maxLS <= 128; maxLS *= 2) {
		for (int32_t tau = 2; tau <= 8; tau *= 2) {
			measure_for_text(filename, tau, maxLS, true);
		}
	}
  
  return 0;
}
