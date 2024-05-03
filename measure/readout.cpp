/*******************************************************************************
 * This file is part of pasta::block_tree
 *
 * Copyright (C) 2023 Florian Kurpicz <florian@kurpicz.org>
 *
 * pasta::block_tree is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * pasta::block_tree is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with pasta::block_tree.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <string>
#include <cstdint>
#include <iostream>
#include <random>
#include <fstream>

#include <pasta/block_tree/construction/block_tree_lpf.hpp>

double calculate_entropy(std::vector<uint8_t> text) {
	std::vector<int32_t> freqs;
	freqs.resize(256);
	for (size_t i = 0; i < 256; i++) {
		freqs[i] = 0;
	}

	for (size_t i = 0; i < text.size(); i++) {
		freqs[text[i]]++;
	}
	/*
	for (size_t i = 0; i < 256; i++) {
		std::cout << (char) i << ": " << freqs[i] << "\n";
	}
	*/

	double entropy = 0;
	for (size_t i = 0; i < 256; i++) {
		if (freqs[i] == 0) continue;
		double freq = (double) freqs[i] / (double) text.size();
		entropy += freq * log2(freq);
	}

	return entropy * -1.0;
}

int32_t main()
{

	std::string filename = "./english.25MB";
	std::ifstream file(filename);
	std::ostringstream ss;
	ss << file.rdbuf();
	const std::string& s = ss.str();
	std::vector<uint8_t> text(s.begin(), s.end());

	// check:
	//std::copy(vec.begin(), vec.end(), std::ostream_iterator<char>(std::cout));

  // Generate some random text
  /*
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint8_t> dist(0, 15);

  size_t const string_length = 100;
  std::vector<uint8_t> text;
  text.resize(string_length);
  for (size_t i = 0; i < text.size(); ++i) {
    text[i] = dist(gen);
  }
  */

  // Build the block tree for the random text
  auto* bt = pasta::make_block_tree_lpf<uint8_t, int32_t>(text, 2, 1, true);

	std::vector<uint8_t> leave_text = {};
	std::vector<uint8_t> remapped_text = {};

	// Use the block tree to access individual characters of the text
  	std::cout << "# Access" << "\n";
  	for (size_t i = 0; i < text.size(); ++i) {
    	std::cout << bt->access(i) << ", ";
		remapped_text.push_back(bt->access(i));
  	}
  	std::cout << "\n";

  	std::cout << "# Leaves" << "\n";
  	for (size_t i = 0; i < bt->compressed_leaves_.size(); ++i) {
		uint8_t letter = bt->compressed_leaves_[i];
		leave_text.push_back(letter);
		std::cout << +leave_text[i] << ", ";
  	}
  	std::cout << "\n";

	std::cout << "Entropy of Leaves: " << calculate_entropy(leave_text) << "\n";
	std::cout << "Entropy of Text: " << calculate_entropy(text) << "\n";
	std::cout << "Entropy of Remapped Text: " << calculate_entropy(remapped_text) << "\n";

  	// Clean-up
  	delete bt;
  
  return 0;
}

/******************************************************************************/
