
#include <iostream>
#include <vector>
#include <bit_vector.hpp>
#include <support/find_l2_flat_with.hpp>
#include <support/flat_rank.hpp>
#include <support/flat_rank_select.hpp>
#include <support/optimized_for.hpp>
#include <support/rank.hpp>
#include <support/rank_select.hpp>
#include <support/wide_rank.hpp>
#include <support/wide_rank_select.hpp>
#include <omp.h>
#include <sdsl/int_vector.hpp>
#ifndef BLOCK_TREE_BV_BLOCKTREE_HPP
#define BLOCK_TREE_BV_BLOCKTREE_HPP

template<typename input_type, typename size_type>
class BV_Block_Tree {
public:
    bool CUT_FIRST_LEVELS = true;
    size_type tau_;
    size_type max_leaf_length_;
    size_type s_ = 1;
    size_type leaf_size = 0;
    size_type amount_of_leaves = 0;
    bool rank_support = false;
    std::vector<pasta::BitVector*> block_tree_types_;
    std::vector<pasta::RankSelect<pasta::OptimizedFor::ONE_QUERIES>*> block_tree_types_rs_;
    std::vector<sdsl::int_vector<>*> block_tree_pointers_;
    std::vector<sdsl::int_vector<>*> block_tree_offsets_;
    std::vector<sdsl::int_vector<>*> block_tree_encoded_;
    std::vector<int64_t> block_size_lvl_;
    std::vector<int64_t> block_per_lvl_;
    std::vector<input_type> leaves_;
    std::unordered_map<input_type, size_type> chars_index_;
    std::vector<input_type> chars_;
    size_type u_chars_;
    std::vector<std::vector<int64_t>> top_level_c_ranks_;
    std::vector<std::vector<sdsl::int_vector<>>> c_ranks_;
    std::vector<std::vector<sdsl::int_vector<>>> pointer_c_ranks_;
    int64_t access_encoded(size_type index) {
        int64_t block_size = block_size_lvl_[0];
        int64_t blk_pointer = index / block_size;
        int64_t off = index % block_size;
        int64_t child = 0;
        for (size_type i = 0; i < block_tree_types_.size(); i++) {
            auto& lvl = *block_tree_types_[i];
            auto& lvl_rs = *block_tree_types_rs_[i];
            auto& lvl_enc = *block_tree_encoded_[i];
            if (lvl[blk_pointer] == 0) {
                size_type blk = lvl_rs.rank0(blk_pointer);
                size_type enc = lvl_enc[blk];
                blk_pointer = enc / block_size;
                off += enc % block_size;
                if (off >= block_size) {
                    blk_pointer++;
                    off -= block_size;
                }
            }
            block_size /= tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer = lvl_rs.rank1(blk_pointer) * tau_ + child;
        }
        return leaves_[blk_pointer * leaf_size + off];
    }
    int64_t access(size_type index) {
        int64_t block_size = block_size_lvl_[0];
        int64_t blk_pointer = index / block_size_lvl_[0];
        int64_t off = index % block_size_lvl_[0];
        int64_t child;
        for (size_type i = 0; i < block_tree_types_.size(); i++) {
            auto& lvl = *block_tree_types_[i];
            auto& lvl_rs = *block_tree_types_rs_[i];
            auto& lvl_ptr = *block_tree_pointers_[i];
            auto& lvl_off = *block_tree_offsets_[i];
            if (lvl[blk_pointer] == 0) {
                size_type blk = lvl_rs.rank0(blk_pointer);
                off = off + lvl_off[blk];
                blk_pointer = lvl_ptr[blk];
                if (off >= block_size) {
                    blk_pointer++;
                    off -= block_size;
                }
            }
            block_size /= tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer = lvl_rs.rank1(blk_pointer) * tau_ + child;
        }
        return leaves_[blk_pointer * leaf_size + off];
    };
    int64_t select(input_type c, size_type j) {

        auto c_index = chars_index_[c];
        auto& top_level = *block_tree_types_[0];

        auto& top_level_rs = *block_tree_types_rs_[0];
        auto& top_level_ptr = *block_tree_pointers_[0];
        auto& top_level_off = *block_tree_offsets_[0];
        size_type current_block = (j - 1) / block_size_lvl_[0];
        size_type end_block = c_ranks_[c_index][0].size() - 1;
        int64_t block_size = block_size_lvl_[0];
        // find first level block containing the jth occurrence of c with a bin search
        while (current_block != end_block) {
            size_type m = current_block + (end_block-current_block)/2;
            size_type f = (m == 0)? 0 : c_ranks_[c_index][0][m - 1];
            if (f < j) {
                if (end_block - current_block == 1) {
                    if (c_ranks_[c_index][0][m] < j) {
                        current_block = m + 1;
                    }
                    break;
                }
                current_block = m;
            } else {
                end_block = m-1;
            }
        }

        // accumulator
        int64_t s = current_block * block_size - 1;
        // index that indicates how many c's are still unaccounted for
        j -= (current_block == 0) ? 0 : c_ranks_[c_index][0][current_block - 1];
        // we translate unmarked blocks on the top level independently as it differs from the other levels
        if (!top_level[current_block]) {
            int64_t blk = top_level_rs.rank0(current_block);
            current_block = top_level_ptr[blk];
            int64_t g = top_level_off[blk];
            int64_t rank_d = (current_block == 0)? c_ranks_[c_index][0][0] : c_ranks_[c_index][0][current_block] - c_ranks_[c_index][0][current_block - 1];
            rank_d -= pointer_c_ranks_[c_index][0][blk];
            if (rank_d < j) {
                j -= rank_d;
                s += (block_size - g);
                current_block++;
            } else {
                j += pointer_c_ranks_[c_index][0][blk];
                s -= g;

            }
        }
        int64_t i = 1;
        while (i < block_tree_types_.size()) {
            auto& current_level = *block_tree_types_[i];
            auto& current_level_rs= *block_tree_types_rs_[i];
            auto& current_level_ptr= *block_tree_pointers_[i];
            auto& current_level_off= *block_tree_offsets_[i];
            auto& prev_level_rs= *block_tree_types_rs_[i - 1];
            current_block = prev_level_rs.rank1(current_block) * tau_;
            block_size /= tau_;
            int64_t k = current_block;
            while ((int64_t)c_ranks_[c_index][i][current_block] < j) {
                current_block++;
            }
            j -= (current_block == k)? 0 : c_ranks_[c_index][i][current_block - 1];
            s += (current_block - k) * block_size;
            if (!current_level[current_block]) {
                int64_t blk = current_level_rs.rank0(current_block);
                current_block = current_level_ptr[blk];
                int64_t g = current_level_off[blk];
                int64_t rank_d = (current_block % tau_ == 0)? c_ranks_[c_index][i][current_block] : c_ranks_[c_index][i][current_block] - c_ranks_[c_index][i][current_block - 1];
                rank_d -= pointer_c_ranks_[c_index][i][blk];
                if (rank_d < j) {
                    j -= rank_d;
                    s += (block_size - g);
                    current_block++;
                } else {
                    j += pointer_c_ranks_[c_index][i][blk];
                    s -= g;

                }
            }
            i++;
        }

        current_block = (*block_tree_types_rs_[i - 1]).rank1(current_block) * tau_;
        int64_t l = 0;
        while (j > 0) {
            if (leaves_[current_block * leaf_size + l] == c) j--;
            l++;
        }
        return s + l;
//        size_type h = current_block * block_size_lvl_[0];
//        size_type first_lvl_prefix = (current_block == 0) ? 0 : c_ranks_[chars_index_[c]][0][current_block - 1];
//        size_type current_length = block_size_lvl_[0];
//        size_type i = 0;
//        while (i < block_tree_types_.size() - 1) {
//            if ((*block_tree_types_[i])[current_block]) {
//                size_type first_child = (*block_tree_types_rs_[i]).rank1(current_block) * tau_;
//                size_type child = first_child;
//                int last_possible_child = (first_child + tau_ - 1 > c_ranks_[chars_index_[c]][i].size() - 1) ?  c_ranks_[chars_index_[c]][i].size() - 1 : first_child + tau_ -1;
//                while (child < last_possible_child && j > c_ranks_[chars_index_[c]][i + 1][current_block - 1]) {
//                    child++;
//                }
//
//            }
//
//        }
//        return 2;
        return current_block;
    }
    int64_t rank_base(input_type c, size_type index) {
        pasta::BitVector& top_level = *block_tree_types_[0];
        auto& top_level_rs = *block_tree_types_rs_[0];
        auto& top_level_ptr = *block_tree_pointers_[0];
        auto& top_level_off = *block_tree_offsets_[0];
        int64_t c_index = chars_index_[c];
        int64_t block_size = block_size_lvl_[0];
        int64_t blk_pointer = index / block_size;
        int64_t off = index % block_size;
        int64_t rank = (blk_pointer == 0) ? 0 : c_ranks_[c_index][0][blk_pointer - 1];
        int64_t child = 0;
        if (top_level[blk_pointer]) {
            block_size /= tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer = top_level_rs.rank1(blk_pointer) * tau_ + child;
        } else {
            size_type blk = top_level_rs.rank0(blk_pointer);
            rank -= pointer_c_ranks_[c_index][0][blk];
            size_type to = off + top_level_off[blk];
            off = off + top_level_off[blk];
            blk_pointer = top_level_ptr[blk];
            child = blk_pointer;
            if (to >= block_size) {
                int64_t adder = (child == 0) ? c_ranks_[c_index][0][blk_pointer] : c_ranks_[c_index][0][blk_pointer] - c_ranks_[c_index][0][blk_pointer - 1];
                rank += adder;
                blk_pointer++;
                off = to - block_size;
            }
            block_size = block_size/tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer = top_level_rs.rank1(blk_pointer) * tau_ + child;
        }
        // we first calculate the
        size_type i = 1;
        while (i < block_tree_types_.size()) {
            rank += (child == 0) ? 0 : c_ranks_[c_index][i][blk_pointer - 1];
            if ((*block_tree_types_[i])[blk_pointer]) {
                size_type rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
                block_size /= tau_;
                child = off / block_size;
                off = off % block_size;
                blk_pointer = rank_blk * tau_ + child;
                i++;
            } else {
                size_type blk = block_tree_types_rs_[i]->rank0(blk_pointer);
                rank -= pointer_c_ranks_[c_index][i][blk];
                size_type ptr_off = (*block_tree_offsets_[i])[blk];
                size_type to = off + ptr_off;
                off = off + ptr_off;
                blk_pointer = (*block_tree_pointers_[i])[blk];
                child = blk_pointer % tau_;

                if (to >= block_size) {
                    auto adder = (child == 0) ? c_ranks_[c_index][i][blk_pointer] : c_ranks_[c_index][i][blk_pointer] - c_ranks_[c_index][i][blk_pointer - 1];
                    rank += adder;
                    blk_pointer++;
                    child = blk_pointer % tau_;
                    off = to - block_size;
                }
                auto remove_prefix = (child == 0) ? 0 : c_ranks_[c_index][i][blk_pointer - 1];
                rank -= remove_prefix;
            }
        }
        size_type prefix_leaves = blk_pointer - child;
        for (int j = 0; j < child * leaf_size; j++) {
            if ((leaves_)[prefix_leaves * leaf_size + j] == c) rank++;
        }
        for (int j = 0; j <= off; j++) {
            if ((leaves_)[blk_pointer * leaf_size + j] == c) rank++;
        }
        return rank;
    }
    int64_t rank(input_type c, size_type index) {
        pasta::BitVector& top_level = *block_tree_types_[0];
        auto& top_level_rs = *block_tree_types_rs_[0];
        auto& top_level_ptr = *block_tree_pointers_[0];
        auto& top_level_off = *block_tree_offsets_[0];
        int64_t c_index = chars_index_[c];
        int64_t block_size = block_size_lvl_[0];
        int64_t blk_pointer = index / block_size;
        int64_t off = index % block_size;
        int64_t rank = (blk_pointer == 0) ? 0 : c_ranks_[c_index][0][blk_pointer - 1];
        int64_t child = 0;
        if (top_level[blk_pointer]) {
            block_size /= tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer = top_level_rs.rank1(blk_pointer) * tau_ + child;
        } else {
            size_type blk = top_level_rs.rank0(blk_pointer);
            rank -= pointer_c_ranks_[c_index][0][blk];
            off = off + top_level_off[blk];
            blk_pointer = top_level_ptr[blk];
            child = blk_pointer;
            if (off >= block_size) {
                rank += (child == 0) ? c_ranks_[c_index][0][blk_pointer] : c_ranks_[c_index][0][blk_pointer] - c_ranks_[c_index][0][blk_pointer - 1];
                blk_pointer++;
                off = off - block_size;
            }
            block_size = block_size/tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer = top_level_rs.rank1(blk_pointer) * tau_ + child;
        }
        // we first calculate the 
        size_type i = 1;
        while (i < block_tree_types_.size()) {
            rank += (child == 0) ? 0 : c_ranks_[c_index][i][blk_pointer - 1];
            if ((*block_tree_types_[i])[blk_pointer]) {
                size_type rank_blk = block_tree_types_rs_[i]->rank1(blk_pointer);
                block_size /= tau_;
                child = off / block_size;
                off = off % block_size;
                blk_pointer = rank_blk * tau_ + child;
                i++;
            } else {
                size_type blk = block_tree_types_rs_[i]->rank0(blk_pointer);
                rank -= pointer_c_ranks_[c_index][i][blk];
                size_type ptr_off = (*block_tree_offsets_[i])[blk];
                off = off + ptr_off;
                blk_pointer = (*block_tree_pointers_[i])[blk];
                child = blk_pointer % tau_;
                if (off >= block_size) {
                    rank += (child == 0) ? c_ranks_[c_index][i][blk_pointer] : c_ranks_[c_index][i][blk_pointer] - c_ranks_[c_index][i][blk_pointer - 1];
                    blk_pointer++;
                    child = blk_pointer % tau_;
                    off = off - block_size;
                }
                auto remove_prefix = (child == 0) ? 0 : c_ranks_[c_index][i][blk_pointer - 1];
                rank -= remove_prefix;
            }
        }
        size_type prefix_leaves = blk_pointer - child;
        for (int j = 0; j < child * leaf_size; j++) {
            if ((leaves_)[prefix_leaves * leaf_size + j] == c) rank++;
        }
        for (int j = 0; j <= off; j++) {
            if ((leaves_)[blk_pointer * leaf_size + j] == c) rank++;
        }
        return rank;
    };
    int64_t print_space_usage() {
        int64_t space_usage = sizeof(tau_) + sizeof(max_leaf_length_) + sizeof(s_) + sizeof(leaf_size);
        for (auto bv: block_tree_types_) {
            space_usage += bv->size()/8;
        }
        for (auto rs: block_tree_types_rs_) {
            space_usage += rs->space_usage();
        }
        for (const auto iv: block_tree_pointers_) {
            space_usage += (int64_t) sdsl::size_in_bytes(*iv);
        }
        for (const auto iv: block_tree_offsets_) {
            space_usage +=  (int64_t) sdsl::size_in_bytes(*iv);
        }
        if (rank_support) {
            for (auto c: chars_) {
                int64_t sum= 0;
                for (auto lvl: pointer_c_ranks_[chars_index_[c]]) {
                    sum += sdsl::size_in_bytes(lvl);
                }
                for (auto lvl: c_ranks_[chars_index_[c]]) {
                    sum += sdsl::size_in_bytes(lvl);
                }
                space_usage += sum;
            }
        }

        for (auto v: block_size_lvl_) {
            space_usage += sizeof(v);
        }
        for (auto v: block_per_lvl_) {
            space_usage += sizeof(v);
        }
        space_usage += leaves_.size() * sizeof(input_type);

        return space_usage;
    };
    int32_t add_rank_support() {
        rank_support = true;
        c_ranks_.resize(chars_.size(), std::vector<sdsl::int_vector<0>>());
        pointer_c_ranks_.resize(chars_.size(), std::vector<sdsl::int_vector<0>>());
        for (int i = 0; i < c_ranks_.size(); i++) {
            c_ranks_[i].resize(block_tree_types_.size(), sdsl::int_vector<0>());
            for (int j = 0; j < c_ranks_[i].size(); j++) {
                c_ranks_[i][j].resize(block_tree_types_[j]->size());
            }
        }
        for (int i = 0; i < pointer_c_ranks_.size(); i++) {
            pointer_c_ranks_[i].resize(block_tree_pointers_.size(), sdsl::int_vector<0>());
            for (int j = 0; j < pointer_c_ranks_[i].size(); j++) {
                pointer_c_ranks_[i][j].resize(block_tree_pointers_[j]->size());
            }
        }
        for (auto c: chars_) {
            for (size_type i = 0; i < block_tree_types_[0]->size(); i++) {
                rank_block(c, 0, i);
            }
            size_type temp1 = 0;
            size_type temp2 = 0;
            size_type max = 0;
            for (size_type i = 1; i < block_tree_types_[0]->size(); i++) {
                c_ranks_[chars_index_[c]][0][i] += c_ranks_[chars_index_[c]][0][i - 1];
                if (c_ranks_[chars_index_[c]][0][i] > max) {
                    max = c_ranks_[chars_index_[c]][0][i];
                }
            }
            for (size_type i = 1; i < block_tree_types_.size(); i++) {
                size_type counter = tau_;
                size_type acc = 0;
                for (size_type j = 0; j < block_tree_types_[i]->size(); j++) {
                    size_type temp = c_ranks_[chars_index_[c]][i][j];
                    c_ranks_[chars_index_[c]][i][j] += acc;
                    acc += temp;
                    counter--;
                    if (counter == 0) {
                        acc = 0;
                        counter = tau_;
                    }
                }
//                size_type j = 0;
//                while (j < block_tree_types_[i]->size() - 1) {
//                    temp1 = c_ranks_[chars_index_[c]][i][j];
//                    temp2 = 0;
//                    for (size_type k = 0; k < tau_ - 1; k++) {
//                        std::cout << k << " " << block_tree_types_[i]->size() << std::endl;
//                        temp1 = c_ranks_[chars_index_[c]][i][j + k];
//                        temp2 = c_ranks_[chars_index_[c]][i][j + k + 1];
//                        if (c == 'e' && i == 8) {
//                            std::cout << "dieses " << temp1 << " " << temp2 << std::endl;
//                        }
//                        c_ranks_[chars_index_[c]][i][j + k + 1] += c_ranks_[chars_index_[c]][i][j + k];
//                        temp1 = temp2;
//                    };
//                    j += tau_;
//                }
//                std::cout << c << " " << "rank dones3" << i << std::endl;
//                temp1 = c_ranks_[chars_index_[c]][i][j];
//                temp2 = 0;
//                while (j < block_tree_types_[i]->size() - 1) {
//                    temp1 = c_ranks_[chars_index_[c]][i][j];
//                    temp2 = c_ranks_[chars_index_[c]][i][j + 1];
//                    c_ranks_[chars_index_[c]][i][j + 1] += c_ranks_[chars_index_[c]][i][j];
//                    temp1 = temp2;
//                    j++;
//                }
            }
                for (size_type i = 0; i < pointer_c_ranks_[chars_index_[c]].size(); i++) {
                    sdsl::util::bit_compress(pointer_c_ranks_[chars_index_[c]][i]);
                }
                for (size_type i = 0; i < c_ranks_[chars_index_[c]].size(); i++) {
                    sdsl::util::bit_compress(c_ranks_[chars_index_[c]][i]);
                }
        }
//        c_ranks_.resize(u_chars_, )
//        for (auto c: chars_) {
//            size_type rank_c = 0;
//            for (size_type i = 0; i < block_tree_types_[0]->size(); i++) {
//                rank_c += rank_block(c,0,i);
//            }
//            std::cout << c << ": " << rank_c << std::endl;
//            for (int i = 0; i < c_ranks_[chars_index_[c]].size(); i++) {
//                for (int j = 0; j < c_ranks_[chars_index_[c]][i].size(); j++) {
//                    std::cout << c_ranks_[chars_index_[c]][i][j] << " ";
//                }
//                std::cout << std::endl;
//            }
//            for (int i = 0; i < pointer_c_ranks_[chars_index_[c]].size(); i++) {
//                for (int j = 0; j < pointer_c_ranks_[chars_index_[c]][i].size(); j++) {
//                    std::cout << pointer_c_ranks_[chars_index_[c]][i][j] << " ";
//                }
//                std::cout << std::endl;
//            }
//        }
//        for (int64_t i = 0; i < amount_of_leaves; i++) {
//            leaf_ranks[chars_index_[c]][i]  = rank_leaf(c,i, leaf_size);
//        }
//        std::vector<std::vector<std::vector<size_type>>> ranks;
//        std::vector<std::vector<std::vector<size_type>>> pointer_ranks;
//        for (int i = block_tree_types_.size() - 1; i >= 0; i++) {
//            for (int j = 0; j < block_tree_types_[i]->size(); j++) {
//                size_type child = j % tau_;
//                std::vector<size_type> accumulator = std::vector<size_type>(u_chars_, 0);
//                for (int k = 0; k < tau_; k++) {
////                    for (auto c: chars_) {
////                        ranks[i][j][chars_index_[c]]
////                    }
//                }
//
//                // ranks as child
//                if (!(*block_tree_types_[i])[j]) {
//                    // ranks internal
//                }
//            }
//        }

        return 0;
    }
    int32_t add_encoded_pointers() {
        auto block_size = block_size_lvl_[0];
        for (size_type i = 0; i < block_tree_pointers_.size(); i++) {
            auto p = new sdsl::int_vector(block_tree_pointers_[i]->size(), 0, 64);
            for (size_type j = 0; j < block_tree_pointers_[i]->size(); j++) {
                (*p)[j] = (*block_tree_pointers_[i])[j] * block_size + (*block_tree_offsets_[i])[j];
            }
            block_size /= tau_;
            sdsl::util::bit_compress(*p);
            this->block_tree_encoded_.push_back(p);
        }
        return 0;
    }
    int32_t add_rank_support_omp(int32_t threads) {
        rank_support = true;
        c_ranks_.resize(chars_.size(), std::vector<sdsl::int_vector<0>>());
        pointer_c_ranks_.resize(chars_.size(), std::vector<sdsl::int_vector<0>>());
        for (int i = 0; i < c_ranks_.size(); i++) {
            c_ranks_[i].resize(block_tree_types_.size(), sdsl::int_vector<0>());
            for (int j = 0; j < c_ranks_[i].size(); j++) {
                c_ranks_[i][j].resize(block_tree_types_[j]->size());
            }
        }
        for (int i = 0; i < pointer_c_ranks_.size(); i++) {
            pointer_c_ranks_[i].resize(block_tree_pointers_.size(), sdsl::int_vector<0>());
            for (int j = 0; j < pointer_c_ranks_[i].size(); j++) {
                pointer_c_ranks_[i][j].resize(block_tree_pointers_[j]->size());
            }
        }
        omp_set_num_threads(threads);

#pragma omp parallel for default(none)
            for (auto c: chars_) {
                for (size_type i = 0; i < block_tree_types_[0]->size(); i++) {
                    rank_block(c, 0, i);
                }
                size_type temp1 = 0;
                size_type temp2 = 0;
                size_type max = 0;
                for (size_type i = 1; i < block_tree_types_[0]->size(); i++) {
                    c_ranks_[chars_index_[c]][0][i] += c_ranks_[chars_index_[c]][0][i - 1];
                    if (c_ranks_[chars_index_[c]][0][i] > max) {
                        max = c_ranks_[chars_index_[c]][0][i];
                    }
                }
                for (size_type i = 1; i < block_tree_types_.size(); i++) {
                    size_type counter = tau_;
                    size_type acc = 0;
                    for (size_type j = 0; j < block_tree_types_[i]->size(); j++) {
                        size_type temp = c_ranks_[chars_index_[c]][i][j];
                        c_ranks_[chars_index_[c]][i][j] += acc;
                        acc += temp;
                        counter--;
                        if (counter == 0) {
                            acc = 0;
                            counter = tau_;
                        }
                    }
//                size_type j = 0;
//                while (j < block_tree_types_[i]->size() - 1) {
//                    temp1 = c_ranks_[chars_index_[c]][i][j];
//                    temp2 = 0;
//                    for (size_type k = 0; k < tau_ - 1; k++) {
//                        std::cout << k << " " << block_tree_types_[i]->size() << std::endl;
//                        temp1 = c_ranks_[chars_index_[c]][i][j + k];
//                        temp2 = c_ranks_[chars_index_[c]][i][j + k + 1];
//                        if (c == 'e' && i == 8) {
//                            std::cout << "dieses " << temp1 << " " << temp2 << std::endl;
//                        }
//                        c_ranks_[chars_index_[c]][i][j + k + 1] += c_ranks_[chars_index_[c]][i][j + k];
//                        temp1 = temp2;
//                    };
//                    j += tau_;
//                }
//                std::cout << c << " " << "rank dones3" << i << std::endl;
//                temp1 = c_ranks_[chars_index_[c]][i][j];
//                temp2 = 0;
//                while (j < block_tree_types_[i]->size() - 1) {
//                    temp1 = c_ranks_[chars_index_[c]][i][j];
//                    temp2 = c_ranks_[chars_index_[c]][i][j + 1];
//                    c_ranks_[chars_index_[c]][i][j + 1] += c_ranks_[chars_index_[c]][i][j];
//                    temp1 = temp2;
//                    j++;
//                }
                }
                for (size_type i = 0; i < pointer_c_ranks_[chars_index_[c]].size(); i++) {
                    sdsl::util::bit_compress(pointer_c_ranks_[chars_index_[c]][i]);
                }
                for (size_type i = 0; i < c_ranks_[chars_index_[c]].size(); i++) {
                    sdsl::util::bit_compress(c_ranks_[chars_index_[c]][i]);
                }
            }


//        c_ranks_.resize(u_chars_, )
//        for (auto c: chars_) {
//            size_type rank_c = 0;
//            for (size_type i = 0; i < block_tree_types_[0]->size(); i++) {
//                rank_c += rank_block(c,0,i);
//            }
//            std::cout << c << ": " << rank_c << std::endl;
//            for (int i = 0; i < c_ranks_[chars_index_[c]].size(); i++) {
//                for (int j = 0; j < c_ranks_[chars_index_[c]][i].size(); j++) {
//                    std::cout << c_ranks_[chars_index_[c]][i][j] << " ";
//                }
//                std::cout << std::endl;
//            }
//            for (int i = 0; i < pointer_c_ranks_[chars_index_[c]].size(); i++) {
//                for (int j = 0; j < pointer_c_ranks_[chars_index_[c]][i].size(); j++) {
//                    std::cout << pointer_c_ranks_[chars_index_[c]][i][j] << " ";
//                }
//                std::cout << std::endl;
//            }
//        }
//        for (int64_t i = 0; i < amount_of_leaves; i++) {
//            leaf_ranks[chars_index_[c]][i]  = rank_leaf(c,i, leaf_size);
//        }
//        std::vector<std::vector<std::vector<size_type>>> ranks;
//        std::vector<std::vector<std::vector<size_type>>> pointer_ranks;
//        for (int i = block_tree_types_.size() - 1; i >= 0; i++) {
//            for (int j = 0; j < block_tree_types_[i]->size(); j++) {
//                size_type child = j % tau_;
//                std::vector<size_type> accumulator = std::vector<size_type>(u_chars_, 0);
//                for (int k = 0; k < tau_; k++) {
////                    for (auto c: chars_) {
////                        ranks[i][j][chars_index_[c]]
////                    }
//                }
//
//                // ranks as child
//                if (!(*block_tree_types_[i])[j]) {
//                    // ranks internal
//                }
//            }
//        }

        return 0;
    }
    inline size_type leading_zeros(int32_t val) {
        return __builtin_clz(static_cast<unsigned int>(val) | 1);
    }
    inline size_type leading_zeros(int64_t val) {
        return __builtin_clzll(static_cast<unsigned long long>(val) | 1);
    }
    void calculate_padding(int64_t &padding, int64_t text_length, int64_t &height, int64_t &blk_size) {
        int64_t tmp_padding = this->s_ * this->tau_;
        int64_t h = 1;
        blk_size = tau_;
        while (tmp_padding < text_length) {
            tmp_padding *= this->tau_;
            blk_size *= this->tau_;
            h++;
        }
        height = h;
        padding = tmp_padding - text_length;
//        std::cout << "Padding: " << padding << " h: " << h << " SIZE: " << tmp_padding << " BLK_SIZE: " <<  blk_size <<   std::endl;
    }
    size_type rank_block(input_type c, size_type i, size_type j) {
        if (j >= block_tree_types_[i]->size()) {
            return 0;
        }
        size_type rank_c = 0;
        if ((*block_tree_types_[i])[j] == 1) {
            if (i != block_tree_types_.size() - 1) {
                size_type rank_blk = block_tree_types_rs_[i]->rank1(j);
                for (size_type k = 0; k < tau_; k++) {
                    rank_c += rank_block(c, i + 1, rank_blk * tau_ + k);
                }
            } else {
                size_type rank_blk = block_tree_types_rs_[i]->rank1(j);
                for (size_type k = 0; k < tau_; k++) {
                    rank_c += rank_leaf(c, rank_blk * tau_ + k, leaf_size);
                }
            }
        } else {
            size_type rank_0 = block_tree_types_rs_[i]->rank0(j);
            size_type ptr = (*block_tree_pointers_[i])[rank_0];
            size_type off = (*block_tree_offsets_[i])[rank_0];
            size_type rank_g = 0;
            rank_c += c_ranks_[chars_index_[c]][i][ptr];
            if (off != 0) {
                rank_g = part_rank_block(c,i, ptr, off);
                size_type rank_2nd = part_rank_block(c,i,ptr + 1, off);
                rank_c -= rank_g;
                rank_c += rank_2nd;
            }
            pointer_c_ranks_[chars_index_[c]][i][rank_0] = rank_g;
        }
        c_ranks_[chars_index_[c]][i][j] = rank_c;
        return rank_c;
    }
    size_type part_rank_block(input_type c, size_type i, size_type j, size_type g) {
        if (j >= block_tree_types_[i]->size()) {
            return 0;
        }
        size_type rank_c = 0;
        size_type blk_size = block_size_lvl_[i];
        if ((*block_tree_types_[i])[j] == 1) {
            if (i != block_tree_types_.size() - 1) {
                size_type rank_blk = block_tree_types_rs_[i]->rank1(j);
                size_type k = 0;
                size_type k_sum = 0;
                for (k = 0; k < tau_ && k_sum + block_size_lvl_[i + 1] <= g; k++) {
                    rank_c += c_ranks_[chars_index_[c]][i + 1][rank_blk * tau_ + k];
                    k_sum+=block_size_lvl_[i + 1];
                }

                if (k_sum != g) {
                    rank_c += part_rank_block(c, i + 1, rank_blk * tau_ + k, g - k_sum);
                }
            } else {
                size_type rank_blk = block_tree_types_rs_[i]->rank1(j);
                size_type k = 0;
                size_type k_sum = 0;
                for (k = 0; k < tau_ && k_sum + leaf_size <= g; k++) {
                    rank_c += rank_leaf(c, rank_blk * tau_ + k, leaf_size);
                    k_sum += leaf_size;
                }

                if (k_sum != g) {
                    rank_c += rank_leaf(c, rank_blk * tau_ + k, g % leaf_size);
                }
            }
        } else {
            size_type rank_0 = block_tree_types_rs_[i]->rank0(j);
            size_type ptr = (*block_tree_pointers_[i])[rank_0];
            size_type off = (*block_tree_offsets_[i])[rank_0];
            if (g + off >= block_size_lvl_[i]) {
                rank_c += c_ranks_[chars_index_[c]][i][ptr] - pointer_c_ranks_[chars_index_[c]][i][rank_0] + part_rank_block(c, i, ptr + 1, g + off - block_size_lvl_[i]);
            } else {
                rank_c += part_rank_block(c, i, ptr, g + off) - pointer_c_ranks_[chars_index_[c]][i][rank_0];
            }
        }
        return rank_c;
    }
    size_type rank_leaf(input_type c, size_type leaf_index,size_type i) {

        if (leaf_index * leaf_size >= leaves_.size()) {
            return 0;
        }
//        size_type x = leaves_.size() - leaf_index * this->tau_;
//        i = std::min(i, x);
        size_type result = 0;
        for (size_type ind = 0; ind < i; ind++) {
            if (leaves_[leaf_index * leaf_size + ind] == c) {
                result++;
            }
        }
        return result;
    }

    size_type map_unique_chars(std::vector<input_type> &text) {
        this->u_chars_ = 0;
        input_type i = 0;
        for (auto a: text) {
            if (chars_index_.find(a) == chars_index_.end()) {
                chars_index_[a] = i;
                i++;
                chars_.push_back(a);
            }
        }
        this->u_chars_ = i;
        return 0;
    };
    size_type find_next_smallest_index_binary_search(size_type i, std::vector<int64_t>& pVector) {
        int64_t l = 0;
        int64_t r = pVector.size();
        while (l < r) {
            int64_t m = std::floor((l+r)/2);
            if (i < pVector[m]) {
                r = m;
            } else {
                l = m + 1;
            }
        }
        return r -1 ;
    };
    int64_t find_next_smallest_index_linear_scan(size_type i, std::vector<size_type>& pVector) {
        int64_t b = 0;
        while (b < pVector.size() && i >= pVector[b]) {
            b++;
        }
        return b - 1;
    };
    size_type find_next_smallest_index_block_tree(size_type index) {
        size_type block_size = this->block_size_lvl_[0];
        size_type blk_pointer = index / block_size;
        size_type off = index % block_size;
        size_type child = 0;
        for (size_type i = 0; i < this->block_tree_types_.size(); i++) {
            if ((*this->block_tree_types_[i])[blk_pointer] == 0) {
                return -1;
            }
            if (off > 0 && (*this->block_tree_types_[i])[blk_pointer + 1] == 0) {
                return -1;
            }
            size_type rank_blk = this->block_tree_types_rs_[i]->rank1(blk_pointer);
            blk_pointer = rank_blk * this->tau_;
            block_size /= this->tau_;
            child = off / block_size;
            off = off % block_size;
            blk_pointer += child;
        }
        return blk_pointer;
    };
};
#endif //BLOCK_TREE_BV_BLOCKTREE_HPP
