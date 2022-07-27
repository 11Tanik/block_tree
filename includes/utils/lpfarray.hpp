
#include "libsais.h"
#include "libsais64.h"
#include <iostream>
#include <vector>
#ifndef BLOCK_TREE_LPFARRAY_H
#define BLOCK_TREE_LPFARRAY_H
int32_t lpf_array64(std::vector<uint8_t> &text, std::vector<int64_t> &lpf, std::vector<int64_t> &lpf_ptr) {
    std::vector<int64_t> sa(text.size());
    std::vector<int64_t> plcp(text.size());
    std::vector<int64_t> lcp(text.size());
    libsais64(text.data(), sa.data() , text.size(), 0, NULL);
    libsais64_plcp(text.data(), sa.data(), plcp.data(), text.size());
    libsais64_lcp(plcp.data(), sa.data(), lcp.data(), text.size());
    std::vector<int64_t> isu(text.size());
    std::vector<int64_t> prev(text.size());
    std::vector<int64_t> next(text.size());
    for(int64_t i = 0; i < text.size(); i++) {
        isu[sa[i]] = i;
    }
    for (int64_t r = 0; r < text.size() - 1; r++) {
        prev[r] = r-1;
        next[r] = r+1;
    }
    for (int64_t i = text.size() - 1; i >= 0; i--) {
        int64_t r = isu[i];
        lpf[i] = std::max(lcp[r], lcp[next[r]]);
        if (lcp[r] <= lcp[next[r]]) {
            lpf_ptr[i] = sa[next[r]];
            lcp[next[r]] = lcp[r];
        } else {
            lpf_ptr[i] = sa[prev[r]];
        }
        if (lpf[i] == 0) {
            lpf_ptr[i] = -1;
        }
        if (prev[r] >= 0) {
            next[prev[r]] = next[r];
        }
        if (next[r] < text.size()) {
            prev[next[r]] = prev[r];
        }
    }
    return 0;
}

int32_t lpf_array(std::vector<uint8_t> &text, std::vector<int32_t> &lpf) {
    std::vector<int32_t> sa(text.size());
    std::vector<int32_t> plcp(text.size());
    std::vector<int32_t> lcp(text.size());
    libsais(text.data(), sa.data() , text.size(), 0, NULL);
    libsais_plcp(text.data(), sa.data(), plcp.data(), text.size());
    libsais_lcp(plcp.data(), sa.data(), lcp.data(), text.size());
    std::vector<int32_t> isu(text.size());
    std::vector<int32_t> prev(text.size());
    std::vector<int32_t> next(text.size());
    for(int i = 0; i < text.size(); i++) {
        isu[sa[i]] = i;
    }
    for (int r = 0; r < text.size(); r++) {
        prev[r] = r-1;
        next[r] = r+1;
    }
    for (int i = text.size() - 1; i >= 0; i--) {
        int r = isu[i];
        lpf[i] = std::max(lcp[r], lcp[next[r]]);
        lcp[next[r]] = std::min(lcp[r], lcp[next[r]]);
        if (prev[r] >= 0) {
            next[prev[r]] = next[r];
        }
        if (next[r] < text.size()) {
            prev[next[r]] = prev[r];
        }
    }
    return 0;
}
int32_t lpf_array(std::string &text, std::vector<int32_t> &lpf, std::vector<int32_t> &lpf_ptr) {
    std::vector<int32_t> sa(text.size());
    std::vector<int32_t> plcp(text.size());
    std::vector<int32_t> lcp(text.size());
    libsais(reinterpret_cast<const uint8_t *>(text.c_str()), sa.data() , text.size(), 0, NULL);
    libsais_plcp(reinterpret_cast<const uint8_t *>(text.c_str()), sa.data(), plcp.data(), text.size());
    libsais_lcp(plcp.data(), sa.data(), lcp.data(), text.size());
    std::vector<int32_t> isu(text.size());
    std::vector<int32_t> prev(text.size());
    std::vector<int32_t> next(text.size());
    for(int i = 0; i < text.size(); i++) {
        isu[sa[i]] = i;
    }
    for (int r = 0; r < text.size() - 1; r++) {
        prev[r] = r-1;
        next[r] = r+1;
    }
    for (int i = text.size() - 1; i >= 0; i--) {
        int r = isu[i];
        lpf[i] = std::max(lcp[r], lcp[next[r]]);
        lcp[next[r]] = std::min(lcp[r], lcp[next[r]]);
        if (prev[r] >= 0) {
            next[prev[r]] = next[r];
        }
        if (next[r] < text.size()) {
            prev[next[r]] = prev[r];
        }
    }
    return 0;
}
#endif //BLOCK_TREE_LPFARRAY_H
