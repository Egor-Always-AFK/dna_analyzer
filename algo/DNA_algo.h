//
// Created by shoko on 09.09.2023.
//

#ifndef DNA_ANALYZER_DNA_ALGO_H
#define DNA_ANALYZER_DNA_ALGO_H

#include <algorithm>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <vector>


namespace DnaAlgorithm {
    static std::map<char, int> alphabet = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};

    struct newSequence {
        int score = 0;
        std::string string1;
        std::string string2;
    };

    bool correctSequence(const std::string &string);

    std::vector<size_t> rabinKarpSearch(const std::string &needle, const std::string &haystack, int q);

    newSequence newSequenceAlignment(const int match_s, const int mismatch_s, int gap_s,
                                     const std::string &str1, const std::string &str2);

    bool regExprChecker(const std::string &string1, const std::string &string2);

    int kSimilar(const std::string &string1, const std::string &string2);

    std::string minWindowString(std::string &s, std::string &t);
} // namespace DnaAlgorithm

#endif //DNA_ANALYZER_DNA_ALGO_H
