//
// Created by shoko on 09.09.2023.
//

#include "DNA_algo.h"

namespace DnaAlgorithm {

    bool matchHere(const std::string &regexp, size_t r_ind, const std::string &text,
                   size_t t_ind);

    bool correctSequence(const std::string &string) {
        for (size_t i = 0; i < string.size(); i++)
            if (!alphabet.count(string[i]))
                return false;
        return true;
    }

    std::vector<size_t> rabinKarpSearch(const std::string &needle, const std::string &haystack, int q) {
        std::vector<size_t> result;
        if (needle.size() > haystack.size() || needle.size() == 0 || haystack.size() == 0)
            return result;
        int hay_hash = 0;
        int needle_hash = 0;
        int h = 1;
        int alphabet_size = (int)alphabet.size();

        for (size_t i = 0; i < needle.size() - 1; i++)
            h = ((h * alphabet_size) % q);

        for (size_t i = 0; i < needle.size(); i++) {
            hay_hash = ((alphabet_size * hay_hash + alphabet[haystack[i]]) % q);
            needle_hash = ((alphabet_size * needle_hash + alphabet[needle[i]]) % q);
        }

        for (size_t i = 0; i <= haystack.size() - needle.size(); i++) {
            if (needle_hash == hay_hash && haystack.compare(i , needle.size(), needle) == 0) {
                result.push_back(i);
            }

            hay_hash = (alphabet_size * (hay_hash - alphabet[haystack[i]] * h) + alphabet[haystack[i + needle.size()]]) % q;

            if (hay_hash < 0)
                hay_hash = (hay_hash + q);
        }

        return result;
    }

    newSequence newSequenceAlignment(const int match_s, const int mismatch_s, int gap_s,
                                     const std::string &str1, const std::string &str2) {

        std::vector<std::vector<char>> backtrack_matrix(
                str1.size() + 1, std::vector<char>(str2.size() + 1, 'n'));

        std::vector<std::vector<int>> matrix(str1.size() + 1,
                                             std::vector<int>(str2.size() + 1));
        matrix[0][0] = 0;
        for (int i = 0; i < (int)str1.size() + 1; i++) {
            matrix[i][0] = i * gap_s;
            backtrack_matrix[i][0] = 'u';
        }
        for (int i = 0; i < (int)str2.size() + 1; i++) {
            matrix[0][i] = i * gap_s;
            backtrack_matrix[0][i] = 'l';
        }

        for (int i = 1; i < (int)str1.size() + 1; i++)
            for (int j = 1; j < (int)str2.size() + 1; j++) {
                int diagonal = str1[i - 1] == str2[j - 1]
                               ? matrix[i - 1][j - 1] + match_s
                               : matrix[i - 1][j - 1] + mismatch_s;
                int up = matrix[i - 1][j] + gap_s;
                int left = matrix[i][j - 1] + gap_s;

                if (up >= left && up >= diagonal) {
                    matrix[i][j] = up;
                    backtrack_matrix[i][j] = 'u';
                } else if (diagonal >= up && diagonal >= left) {
                    matrix[i][j] = diagonal;
                    backtrack_matrix[i][j] = 'd';
                } else if (left >= up && left >= diagonal) {
                    matrix[i][j] = left;
                    backtrack_matrix[i][j] = 'l';
                }
            }

        newSequence nw_res;
        int res_len = std::max((int)str1.size(), (int)str2.size());
        nw_res.string1.reserve(res_len);
        nw_res.string2.reserve(res_len);
        nw_res.score = matrix.back().back();

        int i = (int)str1.size();
        int j = (int)str2.size();
        while (i != 0 || j != 0) {
            char move = backtrack_matrix[i][j];
            switch(move) {
                case 'd':
                    nw_res.string1.insert(0, 1, str1[i - 1]);
                    nw_res.string2.insert(0, 1, str2[j - 1]);
                    i -= 1;
                    j -= 1;
                    break;
                case 'u':
                    nw_res.string1.insert(0, 1, str1[i - 1]);
                    nw_res.string2.insert(0, 1, '-');
                    i -= 1;
                    break;
                case 'l':
                    nw_res.string1.insert(0, 1, '-');
                    nw_res.string2.insert(0, 1, str2[j - 1]);
                    j -= 1;
            }
        }

        return nw_res;
    }

    int matchStar(const std::string &regexp, size_t r_ind, const std::string &text,
                  size_t t_ind) {
        do {
            if (matchHere(regexp, r_ind, text, t_ind)) return 1;
        } while (t_ind++ < text.size());
        return 0;
    }

    int matchPlus(char ch, const std::string &regexp, size_t r_ind,
                  const std::string &text, size_t t_ind) {
        do {
            if (matchHere(regexp, r_ind, text, t_ind))
                return 1;
        } while (t_ind < text.size() && (text[t_ind++] == ch || ch == '.'));
        return 0;
    }

    bool matchHere(const std::string &regexp, size_t r_ind, const std::string &text,
                  size_t t_ind) {
        if (regexp.size() == r_ind)
            return text.size() == t_ind;

        if (r_ind + 1 < regexp.size() && regexp[r_ind + 1] == '+')
            return matchPlus(regexp[r_ind], regexp, r_ind + 2, text, t_ind);

        if (regexp[r_ind] == '?')
            return (matchHere(regexp, r_ind + 1, text, t_ind + 1) ||
                    matchHere(regexp, r_ind + 1, text, t_ind));

        if (regexp[r_ind] == '*')
            return matchStar(regexp, r_ind + 1, text, t_ind);

        if (regexp[r_ind] == '.' || regexp[r_ind] == text[t_ind])
            return matchHere(regexp, r_ind + 1, text, t_ind + 1);

        return 0;
    }

    bool regExprChecker(const std::string &string1, const std::string &string2) {
        return matchHere(string1, 0, string2, 0);
    }

    int kSimilar(const std::string &string1, const std::string &string2) {
        std::queue<std::string> q;
        std::set<std::string> visited;
        int count = 0;

        q.push(string1);
        visited.insert(string1);

        while (!q.empty()) {
            size_t q_size = q.size();

            for (size_t i = 0; i < q_size; ++i) {
                std::string s = q.front();
                q.pop();

                size_t p = 0;
                while (p < s.length() && s[p] == string2[p])
                    ++p;

                if (p == s.length())
                    return count;

                for (size_t j = p + 1; j < s.length(); ++j) {
                    if (s[j] != string2[j] && s[j] == string2[p]) {
                        std::string tmp = s;
                        std::swap(tmp[j], tmp[p]);
                        if (!visited.count(tmp)) {
                            q.push(tmp);
                            visited.insert(tmp);
                        }
                    }
                }
            }

            ++count;
        }

        return count;
    }

    bool contains(const std::unordered_map<char, int> &counts, char ch) {
        return counts.find(ch) != counts.end();
    }

    std::string minWindowString(std::string &s, std::string &t) {
        const int m = (int)s.length();
        const int n = (int)t.length();
        if (m < n || m == 0 || n == 0)
            return "";

        std::unordered_map<char, int> counts;
        int balance = 0;
        for (const auto &ch : t) {
            ++counts[ch];
            ++balance;
        }

        int left = 0, right = -1, start = -1, minLen = m;
        while (right < m) {
            if (balance == 0) {
                int len = right - left + 1;
                if (len < minLen) {
                    minLen = len;
                    start = left;
                }
                if (contains(counts, s[left]) && (++counts[s[left]] >= 1))
                    ++balance;

                ++left;
            } else {
                if (++right >= m)
                    break;

                if (contains(counts, s[right]) && (--counts[s[right]] >= 0))
                    --balance;

                if (start > -1) {
                    if (contains(counts, s[left]) && (++counts[s[left]] >= 1))
                        ++balance;

                    ++left;
                }
            }
        }

        return (start > -1) ? s.substr(start, minLen) : "";
    }

} // DnaAlgorithm