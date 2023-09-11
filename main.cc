#include <fstream>

#include "algo/DNA_algo.h"
#include "menu.h"

void MenuItemDnaSearch() {
  std::string haystack_path =
      Menu::GetUserInput("Enter path to text file(Haystack(for example <tests/part-01-haystack.txt>)):\n");
  std::string needle_path =
      Menu::GetUserInput("Enter path to text file(Needle(for example <tests/part-01-needle.txt>)):\n");

  std::ifstream haystack_input(haystack_path);
  std::ifstream needle_input(needle_path);

  std::string needle, haystack;
  haystack_input >> haystack;
  needle_input >> needle;

  if (!haystack_input.is_open() || haystack_input.fail())
    throw std::invalid_argument("Can't read haystack");

  if (!needle_input.is_open() || needle_input.fail())
    throw std::invalid_argument("Can't read needle");

  if (!DnaAlgorithm::correctSequence(needle))
    throw std::invalid_argument(
        "Needle has wrong dna format. Only use alphabet {A, C, G, T } "
        "uppercase");

  if (!DnaAlgorithm::correctSequence(haystack))
    throw std::invalid_argument(
        "Haystack has wrong dna format. Only use alphabet {A, C, G, T } "
        "uppercase");

  int q = 13;
  std::vector<size_t> res = DnaAlgorithm::rabinKarpSearch(needle, haystack, q);

  std::cout << "Occurrences in a string:" << std::endl;
  for (auto r : res) std::cout << r << " ";
}

void MenuItemNwSequenceAlignment() {
  std::string input_file_path =
      Menu::GetUserInput("Enter path to text input file:\n");

  std::ifstream input(input_file_path);

  int match_s, mismatch_s, gap_s;
  std::string str1, str2;

  input >> match_s >> mismatch_s >> gap_s >> str1 >> str2;

  if (!input.is_open() || input.fail())
    throw std::invalid_argument("Can't read input file");

  if (!DnaAlgorithm::correctSequence(str1) ||
      !DnaAlgorithm::correctSequence(str2))
    throw std::invalid_argument(
        "One of strings has wrong dna format. Only use alphabet {A, C, G, T } "
        "uppercase");

  DnaAlgorithm::newSequence res =
      DnaAlgorithm::newSequenceAlignment(match_s, mismatch_s, gap_s, str1, str2);

  std::string lines;
  for (size_t i = 0; i < res.string1.size(); i++)
    lines.push_back(res.string1.at(i) == res.string2.at(i) ? '|' : ' ');

  std::cout << "Score and aligned strings: " << std::endl
            << res.score << std::endl
            << res.string1 << std::endl
            << lines << std::endl
            << res.string2;
}

void MenuItemRegExpr() {
  std::string input_file_path =
      Menu::GetUserInput("Enter path to text input file:\n");

  std::string regexp, text;

  std::ifstream input(input_file_path);
  input >> text >> regexp;

  if (!input.is_open() || input.fail())
    throw std::invalid_argument("Can't read input file");

  if (!DnaAlgorithm::correctSequence(text))
    throw std::invalid_argument(
        "Text has wrong dna format. Only use alphabet {A, C, G, T } uppercase");

  bool match = DnaAlgorithm::regExprChecker(regexp, text);

  std::cout << "Match" << std::endl << (match ? "True" : "False");
}

void MenuItemKSimilar() {
  std::string input_file_path =
      Menu::GetUserInput("Enter path to text input file:\n");

  std::string str1, str2;

  std::ifstream input(input_file_path);
  input >> str1 >> str2;

  if (!input.is_open() || input.fail())
    throw std::invalid_argument("Can't read input file");

  if (!DnaAlgorithm::correctSequence(str1) ||
      !DnaAlgorithm::correctSequence(str2))
    throw std::invalid_argument(
        "One of strings has wrong dna format. Only use alphabet {A, "
        "C, G, T } "
        "uppercase");

  std::cout << "K = " << std::endl << DnaAlgorithm::kSimilar(str1, str2);
}

void MenuItemMinimumWindow() {
  std::string input_file_path =
      Menu::GetUserInput("Enter path to text input file:\n");

  std::string str1, str2;

  std::ifstream input(input_file_path);
  input >> str1 >> str2;

  if (!input.is_open() || input.fail())
    throw std::invalid_argument("Can't read input file");

  if (!DnaAlgorithm::correctSequence(str1) ||
      !DnaAlgorithm::correctSequence(str2))
    throw std::invalid_argument(
        "One of strings has wrong dna format. Only use alphabet {A, "
        "C, G, T } "
        "uppercase");

  std::cout << "Minimum window substr" << std::endl
            << DnaAlgorithm::minWindowString(str1, str2);
}

int main() {
  Menu menu;
  menu.AddMenuItem({"Exact DNA search", MenuItemDnaSearch});
  menu.AddMenuItem(
      {"NW sequence alignment project", MenuItemNwSequenceAlignment});
  menu.AddMenuItem({"Matching regular expressions", MenuItemRegExpr});
  menu.AddMenuItem({"K-similar strings", MenuItemKSimilar});
  menu.AddMenuItem({"Minimum window substring", MenuItemMinimumWindow});

  menu.Start();

  return 0;
}
