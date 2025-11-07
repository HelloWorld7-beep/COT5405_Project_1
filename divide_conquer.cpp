#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <chrono>
#include <random>

using namespace std;

// function for testing to generate random DNA sequeneces
string generateRandomSequence(size_t length) {
    static const char nucleotides[] = {'A','T','G','C'};
    static random_device rd;
    static mt19937 gen(rd());
    static uniform_int_distribution<> dis(0, 3);

    string seq;
    for (size_t i = 0; i < length; ++i) {
        seq += nucleotides[dis(gen)];
    }
    return seq;
}

// structure to hold sequence information
struct SequenceInfo {
    string id;
    string sequence;
    string description;

    SequenceInfo(const string& _id, const string& _seq, const string& _desc = "")
        : id(_id), sequence(_seq), description(_desc) {}
};

// alignment using Needleman-Wunsch algorithm
class PairwiseAlignment {
public:
    int match_score;
    int mismatch_score;
    int gap_penalty;

    PairwiseAlignment(int match = 2, int mismatch = -1, int gap = -2)
        : match_score(match), mismatch_score(mismatch), gap_penalty(gap) {}

    tuple<string, string, int> needlemanWunsch(const string& seq1, const string& seq2) {
        size_t m = seq1.size();
        size_t n = seq2.size();
        vector<vector<int>> score_matrix(m + 1, vector<int>(n + 1, 0));

        for (size_t i = 0; i <= m; ++i) {
            score_matrix[i][0] = i * gap_penalty;
        }
        for (size_t j = 0; j <= n; ++j) {
            score_matrix[0][j] = j * gap_penalty;
        }

        for (size_t i = 1; i <= m; ++i) {
            for (size_t j = 1; j <= n; ++j) {
                int matchMismatch = score_matrix[i - 1][j - 1] + getMatchScore(seq1[i - 1], seq2[j - 1]);
                int del = score_matrix[i - 1][j] + gap_penalty;
                int ins = score_matrix[i][j - 1] + gap_penalty;
                score_matrix[i][j] = std::max(matchMismatch, std::max(del, ins));
            }
        }

        string aligned1, aligned2;
        traceback(score_matrix, seq1, seq2, aligned1, aligned2);

        return make_tuple(aligned1, aligned2, score_matrix[m][n]);
    }

private:
    int getMatchScore(char a, char b) {
        return (toupper(a) == toupper(b)) ? match_score : mismatch_score;
    }

    // perform traceback to get aligned sequences
    void traceback(const vector<vector<int>>& score_matrix, const string& seq1,
                   const string& seq2, string& aligned1, string& aligned2) {
        size_t i = seq1.size();
        size_t j = seq2.size();

        while (i > 0 || j > 0) {
            int score = score_matrix[i][j];

            if (i > 0 && j > 0 && score == score_matrix[i - 1][j - 1] + getMatchScore(seq1[i - 1], seq2[j - 1])) {
                aligned1 = seq1[i - 1] + aligned1;
                aligned2 = seq2[j - 1] + aligned2;
                --i; --j;
            } else if (i > 0 && score == score_matrix[i - 1][j] + gap_penalty) {
                aligned1 = seq1[i - 1] + aligned1;
                aligned2 = '-' + aligned2;
                --i;
            } else if (j > 0 && score == score_matrix[i][j - 1] + gap_penalty) {
                aligned1 = '-' + aligned1;
                aligned2 = seq2[j - 1] + aligned2;
                --j;
            } else {
                if (i > 0) { aligned1 = seq1[i - 1] + aligned1; aligned2 = '-' + aligned2; --i; }
                else if (j > 0) { aligned1 = '-' + aligned1; aligned2 = seq2[j - 1] + aligned2; --j; }
            }
        }
    }
};

// Multiple Sequence Alignment using Divide and Conquer
class MultipleSequenceAlignment {
public:
    int match_score;
    int mismatch_score;
    int gap_penalty;

    MultipleSequenceAlignment(int match = 2, int mismatch = -1, int gap = -2)
        : match_score(match), mismatch_score(mismatch), gap_penalty(gap) {}

    // Divide and conquer MSA
    vector<string> divideAndConquerAlign(const vector<SequenceInfo>& sequences, size_t minGroupSize = 2) {
        if (sequences.size() <= minGroupSize) return alignSmallGroup(sequences);

        size_t mid = sequences.size() / 2;
        vector<SequenceInfo> left(sequences.begin(), sequences.begin() + mid);
        vector<SequenceInfo> right(sequences.begin() + mid, sequences.end());

        vector<string> leftAligned = divideAndConquerAlign(left, minGroupSize);
        vector<string> rightAligned = divideAndConquerAlign(right, minGroupSize);

        return mergeAlignments(leftAligned, rightAligned);
    }

    // Generate consensus sequence from alignment
    string getConsensus(const vector<string>& alignment) {
        if (alignment.empty()) return "";
        string consensus;
        size_t len = alignment[0].size();

        for (size_t pos = 0; pos < len; ++pos) {
            map<char, int> counts{{'A',0},{'T',0},{'G',0},{'C',0},{'U',0},{'-',0}};
            for (const auto& seq : alignment) {
                if (pos < seq.size()) {
                    char c = toupper(seq[pos]);
                    if (counts.find(c) != counts.end()) counts[c]++;
                }
            }

            int maxCount = -1;
            char maxChar = '-';
            for (auto it = counts.begin(); it != counts.end(); ++it) {
                char c = it->first;
                int count = it->second;
                if (count > maxCount && c != '-') { maxCount = count; maxChar = c; }
            }
            consensus += maxChar;
        }
        return consensus;
    }

private:
    // Align small group of sequences directly
    vector<string> alignSmallGroup(const vector<SequenceInfo>& sequences) {
        if (sequences.size() == 1) {
            return { sequences[0].sequence };
        }
        if (sequences.size() == 2) {
            PairwiseAlignment pa(match_score, mismatch_score, gap_penalty);
            string a1, a2;
            int score;
            tie(a1, a2, score) = pa.needlemanWunsch(sequences[0].sequence, sequences[1].sequence);
            return { a1, a2 };
        }

        vector<string> result = alignSmallGroup({ sequences[0], sequences[1] });
        for (size_t i = 2; i < sequences.size(); ++i) {
            string consensus = getConsensus(result);
            PairwiseAlignment pa(match_score, mismatch_score, gap_penalty);
            string alignedConsensus, alignedNext;
            int score;
            tie(alignedConsensus, alignedNext, score) = pa.needlemanWunsch(consensus, sequences[i].sequence);

            // Simplified update
            result.push_back(alignedNext);
        }
        return result;
    }

    // Merge two alignments
    vector<string> mergeAlignments(const vector<string>& left, const vector<string>& right) {
        string leftCons = getConsensus(left);
        string rightCons = getConsensus(right);
        PairwiseAlignment pa(match_score, mismatch_score, gap_penalty);
        string alignedLeft, alignedRight;
        int score;
        tie(alignedLeft, alignedRight, score) = pa.needlemanWunsch(leftCons, rightCons);

        // Simplified merge: just concatenate
        vector<string> merged = left;
        merged.insert(merged.end(), right.begin(), right.end());
        return merged;
    }
};

int main() {
    vector<int> numSequences = {2, 8, 32, 128, 512, 2048, 8192}; // number of sequences to test
    size_t seqLength = 50; // fixed length for all sequences
    MultipleSequenceAlignment msa;

    cout << "n,runtime_seconds\n"; // CSV header

    for (int n : numSequences) {
        vector<SequenceInfo> sequences;
        for (int i = 0; i < n; ++i) {
            sequences.push_back({ "seq" + to_string(i+1), generateRandomSequence(seqLength) });
        }

        auto start = chrono::high_resolution_clock::now();
        vector<string> aligned = msa.divideAndConquerAlign(sequences);
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> duration = end - start;

        cout << n << "," << duration.count() << "\n";

        if (sequences.size() <= 10) { // print only for small sets
            cout << "Aligned Sequences:\n";
            for (size_t i = 0; i < aligned.size(); ++i) {
                cout << sequences[i].id << ": " << aligned[i] << "\n";
            }
        }
    }

    return 0;
}
