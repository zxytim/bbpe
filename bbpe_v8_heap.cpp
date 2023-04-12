#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// from
// https://stackoverflow.com/questions/2193544/how-to-print-additional-information-when-assert-fails
#define ASSERT(left, operator, right)                                                                  \
    do {                                                                                               \
        if (!((left) operator(right))) {                                                               \
            std::ostringstream ss;                                                                     \
            ss << "ASSERT FAILED: " << #left                                                           \
               << #                                                                                    \
                  operator<< #right << " @ " << __FILE__ << " "                                        \
                                                            "(" << __LINE__ << ")."                    \
                                                                               " " << #left << "=" <<( \
                                                                                   left)               \
               << "; " << #right << "=" << (right) << std::endl;                                       \
            throw std::runtime_error(ss.str());                                                        \
        }                                                                                              \
    } while (0)

std::vector<uint32_t> read_binary_file(const std::string &filename) {
    // Open the binary file in input mode and binary mode.
    std::ifstream file(filename, std::ios::in | std::ios::binary);

    // Check if the file was opened successfully.
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file");
    }

    // Seek to the end of the file to get its size.
    file.seekg(0, std::ios::end);
    std::streamsize file_size = file.tellg();
    file.seekg(0, std::ios::beg);

    // Create a vector to hold the unsigned integers.
    std::vector<uint32_t> data(file_size);

    // Read the bytes from the file and store them as unsigned integers.
    for (std::size_t i = 0; i < file_size; ++i) {
        unsigned char byte;
        if (!file.read(reinterpret_cast<char *>(&byte), 1)) {
            throw std::runtime_error("Failed to read from file");
        }
        data[i] = static_cast<uint32_t>(byte);
    }

    // Check every value in the vector is less than 256.
    for (std::size_t i = 0; i < file_size; ++i) {
        ASSERT(data[i], <, 256);
    }

    // Close the file.
    file.close();

    // Return the vector of unsigned integers.
    return data;
}

// Timer class t
class Timer {
  public:
    Timer() {
        start_time = std::chrono::high_resolution_clock::now();
        last_tick_time = start_time;
    }

    // Print time elapsed since the last call to print_time_elapsed.
    // return duration
    double tick(const std::string &message) {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_time - last_tick_time);
        std::cout << message << " took " << duration.count() << " ms"
                  << std::endl;
        last_tick_time = end_time;

        return duration.count() / 1000.0;
    }

  private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point last_tick_time;
};

// 0xe08c1d668b756f82;

const uint64_t GOLDEN_RATIO = 0x9e3779b97f4a7c15;   // golden ratio
const uint64_t LARGE_PRIME = 18446744073709551557u; // 2^64 - 59

uint64_t hash_scramble(uint64_t x) {
    // prevent 0 from being hashed to 0
    x *= GOLDEN_RATIO;
    x += LARGE_PRIME;
    x ^= x >> 32;
    return x;
}

/*
I need a hash merge function over the uint64_t domain that:
1. Associative: hash_merge(hash_merge(a, b), c) == hash_merge(a, hash_merge(b,
c))
2. Non-zero self-merging: hash_merge(a, a) != 0, if a != 0
3. No fixed point: hash_merge(a, b) not in {a, b}, no matter a or b is zero or
not.

It is vital that associativity holds even if a, b, c are the same.

Please analyze the problem thoroughly,  come up with an idea, provide a C++
implementation and a formal proof that these condition holds.

Note: GPT 3.5 does not come up with a good solution. I am writing one myself.
*/

uint64_t hash_merge(uint64_t a, uint64_t b) {
    const uint64_t GOLDEN_RATIO = 0x9e3779b97f4a7c15;
    return (a ^ b ^ GOLDEN_RATIO);
}

// Code as in "codebook".
class Code {
  public:
    Code *left, *right; // How this code is merged from two other codes.
    std::vector<uint8_t>
        chars; // the length of this code is the length of this vector

    inline size_t length() const { return chars.size(); }
    std::string to_string() const {
        std::string s;
        for (size_t i = 0; i < chars.size(); ++i) {
            s += chars[i];
        }
        return s;
    }

    uint64_t hash;
    uint64_t idx;                 // index in the codebook
    std::set<uint64_t> positions; // start positions of this code in the input

    // frequency is just the length of positions
    inline size_t frequency() const { return positions.size(); }
};

bool code_compare_func(Code *a, Code *b) {
    // Compare two codes.
    // Return true if a should be popped before b.
    if (a->frequency() != b->frequency()) {
        return a->frequency() > b->frequency();
    } else {
        return a->hash < b->hash;
    }
}

class CandidatePriorityQueue {
    // A priority queue of candidates.
    // It is a normal priority queue plus support for fast removal of elements.
  public:
    // The priority queue storage.
    std::vector<Code *> data;

    // The index of each code in the priority queue.
    std::unordered_map<Code *, size_t> idx_map;

    void bubble_up(size_t idx) {
        while (idx > 0) {
            size_t parent_idx = (idx - 1) / 2;
            if (code_compare_func(data[idx], data[parent_idx])) {
                std::swap(data[idx], data[parent_idx]);
                idx_map[data[idx]] = idx;
                idx_map[data[parent_idx]] = parent_idx;
                idx = parent_idx;
            } else {
                break;
            }
        }
    }

    void bubble_up(Code *code) {
        auto it = idx_map.find(code);
        assert(it != idx_map.end());
        size_t idx = it->second;
        bubble_up(idx);
    }

    void push(Code *code) {
        data.push_back(code);
        size_t cur_idx = data.size() - 1;
        idx_map[code] = cur_idx;

        // Bubble up.
        bubble_up(cur_idx);
    }

    void bubble_down(size_t idx) {
        while (idx < data.size()) {
            size_t left_idx = idx * 2 + 1;
            size_t right_idx = idx * 2 + 2;
            size_t min_idx = idx;
            if (left_idx < data.size() &&
                code_compare_func(data[left_idx], data[min_idx])) {
                min_idx = left_idx;
            }
            if (right_idx < data.size() &&
                code_compare_func(data[right_idx], data[min_idx])) {
                min_idx = right_idx;
            }
            if (min_idx != idx) {
                std::swap(data[idx], data[min_idx]);
                idx_map[data[idx]] = idx;
                idx_map[data[min_idx]] = min_idx;
                idx = min_idx;
            } else {
                break;
            }
        }
    }

    void bubble_down(Code *code) {
        auto it = idx_map.find(code);
        assert(it != idx_map.end());
        size_t idx = it->second;
        bubble_down(idx);
    }

    void remove(Code *code) {
        // Remove the code from the priority queue.
        // This is done by swapping the code with the last element,
        // and then pop the last element.
        auto it = idx_map.find(code);
        if (it == idx_map.end()) {
            return;
        }
        size_t idx = it->second;
        idx_map.erase(it);

        // Swap with the last element.
        size_t last_idx = data.size() - 1;
        std::swap(data[idx], data[last_idx]);
        idx_map[data[idx]] = idx;
        data.pop_back();

        // Bubble down.
        bubble_down(idx);
    }

    Code *top() { return data[0]; }

    void pop() { remove(data[0]); }

    size_t size() { return data.size(); }
};

// map from (code_a, code_b) -> hash(concat(code_a.chars, code_b.chars))
std::map<std::pair<Code *, Code *>, uint64_t> code_hash_cache;

std::unordered_map<size_t, std::vector<std::pair<Code *, Code *>>>
    hash2code_pairs;
std::unordered_map<std::string, Code *> str2code;

uint64_t compute_hash(const std::vector<uint8_t> &chars) {
    uint64_t hash = hash_scramble(chars[0]);
    for (size_t i = 1; i < chars.size(); ++i) {
        uint64_t x = hash_scramble(chars[i]);
        hash = hash * GOLDEN_RATIO + x;
    }
    return hash;
}

uint64_t code_hash_merge(Code *a, Code *b) {
    // Check if we have already computed the hash.
    auto it = code_hash_cache.find(std::make_pair(a, b));
    if (it != code_hash_cache.end()) {
        return it->second;
    }

    // We want to use a simple function such as
    //      hash_merge(a->hash, b->hash);
    // However, since did not find a good hash function,
    // we hash the merged chars instead.
    // This is not ideal, but it works.

    std::vector<uint8_t> chars;
    chars.insert(chars.end(), a->chars.begin(), a->chars.end());
    chars.insert(chars.end(), b->chars.begin(), b->chars.end());

    uint64_t hash = compute_hash(chars);

    // Cache the hash.
    assert(code_hash_cache.find(std::make_pair(a, b)) == code_hash_cache.end());
    hash2code_pairs[hash].push_back(std::make_pair(a, b));

    code_hash_cache[std::make_pair(a, b)] = hash;

    return hash;
}

int main() {
    Timer timer;
    // Hyperparameters.
    const size_t num_threads = 2;
    const size_t num_code = 256 + 5000; // The vocabulary size

    // Read the data from the file.
    const std::string fname = "./enwik9.10M";
    // auto data = read_binary_file("./bpe_data.bin");
    // auto data = read_binary_file("./bpe_data/enwik9.100M");

    std::cout << "Reading data from " << fname << std::endl;
    auto data = read_binary_file(fname);
    // auto data = read_binary_file("./bpe_data.ascii.1024");
    //  auto data = read_binary_file("./bpe_data.ascii.8");
    // auto data = read_binary_file("./test_case.0.in");
    //  data.resize(1000000);

    timer.tick("Read data");

    // Initialize 0-255 codes.
    std::vector<Code *> codes;

    // Map from hash to code.
    std::unordered_map<uint64_t, Code *> code_map;

    for (size_t i = 0; i < 256; ++i) {
        Code *code = new Code();
        code->idx = i;
        code->left = nullptr;
        code->right = nullptr;
        code->chars.push_back(i);
        code->hash = hash_scramble(i);
        code_map[code->hash] = code;
        codes.push_back(code);
        str2code[code->to_string()] = code;
    }
    timer.tick("Initialize 0-255 codes");

    // Initialize frequencies of first codes.
    for (size_t i = 0; i < data.size(); ++i) {
        codes[data[i]]->positions.insert(i);
    }
    timer.tick("Initialize frequencies of first codes");

    // Pool of candidates
    // std::set<Code *> candidates;
    CandidatePriorityQueue candidates;

    // Fill all codes to corresponding position.
    std::vector<Code *> done(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        done[i] = codes[data[i]];
    }
    timer.tick("Fill all codes to corresponding position");

    auto add_cand = [&](size_t pos, Code *code_a, Code *code_b) {
        ASSERT(pos, <, data.size());
        // uint64_t hash = hash_merge(code_a->hash, code_b->hash);
        uint64_t hash = code_hash_merge(code_a, code_b);

        Code *ret = nullptr;
        if (code_map.find(hash) == code_map.end()) {
            Code *code = new Code();
            code->left = code_a;
            code->right = code_b;
            code->chars.insert(code->chars.end(), code_a->chars.begin(),
                               code_a->chars.end());
            code->chars.insert(code->chars.end(), code_b->chars.begin(),
                               code_b->chars.end());
            code->hash = hash;
            code->positions.insert(pos);
            code_map[hash] = code;

            auto s = code->to_string();
            assert(str2code.find(s) == str2code.end());
            str2code[code->to_string()] = code;
            // candidates.insert(code);
            candidates.push(code);
            ret = code;
        } else {
            Code *code = code_map[hash];
            code->positions.insert(pos);
            candidates.bubble_up(code);
            ret = code;
        }
        return ret;
    };

    // Generate all consecutive pairs of codes as active codes.
    for (size_t i = 0; i < done.size() - 1; ++i) {
        Code *code_a = done[i];
        Code *code_b = done[i + 1];

        add_cand(i, code_a, code_b);
    }

    std::cout << "#candidates: " << candidates.size() << std::endl;
    timer.tick("Generate all pairs of codes as active codes");

    auto remove_cand = [&](Code *code, size_t position) {
        // If postition not in code->positions, do nothing.
        if (code->positions.find(position) == code->positions.end()) {
            return;
        }
        code->positions.erase(position);

        // Single character code has been added to the vocabulary mandatorily.
        // Such that it is not in the candidate pool.
        if (code->chars.size() == 1) {
            return;
        }
        // If is has been removed before, skip
        if (candidates.idx_map.find(code) == candidates.idx_map.end()) {
            return;
        }

        candidates.bubble_down(code);

        if (code->positions.size() == 0) {
            // candidates.erase(code);
            // delete code;
        }
    };

    // If vocab.txt exists, move it to vocab.txt.old
    std::ifstream infile("vocab.txt");
    if (infile.good()) {
        std::rename("vocab.txt", "vocab.txt.old");
    }

    std::ifstream fvocabstr("vocab.string.txt");
    if (fvocabstr.good()) {
        std::rename("vocab.string.txt", "vocab.string.txt.old");
    }

    // Write the initial vocabulary to a file.
    std::ofstream vocab_file("vocab.txt");
    for (size_t i = 0; i < codes.size(); ++i) {
        Code *c = codes[i];
        vocab_file << i << " " << -1 << " " << -1 << " " << c->frequency()
                   << std::endl;
    }

    const size_t print_interval = 20;
    // Main loop
    while (codes.size() < num_code) {
        if (codes.size() % print_interval == 0) {
            std::cout << "============================" << std::endl;
            std::cout << "Progress: " << codes.size() << "/" << num_code
                      << std::endl;
        }

        if (candidates.size() == 0) {
            std::cout << "No candidate left" << std::endl;
            break;
        }
        // Find the most frequent active code.
        /*
        Code *code = *candidates.begin();
        for (auto it = candidates.begin(); it != candidates.end(); ++it) {
            if ((*it)->frequency() > code->frequency()) {
                // TODO: conside length and alphabet
                code = *it;
            }
        }
        */
        Code *code = candidates.top();

#ifdef DEBUG_TIMING
        timer.tick("Find most frequent active code");
#endif

        // Add the code to the vocabulary.
        codes.push_back(code);
        code->idx = codes.size() - 1;

        // Remove the code from the active codes.
        // candidates.erase(code);
        candidates.pop();

        // Update the sequence.
        auto positions = code->positions;
        for (auto &pos : positions) {
            // If this specific incarnation of code is not in the sequence,
            // it is been removed by previous iterations.
            if (code->positions.find(pos) == code->positions.end())
                continue;

            size_t l = pos,
                   r = pos + code->length() - 1; // inclusive

            Code *g = code;
            size_t g_start = l;

            Code *a = done[l], *b = done[r];
            size_t a_start = l, b_start = l + a->length();

            Code *c = nullptr, *e = nullptr, *j = nullptr;
            size_t c_start = -1, e_start = -1, j_start = -1;
            if (l > 0) {
                c = done[l - 1];
                c_start = l - c->length();

                // e is the merge of c and a
                // auto hash_e = hash_merge(c->hash, a->hash);
                auto hash_e = code_hash_merge(c, a);
                // assert hash_e in code_map
                assert(code_map.find(hash_e) != code_map.end());
                e = code_map[hash_e];
                e_start = c_start;

                // j is the merge of c and g
                // auto hash_j = hash_merge(c->hash, g->hash);
                auto hash_j = code_hash_merge(c, g);
                // j = code_map[hash_j];
                j_start = c_start;
            }

            Code *d = nullptr, *f = nullptr, *k = nullptr;
            size_t d_start = -1, f_start = -1, k_start = -1;
            if (r < done.size() - 1) {
                d = done[r + 1];
                d_start = r + 1;

                // f is the merge of b and d
                // auto hash_f = hash_merge(b->hash, d->hash);
                auto hash_f = code_hash_merge(b, d);
                assert(code_map.find(hash_f) != code_map.end());
                f = code_map[hash_f];
                f_start = b_start;

                // k is the merge of g and d
                // auto hash_k = hash_merge(g->hash, d->hash);
                auto hash_k = code_hash_merge(g, d);
                // k = code_map[hash_k];
                k_start = g_start;
            }

            remove_cand(a, a_start);

            remove_cand(b, b_start);

            if (e)
                remove_cand(e, e_start);
            if (f)
                remove_cand(f, f_start);
            if (l > 0)
                add_cand(j_start, c, g);
            if (r < done.size() - 1)
                add_cand(k_start, g, d);

            done[l] = done[r] = g;
        }
#ifdef DEBUG_TIMING
        timer.tick("Update sequence");
#endif

#ifdef XDEBUG
        // Check all tokens cover the sequence exactly once.
        size_t now = 0;
        while (now < done.size()) {
            Code *code = done[now];
            now += code->length();
        }
        ASSERT(now, ==, done.size());

        timer.tick("Check all tokens cover the sequence exactly once");
#endif

        // output the code

        if (codes.size() % print_interval == 0) {
            std::cout << "Code: " << '"' << code->to_string() << '"'
                      << std::endl;
            std::cout << "Code length: " << code->length() << std::endl;
            std::cout << "Code frequency: " << code->frequency() << std::endl;
        }
        std::ofstream fvocabstr2;
        fvocabstr2.open("vocab.string.txt", std::ios::app);
        fvocabstr2 << code->to_string() << std::endl;

        // Write the code to a file.
        std::ofstream vocab_file2;
        vocab_file2.open("vocab.txt", std::ios::app);
        vocab_file2 << code->idx << " " << code->left->idx << " "
                    << code->right->idx << " "
                    << code->frequency() //<< " " << code->to_string()
                    << std::endl;
    }
    return 0;
}