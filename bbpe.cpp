#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <thread>

typedef unsigned int uint;
typedef unsigned long long ull;

// Typedef vocab to be a vector of pairs of ints.
typedef std::vector<std::pair<uint, uint>> vocab_t;

// Typedef frequencies to be an unordered map of pairs of ints (compressed to a
// single int) to ints.
typedef std::unordered_map<ull, ull> freq_count_t;

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

std::vector<uint> read_binary_file(const std::string &filename) {
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
    std::vector<uint> data(file_size);

    // Read the bytes from the file and store them as unsigned integers.
    for (std::size_t i = 0; i < file_size; ++i) {
        unsigned char byte;
        if (!file.read(reinterpret_cast<char *>(&byte), 1)) {
            throw std::runtime_error("Failed to read from file");
        }
        data[i] = static_cast<uint>(byte);
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

int main() {
    Timer timer;
    auto tokens = read_binary_file("bpe_data.bin");
    size_t num_tokens = tokens.size();

    timer.tick("Read file");

    // print number of bytes
    std::cout << "Number of bytes: " << tokens.size() << std::endl;

    // Initialize
    vocab_t vocab;
    for (uint i = 0; i < 256; i++) {
        vocab.push_back(std::make_pair(i, -1));
    }

    std::vector<ull> frequencies(256);
    for (size_t i = 0; i < num_tokens; i++) {
        frequencies[tokens[i]]++;
    }

    // Target vocabulary size
    uint vocab_size_target = 265;

    // If vocab.txt exists, move it to vocab.txt.old
    std::ifstream infile("vocab.txt");
    if (infile.good()) {
        std::rename("vocab.txt", "vocab.txt.old");
    }

    // Write the initial vocabulary to a file.
    std::ofstream vocab_file("vocab.txt");
    for (uint i = 0; i < 256; ++i) {
        vocab_file << vocab[i].first << " " << vocab[i].second << " "
                   << frequencies[i] << std::endl;
    }

    timer.tick("Initialize");

    // Record the start time
    auto start_time = std::chrono::high_resolution_clock::now();

    while (vocab.size() < vocab_size_target) {
        // Start time of current iteration
        auto start_iter = std::chrono::high_resolution_clock::now();

        // (HEAVY LIFTING) Find the most frequent pair of tokens
        /* 
        freq_count_t freq_count;
        ull vocab_base = vocab.size() + 1;
        for (size_t i = 0; i < num_tokens - 1; i++) {
            uint token1 = tokens[i];
            uint token2 = tokens[i + 1];
            ull pair_idx = token1 * vocab_base + token2;
            freq_count[pair_idx]++;
        }
         */
        // (HEAVY LIFTING) Find the most frequent pair of tokens, using thread
        const size_t num_threads = 4;
        std::vector<std::thread> threads;
        std::vector<freq_count_t> freq_counts(num_threads);
        ull vocab_base = vocab.size() + 1;
        for (size_t i = 0; i < num_threads; i++) {
            threads.emplace_back([&tokens, &freq_counts, i, num_threads,
                                  vocab_base, num_tokens]() {
                for (size_t j = i; j < num_tokens - 1; j += num_threads) {
                    uint token1 = tokens[j];
                    uint token2 = tokens[j + 1];
                    ull pair_idx = token1 * vocab_base + token2;
                    freq_counts[i][pair_idx]++;
                }
            });
        }
        for (auto &t : threads) {
            t.join();
        }
        freq_count_t freq_count;
        for (auto &freq_count_i : freq_counts) {
            for (auto &pair : freq_count_i) {
                freq_count[pair.first] += pair.second;
            }
        }
        timer.tick("Count pairs");

        // Find the most frequent pair
        auto most_common =
            std::max_element(freq_count.begin(), freq_count.end(),
                             [](const auto &p1, const auto &p2) {
                                 return p1.second < p2.second;
                             });

        const auto most_common_pair = std::make_pair(
            most_common->first / vocab_base, most_common->first % vocab_base);
        timer.tick("Find most common pair");

        // Add the new token to the vocabulary
        vocab.push_back(most_common_pair);
        frequencies.push_back(most_common->second);

        // Write to vocab.txt incrementally
        std::ofstream vocab_file2;
        vocab_file2.open("vocab.txt", std::ios::app);
        vocab_file2 << most_common_pair.first << " " << most_common_pair.second
                    << " " << most_common->second << std::endl;

        // (HEAVY LIFTING) Replace the most frequent pair with the new token
        /* 
        size_t j = 0;
        uint new_token_id = vocab.size() - 1;
        for (size_t i = 0; i < num_tokens - 1; i++) {
            ull pair_idx = tokens[i] * vocab_base + tokens[i + 1];
            if (pair_idx == most_common->first) {
                tokens_new[j] = new_token_id;
                i++;
            } else {
                tokens_new[j] = tokens[i];
            }
            if (i == num_tokens - 2) {
                tokens_new[j + 1] = tokens[i + 1];
                j++;
            }
            j++;
        } */

        // (HEAVY LIFTING) Replace the most frequent pair with the new token, using thread
        uint new_token_id = vocab.size() - 1;
        std::vector<std::thread> threads2;
        std::vector<std::vector<uint>> tokens_new_threaded(num_threads);
        
        for (size_t i = 0; i < num_threads; i++) {
            threads2.emplace_back([&tokens, &tokens_new_threaded, i, num_threads,
                                   vocab_base, num_tokens, most_common,
                                   new_token_id]() {
                // thread i responsible for 
                // tokens[ceil(num_tokens/num_threads)*i:ceil(num_tokens/num_threads)*(i+1)]
                size_t block_size = std::ceil(num_tokens / num_threads);
                size_t start = block_size * i;
                size_t end = std::min(block_size * (i + 1), num_tokens);

                for (size_t j = start; j < end; j++) {
                    ull pair_idx = tokens[j] * vocab_base + tokens[j + 1];
                    if (pair_idx == most_common->first) {
                        tokens_new_threaded[i].push_back(new_token_id);
                        j++;
                    } else {
                        tokens_new_threaded[i].push_back(tokens[j]);
                    }
                    if (j == tokens.size() - 2) {
                        tokens_new_threaded[i].push_back(tokens[j + 1]);
                    }
                }
            });
        }
        for (auto &t : threads2) {
            t.join();
        }
        // Write back to tokens_new
        size_t j = 0;
        for (auto &tokens_new_thread : tokens_new_threaded) {
            std::copy(tokens_new_thread.begin(), tokens_new_thread.end(),
                      tokens.begin() + j);
            j += tokens_new_thread.size();
        }
        auto old_num_tokens = num_tokens;
        num_tokens = j;

        timer.tick("Replace pairs");

        // Print log
        std::cout << "==============================" << std::endl;
        std::cout << "Iteration " << vocab.size() - 256 << " finished. "
                  << std::endl;

        // Print the vocab merge
        std::cout << "Merge: " << most_common_pair.first << " "
                  << most_common_pair.second << " " << most_common->second
                  << std::endl;
        std::cout << "#tokens: " << old_num_tokens << "->" << num_tokens
                  << std::endl;

        // Duration of current iteration
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_iter = now - start_iter;
        std::cout << "Duration of current iteration: " << elapsed_iter.count()
                  << " s" << std::endl;

        // Time elased and ETA
        std::chrono::duration<double> elapsed_time = now - start_time;
        std::cout << "Time elapsed: " << elapsed_time.count() << " s"
                  << std::endl;
        double eta = elapsed_time.count() / (vocab.size() - 256) *
                     (vocab_size_target - 256);
        std::cout << "ETA: " << eta << " s" << std::endl;
    }
}