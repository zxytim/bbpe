#include <fstream>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <chrono>


typedef unsigned int uint;
typedef unsigned long long ull;

// from https://stackoverflow.com/questions/2193544/how-to-print-additional-information-when-assert-fails
#define ASSERT(left, operator, right) \
do { \
    if (!((left) operator (right))) { \
        std::ostringstream ss; \
        ss << "ASSERT FAILED: " << #left << #operator << #right << " @ " << __FILE__ << " (" << __LINE__ << "). " << #left << "=" << (left) << "; " << #right << "=" << (right) << std::endl; \
        throw std::runtime_error(ss.str()); \
    } \
} while(0)


std::vector<uint> read_binary_file(const std::string& filename) {
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
        if (!file.read(reinterpret_cast<char*>(&byte), 1)) {
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

// read all files in a directory and return a concatenated vector of their contents, using std::filesystem routines
std::vector<uint> read_binary_files(const std::string& directory) {
    std::vector<uint> data;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        std::vector<uint> file_data = read_binary_file(entry.path());
        data.insert(data.end(), file_data.begin(), file_data.end());
    }
    return data;
}



// Typedef vocab to be a vector of pairs of ints.
typedef std::vector<std::pair<uint, uint>> vocab_t;

// Typedef frequencies to be an unordered map of pairs of ints (compressed to a single int) to ints.
typedef std::unordered_map<uint, uint> freqs_t;


int main()
{
    std::vector<uint> data = read_binary_file("bpe_data.bin");
    //std::vector<uint> data = read_binary_files("bpe_data");

    // print number of bytes
    std::cout << "Total number of bytes: " << data.size() << std::endl;

    // Create the initial vocabulary of 0-255.
    vocab_t vocab;
    for (uint i = 0; i < 256; ++i)
    {
        vocab.push_back(std::make_pair(i, -1));
    }

    std::vector<ull> frequencies;
    for (uint i = 0; i < 256; ++i)
    {
        frequencies.push_back(0);
    }

    for (size_t i = 0; i < data.size() - 1; ++i)
    {
        frequencies[data[i]]++;
    }

    /*
    std::vector<std::string> vocab_strings;
    for (uint i = 0; i < 256; ++i)
    {
        vocab_strings.push_back(std::string(1, static_cast<char>(i)));
    }
    */

    // Set the vocabulary size.
    uint vocab_size = 265;

    // Record the start time
    auto start = std::chrono::high_resolution_clock::now();

    // If vocab.txt exists, move it to vocab.txt.old
    std::ifstream infile("vocab.txt");
    if (infile.good()) {
        std::rename("vocab.txt", "vocab.txt.old");
    }
    
    // Write the initial vocabulary to a file.
    std::ofstream vocab_file("vocab.txt");
    for (uint i = 0; i < 256; ++i)
    {
        vocab_file << vocab[i].first << " " << vocab[i].second << " " << frequencies[i] << std::endl;
    }

    // While the vocabulary size is less than the desired size.
    while (vocab.size() < vocab_size)
    {
        // Start time of current iteration
        auto start_iter = std::chrono::high_resolution_clock::now();

        // Print the current vocabulary size.
        std::cout << "Vocab size: " << vocab.size() << std::endl;

        // Count the frequencies of bigrams.
        freqs_t freqs;
        uint vocab_base = vocab.size() + 1;
        for (size_t i = 0; i < data.size() - 1; ++i)
        {
            auto pair_idx = data[i] * vocab_base + data[i + 1];
            // assert data[i] < vocab_base and data[i + 1] < vocab_base
            ASSERT(data[i], <, vocab_base);
            ASSERT(data[i + 1], <, vocab_base);

            freqs[pair_idx]++;
        }

        // Find the most frequent bigram.
        auto most_common = std::max_element(freqs.begin(), freqs.end(),
                                            [](const auto &a, const auto &b)
                                            { return a.second < b.second; });

        // Print the most frequent bigram.
        const std::pair<uint, uint> most_common_pair = std::make_pair(most_common->first / vocab_base, most_common->first % vocab_base);
        std::cout << "Most frequent bigram: (" << most_common_pair.first << ", " << most_common_pair.second << ")" << std::endl;
        std::cout << "Frequency: " << most_common->second << std::endl;

        // Add the most frequent bigram to the vocabulary.
        vocab.push_back(most_common_pair);
        frequencies.push_back(most_common->second);

        // Write to vocab.txt incrementally
        std::ofstream vocab_file2;
        vocab_file2.open("vocab.txt", std::ios::app);
        vocab_file2 << most_common_pair.first << " " << most_common_pair.second << " " << most_common->second << std::endl;
        
        /* vocab_strings.push_back(vocab_strings[most_common_pair.first] + vocab_strings[most_common_pair.second]);
        std::cout << "New token: " << vocab_strings[vocab_strings.size() - 1] << std::endl; */

        // Replace the most frequent bigram with the new token.
        std::vector<uint> new_data;
        for (size_t i = 0; i < data.size() - 1; ++i)
        {
            auto pair_idx = data[i] * vocab_base + data[i + 1];
            if (pair_idx == most_common->first)
            {
                new_data.push_back(vocab.size() - 1);
                i++;
            }
            else
            {
                new_data.push_back(data[i]);
            }

            if (i == data.size() - 2)
            {
                new_data.push_back(data[i + 1]);
            }
        }

        data = new_data;

        // Print the "Elapsed > ETA" in seconds with 1 decimal place.
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        auto eta = (elapsed / (vocab.size() - 256)) * (vocab_size - vocab.size());
        std::cout << "Duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start_iter).count() / 1000.0 << "s" << std::endl;
        std::cout << "Elapsed: " << elapsed / 1000.0 << "s > ETA: " << eta / 1000.0 << "s" << std::endl;       
    }


    // Save the vocabulary strings to a file.
    /*
    std::ofstream vocab_strings_file("vocab_strings.txt");
    for (const auto &str : vocab_strings)
    {
        vocab_strings_file << str << std::endl;
    }
    */

    return 0;
}
