huffman-codes
=============

A simple implementation of Huffman Codes in C++

The core algorithm is taken from the CLR book (Introduction of Algorithms)
Chapter 16.3, and directly used to implement the 'build_tree()' routine.

After the tree is built, a code table that maps a character to a binary
code is built from the tree, and used for encoding text. Decoding is done
by traversing the Huffman tree, as prescribed by the algorithm.

Binary codes are represented by std::vector<bool>, which is a specialized
vector that optimizes space.
