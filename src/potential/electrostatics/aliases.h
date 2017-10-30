#ifndef ALIASES_H
#define ALIASES_H


#include <tuple>
#include <vector>


// Type aliases, unless you like scope resolution operators
using Id = std::string;

template <typename T, typename U>
using Pair = std::pair<T,U>;

template <typename T>
using Vec = std::vector<T>;

template <typename K, typename V>
using AssocVec = Vec<Pair<K,V>>;

using IndexPair = Pair<size_t, size_t>;


#endif // ALIASES_H
