#include <options.h>

EuropeanOption::EuropeanOption(const Underlying &underlying)
 : Option(underlying) {
    uint j = underlying.idx;

    if (options.size() <= j) {
        options.emplace_back(std::vector<uint>());
    };

    unsigned i = options[j].size();
    options[j].emplace_back(j);
    idx = std::make_tuple(i, j);
}