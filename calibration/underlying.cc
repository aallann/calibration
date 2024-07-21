#include <underlying.h>

Underlying::Underlying(double S) : S(S) {
    idx = underlyings.size();
    underlyings.emplace_back(idx);
}