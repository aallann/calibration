#ifndef UNDERLYING_H
#define UNDERLYING_H

#include <vector>

class Underlying final {
   public:
    explicit Underlying(double S);
    ~Underlying() = default;

    inline double getSpotPrice() const {
        return this->S;
    }

    inline unsigned getIndex() const {
        return this->idx;
    }

    inline static std::vector<unsigned> getUnderlyings() {
        return underlyings;
    }

    friend class EuropeanOption;

   private:
    double S;
    uint idx;
    inline static std::vector<uint> underlyings = {};
};

#endif  // UNDERLYING_H