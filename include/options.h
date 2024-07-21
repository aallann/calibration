#ifndef OPTIONS_H
#define OPTIONS_H

#include <underlying.h>
#include <utils/__types__.h>

#include <memory>
#include <tuple>
#include <vector>

class Option {
   public:
    Option(const Underlying &underlying) : underlying(underlying) {};
    virtual ~Option() = default;

   protected:
    Underlying underlying;
};

class EuropeanOption final : public Option {
   public:
    EuropeanOption(const Underlying &underlying);
    ~EuropeanOption() override = default;

    inline void static setStrikes(std::shared_ptr<array> strikes) {
        EuropeanOption::strikes = strikes;
    };

    inline void static setTenors(std::shared_ptr<array> tau) {
        EuropeanOption::tau = tau;
    };

    inline std::tuple<unsigned, unsigned> getIndex() const {
        return idx;
    };

    inline double getStrike() const {
        auto [i, j] = this->idx;
        return (*strikes)(i, j);
    };

    inline double getTenor() const {
        auto [i, j] = this->idx;
        return (*tau)(i, j);
    };

    inline static std::vector<std::vector<unsigned>> getOptions() {
        return options;
    };

    inline static std::shared_ptr<array> strikes = nullptr;
    inline static std::shared_ptr<array> tau = nullptr;
    inline static std::vector<std::vector<unsigned>> options = { {} };

    friend class BlackScholes;

   private:
    std::tuple<unsigned, unsigned> idx;
};

#endif  // OPTIONS_H