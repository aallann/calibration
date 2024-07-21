#ifndef OPTIONS_H
#define OPTIONS_H

#include <underlying.h>
#include <utils/__types__.h>

#include <memory>
#include <tuple>
#include <vector>

class Option {
   public:
    explicit Option(const Underlying &underlying)
    : underlying(underlying) {};
    virtual ~Option() = default;

   protected:
    Underlying underlying;
};

class EuropeanOption final : public Option {
   public:
    explicit EuropeanOption(const Underlying &underlying);
    ~EuropeanOption() override = default;

    inline static void setStrikes(std::shared_ptr<array> strikes) {
        EuropeanOption::strikes = strikes;
    };

    inline static void setTenors(std::shared_ptr<array> tau) {
        EuropeanOption::tau = tau;
    };

    inline std::tuple<uint, uint> getIndex() const {
        return this->idx;
    };

    inline double getStrike() const {
        auto [i, j] = this->idx;
        return (*strikes)(i, j);
    };

    inline double getTenor() const {
        auto [i, j] = this->idx;
        return (*tau)(i, j);
    };

    inline static std::vector<std::vector<uint>> getOptions() {
        return options;
    };

    friend class BlackScholes;

   private:
    std::tuple<uint, uint> idx;
    inline static std::shared_ptr<array> strikes = nullptr;
    inline static std::shared_ptr<array> tau = nullptr;
    inline static std::vector<std::vector<uint>> options = { {} };
};

#endif  // OPTIONS_H