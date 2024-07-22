#ifndef BLACK_SCHOLES_H
#define BLACK_SCHOLES_H

#include <model.h>
#include <options.h>
#include <utils/standard_gaussian.h>

#include <cmath>

class BlackScholes : public Model {
   public:
    using Model::Model;
    ~BlackScholes() override = default;

    array d1(std::vector<EuropeanOption> &options) const;

    array d2(std::vector<EuropeanOption> &options) const;

    std::shared_ptr<array> strikeFromDelta();

    array deltaFromStrike();

    array price(std::vector<EuropeanOption> &options);

    array impVol(std::vector<EuropeanOption> &options);

    inline static void setStrikes(std::shared_ptr<array> strikes) {
        BlackScholes::strikes = strikes;
    }

    inline static void setPrices(array prices) {
        BlackScholes::prices = prices;
    }

    inline static array getPrices(){
        return BlackScholes::prices;
    }

   private:
    inline static StandardGaussian standardGaussian;
    inline static std::shared_ptr<array> strikes = nullptr;
    static array prices;
};

#endif  // BLACK_SCHOLES_H