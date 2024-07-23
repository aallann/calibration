#ifndef HESTON_H
#define HESTON_H

#include <model.h>

#include <cmath>

class Heston : public Model{
    public:
     Heston(
        double S,
        uint nDims,
        std::shared_ptr<array> vol = nullptr,
        std::shared_ptr<array> tau = nullptr,
        std::shared_ptr<array> strikes = nullptr,
        std::shared_ptr<array> rfb = nullptr,
        std::shared_ptr<array> rft = nullptr,
        std::shared_ptr<array> deltas = nullptr,
        std::shared_ptr<array> omegas = nullptr
     );

     ~Heston() override = default;

     // auxiliary functions 
     complex_array c(complex_matrix &_xi_) const;
     complex_array h(complex_matrix &_xi_) const;
     complex_array s(complex_matrix &_xi_, complex_array &c, complex_array &h) const;
     complex_array g(complex_array &c, complex_array &d) const;
     complex_array A(complex_array &g, complex_array &s, int i) const;
     complex_array B(complex_array &g, complex_array &s, int i) const;

     complex_array phi(complex_matrix &_xi_, int i) const;
     array price(const array &p);

    inline static void setStrikes(std::shared_ptr<array> strikes_ptr) {
        Heston::strikes = strikes_ptr;
    };

    inline static void setPrices(std::shared_ptr<array> prices_ptr) {
        Heston::prices = prices_ptr;
    };

    inline static array getPrices(){
        return (*Heston::prices);
    }

    private:
     inline void _updateState(const array &p);
     inline void _discretiseSpace();

     inline static std::shared_ptr<array> prices = nullptr;
     inline static std::shared_ptr<array> strikes = nullptr;

     array kappa, vbar, sigma, rho, v0, am, an, rm, rn;
     ndarray<double, 1, 1> rm_init, rn_init;
     ndarray<double, 1, 1024> xi;
     double N2, dxi, upperBound, nGrid;
     uint nDims, xWidth;
};

#endif  // HESTON_H