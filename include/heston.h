#ifndef HESTON_H
#define HESTON_H

#include <model.h>

#include <cmath>
#include <variant>

class Heston : public Model{
    public:
     Heston(
        double S,
        uint nDims,
        std::shared_ptr<array> vol,
        std::shared_ptr<array> tau,
        std::shared_ptr<array> strikes,
        std::shared_ptr<array> rfb,
        std::shared_ptr<array> rft,
        std::shared_ptr<array> deltas,
        std::shared_ptr<array> omegas
     );

     ~Heston() override = default;

     // auxiliary equations
     complex_array c(complex_matrix &_xi_) const;

     complex_array h(complex_matrix &_xi_) const;

     complex_array s(
         complex_matrix &_xi_,
         complex_array &c, 
         complex_array &h
      ) const;

     complex_array g(
         complex_array &c,
         complex_array &d
      ) const;

     complex_array A(
         complex_array &g,
         complex_array &s,
         int i
      ) const;

     complex_array B(
         complex_array &g,
         complex_array &s,
         int i
      ) const;

     // auxiliary partial derivatives
     complex_array dg(
         std::variant<int, complex_array> &du, 
         complex_array &dv, 
         complex_array &c, 
         complex_array &s
      ) const;

     complex_array dA(
         complex_array &du, 
         complex_array &dv, 
         complex_array &s,
         complex_array &g,
         int i
      ) const;

     complex_array dB(
         complex_array &du,
         complex_array &dv, 
         complex_array &s,
         complex_array &g,
         int i
      ) const;

     // charf sensitivities to Heston parameters
       complex_array Kappa(
         complex_array &c,
         complex_array &g,
         complex_array &s,
         complex_array &A,
         complex_array &B,
         int i
      ) const;

     complex_array Vbar(
         complex_array &c,
         complex_array &s,
         complex_array &A,
         int i
      ) const;

     complex_array Sigma(
         complex_matrix &_xi_,
         complex_array &c,
         complex_array &h,
         complex_array &s,
         complex_array &g,
         complex_array &A,
         complex_array &B,
         int i
      ) const;
   
      complex_array Rho(
         complex_matrix &_xi_,
         complex_array &c,
         complex_array &s,
         complex_array &g,
         complex_array &A,
         complex_array &B,
         int i
      ) const;

      complex_array V0(
         complex_array &c,
         complex_array &s,
         complex_array &B,
         int i
      ) const;

      complex_array Am(
         complex_matrix &_xi_,
         complex_array &c,
         complex_array &s,
         complex_array &g,
         complex_array &A,
         complex_array &B,
         int i
      ) const;

      complex_array An(
         complex_matrix &_xi_,
         complex_array &c,
         complex_array &s,
         complex_array &g,
         complex_array &A,
         complex_array &B,
         int i
      ) const;

      complex_array Rm(
         complex_matrix &_xi_,
         complex_array &c,
         complex_array &s,
         complex_array &g,
         complex_array &A,
         complex_array &B,
         int i
      ) const;

      complex_array Rn(
         complex_matrix &_xi_,
         complex_array &c,
         complex_array &s,
         complex_array &g,
         complex_array &A,
         complex_array &B,
         int i
      ) const;

      complex_array Rn_init(
         complex_matrix &_xi_,
         int i
      ) const;

      complex_array Rm_init(
         complex_matrix &_xi_,
         int i
      ) const;
   
      complex_array phi(
         complex_matrix &_xi_,
         complex_array &c,
         complex_array &s,
         complex_array &A,
         complex_array &B,
         int i
      ) const;

     array price(const array &p);

     array gradient(const array &p);

     inline static void setStrikes(std::shared_ptr<array> strikes_ptr) {
        Heston::strikes = strikes_ptr;
     }

     inline static void setPrices(std::shared_ptr<array> prices_ptr) {
        Heston::prices = prices_ptr;
     }

     inline static array getPrices(){
        return (*Heston::prices);
     }

    private:
     void _updateParams(const array &p);
     void _discretiseSpace();

     inline static std::shared_ptr<array> prices = nullptr;
     inline static std::shared_ptr<array> strikes = nullptr;

     array kappa, vbar, sigma, rho, v0, am, an, rm, rn;
     ndarray<double, 1, 1> rm_init, rn_init;
     ndarray<double, 1, 1024> xi;
     double N2, dxi, upperBound, nGrid;
     uint nDims, xWidth;
};

#endif  // HESTON_H