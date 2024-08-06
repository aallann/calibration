#ifndef YIELD_CURVE_H
#define YIELD_CURVE_H

#include <utils/__types__.h>

#include <memory>
#include <stdexcept>
#include <string>

class YieldCurve {
    public:
     YieldCurve(std::string currency, unsigned nDims);
     ~YieldCurve() = default;

     array s(array &h);

     arrat g(array &s);

     array A(
         array &g,
         array &s,
         int i
      ) const;

     array B(
         array &g,
         array &s,
         int i
      ) const;

     // auxiliary partial derivatives
     array dg(
         std::variant<int, array> &du, 
         array &dv, 
         array &c, 
         array &s
      ) const;

     array dA(
         array &du, 
         array &dv, 
         array &s,
         array &g,
         int i
      ) const;

     array dB(
         array &du,
         array &dv, 
         array &s,
         array &g,
         int i
      ) const;

     // yield curve sensitivities to Heston params
     array Kappa(
        array &s,
        array &g,
        array &A,
        array &B,
        int i
     ) const;

    array Vbar(
        array &s,
        array &g, 
        array &A,
        array &B,
        int i
     ) const;

     array Sigma(
        array &s, 
        array &g,
        array &A,
        array &B,
        int i
     ) const;

     array V0(
        array &s,
        array &g, 
        array &A,
        array &B,
        int i
     ) const;

     array R(
        array &s,
        array &g, 
        array &A,
        array &B,
        int i
     ) const;

     array (
        array &s,
        array &g, 
        array &A,
        array &B,
        int i
     ) const;

     array operator()(const array &p);
     
     array gradient(const array &p);

     inline static void setTenors(std::shared_ptr<array> tau) {
        YieldCurve::tau = tau;
     }

     inline static std::shared_ptr<array> tau = nullptr;

    private:
     void _updateParams(const array &p);

     array kappa, vbar, sigma, rho, v0, am, an, rm, rn;
     ndarray<double, 1, 1> rm_init, rn_init;

     std::string currency;
     unsigned nDims;  
}
#endif // YIELD_CURVE_H