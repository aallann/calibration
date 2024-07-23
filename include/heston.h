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
        std::shared_ptr<array> rfb = nullptr,
        std::shared_ptr<array> rft = nullptr,
        std::shared_ptr<array> deltas = nullptr,
        std::shared_ptr<array> omegas = nullptr
     );

     ~Heston() override = default;

     array price(const array &p) const;

     inline void _updateState(const array &p) {
        kappa = p.block(0, 0, nDims, 1);
        vbar = p.block(nDims, 0, nDims, 1);
        sigma = p.block(2 * nDims, 0, nDims, 1);
        rho = p.block(3 * nDims, 0, nDims, 1);
        v0 = p.block(4 * nDims, 0, nDims, 1);
        am = p.block(5 * nDims, 0, nDims, 1);
        an = p.block(6 * nDims, 0, nDims, 1);
        rm = p.block(7 * nDims, 0, nDims, 1);
        rn = p.block(8 * nDims, 0, nDims, 1);
        rm_init = (array(1, 1) << p(9 * nDims, 0)).finished();
        rn_init = (array(1, 1) << p(p.rows() - 1, 0)).finished();
     }

     inline void _discretiseSpace(){
        xWidth = 20;
        nGrid = std::pow(2, 10);
        N2 = 0.5 * static_cast<double>(nGrid);
        upperBound = static_cast<double>(xWidth) / 2.0;
        dxi = M_PI / upperBound;
        xi.row(0) = dxi * xi.row(0).setLinSpaced(nGrid, -N2, N2);
     }
     
     array kappa, vbar, sigma, rho, v0, am, an, rm, rn;
     ndarray<double, 1, 1> rm_init, rn_init;

     ndarray<double, 1, 1024> xi;
     double N2, dxi, upperBound, nGrid;
     uint nDims, xWidth;
};

#endif  // HESTON_H