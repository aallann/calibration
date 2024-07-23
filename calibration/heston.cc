#include <heston.h>

Heston::Heston(
    double S,
    uint nDims,
    std::shared_ptr<array> vol,
    std::shared_ptr<array> tau,
    std::shared_ptr<array> rfb,
    std::shared_ptr<array> rft,
    std::shared_ptr<array> deltas,
    std::shared_ptr<array> omegas
) : nDims(nDims),
    Model(S, vol, tau, rfb, rft, deltas, omegas) {
    this->_discretiseSpace();
};