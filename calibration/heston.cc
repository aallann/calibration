#include <heston.h>

static std::complex<double> _i_(0, 1);

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
}

complex_array Heston::c(complex_matrix &_xi_) const {
    return (-(rho * sigma * (an - am) * _i_).matrix() * _xi_).array().colwise() 
            + kappa;
}

complex_array Heston::h(complex_matrix &_xi_) const {
    return (((an - am).square().matrix() * _xi_.cwiseProduct(_xi_)).array() +
        ((_i_ * (an.square() - am.square())).matrix() * _xi_).array() -
            ((2. * _i_ * (rn - rm)).matrix() * _xi_).array()).colwise() + (2 * rn);
}

complex_array Heston::s(
    complex_matrix &_xi_, complex_array &c, complex_array &h
) const {
    return sqrt(c.square().colwise() + sigma.square() * h);
}

complex_array Heston::g(complex_array &c, complex_array &s) const {
    return (c - s) / (c + s);
}

complex_array Heston::A(complex_array &g, complex_array &s, int i) const {
    return (1 - g * exp(-s * (*tau)(i, 0))) / (1 - g);
};

complex_array Heston::B(complex_array &g, complex_array &s, int i) const {
    return (1 - exp(-s * (*tau)(i, 0))) / (1 - g * exp(-s * (*tau)(i, 0)));
};

complex_array Heston::phi(complex_matrix &_xi_, int i) const {
    complex_array c = this->c(_xi_);
    complex_array h = this->h(_xi_);
    complex_array s = this->s(_xi_, c, h);
    complex_array g = this->g(c, s);

    complex_array A = this->A(g, s, i);
    complex_array B = this->B(g, s, i);

    return exp(
        (-rn_init.replicate(1, _xi_.cols()) + 
            ((_i_ * (rn_init - rm_init)).matrix() * _xi_).array() * (*tau)(i, 0) +
        ((kappa * vbar / sigma.square()).matrix().transpose() *
            ((c - s) * (*tau)(i, 0) - 2 * log(A)).matrix())
                 .array() +
             ((v0 / sigma.square()).matrix().transpose() * ((c - s) * B).matrix())
                 .array()));
}

array Heston::price(const array &p) {
    this->_updateState(p);

    array prices((*tau).rows(), (*tau).cols());
    for (int i = 0; i < (*tau).rows(); i++) {
        double alpha = -4 * (*omegas)(i);  // damping
        complex_matrix _xi_ = xi + _i_ * alpha;

        complex_array _phi_ = this->phi(_xi_, i);

        array moneyness = (array(1, 1) << (*strikes)(i, 0) / S).finished();
        complex_array payoff = S * ((moneyness.replicate(1, xi.cols())).pow(_i_ * xi + alpha + 1) /
                                    ((_i_ * xi + alpha) * (_i_ * xi + alpha + 1)));

        array integrand = real((conj(payoff) * _phi_));
        prices.row(i) = integrand.sum() * dxi / (2 * M_PI);
    }

    return array(1, 1);
}