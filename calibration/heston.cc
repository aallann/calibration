#include <heston.h>

static std::complex<double> _i_(0, 1);

Heston::Heston(
    double S,
    uint nDims,
    std::shared_ptr<array> vol,
    std::shared_ptr<array> tau,
    std::shared_ptr<array> strikes,
    std::shared_ptr<array> rfb,
    std::shared_ptr<array> rft,
    std::shared_ptr<array> deltas,
    std::shared_ptr<array> omegas
) : nDims(nDims), 
    Model(S, vol, tau, rfb, rft, deltas, omegas) {
    this->_discretiseSpace();
    Heston::strikes = strikes;
}

// auxiliary equations
complex_array Heston::c(complex_matrix &_xi_) const {
    return kappa.replicate(1, _xi_.cols()) - 
        ((rho * sigma * (an - am) * _i_).matrix() * _xi_).array();
}

complex_array Heston::h(complex_matrix &_xi_) const {
    return ((an - am).square().matrix() * _xi_.cwiseProduct(_xi_)).array() +
                ((_i_ * (an.square() - am.square())).matrix() * _xi_).array() -
            ((2. * _i_ * (rn - rm)).matrix() * _xi_).array() + 
                    (2 * rn).replicate(1, _xi_.cols()); 
}

complex_array Heston::s(
    complex_matrix &_xi_,
    complex_array &c,
    complex_array &h
) const {
    return sqrt(c.square() 
        + sigma.square().replicate(1, _xi_.cols()) * h);
}

complex_array Heston::g(
    complex_array &c,
    complex_array &s
) const {
    return (c - s) / (c + s);
}

complex_array Heston::A(
    complex_array &g,
    complex_array &s,
    int i
) const {
    return (1 - g * exp(-s * (*tau)(i, 0))) 
        / (1 - g);
}

complex_array Heston::B(
    complex_array &g,
    complex_array &s,
    int i
) const {
    return (1 - exp(-s * (*tau)(i, 0))) / 
        (1 - g * exp(-s * (*tau)(i, 0)));
}

// auxiliary partials
complex_array Heston::dg(
    std::variant<int, complex_array> &du,
    complex_array &dv, 
    complex_array &c, 
    complex_array &s
) const {
    return std::visit([&dv, &c, &s](auto&& arg) -> complex_array {
        return ((arg - dv) * (c + s) - (c - s) * (arg + dv)) 
            / (c + s).square();
    }, du);
}

complex_array Heston::dA(
    complex_array &du, 
    complex_array &dv, 
    complex_array &s,
    complex_array &g,
    int i
) const {
    return (*tau)(i, 0) * exp(-s * (*tau)(i, 0)) * g / 
        (1 - g) * du + (1 - exp(-s * (*tau)(i, 0))) / 
            (1 - g).square() * dv;
}

complex_array Heston::dB(
    complex_array &du,
    complex_array &dv, 
    complex_array &s,
    complex_array &g,
    int i
) const {
    return (*tau)(i, 0) * exp(-s * (*tau)(i, 0)) * (1 - g) / 
        (1 - g * exp(-s * (*tau)(i, 0))).square() * du + 
            (1 - exp(-s * (*tau)(i, 0))) * exp(-s * (*tau)(i, 0)) / 
                (1 - g * exp(-s * (*tau)(i, 0))).square() * dv;
}


// charf sensitivities to Heston parameters
complex_array Heston::Kappa(
    complex_array &c,
    complex_array &g,
    complex_array &s,
    complex_array &A,
    complex_array &B,
    int i
) const{
    std::variant<int, complex_array> dcdk = 1;
    complex_array dhdk = c / s;
    complex_array dgdk = this->dg(dcdk, dhdk, c, s);
    complex_array dAdk = this->dA(dhdk, dgdk, s, g, i);
    complex_array dBdk = this->dB(dhdk, dgdk, s, g, i);

    return ((vbar / sigma.square()).replicate(1, c.cols()) *
                ((c - s) * (*tau)(i, 0) - 2 * log(A) +
                    kappa.replicate(1, c.cols()) * (1 - dhdk) *
                        (*tau)(i, 0) -
                    2 * kappa.replicate(1, c.cols()) * dAdk / A) +
                (v0 / sigma.square()).replicate(1, c.cols()) *
                    ((1 - dhdk) * B + (c - s) * dBdk));
}

complex_array Heston::Vbar(
    complex_array &c,
    complex_array &s,
    complex_array &A,
    int i
) const {
    return (kappa / sigma.square()).replicate(1, c.cols()) *
        ((c - s) * (*tau)(i, 0) - 2 * log(A));
}

complex_array Heston::Sigma(
    complex_matrix &_xi_,
    complex_array &c,
    complex_array &h,
    complex_array &s,
    complex_array &g,
    complex_array &A,
    complex_array &B,
    int i
) const {
    complex_array dcdsigma = -_i_ * _xi_.array() * 
            ((an - am) * rho).replicate(1, c.cols());
    complex_array dhdsigma = (c * dcdsigma + 
            sigma.replicate(1, _xi_.cols()) * h) / (s);
    std::variant<int, complex_array> dc = dcdsigma;
    complex_array dgdsigma = this->dg(dc, dhdsigma, c, s);
    complex_array dAdsigma = this->dA(dhdsigma, dgdsigma, s, g, i);
    complex_array dBdsigma = this->dB(dhdsigma, dgdsigma, s, g, i);

    return (((kappa * vbar * (*tau)(i, 0))).replicate(1, c.cols()) *
        ((dcdsigma - dhdsigma) * sigma.square().replicate(1, c.cols()) -
            2 * (c - s) * sigma.replicate(1, c.cols())) /
        sigma.pow(4).replicate(1, c.cols()) -
            2 * (kappa * vbar).replicate(1, c.cols()) *
        (sigma.square().replicate(1, c.cols()) / A * dAdsigma -
            2 * sigma.replicate(1, c.cols()) * log(A)) /
        sigma.pow(4).replicate(1, c.cols()) +
            v0.replicate(1, c.cols()) *
        (((dcdsigma - dhdsigma) * B + (c - s) * dBdsigma) *
            sigma.square().replicate(1, c.cols()) -
        2 * (c - s) * B * sigma.replicate(1, c.cols())) /
            sigma.pow(4).replicate(1, c.cols()));
};
   
complex_array Heston::Rho(
    complex_matrix &_xi_, 
    complex_array &c,
    complex_array &s,
    complex_array &g,
    complex_array &A,
    complex_array &B,
    int i
) const {
    complex_array dcdrho = -_i_ * _xi_.array() * 
        ((an - am) * sigma).replicate(1, _xi_.cols());
    complex_array dhdrho = c / s * dcdrho;
    std::variant<int, complex_array> dc = dcdrho;
    complex_array dgdrho = this->dg(dc, dhdrho, c, s);
    complex_array dAdrho = this->dA(dhdrho, dgdrho, s, g, i);
    complex_array dBdrho = this->dB(dhdrho, dgdrho, s, g, i);

    return  (((kappa * vbar / sigma.square())).replicate(1, c.cols()) *
        ((dcdrho - dhdrho) * (*tau)(i, 0) - 2 * dAdrho / A) +
            (v0 / sigma.square()).replicate(1, c.cols()) * 
        ((dcdrho - dhdrho) * B + (c - s) * dBdrho));
}

complex_array Heston::V0(
    complex_array &c,
    complex_array &s,
    complex_array &B,
    int i
) const {
    return (c - s) * B / sigma.square().replicate(1, c.cols());
}

complex_array Heston::Am(
    complex_matrix &_xi_, 
    complex_array &c,
    complex_array &s,
    complex_array &g,
    complex_array &A,
    complex_array &B,
    int i
) const {
    complex_array dcdam = _i_ * _xi_.array() * 
        (rho * sigma).replicate(1, c.cols());
    complex_array dzdam = -2 * ((an - am).matrix() * _xi_.cwiseProduct(_xi_)).array() -
        2.0 * _i_ * c.array() * am.replicate(1, _xi_.cols());
    complex_array dhdam = (2 * c * dcdam + 
        sigma.replicate(1, _xi_.cols()).square() * dzdam) / (2 * s);
    std::variant<int, complex_array> dc = dcdam;
    complex_array dgdam = this->dg(dc, dhdam, c, s);
    complex_array dAdam = this->dA(dhdam, dgdam, s, g, i);
    complex_array dBdam = this->dB(dhdam, dgdam, s, g, i);

    return (((kappa * vbar / sigma.square())).replicate(1, c.cols()) *
        ((dcdam - dhdam) * (*tau)(i, 0) - 2 * dAdam / A) +
            (v0 / sigma.square()).replicate(1, c.cols()) *
                ((dcdam - dhdam) * B + (c - s) * dBdam));
};

complex_array Heston::An(
    complex_matrix &_xi_,
    complex_array &c,
    complex_array &s,
    complex_array &g,
    complex_array &A,
    complex_array &B,
    int i
) const { 
    complex_array dcdan = -_i_ * _xi_.array() * 
        (rho * sigma).replicate(1, _xi_.cols());
    complex_array dzdan = 2 * ((an - am).matrix() * _xi_.cwiseProduct(_xi_)).array() +
        2.0 * _i_ * _xi_.array() * an.replicate(1, _xi_.cols());
    complex_array dhdan = (2 * c * dcdan + 
        sigma.replicate(1, _xi_.cols()).square() * dzdan) / (2 * s);
    std::variant<int, complex_array> dc = dcdan;
    complex_array dgdan = this->dg(dc, dhdan, c, s);
    complex_array dAdan = this->dA(dhdan, dgdan, s, g, i);
    complex_array dBdan = this->dB(dhdan, dgdan, s, g, i);

    return (((kappa * vbar / sigma.square())).replicate(1, c.cols()) *
        ((dcdan - dhdan) * (*tau)(i, 0) - 2 * dAdan / A) +
            (v0 / sigma.square()).replicate(1, c.cols()) *
                ((dcdan - dhdan) * B + (c - s) * dBdan));
};

complex_array Heston::Rm(
    complex_matrix &_xi_, 
    complex_array &c,
    complex_array &s,
    complex_array &g,
    complex_array &A,
    complex_array &B,
    int i
) const { 
    complex_array dhdrm = _i_ * _xi_.array() * 
        sigma.replicate(1, _xi_.cols()).square() / s;
    std::variant<int, complex_array> dc = 0;
    complex_array dgdrm = this->dg(dc, dhdrm, c, s);
    complex_array dAdrm = this->dA(dhdrm, dgdrm, s, g, i);
    complex_array dBdrm = this->dB(dhdrm, dgdrm, s, g, i);

    return (((kappa * vbar / sigma.square())).replicate(1, c.cols()) *
        (-dhdrm * (*tau)(i, 0) - 2 * dAdrm / A) +
            (v0 / sigma.square()).replicate(1, c.cols()) *
                (-dhdrm * B + (c - s) * dBdrm));
};

complex_array Heston::Rn(
    complex_matrix &_xi_,
    complex_array &c,
    complex_array &s,
    complex_array &g,
    complex_array &A,
    complex_array &B,
    int i
) const {
    complex_array dhdrn = (1 - _i_ * _xi_.array()) * 
    sigma.replicate(1, _xi_.cols()).square() / s;
    std::variant<int, complex_array> dc = 0;
    complex_array dgdrn = this->dg(dc, dhdrn, c, s);
    complex_array dAdrn = this->dA(dhdrn, dgdrn, s, g, i);
    complex_array dBdrn = this->dB(dhdrn, dgdrn, s, g, i);

    return (((kappa * vbar / sigma.square())).replicate(1, c.cols()) *
        (-dhdrn * (*tau)(i, 0) - 2 * dAdrn / A) + (v0 / sigma.square())
            .replicate(1, c.cols()) *
                (-dhdrn * B + (c - s) * dBdrn));
};

complex_array Heston::Rm_init(
    complex_matrix &_xi_,
    int i
) const { return -_i_ * (*tau)(i, 0) * _xi_.array(); };

complex_array Heston::Rn_init(
    complex_matrix &_xi_, 
    int i
    ) const { return (-this->Rm_init(_xi_, i)) - 1; };

complex_array Heston::phi(
    complex_matrix &_xi_,
    complex_array &c,
    complex_array &s,
    complex_array &A,
    complex_array &B,
    int i
) const {
    return exp(
        (-rn_init.replicate(1, _xi_.cols()) + 
            ((_i_ * (rn_init - rm_init)).matrix() * _xi_).array() * (*tau)(i, 0) +
        ((kappa * vbar / sigma.square()).matrix().transpose() *
            ((c - s) * (*tau)(i, 0) - 2 * log(A)).matrix()).array() +
        ((v0 / sigma.square()).matrix().transpose() * ((c - s) * B).matrix())
            .array()));
}

array Heston::price(const array &p) {
    this->_updateParams(p);

    array prices((*tau).rows(), (*tau).cols());
    for (int i = 0; i < (*tau).rows(); i++) {
        double alpha = -4 * (*omegas)(i);       // damping
        complex_matrix _xi_ = xi + _i_ * alpha; // shifted xi

        complex_array c = this->c(_xi_);
        complex_array h = this->h(_xi_);
        complex_array s = this->s(_xi_, c, h);
        complex_array g = this->g(c, s);

        complex_array A = this->A(g, s, i);
        complex_array B = this->B(g, s, i);

        complex_array _phi_ = this->phi(_xi_, c, s, A, B, i);

        array moneyness = (array(1, 1) << (*strikes)(i, 0) / S).finished();
        complex_array payoff = S * 
            ((moneyness.replicate(1, xi.cols())).pow(_i_ * xi + alpha + 1) /
                ((_i_ * xi + alpha) * (_i_ * xi + alpha + 1)));

        array integrand = real((conj(payoff) * _phi_));
        prices.row(i) = integrand.sum() * dxi / (2 * M_PI);
    }

    return prices;
}

array Heston::gradient(const array &p) {
    this->_updateParams(p);

    array priceGradient(p.rows(), (*tau).rows());
    for (int i = 0; i < (*tau).rows(); i++) {
        double alpha = -4 * (*omegas)(i); // damping
        complex_matrix _xi_ = xi + _i_ * alpha;

        complex_array c = this->c(_xi_);
        complex_array h = this->h(_xi_);
        complex_array s = this->s(_xi_, c, h);
        complex_array g = this->g(c, s);

        complex_array A = this->A(g, s, i);
        complex_array B = this->B(g, s, i);

        complex_array grad = array(p.rows(), _xi_.cols());

        grad.block(0, 0, nDims, _xi_.cols()) = this->Kappa(c, g, s, A, B, i);
        grad.block(nDims, 0, nDims, _xi_.cols()) = this->Vbar(c, s, A, i);
        grad.block(2 * nDims, 0, nDims, _xi_.cols()) = this->Sigma(_xi_, c, h, s, g, A, B, i);
        grad.block(3 * nDims, 0, nDims, _xi_.cols()) = this->Rho(_xi_, c, s, g, A, B, i);
        grad.block(4 * nDims, 0, nDims, _xi_.cols()) = this->V0(c, s, B, i);
        grad.block(5 * nDims, 0, nDims, _xi_.cols()) = this->Am(_xi_, c, s, g, A, B, i);
        grad.block(6 * nDims, 0, nDims, _xi_.cols()) = this->An(_xi_, c, s, g, A, B, i);
        grad.block(7 * nDims, 0, nDims, _xi_.cols()) = this->Rm(_xi_, c, s, g, A, B, i);
        grad.block(8 * nDims, 0, nDims, _xi_.cols()) = this->Rn(_xi_, c, s, g, A, B, i);
        grad.block(9 * nDims, 0, 1, _xi_.cols()) = this->Rm_init(_xi_, i);
        grad.block(9 * nDims + 1, 0, 1, _xi_.cols()) = this->Rn_init(_xi_, i); 

        complex_array _phi_ = this->phi(_xi_, c, s, A, B, i);
        
        array moneyness = (array(1, 1) << (*strikes)(i, 0) / S).finished();
        complex_array payoff =
            S * ((moneyness.replicate(1, xi.cols())).pow(_i_ * xi + alpha + 1) /
                    ((_i_ * xi + alpha) * (_i_ * xi + alpha + 1)));

        array integrand = real(conj(payoff)
            .replicate(grad.rows(), 1) * (grad * _phi_.replicate(grad.rows(), 1)));
        priceGradient.col(i) = (integrand.rowwise().sum() * dxi / (2 * M_PI));
    }

    return priceGradient.transpose();
}

// internal state 
void Heston::_updateParams(const array &p) {
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

// internal state 
void Heston::_discretiseSpace(){
    xWidth = 20;
    nGrid = std::pow(2, 10);
    N2 = 0.5 * static_cast<double>(nGrid);
    upperBound = static_cast<double>(xWidth) / 2.0;
    dxi = M_PI / upperBound;
    xi.row(0) = dxi * xi.row(0).setLinSpaced(nGrid, -N2, N2);
}