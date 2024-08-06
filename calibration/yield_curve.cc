#include<yield_curve.h>

YieldCurve::YieldCurve(std::string currency, unsigned nDims) : currency(currency), nDims(nDims) {
    if (currency != "base" && currency != "term") {
            throw std::invalid_argument("Currency may only be 'base' or 'term'");
    }
}

array YieldCurve::s(array &h) {
    return sqrt(kappa.square() + sigma.square() * h);
}

array YieldCurve::g(array &s) {
    return (kappa - s) / (kappa + s);
}

array YieldCurve::A(
    array &s,
    array &g,
    int i
) const {
    return (1 - g * exp(-s * (*tau)(i, 0))) 
        / (1 - g);
}

array YieldCurve::B(
    array &s,
    array &g,
    int i
) const {
    return (1 - exp(-s * (*tau)(i, 0))) / 
        (1 - g * exp(-s * (*tau)(i, 0)));
}

// auxiliary partials
array YieldCurve::dg(
    std::variant<int, array> &du,
    array &dv, 
    array &c, 
    array &s
) const {
    return std::visit([&dv, &c, &s](auto&& arg) -> array {
        return ((arg - dv) * (c + s) - (c - s) * (arg + dv)) 
            / (c + s).square();
    }, du);
}

array YieldCurve::dA(
    array &du, 
    array &dv, 
    array &s,
    array &g,
    int i
) const {
    return (*tau)(i, 0) * exp(-s * (*tau)(i, 0)) * g / 
        (1 - g) * du + (1 - exp(-s * (*tau)(i, 0))) / 
            (1 - g).square() * dv;
}

array YieldCurve::dB(
    array &du,
    array &dv, 
    array &s,
    array &g,
    int i
) const {
    return (*tau)(i, 0) * exp(-s * (*tau)(i, 0)) * (1 - g) / 
        (1 - g * exp(-s * (*tau)(i, 0))).square() * du + 
            (1 - exp(-s * (*tau)(i, 0))) * exp(-s * (*tau)(i, 0)) / 
        (1 - g * exp(-s * (*tau)(i, 0))).square() * dv;
}

     // yield curve sensitivities to Heston params
array YieldCurve::Kappa(
   array &s,
   array &g,
   array &A,
   array &B,
   int i
) const {
    array dsdk = kappa / s;
    array dgdk = this->dg(1, dsdk, kappa, s);
    array dAdk = this->dA(dsdk, dgdk, s, g, i);
    array dBdk = this->dB(dsdk, dgdk, s, g, i);
    
    return -vbar / sigma.square() * (
        (kappa - s) + kappa * (1 - dsdk) - 2 * 
            (kappa * dAdk / A + log(A)) / (*tau)(i, 0)
        ) - v0 / sigma.square() * (
            (1 - dsdk) * B + (kappa - s) * dBdk
        ) / (*tau)(i, 0);
}

array YieldCurve::Sigma(
    array &s, 
    array &g,
    array &A,
    array &B,
    int i
) const {
    array dsdsigma = sigma / s * h;
    array dgdsigma = this->dg(0, dsdsigma, s);
    array dAdsigma = this->dA(dsdsigma, dgdsigma, s, g, i);
    array dBdsigma = this->dB(dsdsigma, dgdsigma, s, g, i); 
    
    return kappa * vbar * (*tau)(i, 0) * 
        (dsdsigma * sigma.square() + 2 * (kappa - s) * sigma) / 
            sigma.pow(4) / (*tau)(i, 0) + 2 * kappa * vbar * 
        (sigma.square() / A * dAdsigma - 2 * sigma * log(A)) / 
            sigma.pow(4) / (*tau)(i, 0) + v0 * (
        (dsdsigma * B - (kappa - s) * dBdsigma) * 
            sigma.square() + 2 * (kappa - s) * B * sigma
        ) / sigma.pow(4) / (*tau).row(i);
}

array YieldCurve::R(
    array &s,
    array &g, 
    array &A,
    array &B,
    int i
) const {
    array dsdr = sigma.square() / s;
    array dgdr = this->dg(0, dsdr, s);
    array dAdr = this->dA(dsdr, dgdr, s, g, i);
    array dBdr = this->dB(dsdr, dgdr, s, g, i);

    return kappa * vbar / sigma.square() * 
        (dsdr + 2 * dAdr / A / (*tau).row(i)) +
            v0 / sigma.square() * 
        (dsdr * B - (kappa - s) * dBdr / (*tau).row(i));
}


YieldCurve::operator()(const array &p) {
    this->_updateParams(p);

    array h(1, 1);
    double r;
    if (this->currency == "base") {
        r = rm_init;
        h << 2 * rm;
    } else if (this->currency == "term") {
        r = rn_init;
        h << 2 * rn;
    };

    array s = this->s(h);
    array g = this->g(s);

    array yield = array(tau->rows(), tau->cols());
    for (int i = 0; i <= (*tau)->rows()) {
        array A = this->A(s, g, i);
        array B = this->B(s, g, i);

        yield.row(i) = r - (
                (kappa * vbar / sigma.square()).matrix().transpose() *
                    (kappa - s - 2 * log(A)) / (*tau)(i, 0).matrix()
            ).array() - 
                ((v0 / sigma.square()).matrix().transpose() *
                    (((kappa - s) * (B) / (*tau)(i, 0))
                ).matrix()
            ).array();
    }

    return yield;
}

array YieldCurve::gradient(const array &p) {
    this->_updateParams(p);

    array h(1, 1);
    double r;
    if (this->currency == "base") {
        r = rm_init;
        h << 2 * rm;
    } else if (this->currency == "term") {
        r = rn_init;
        h << 2 * rn;
    }

    array s = this->s(h);
    array g = this->g(s);
    
    array yieldGradient(p.rows(), tau->rows());
    for (int i = 0; i <= (*tau).rows(); i++) {
        array A = this->A(s, g, i);
        array B = this->B(s, g, i);

    }
    return yieldGradient.transpose()
}

// internal state 
void YieldCurve::_updateParams(const array &p) {
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