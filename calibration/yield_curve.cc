#include<yield_curve.h>
#include <iostream>

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
    return (1 - g * exp(-s * (*tau).row(i))) 
        / (1 - g);
}

array YieldCurve::B(
    array &s,
    array &g,
    int i
) const {
    return (1 - exp(-s * (*tau).row(i))) / 
        (1 - g * exp(-s * (*tau).row(i)));
}

// auxiliary partials
array YieldCurve::dg(
    int du,
    array &dv, 
    array &s
) const { 
    return ((du - dv) * (kappa + s) - (kappa - s) * (du + dv)) / (kappa + s).square(); 
}

array YieldCurve::dA(
    array &du, 
    array &dv, 
    array &s,
    array &g,
    int i
) const {
    return (*tau)(i, 0) * exp(-s * (*tau).row(i)) * g / 
        (1 - g) * du + (1 - exp(-s * (*tau).row(i))) / 
            (1 - g).square() * dv;
}

array YieldCurve::dB(
    array &du,
    array &dv, 
    array &s,
    array &g,
    int i
) const {
    return (*tau).row(i) * exp(-s * (*tau).row(i)) * (1 - g) / 
        (1 - g * exp(-s * (*tau)(i, 0))).square() * du + 
            (1 - exp(-s * (*tau)(i, 0))) * exp(-s * (*tau).row(i)) / 
        (1 - g * exp(-s * (*tau)(i, 0))).square() * dv;
}

// yield curve sensitivities to Heston params
void YieldCurve::Kappa(
   array &grad,
   array &s,
   array &g,
   array &A,
   array &B,
   int i
) const {
    array dsdk = kappa / s;
    array dgdk = this->dg(1, dsdk, s);
    array dAdk = this->dA(dsdk, dgdk, s, g, i);
    array dBdk = this->dB(dsdk, dgdk, s, g, i);
    
    grad.block(0, i, nDims, 1) = 
        -vbar / sigma.square() * (
            (kappa - s) + kappa * (1 - dsdk) - 2 * 
                (kappa * dAdk / A + log(A)) / (*tau).row(i)
        ) - v0 / sigma.square() * (
            (1 - dsdk) * B + (kappa - s) * dBdk
            ) / (*tau).row(i);
}

void YieldCurve::Vbar(
    array &grad,
    array &s,
    array &g,
    array &A,
    array &B,
    int i
) const {
    grad.block(nDims, i, nDims, 1) =
        -kappa / sigma.square() * ((kappa - s) - 2 * 
            log(A) / (*tau).row(i));
}

void YieldCurve::Sigma(
    array &grad,
    array &s,
    array &h, 
    array &g,
    array &A,
    array &B,
    int i
) const {
    array dsdsigma = sigma / s * h;
    array dgdsigma = this->dg(0, dsdsigma, s);
    array dAdsigma = this->dA(dsdsigma, dgdsigma, s, g, i);
    array dBdsigma = this->dB(dsdsigma, dgdsigma, s, g, i); 
    
    grad.block(2 * nDims, i, nDims, 1) =
        kappa * vbar * (*tau).row(i) * (dsdsigma * sigma.square() + 
            2 * (kappa - s) * sigma) / sigma.pow(4) /
                (*tau).row(i) +
            2 * kappa * vbar * (sigma.square() / 
                A * dAdsigma - 
            2 * sigma * log(A)) / sigma.pow(4) / 
                (*tau).row(i) + v0 * ((dsdsigma * B - (kappa - s) * 
            dBdsigma) * sigma.square() + 
                2 * (kappa - s) * B * sigma) /
            sigma.pow(4) / (*tau).row(i);
}

void YieldCurve::V0(
    array &grad,
    array &s,
    array &g,
    array &A,
    array &B,
    int i
) const {
    grad.block(4 * nDims, i, nDims, 1) = 
        -(kappa - s) * B / sigma.square() / (*tau).row(i);
}

void YieldCurve::R(
    array &grad,
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

    grad.block(7 * nDims, i, nDims, 1) = 
        kappa * vbar / sigma.square() * 
            (dsdr + 2 * dAdr / A / (*tau).row(i)) +
        v0 / sigma.square() * 
            (dsdr * B - (kappa - s) * dBdr / (*tau).row(i));
}


array YieldCurve::operator()(const array &p) {
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
    for (int i = 0; i <= tau->rows(); i++) {
        array A = this->A(s, g, i);
        array B = this->B(s, g, i);

        yield.row(i) = r - (
            (kappa * vbar / sigma.square()).matrix().transpose() 
                * (kappa - s - 2 * log((1 - g * exp(-s * (*tau).row(i))) / 
            (1 - g)) / (*tau).row(i)).matrix()
                ).array() -
                ((v0 / sigma.square()).matrix().transpose() *
             (((kappa - s) * (1 - exp(-s * (*tau).row(i))) / 
                (1 - g * exp(-s * (*tau).row(i)))) / (*tau).row(i)).matrix())
                .array();
    }

    return yield;
}

array YieldCurve::gradient(const array &p) {
    this->_updateParams(p);

    array h(1, 1);
    double r;
    if (currency == "base") {
        r = rm_init;
        h << 2 * rm;
    } else if (currency == "term") {
        r = rn_init;
        h << 2 * rn;
    }

    array s = this->s(h);
    array g = this->g(s);
    
    array yieldGradient(p.rows(), tau->rows());
    for (int i = 0; i <= tau->rows(); i++) {
        array A = this->A(s, g, i);
        array B = this->B(s, g, i);

        std::cout << "yes" << std::endl;

        Kappa(yieldGradient, s, g, A, B, i);
        Vbar(yieldGradient, s, g, A, B, i);
        Sigma(yieldGradient, s, h, g, A, B, i);
        V0(yieldGradient, s, g, A, B, i);
        R(yieldGradient, s, g, A, B, i);

        std::cout << "yes" << std::endl;

        yieldGradient.block(8 * nDims, i, nDims, 1) = array::Zero(nDims, 1);
        std::cout << "yes" << std::endl;
        yieldGradient.block(9 * nDims, i, 1, 1) = 1;
        yieldGradient.block(9 * nDims + 1, i, 1, 1) = 0;
        std::cout << "yes" << std::endl;
        if (this->currency == "term") {
            // if curremcy is not base, simply swap rows; this operation is efficient due to how the
            // instructions are handled by the compiler... under the hood, it is just a pointer change
            // with respect to array indices, rather than relocations or memory copy operations.
            yieldGradient.block(7 * nDims, i, nDims, 1).swap(yieldGradient.block(8 * nDims, i, nDims, 1));
            yieldGradient.block(9 * nDims, i, nDims, 1).swap(yieldGradient.block(9 * nDims + 1, i, nDims, 1));
        }
    } 
    
    return yieldGradient.transpose();
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
    rm_init = p(9 * nDims, 0);
    rn_init = p(p.rows() - 1, 0);
}