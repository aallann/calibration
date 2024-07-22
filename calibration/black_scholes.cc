#include <black_scholes.h>

BlackScholes::BlackScholes(double S) : Model(S) {};

array BlackScholes::d1(std::vector<EuropeanOption> &options) const {
    return (log(S / (*strikes)) + 
        ((*rft) - (*rfb) + 0.5 * (*vol).square()) * (*tau)) /
           ((*vol) * sqrt(*tau));
}

array BlackScholes::d2(std::vector<EuropeanOption> &options) const {
    return d1(options) - (*vol) * sqrt(*tau);
}

std::shared_ptr<array> BlackScholes::strikeFromDelta() {
    uint n = vol->rows();
    uint m = vol->cols();

    array fwdDeltas(n, m) = exp((*rfb) * (*tau)) * (*deltas);

    array normInvFwdDeltas(n, m);
    for (uint i = 0; i < n; i++) {
        for (uint j = 0; j < m; j++) {
            normInvFwdDeltas(i, j) =
                standardGaussian.quantile((*omegas)(i, j) * fwdDeltas(i, j));
        }
    }

    array _strikes = S * exp((*rfb) +
        (-(*omegas) * normInvFwdDeltas) * (*vol) * sqrt(*tau) +
        ((*rft) - (*rfb) + 0.5 * (*vol).square()) * (*tau));

    strikes = std::make_shared<array>(_strikes);

    setStrikes(strikes);

    return strikes;
}

array BlackScholes::getDeltaFromStrike() {
    //TODO
    return array();
}

array BlackScholes::price(std::vector<EuropeanOption> &options) {
    array d1 = this->d1(options);
    array d2 = this->d2(options);

    array prices((*vol).rows(), (*vol).cols());
    for (const auto &option : options) {
        auto [i, j] = option.idx;
        double omega = (*omegas)(i, j);

        prices(i, j) = omega * (
            (S * std::exp(-(*rfb)(i, j) * (*tau)(i, j))
                * standardGaussian.cdf(omega * d1(i, j)))
            - ((*option.strikes)(i, j) * std::exp(-(*rft)(i, j) * (*tau)(i, j))
                * standardGaussian.cdf(omega * d2(i, j)))
        );
    }

    return prices;
}

array BlackScholes::impVol(std::vector<EuropeanOption> &options) {
        //TODO
    return array();
}