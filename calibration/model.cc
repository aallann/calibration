#include <model.h>


Model::Model(
    double S,
    std::shared_ptr<array> vol,
    std::shared_ptr<array> tau,
    std::shared_ptr<array> rfb,
    std::shared_ptr<array> rft,
    std::shared_ptr<array> deltas,
    std::shared_ptr<array> omegas
) : S(S) {
    Model::vol = vol;
    Model::tau = tau;
    Model::rfb = rfb;
    Model::rft = rft;
    Model::deltas = deltas;
    Model::omegas = omegas;
};

void Model::setVol(std::shared_ptr<array> vol) {
    Model::vol = vol;
};

void Model::setTenors(std::shared_ptr<array> tau) {
    Model::tau = tau;
};

void Model::setBaseYieldRates(std::shared_ptr<array> rfb) {
    Model::rfb = rfb;
};

void Model::setTermYieldRates(std::shared_ptr<array> rft) {
    Model::rft = rft;
};

void Model::setDeltas(std::shared_ptr<array> deltas) {
    Model::deltas = deltas;
};

void Model::setOmegas(std::shared_ptr<array> omegas) {
    Model::omegas = omegas;
};