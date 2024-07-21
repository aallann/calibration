#include<model.h>


Model::Model(
        double S,
        std::shared_ptr<array> vol,
        std::shared_ptr<array> tau,
        std::shared_ptr<array> rfb,
        std::shared_ptr<array> rft,
        std::shared_ptr<array> deltas,
        std::shared_ptr<array> omegas
    ):
        S(S),
        vol(vol),
        tau(tau),
        rfb(rfb),
        rft(rft),
        deltas(deltas),
        omegas(omegas) {};

static void Model::setVol(std::shared_ptr<array> vol) {
    Model::vol = vol;
};

static void Model::setTau(std::shared_ptr<array> tau) {
    Model::tau = tau;
};

static void Model::setBaseYieldRates(std::shared_ptr<array> rfb) {
    Model::rfb = rfb;
};

static void Model::setTermYieldRates(std::shared_ptr<array> rft) {
    Model::rft = rft;
};

static void Model::setDeltas(std::shared_ptr<array> deltas) {
    Model::deltas = deltas;
};

static void Model::setOmegas(std::shared_ptr<array> omegas) {
    Model::omegas = omegas;
};