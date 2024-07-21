#ifndef MODEL_H
#define MODEL_H

#include <utils/__all__.h>

class Model {
   public:
    Model(
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

    virtual ~Model() = default;

    inline static void setVol(std::shared_ptr<array> vol) {
        Model::vol = vol;
    };

    inline static void setTau(std::shared_ptr<array> tau) {
        Model::tau = tau;
    };

    inline static void setBaseYieldRates(std::shared_ptr<array> rfb) {
        Model::rfb = rfb;
    };

    inline static void setTermYieldRates(std::shared_ptr<array> rft) {
        Model::rft = rft;
    };

    inline static void setDeltas(std::shared_ptr<array> deltas) {
        Model::deltas = deltas;
    };

    inline static void setOmegas(std::shared_ptr<array> omegas) {
        Model::omegas = omegas;
    };

   protected:
    double S;
    inline static std::shared_ptr<array> vol = nullptr;
    inline static std::shared_ptr<array> tau = nullptr;
    inline static std::shared_ptr<array> rfb = nullptr;
    inline static std::shared_ptr<array> rft = nullptr;
    inline static std::shared_ptr<array> deltas = nullptr;
    inline static std::shared_ptr<array> omegas = nullptr;
};

#endif  // MODEL_H