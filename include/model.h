#ifndef MODEL_H
#define MODEL_H

#include <utils/__types__.h>

#include <memory>

class Model {
   public:
    Model(
        double S,
        std::shared_ptr<array> vol = nullptr,
        std::shared_ptr<array> tau = nullptr,
        std::shared_ptr<array> rfb = nullptr,
        std::shared_ptr<array> rft = nullptr,
        std::shared_ptr<array> deltas = nullptr,
        std::shared_ptr<array> omegas = nullptr
    );

    virtual ~Model() = default;

    static void setVol(std::shared_ptr<array> vol);

    static void setTau(std::shared_ptr<array> tau);

    static void setBaseYieldRates(std::shared_ptr<array> rfb);

    static void setTermYieldRates(std::shared_ptr<array> rft);

    static void setDeltas(std::shared_ptr<array> deltas);

    static void setOmegas(std::shared_ptr<array> omegas);

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