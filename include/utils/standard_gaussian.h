#ifndef STANDARD_GAUSSIAN_H
#define STANDARD_GAUSSIAN_H

#include <boost/math/distributions.hpp>

class StandardGaussian {
   public:
    StandardGaussian() : dist(0.0, 1.0) {};

    inline double quantile(const double &p) const {
        return boost::math::quantile(dist, p);
    };

    inline double cdf(const double &x) const {
        return boost::math::cdf(dist, x);
    };

    inline double pdf(const double &x) const {
        return boost::math::pdf(dist, x);
    };

   private:
    boost::math::normal_distribution<double> dist;
};

#endif  // STANDARD_GAUSSIAN_H