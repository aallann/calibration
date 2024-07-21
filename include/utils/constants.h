#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <utils/__types__.h>

namespace input {
    // volatility
    static const array _VOLATILITY_ = (array(165, 1) << 8.414,
                                          8.268,
                                          8.179,
                                          8.053,
                                          7.982,
                                          7.94,
                                          7.982,
                                          8.053,
                                          8.161,
                                          8.247,
                                          8.388,
                                          8.606,
                                          8.434,
                                          8.325,
                                          8.17,
                                          8.079,
                                          8.018,
                                          8.036,
                                          8.1,
                                          8.2,
                                          8.286,
                                          8.424,
                                          8.721,
                                          8.539,
                                          8.423,
                                          8.259,
                                          8.149,
                                          8.06,
                                          8.051,
                                          8.086,
                                          8.168,
                                          8.236,
                                          8.35,
                                          9.163,
                                          8.844,
                                          8.631,
                                          8.348,
                                          8.154,
                                          7.99,
                                          7.961,
                                          8.002,
                                          8.109,
                                          8.206,
                                          8.375,
                                          9.477,
                                          9.088,
                                          8.82,
                                          8.458,
                                          8.204,
                                          7.978,
                                          7.891,
                                          7.903,
                                          7.975,
                                          8.057,
                                          8.205,
                                          9.822,
                                          9.365,
                                          9.053,
                                          8.611,
                                          8.31,
                                          8.038,
                                          7.92,
                                          7.914,
                                          7.983,
                                          8.065,
                                          8.215,
                                          10.123,
                                          9.618,
                                          9.253,
                                          8.752,
                                          8.414,
                                          8.102,
                                          7.959,
                                          7.937,
                                          8.002,
                                          8.088,
                                          8.242,
                                          10.341,
                                          9.788,
                                          9.4,
                                          8.858,
                                          8.489,
                                          8.148,
                                          7.986,
                                          7.953,
                                          8.02,
                                          8.108,
                                          8.262,
                                          10.999,
                                          10.348,
                                          9.874,
                                          9.219,
                                          8.784,
                                          8.388,
                                          8.186,
                                          8.136,
                                          8.191,
                                          8.283,
                                          8.441,
                                          11.559,
                                          10.822,
                                          10.275,
                                          9.527,
                                          9.039,
                                          8.6,
                                          8.371,
                                          8.308,
                                          8.365,
                                          8.457,
                                          8.624,
                                          12.045,
                                          11.253,
                                          10.674,
                                          9.891,
                                          9.378,
                                          8.918,
                                          8.683,
                                          8.624,
                                          8.696,
                                          8.808,
                                          8.995,
                                          12.158,
                                          11.378,
                                          10.808,
                                          10.038,
                                          9.53,
                                          9.07,
                                          8.815,
                                          8.732,
                                          8.792,
                                          8.892,
                                          9.051,
                                          12.599,
                                          11.808,
                                          11.244,
                                          10.478,
                                          9.955,
                                          9.465,
                                          9.18,
                                          9.073,
                                          9.111,
                                          9.192,
                                          9.318,
                                          12.808,
                                          12.022,
                                          11.456,
                                          10.696,
                                          10.168,
                                          9.67,
                                          9.363,
                                          9.239,
                                          9.234,
                                          9.287,
                                          9.373,
                                          12.971,
                                          12.195,
                                          11.637,
                                          10.889,
                                          10.369,
                                          9.865,
                                          9.546,
                                          9.406,
                                          9.382,
                                          9.42,
                                          9.481)
                                          .finished() *
                                      0.01;

    // maturities
    static const array _TAU_ = (array(15, 1) << 14,
                                   21,
                                   30,
                                   60,
                                   90,
                                   120,
                                   150,
                                   180,
                                   270,
                                   360,
                                   540,
                                   720,
                                   1080,
                                   1440,
                                   1800)
                                   .finished() /
                               360;

    // deltas
    static const array _DELTAS_ = (array(1, 11) << -0.05,
        -0.1,
        -0.15,
        -0.25,
        -0.35,
        -0.5,
        0.35,
        0.25,
        0.15,
        0.1,
        0.05)
                                      .finished();

    // binary call/put : +1./-1. indicator
    static const array _OMEGA_ = (_DELTAS_ / _DELTAS_.abs());

    // term yield curve
    static const array _TERM_YIELD_ = (array(15, 1) << 0.023,
                                          0.023,
                                          0.023,
                                          0.023,
                                          0.023,
                                          0.022,
                                          0.021,
                                          0.020,
                                          0.029,
                                          0.038,
                                          0.082,
                                          0.126,
                                          0.315,
                                          0.550,
                                          0.779)
                                          .finished() /
                                      100;

    // term yield curve
    static const array _BASE_YIELD_ = (array(15, 1) << 0.053,
                                          0.053,
                                          0.053,
                                          0.068,
                                          0.081,
                                          0.090,
                                          0.099,
                                          0.109,
                                          0.125,
                                          0.140,
                                          0.197,
                                          0.247,
                                          0.363,
                                          0.560,
                                          0.749)
                                          .finished() /
                                      100;

    static double _SPOT_ = 1.3306;

};  // namespace input

namespace validation {
    static const array _BS_PRICES_ = (array(165, 1) << 0.0005,
        0.0006,
        0.0007,
        0.0011,
        0.0013,
        0.0016,
        0.0019,
        0.0021,
        0.0027,
        0.0034,
        0.0043,
        0.0051,
        0.0066,
        0.0079,
        0.0091,
        0.0010,
        0.0013,
        0.0016,
        0.0023,
        0.0029,
        0.0035,
        0.0040,
        0.0045,
        0.0059,
        0.0071,
        0.0092,
        0.0109,
        0.0141,
        0.0169,
        0.0195,
        0.0017,
        0.0021,
        0.0025,
        0.0037,
        0.0047,
        0.0055,
        0.0063,
        0.0071,
        0.0092,
        0.0111,
        0.0144,
        0.0170,
        0.0221,
        0.0265,
        0.0307,
        0.0032,
        0.0040,
        0.0048,
        0.0069,
        0.0086,
        0.0101,
        0.0115,
        0.0128,
        0.0165,
        0.0199,
        0.0256,
        0.0303,
        0.0397,
        0.0478,
        0.0557,
        0.0050,
        0.0062,
        0.0075,
        0.0106,
        0.0132,
        0.0154,
        0.0175,
        0.0195,
        0.0249,
        0.0298,
        0.0384,
        0.0456,
        0.0599,
        0.0723,
        0.0846,
        0.0084,
        0.0104,
        0.0125,
        0.0177,
        0.0217,
        0.0254,
        0.0287,
        0.0317,
        0.0404,
        0.0482,
        0.0622,
        0.0741,
        0.0974,
        0.1182,
        0.1389,
        0.0049,
        0.0060,
        0.0072,
        0.0100,
        0.0121,
        0.0140,
        0.0157,
        0.0172,
        0.0214,
        0.0251,
        0.0315,
        0.0366,
        0.0459,
        0.0536,
        0.0607,
        0.0031,
        0.0038,
        0.0046,
        0.0064,
        0.0077,
        0.0089,
        0.0099,
        0.0109,
        0.0135,
        0.0158,
        0.0199,
        0.0231,
        0.0290,
        0.0338,
        0.0382,
        0.0017,
        0.0020,
        0.0024,
        0.0034,
        0.0040,
        0.0047,
        0.0052,
        0.0057,
        0.0071,
        0.0083,
        0.0105,
        0.0122,
        0.0153,
        0.0177,
        0.0200,
        0.0010,
        0.0013,
        0.0015,
        0.0021,
        0.0025,
        0.0029,
        0.0032,
        0.0035,
        0.0044,
        0.0051,
        0.0065,
        0.0075,
        0.0094,
        0.0109,
        0.0123,
        0.0005,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0013,
        0.0015,
        0.0016,
        0.0020,
        0.0023,
        0.0029,
        0.0034,
        0.0042,
        0.0049,
        0.0055)
                                         .finished();

    static const array _BS_STRIKES_ = (array(165, 1) << 1.2950,
        1.2862,
        1.2771,
        1.2522,
        1.2324,
        1.2143,
        1.1980,
        1.1835,
        1.1436,
        1.1087,
        1.0573,
        1.0204,
        0.9540,
        0.9050,
        0.8647,
        1.3033,
        1.2966,
        1.2896,
        1.2714,
        1.2568,
        1.2436,
        1.2316,
        1.2211,
        1.1919,
        1.1663,
        1.1279,
        1.0997,
        1.0481,
        1.0099,
        0.9782,
        1.3087,
        1.3034,
        1.2979,
        1.2838,
        1.2726,
        1.2624,
        1.2534,
        1.2453,
        1.2231,
        1.2038,
        1.1741,
        1.1521,
        1.1111,
        1.0812,
        1.0564,
        1.3166,
        1.3133,
        1.3098,
        1.3012,
        1.2945,
        1.2886,
        1.2833,
        1.2786,
        1.2658,
        1.2548,
        1.2375,
        1.2247,
        1.2002,
        1.1833,
        1.1701,
        1.3227,
        1.3209,
        1.3190,
        1.3145,
        1.3110,
        1.3080,
        1.3054,
        1.3031,
        1.2969,
        1.2917,
        1.2838,
        1.2783,
        1.2675,
        1.2619,
        1.2590,
        1.3308,
        1.3309,
        1.3310,
        1.3314,
        1.3319,
        1.3323,
        1.3329,
        1.3334,
        1.3351,
        1.3369,
        1.3411,
        1.3454,
        1.3532,
        1.3636,
        1.3764,
        1.3389,
        1.3409,
        1.3430,
        1.3482,
        1.3522,
        1.3560,
        1.3594,
        1.3626,
        1.3719,
        1.3804,
        1.3964,
        1.4100,
        1.4327,
        1.4503,
        1.4640,
        1.3451,
        1.3486,
        1.3521,
        1.3611,
        1.3678,
        1.3740,
        1.3796,
        1.3848,
        1.3998,
        1.4136,
        1.4392,
        1.4604,
        1.4982,
        1.5289,
        1.5552,
        1.3532,
        1.3585,
        1.3639,
        1.3779,
        1.3880,
        1.3975,
        1.4061,
        1.4140,
        1.4367,
        1.4577,
        1.4968,
        1.5287,
        1.5874,
        1.6354,
        1.6791,
        1.3588,
        1.3655,
        1.3722,
        1.3899,
        1.4024,
        1.4143,
        1.4251,
        1.4350,
        1.4635,
        1.4897,
        1.5393,
        1.5792,
        1.6535,
        1.7143,
        1.7710,
        1.3675,
        1.3762,
        1.3849,
        1.4085,
        1.4249,
        1.4405,
        1.4548,
        1.4680,
        1.5057,
        1.5407,
        1.6073,
        1.6597,
        1.7591,
        1.8403,
        1.9180)
                                          .finished();

    static const array _HESTON_PRICES_ = (array(165, 1) << 0.0006,
        0.0007,
        0.0008,
        0.0011,
        0.0014,
        0.0015,
        0.0016,
        0.0018,
        0.0021,
        0.0023,
        0.0029,
        0.0036,
        0.0044,
        0.0051,
        0.0056,
        0.0013,
        0.0015,
        0.0018,
        0.0026,
        0.0031,
        0.0035,
        0.0039,
        0.0042,
        0.0050,
        0.0056,
        0.0069,
        0.0084,
        0.0101,
        0.0117,
        0.0130,
        0.0020,
        0.0024,
        0.0029,
        0.0041,
        0.0051,
        0.0058,
        0.0064,
        0.0070,
        0.0084,
        0.0096,
        0.0118,
        0.0141,
        0.0170,
        0.0197,
        0.0219,
        0.0037,
        0.0045,
        0.0053,
        0.0077,
        0.0095,
        0.0109,
        0.0122,
        0.0134,
        0.0163,
        0.0188,
        0.0232,
        0.0273,
        0.0331,
        0.0386,
        0.0434,
        0.0056,
        0.0068,
        0.0081,
        0.0117,
        0.0144,
        0.0167,
        0.0187,
        0.0206,
        0.0253,
        0.0294,
        0.0364,
        0.0427,
        0.0525,
        0.0616,
        0.0699,
        0.0090,
        0.0110,
        0.0132,
        0.0188,
        0.0232,
        0.0269,
        0.0302,
        0.0333,
        0.0412,
        0.0482,
        0.0604,
        0.0710,
        0.0889,
        0.1056,
        0.1216,
        0.0053,
        0.0065,
        0.0077,
        0.0109,
        0.0133,
        0.0152,
        0.0169,
        0.0184,
        0.0221,
        0.0251,
        0.0304,
        0.0355,
        0.0456,
        0.0561,
        0.0671,
        0.0034,
        0.0041,
        0.0049,
        0.0069,
        0.0084,
        0.0096,
        0.0107,
        0.0116,
        0.0138,
        0.0156,
        0.0187,
        0.0220,
        0.0284,
        0.0355,
        0.0428,
        0.0017,
        0.0021,
        0.0025,
        0.0035,
        0.0043,
        0.0049,
        0.0054,
        0.0059,
        0.0069,
        0.0078,
        0.0092,
        0.0110,
        0.0145,
        0.0186,
        0.0228,
        0.0010,
        0.0012,
        0.0015,
        0.0020,
        0.0025,
        0.0028,
        0.0031,
        0.0034,
        0.0040,
        0.0045,
        0.0053,
        0.0065,
        0.0088,
        0.0116,
        0.0144,
        0.0004,
        0.0005,
        0.0006,
        0.0008,
        0.0010,
        0.0011,
        0.0012,
        0.0014,
        0.0016,
        0.0018,
        0.0022,
        0.0028,
        0.0040,
        0.0055,
        0.0070)
                                             .finished();

    static const array _BASE_YIELD_CURVE_ = (array(165, 1) << 0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0003,
        0.0003,
        0.0004,
        0.0004,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0012,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0003,
        0.0003,
        0.0004,
        0.0004,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0012,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0003,
        0.0003,
        0.0004,
        0.0004,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0012,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0003,
        0.0003,
        0.0004,
        0.0004,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0012,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0003,
        0.0003,
        0.0004,
        0.0004,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0012,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0003,
        0.0003,
        0.0004,
        0.0004,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0012,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0003,
        0.0003,
        0.0004,
        0.0004,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0012,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0003,
        0.0003,
        0.0004,
        0.0004,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0012,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0003,
        0.0003,
        0.0004,
        0.0004,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0012,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0003,
        0.0003,
        0.0004,
        0.0004,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0012,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0002,
        0.0003,
        0.0003,
        0.0004,
        0.0004,
        0.0006,
        0.0007,
        0.0009,
        0.0011,
        0.0012)
                                                .finished();

    static const array _TERM_YIELD_CURVE_ = (array(165, 1) << 0.0005,
        0.0005,
        0.0006,
        0.0007,
        0.0008,
        0.0009,
        0.0010,
        0.0011,
        0.0014,
        0.0016,
        0.0021,
        0.0026,
        0.0034,
        0.0041,
        0.0047,
        0.0005,
        0.0005,
        0.0006,
        0.0007,
        0.0008,
        0.0009,
        0.0010,
        0.0011,
        0.0014,
        0.0016,
        0.0021,
        0.0026,
        0.0034,
        0.0041,
        0.0047,
        0.0005,
        0.0005,
        0.0006,
        0.0007,
        0.0008,
        0.0009,
        0.0010,
        0.0011,
        0.0014,
        0.0016,
        0.0021,
        0.0026,
        0.0034,
        0.0041,
        0.0047,
        0.0005,
        0.0005,
        0.0006,
        0.0007,
        0.0008,
        0.0009,
        0.0010,
        0.0011,
        0.0014,
        0.0016,
        0.0021,
        0.0026,
        0.0034,
        0.0041,
        0.0047,
        0.0005,
        0.0005,
        0.0006,
        0.0007,
        0.0008,
        0.0009,
        0.0010,
        0.0011,
        0.0014,
        0.0016,
        0.0021,
        0.0026,
        0.0034,
        0.0041,
        0.0047,
        0.0005,
        0.0005,
        0.0006,
        0.0007,
        0.0008,
        0.0009,
        0.0010,
        0.0011,
        0.0014,
        0.0016,
        0.0021,
        0.0026,
        0.0034,
        0.0041,
        0.0047,
        0.0005,
        0.0005,
        0.0006,
        0.0007,
        0.0008,
        0.0009,
        0.0010,
        0.0011,
        0.0014,
        0.0016,
        0.0021,
        0.0026,
        0.0034,
        0.0041,
        0.0047,
        0.0005,
        0.0005,
        0.0006,
        0.0007,
        0.0008,
        0.0009,
        0.0010,
        0.0011,
        0.0014,
        0.0016,
        0.0021,
        0.0026,
        0.0034,
        0.0041,
        0.0047,
        0.0005,
        0.0005,
        0.0006,
        0.0007,
        0.0008,
        0.0009,
        0.0010,
        0.0011,
        0.0014,
        0.0016,
        0.0021,
        0.0026,
        0.0034,
        0.0041,
        0.0047,
        0.0005,
        0.0005,
        0.0006,
        0.0007,
        0.0008,
        0.0009,
        0.0010,
        0.0011,
        0.0014,
        0.0016,
        0.0021,
        0.0026,
        0.0034,
        0.0041,
        0.0047,
        0.0005,
        0.0005,
        0.0006,
        0.0007,
        0.0008,
        0.0009,
        0.0010,
        0.0011,
        0.0014,
        0.0016,
        0.0021,
        0.0026,
        0.0034,
        0.0041,
        0.0047)
                                                .finished();

    static const array _SCALAR_HESTON_INPUT_P_ = (array(11, 1) << 0.254240,
        0.016825,
        0.095551,
        -0.428492,
        0.007233,
        0.000000,
        1.000000,
        0.257204,
        1.062127,
        -0.001714,
        -0.007229)
                                                     .finished();
};  // namespace validation

#endif  //,CONSTANTS_H;