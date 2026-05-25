#ifndef _UTILS_GEOM_PROFILES_HPP_
#define _UTILS_GEOM_PROFILES_HPP_

#include <tuple>
#include <csbq.hpp>

namespace geom_profiles {

template <class Real>
inline Real conv_div_radius(const Real x, const Real r1, const Real r2) {
    const Real cutoff1 = 0.25;
    const Real cutoff2 = 0.35;
    const Real cutoff3 = 0.65;
    const Real cutoff4 = 0.75;

    if (x < cutoff1) {
        return 2 * r1 + r2;
    } else if (x < cutoff2) {
        return r1 * sctl::cos<Real>(sctl::const_pi<Real>() * (x - cutoff1) / (cutoff2 - cutoff1)) + r1 + r2;
    } else if (x < cutoff3) {
        return r2;
    } else if (x < cutoff4) {
        return r1 + r2 - r1 * sctl::cos<Real>(sctl::const_pi<Real>() * (x - cutoff3) / (cutoff4 - cutoff3));
    }

    return 2 * r1 + r2;
}

template <class Real>
inline Real trefoil_radius(const Real x, const Real r_min, const Real r_max) {
    const Real angle = sctl::const_pi<Real>() * (16. * x - 28. / 3.);
    return r_min + (r_max - r_min) * (0.5 * sctl::sin<Real>(angle) + 0.5);
}

template <class Real>
inline std::tuple<Real, Real, Real> trefoil_xyz(const Real x) {
    const Real xminus = x - 0.5;
    const Real x4pi = 4. * sctl::const_pi<Real>() * xminus;
    const Real x8pi = 2. * x4pi;
    const Real xminus2 = xminus * xminus;
    const Real xminus5 = xminus2 * xminus2 * xminus;
    Real xcoeff = xminus2 * 4. - 1.;
    xcoeff = xcoeff / 5.;

    const Real x_ = 0.5 * xminus * sctl::cos<Real>(x4pi) + 8. * xminus5 + 0.5;
    const Real y_ = sctl::sin<Real>(x4pi) * xcoeff + 0.5;
    const Real z_ = sctl::sin<Real>(x8pi) * xcoeff + 0.5;

    return std::make_tuple(x_, y_, z_);
}

}  // namespace geom_profiles

#endif
