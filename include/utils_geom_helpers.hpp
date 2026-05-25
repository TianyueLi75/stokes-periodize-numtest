#ifndef _UTILS_GEOM_HELPERS_HPP_
#define _UTILS_GEOM_HELPERS_HPP_

#include <tuple>
#include <csbq.hpp>

template <class Real> class PeriodicGeom;

template <class Real>
inline sctl::Vector<Real> make_normal_orient(const sctl::Long channel_panels,
                                              const sctl::Long particle_panels,
                                              const sctl::Long ElemOrder,
                                              const sctl::Long FourierOrder,
                                              const Real channel_sign = 1.,
                                              const Real particle_sign = -1.) {
    sctl::Vector<Real> normal_orient;
    const sctl::Long channel_nodes = channel_panels * ElemOrder * FourierOrder;
    const sctl::Long particle_nodes = particle_panels * ElemOrder * FourierOrder;

    for (sctl::Long i = 0; i < channel_nodes * 3; i++) {
        normal_orient.PushBack(channel_sign);
    }
    for (sctl::Long i = 0; i < particle_nodes * 3; i++) {
        normal_orient.PushBack(particle_sign);
    }

    return normal_orient;
}

template <class Real>
inline std::tuple<sctl::SlenderElemList<Real>, sctl::Vector<Real>> build_elem_list(
    PeriodicGeom<Real>& geom,
    const sctl::Vector<sctl::Long>& ElemOrderVec,
    const sctl::Vector<sctl::Long>& FourierOrderVec,
    const sctl::Vector<Real>& Xc,
    const sctl::Vector<Real>& eps,
    const sctl::Vector<Real>& orient,
    const sctl::Vector<Real>& normal_orient,
    const bool use_orient = false) {
    sctl::SlenderElemList<Real> elem_lst;
    const auto normal_orient_ = geom.InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, normal_orient, use_orient);
    return std::make_tuple(elem_lst, normal_orient_);
}

#endif
