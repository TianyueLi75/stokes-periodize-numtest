#ifndef _UTILS_GEOM_HPP_
#define _UTILS_GEOM_HPP_

#include <csbq.hpp>
#include <tuple>
#include "sctl/fmm-wrapper.hpp"  // for ParticleFMM

/**
 * Geometry set up for different test channels with or without particles inside.
 */
template <class Real> class PeriodicGeom {
  static constexpr sctl::Integer COORD_DIM = 3;
  
  public: 
    /**
     * Create SlenderElem object for a straight channel, MPI enabled
     *
     * @param[in] Nelem number of panel quadrature on channel.
     * @param[in] ElemOrder number of Cheb nodes on each panel.
     * @param[in] FourierOder number of trapezoid nodes in azimuthal direction.
     * @param[in] r radius of straight channel.
     * @param[in] ptcls list of Nelem for each particle object to be initiated inside channel.
     * @param[in] ptcls_rs vector location to store radii of particles created by many_sphs().
     * @param[in] ptcls_Xcs vector location to store centerline locations for particles created by many_sphs().
     * @param[in] geom_mode =0 for spheres, =1 for spheroids, =2 for loops
     */
     std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_straight(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const Real r, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    // r1 = 1/2 radius difference between larger and smaller radii; r2 = narrow radius. 2r1+r2 = larger radius.
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build_conv_div(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const Real r1, const Real r2, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_u0s, sctl::Vector<sctl::Long>& ptcls_ifprolate, const sctl::Long ptcl_ord, const sctl::Long N=2560);

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> many_ptcls1(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    // 3 spheres in unit cube
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> many_ptcls3(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> many_ptcls2(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, const sctl::Long Nptcl, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> many_loops2(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, const sctl::Long Nptcl, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_major_rs, sctl::Vector<Real>& ptcls_minor_rs,sctl::Vector<Real>& ptcls_thetas, sctl::Vector<Real>& ptcls_phis);

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> many_spheroids3(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs);

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> many_loops3(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs);

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> loops_system(const sctl::Vector<sctl::Long> ptcls, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Vector<Real>& ptcls_Xcs, const sctl::Vector<Real>& ptcls_major_rs, const sctl::Vector<Real>& ptcls_minor_rs, const sctl::Vector<Real>& ptcls_thetas, const sctl::Vector<Real>& ptcls_phis, sctl::Comm& comm);

    sctl::SlenderElemList<Real> spheroid_system(                                                           
                                                                const sctl::Long Nelem_ptcl, 
                                                                const sctl::Long ElemOrder, 
                                                                const sctl::Long FourierOrder, 
                                                                const sctl::Vector<Real> Xcenter_lst, 
                                                                const sctl::Vector<sctl::Long> if_prolate_lst, // true if prolate
                                                                const sctl::Vector<Real> u0_lst, 
                                                                const sctl::Vector<Real> r_lst, 
                                                                const sctl::Vector<Real> theta_lst, 
                                                                const sctl::Vector<Real> phi_lst, 
                                                                const sctl::Comm comm);

    void add_particles_rotated(sctl::Vector<sctl::Long>& ElemOrderVec, sctl::Vector<sctl::Long>& FourierOrderVec, sctl::Vector<Real>& Xc, sctl::Vector<Real>& eps, sctl::Vector<Real>& orient, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& ptcls_Xcs, const int geom_mode, const sctl::Vector<Real> ptcls_thetas, const sctl::Vector<Real> ptcls_phis);

    sctl::Vector<Real> InitElemList(sctl::SlenderElemList<Real>& elem_lst, const sctl::Vector<sctl::Long>& ElemOrder, const sctl::Vector<sctl::Long>& FourierOrder, const sctl::Vector<Real>& X, const sctl::Vector<Real>& R, const sctl::Vector<Real>& OrientVec, const sctl::Vector<Real>& NormalOrient, bool use_orient = false);

    std::tuple<sctl::Long,sctl::Long> GetGlobalIdx(const sctl::Vector<sctl::Long>& ElemOrder, const sctl::Vector<sctl::Long>& FourierOrder, const sctl::Comm& comm);

    sctl::Vector<Real> vec_nbr_copy(const sctl::Vector<Real> X, const sctl::Integer nbr_range, const sctl::Integer peri_mode);

    sctl::Vector<sctl::Long> vec_nbr_copy(const sctl::Vector<sctl::Long> X, const sctl::Integer nbr_range, const sctl::Integer peri_mode);

    sctl::Vector<Real> X_nbr_copy(const sctl::Vector<Real> X, const sctl::Integer nbr_range, const sctl::Integer peri_mode);

    /**
     * Shared public geometry predicates for filtering targets inside particles.
     * The implementation in utils_geom.cpp routes the rotated and non-rotated
     * paths through the same internal helper logic so the public API stays thin.
     *
     * Given a list of target X, and particle location information, filter out targets inside particles.
     * Returns a vector of exterior targets, and a list of booleans indicating whether original target was inside a particle.
     *
     * @param[in] X vector containing original target points.
     * @param[in] ptcls list of Nelem for each particle object to be initiated inside channel.
     * @param[in] ptcls_rs vector of radii of particles.
     * @param[in] ptcls_Xcs vector of particle centers.
     * @param[in] geom_mode =0 for spheres, =1 for spheroids, =2 for loops
     */
    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> filter_target(const sctl::Vector<Real> X, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real> ptcls_rs, const sctl::Vector<Real> ptcls_Xcs, const int geom_mode);

    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> filter_target_rotated(const sctl::Vector<Real> X, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real> ptcls_rs, const sctl::Vector<Real> ptcls_Xcs, const int geom_mode, const sctl::Vector<Real> ptcls_thetas, const sctl::Vector<Real> ptcls_phis);

    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> filter_spheroids(
                                                                                                const sctl::Vector<Real> X, 
                                                                                                const sctl::Vector<Real> r_all, 
                                                                                                const sctl::Vector<Real> u0_all, 
                                                                                                const sctl::Vector<Real> Xcenter_all, 
                                                                                                const sctl::Vector<sctl::Long> if_prolate_all) ;

    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> filter_spheroids_rotated(
                                                                                                const sctl::Vector<Real> X, 
                                                                                                const sctl::Vector<Real> r_all, 
                                                                                                const sctl::Vector<Real> u0_all, 
                                                                                                const sctl::Vector<Real> Xcenter_all, 
                                                                                                const sctl::Vector<sctl::Long> if_prolate_all,
                                                                                                const sctl::Vector<Real> theta_all, 
                                                                                                const sctl::Vector<Real> phi_all) ;

    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> filter_loops_rotated(
                                                                                                const sctl::Vector<Real> X, 
                                                                                                const sctl::Vector<Real> major_r_all, 
                                                                                                const sctl::Vector<Real> minor_r_all, 
                                                                                                const sctl::Vector<Real> Xcenter_all, 
                                                                                                const sctl::Vector<Real> theta_all, 
                                                                                                const sctl::Vector<Real> phi_all) ;

    /**
     * Given channel centerline location and radius, create a given number of spheres on the interior with some randomness.
     *
     * @param[in] ptcls_Xcs vector location to store centerline locations for particles created by many_sphs().
     * @param[in] ptcls_rs vector location to store radii of particles created by many_sphs().
     * @param[in] Channel_Xc vector of centerline locations for the channel.
     * @param[in] Channel_r radius of channel.
     * @param[in] Nobj number of particles to create.
     */
    void many_sphs(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& Channel_Xc, const sctl::Long Channel_mode, const sctl::Vector<Real>& Channel_eps, const sctl::Long Nobj);

    void add_particles(sctl::Vector<sctl::Long>& ElemOrderVec, sctl::Vector<sctl::Long>& FourierOrderVec, sctl::Vector<Real>& Xc, sctl::Vector<Real>& eps, sctl::Vector<Real>& orient, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    void add_spheroids(sctl::Vector<sctl::Long>& ElemOrderVec, sctl::Vector<sctl::Long>& FourierOrderVec, sctl::Vector<Real>& Xc, sctl::Vector<Real>& eps, sctl::Vector<Real>& orient, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& ptcls_Xcs, const sctl::Vector<Real>& ptcls_u0s, const sctl::Vector<sctl::Long>& ptcls_ifprolate);

    void add_spheroids_rotated(sctl::Vector<sctl::Long>& ElemOrderVec, sctl::Vector<sctl::Long>& FourierOrderVec, sctl::Vector<Real>& Xc, sctl::Vector<Real>& eps, sctl::Vector<Real>& orient, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& ptcls_Xcs, const sctl::Vector<Real>& ptcls_u0s, const sctl::Vector<Real>& ptcls_thetas, const sctl::Vector<Real>& ptcls_phis, const sctl::Vector<sctl::Long>& ptcls_ifprolate);

    void packed_spheroids_conv_div(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_u0s, sctl::Vector<Real>& ptcls_thetas, sctl::Vector<Real>& ptcls_phis, sctl::Vector<sctl::Long>& ptcls_ifprolate, const Real r1, const Real r2, const sctl::Long N = 2560);

    void sphere_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad);
    void spheroid_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad);
    void loop_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real major_r, const Real minor_r); // major_r gives size of overall ring, minor_r gives radius of cross section

    bool outside_sphere(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real pr);
    bool outside_spheroid(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real a, const Real u0, const int if_prolate);
    bool outside_loop(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real major_r, const Real minor_r);
    bool outside_spheroid_rotated(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real a, const Real u0, const int if_prolate, const Real ptheta, const Real pphi);
    bool outside_loop_rotated(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real major_r, const Real minor_r, const Real ptheta, const Real pphi);

    sctl::Vector<Real> exact_field_fmm(const sctl::Vector<Real>& Xtrg, const sctl::Vector<Real>& Xsrc, const sctl::Vector<Real>& sigma, const sctl::Long Ncopy, const sctl::Integer peri_mode);

  private:
    sctl::Comm comm_;
  //   sctl::Long loc_elem_cnt;
  //   sctl::Long loc_elem_dsp;
  //   sctl::Vector<sctl::Long> node_dsp;
    mutable sctl::ParticleFMM<Real,COORD_DIM> fmm;

    sctl::Vector<Real> make_normal_orient(const sctl::Long channel_panels,
                                                const sctl::Long particle_panels,
                                                const sctl::Long ElemOrder,
                                                const sctl::Long FourierOrder,
                                                const Real channel_sign = 1.,
                                                const Real particle_sign = -1.);

    std::tuple<sctl::SlenderElemList<Real>, sctl::Vector<Real>> build_elem_list(
                                                                  PeriodicGeom<Real>& geom,
                                                                  const sctl::Vector<sctl::Long>& ElemOrderVec,
                                                                  const sctl::Vector<sctl::Long>& FourierOrderVec,
                                                                  const sctl::Vector<Real>& Xc,
                                                                  const sctl::Vector<Real>& eps,
                                                                  const sctl::Vector<Real>& orient,
                                                                  const sctl::Vector<Real>& normal_orient,
                                                                  const bool use_orient = false);

    Real conv_div_radius(const Real x, const Real r1, const Real r2);

};


#include <utils_geom.cpp>

#endif 