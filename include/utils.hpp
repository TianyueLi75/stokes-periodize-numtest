#ifndef _PERIODIZE_UTILS_HPP_
#define _PERIODIZE_UTILS_HPP_

#include <csbq.hpp>
#include <tuple>

/**
 * Visualize volume inside SlenderElemList.
 */
template <class Real> class VolumeVis {
    static constexpr sctl::Integer COORD_DIM = 3;
    // sparse target points in channel for example_1_self_convergence.
    static constexpr sctl::Integer s_order = 4;
    static constexpr sctl::Integer t_order = 10;
    static constexpr sctl::Integer r_order = 3; 
    // static constexpr sctl::Integer s_order = 20;
    // static constexpr sctl::Integer t_order = 60;
    // static constexpr sctl::Integer r_order = 12;
    // static constexpr sctl::Integer s_order = 10;
    // static constexpr sctl::Integer t_order = 30;
    // static constexpr sctl::Integer r_order = 12;
  public:

    VolumeVis() = default;

    /**
     * @brief Construct a new VolumeVis object.
     *
     * @param elem_lst the geometry.
     * @param comm MPI communicator.
     * @param shortened whether to restrict to inner target points.
     */
    VolumeVis(const sctl::SlenderElemList<Real>& elem_lst, const sctl::Comm& comm = sctl::Comm::Self(), const bool shortened = false);

    /**
     * @brief Get the coordinates of the discretization points.
     *
     * @return const Vector<Real>& Vector containing the coordinates.
     */
    const sctl::Vector<Real>& GetCoord() const;

    /**
     * @brief Write the volume to a VTK file.
     *
     * @param fname File name.
     * @param F Data associated with the discretization points.
     */
    void WriteVTK(const std::string& fname, const sctl::Vector<Real>& F) const;

    /**
     * @brief Get VTU data.
     *
     * @param vtu_data VTU data object.
     * @param F Data associated with the discretization points.
     */
    void GetVTUData(sctl::VTUData& vtu_data, const sctl::Vector<Real>& F) const;

  private:

    sctl::Comm comm_;
    sctl::Long Nelem;
    sctl::Vector<Real> coord;
};

/**
 * @brief Represents a uniformly discretized cube volume, shifted to have corner at origin.
 *
 * @tparam Real Data type for real numbers.
 */
template <class Real> class CubeVolumeVisShifted {
    static constexpr sctl::Integer COORD_DIM = 3;
  public:

    CubeVolumeVisShifted() = default;

    /**
     * @brief Construct a new CubeVolumeVis object.
     *
     * @param N_ Number of discretization points along one edge of the cube.
     * @param L Length of one edge of the cube.
     * @param comm MPI communicator.
     */
    CubeVolumeVisShifted(const sctl::Long N_, Real L, const sctl::Comm& comm = sctl::Comm::Self());

    /**
     * @brief Get the coordinates of the discretization points.
     *
     * @return const Vector<Real>& Vector containing the coordinates.
     */
    const sctl::Vector<Real>& GetCoord() const;

    /**
     * @brief Write the cube volume to a VTK file.
     *
     * @param fname File name.
     * @param F Data associated with the discretization points.
     */
    void WriteVTK(const std::string& fname, const sctl::Vector<Real>& F) const;

    /**
     * @brief Get VTU data.
     *
     * @param vtu_data VTU data object.
     * @param F Data associated with the discretization points.
     */
    void GetVTUData(sctl::VTUData& vtu_data, const sctl::Vector<Real>& F) const;

  private:

    sctl::Long N, N0;
    sctl::Comm comm;
    sctl::Vector<Real> coord;
};

/**
 * PVFMM cannot handle combined field kernel. Compute SL and DL separately and add them.
 */
template <class Real> class StokesBIO {
  public:

    StokesBIO() = delete;
    StokesBIO(const StokesBIO&) = delete;
    StokesBIO& operator= (const StokesBIO&) = delete;

    StokesBIO(const Real SL_scal, const Real DL_scal, const sctl::Comm comm);

    /**
     * Specify quadrature accuracy tolerance.
     *
     * @param[in] tol quadrature accuracy.
     */
    void SetAccuracy(Real tol);

    /**
     * Add an element-list.
     *
     * @param[in] elem_lst an object (of type ElemLstType, derived from the
     * base class ElementListBase) that contains the description of a list of
     * elements.
     *
     * @param[in] name a string name for this element list.
     */
    template <class ElemLstType> void AddElemList(const ElemLstType& elem_lst, const std::string& name = std::to_string(typeid(ElemLstType).hash_code()));

    /**
     * Get const reference to an element-list.
     *
     * @param[in] name name of the element-list to return.
     *
     * @return const reference to the element-list.
     */
    template <class ElemLstType> const ElemLstType& GetElemList(const std::string& name = std::to_string(typeid(ElemLstType).hash_code())) const;

    /**
     * Delete an element-list.
     *
     * @param[in] name name of the element-list to return.
     */
    void DeleteElemList(const std::string& name);

    /**
     * Delete an element-list.
     */
    template <class ElemLstType> void DeleteElemList();

    /**
     * Set target point coordinates.
     *
     * @param[in] Xtrg the coordinates of target points in array-of-struct
     * order: {x_1, y_1, z_1, x_2, ..., x_n, y_n, z_n}
     */
    void SetTargetCoord(const sctl::Vector<Real>& Xtrg);

    /**
     * Set target point normals.
     *
     * @param[in] Xn_trg the coordinates of target points in array-of-struct
     * order: {nx_1, ny_1, nz_1, nx_2, ..., nx_n, ny_n, nz_n}
     */
    void SetTargetNormal(const sctl::Vector<Real>& Xn_trg);

    /**
     * Get local dimension of the boundary integral operator. Dim(0) is the
     * input dimension and Dim(1) is the output dimension.
     */
    sctl::Long Dim(sctl::Integer k) const;

    /**
     * Setup the boundary integral operator.
     */
    void Setup() const;

    /**
     * Clear setup data.
     */
    void ClearSetup() const;

    /**
     * Evaluate the boundary integral operator.
     *
     * @param[out] U the potential computed at each target point in
     * array-of-struct order.
     *
     * @param[in] F the charge density at each surface discretization node in
     * array-of-struct order.
     */
    void ComputePotential(sctl::Vector<Real>& U, const sctl::Vector<Real>& F) const;

    /**
     * Scale input vector by sqrt of the area of the element.
     * TODO: replace by sqrt of surface quadrature weights (not sure if it makes a difference though)
     */
    void SqrtScaling(sctl::Vector<Real>& U) const;

    /**
     * Scale input vector by inv-sqrt of the area of the element.
     * TODO: replace by inv-sqrt of surface quadrature weights (not sure if it makes a difference though)
     */
    void InvSqrtScaling(sctl::Vector<Real>& U) const;


  private:

    const sctl::Stokes3D_FxU ker_FxU;
    const sctl::Stokes3D_DxU ker_DxU;
    const sctl::Stokes3D_FxUP ker_FxUP;
    const sctl::Stokes3D_FSxU ker_FSxU;

    const sctl::Comm comm_;
    const Real SL_scal_, DL_scal_;
    sctl::BoundaryIntegralOp<Real, sctl::Stokes3D_FxU> LayerPotenSL;
    sctl::BoundaryIntegralOp<Real, sctl::Stokes3D_DxU> LayerPotenDL;
};

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
     * @param[in] nbr_range number of copies in x direction that are not covered by proxy periodization.
     * @param[in] r radius of straight channel.
     * @param[in] ptcls list of Nelem for each particle object to be initiated inside channel.
     * @param[in] ptcls_rs vector location to store radii of particles created by many_sphs().
     * @param[in] ptcls_Xcs vector location to store centerline locations for particles created by many_sphs().
     * @param[in] geom_mode =0 for spheres, =1 for spheroids, =2 for bacteria, =3 for loops
     */
     std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_straight(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real r, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    /**
     * Create SlenderElem object for a sinusoidal channel.
     *
     * @param[in] Nelem number of panel quadrature on channel.
     * @param[in] ElemOrder number of Cheb nodes on each panel.
     * @param[in] FourierOder number of trapezoid nodes in azimuthal direction.
     * @param[in] nbr_range number of copies in x direction that are not covered by proxy periodization.
     * @param[in] r radius of sinusoidal channel.
     * @param[in] mag amplitude magnification to adjust oscillation; set to 0.1 or 0.3 in examples.
     * @param[in] ptcls list of Nelem for each particle object to be initiated inside channel.
     * @param[in] ptcls_rs vector location to store radii of particles created by many_sphs().
     * @param[in] ptcls_Xcs vector location to store centerline locations for particles created by many_sphs().
     * @param[in] geom_mode =0 for spheres, =1 for spheroids, =2 for bacteria, =3 for loops
     */
    sctl::SlenderElemList<Real> build_sinusoidal_serial(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real r, const Real mag, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    /**
     * @brief Sinusoidal pipe with neighbor support distributed among nodes, MPI enabled
     * 
     * @param Nelem 
     * @param ElemOrder 
     * @param FourierOrder 
     * @param nbr_range 
     * @param r 
     * @param mag 
     * @param comm 
     * @return sctl::SlenderElemList<Real> 
     */
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_sinusoidal(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real r, const Real mag, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);
    
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_conv_div(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real r1, const Real r2, const sctl::Comm& comm, sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const sctl::Long ptcl_ord, const int geom_mode);

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_spiral(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real r_spiral, const Real r_channel, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trefoil(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const sctl::Comm& comm, sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const sctl::Long ptcl_ord, const int geom_mode);
    
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> many_ptcls1(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> many_ptcls2(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const sctl::Comm& comm, const sctl::Long Nptcl, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    sctl::SlenderElemList<Real> free_ptcls(const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs);

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_only_ptcls(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real box_sidelen, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    sctl::Vector<Real> InitElemList(sctl::SlenderElemList<Real>& elem_lst, const sctl::Vector<sctl::Long>& ElemOrder, const sctl::Vector<sctl::Long>& FourierOrder, const sctl::Vector<Real>& X, const sctl::Vector<Real>& R, const sctl::Vector<Real>& OrientVec, const sctl::Vector<Real>& NormalOrient, const sctl::Integer nbr_range, const sctl::Integer peri_mode);

    std::tuple<sctl::Long,sctl::Long> GetGlobalIdx(const sctl::Vector<sctl::Long>& ElemOrder, const sctl::Vector<sctl::Long>& FourierOrder, const sctl::Comm& comm);

    sctl::Vector<Real> vec_nbr_copy(const sctl::Vector<Real> X, const sctl::Integer nbr_range, const sctl::Integer peri_mode);

    sctl::Vector<sctl::Long> vec_nbr_copy(const sctl::Vector<sctl::Long> X, const sctl::Integer nbr_range, const sctl::Integer peri_mode);

    sctl::Vector<Real> X_nbr_copy(const sctl::Vector<Real> X, const sctl::Integer nbr_range, const sctl::Integer peri_mode);

    sctl::Vector<Real> form_targets(const sctl::Long r_ord, const sctl::Long azi_ord, const sctl::SlenderElemList<Real>& elem_lst, const sctl::Comm& comm);

    /**
     * Given a list of target X, and particle location information, filter out targets inside particles. 
     * Returns a vector of exterior targets, and a list of booleans indicating whether original target was inside a particle.
     *
     * @param[in] X vector containing original target points.
     * @param[in] ptcls list of Nelem for each particle object to be initiated inside channel.
     * @param[in] ptcls_rs vector of radii of particles.
     * @param[in] ptcls_Xcs vector of centerline locations for particles.
     * @param[in] geom_mode =0 for spheres, =1 for spheroids, =2 for bacteria, =3 for loops
     */
    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> filter_target(const sctl::Vector<Real> X, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real> ptcls_rs, const sctl::Vector<Real> ptcls_Xcs, const int geom_mode);

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

    void packed_sphs_conv_div(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, const Real r1, const Real r2);

    void packed_sphs_trefoil(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, const Real r_min, const Real r_max);

    void sphere_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad);
    void spheroid_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad);
    void bacteria_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad);
    void loop_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad);

    bool outside_sphere(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real pr);
    bool outside_spheroid(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real a, const Real u0);
    bool outside_bacteria(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real pr);
    bool outside_loop(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real pr);
    // void GetDistribution(sctl::Long& elem_cnt_, sctl::Long& elem_dsp_, sctl::Vector<sctl::Long> node_dsp_);

  private:
    sctl::Comm comm_;
  //   sctl::Long loc_elem_cnt;
  //   sctl::Long loc_elem_dsp;
  //   sctl::Vector<sctl::Long> node_dsp;

  };


#include <utils.cpp>

#endif // _PERIODIZE_UTILS_HPP_
