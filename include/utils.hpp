#ifndef _PERIODIZE_UTILS_HPP_
#define _PERIODIZE_UTILS_HPP_

#include <csbq.hpp>
#include <tuple>

/**
 * Visualize volume inside SlenderElemList.
 */
template <class Real> class VolumeVis {
    static constexpr sctl::Integer COORD_DIM = 3;
    static constexpr sctl::Integer s_order = 20;
    static constexpr sctl::Integer t_order = 60;
    static constexpr sctl::Integer r_order = 12;
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
  public: 
    /**
     * Create SlenderElem object for a straight channel.
     *
     * @param[in] Nelem number of panel quadrature on channel.
     * @param[in] ElemOrder number of Cheb nodes on each panel.
     * @param[in] FourierOder number of trapezoid nodes in azimuthal direction.
     * @param[in] nbr_range number of copies in x direction that are not covered by proxy periodization.
     * @param[in] r radius of straight channel.
     * @param[in] ptcls list of Nelem for each particle object to be initiated inside channel.
     * @param[in] ptcls_rs vector location to store radii of particles created by many_sphs().
     * @param[in] ptcls_Xcs vector location to store centerline locations for particles created by many_sphs().
     * @param[in] geom_mode =1 for spheres, =2 for spheroids, =3 for bacteria, =4 for loops
     */
    sctl::SlenderElemList<Real> build_straight(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const Real r, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

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
     * @param[in] geom_mode =1 for spheres, =2 for spheroids, =3 for bacteria, =4 for loops
     */
    sctl::SlenderElemList<Real> build_sinusoidal(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const Real r, const Real mag, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode);

    /**
     * Given a list of target X, and particle location information, filter out targets inside particles. 
     * Returns a vector of exterior targets, and a list of booleans indicating whether original target was inside a particle.
     *
     * @param[in] X vector containing original target points.
     * @param[in] ptcls list of Nelem for each particle object to be initiated inside channel.
     * @param[in] ptcls_rs vector of radii of particles.
     * @param[in] ptcls_Xcs vector of centerline locations for particles.
     * @param[in] geom_mode =1 for spheres, =2 for spheroids, =3 for bacteria, =4 for loops
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
    void many_sphs(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& Channel_Xc, const Real Channel_r, const sctl::Long Nobj);

    void sphere_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad);
    void spheroid_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad);
    void bacteria_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad);
    void loop_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad);

    bool outside_ptcl(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real pr, const int geom_mode);
};


#include <utils.cpp>

#endif // _PERIODIZE_UTILS_HPP_
