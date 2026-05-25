#ifndef _UTILS_VIS_HPP_
#define _UTILS_VIS_HPP_

#include <csbq.hpp>
#include <tuple>

/**
 * Visualize volume inside SlenderElemList.
 */
template <class Real> class VolumeVis {
    static constexpr sctl::Integer COORD_DIM = 3;
    // sparse target points in channel for example_1_self_convergence.
    // static constexpr sctl::Integer s_order = 4;
    // static constexpr sctl::Integer t_order = 10;
    // static constexpr sctl::Integer r_order = 3; 
    // static constexpr sctl::Integer s_order = 20;
    // static constexpr sctl::Integer t_order = 60;
    // static constexpr sctl::Integer r_order = 12;
    static constexpr sctl::Integer s_order = 10;
    static constexpr sctl::Integer t_order = 30;
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
    VolumeVis(const sctl::SlenderElemList<Real>& elem_lst, const sctl::Comm& comm = sctl::Comm::Self());

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

template <class Real> class XsectionVis {
    static constexpr sctl::Integer COORD_DIM = 3;
    static constexpr sctl::Integer s_order = 4;
    static constexpr sctl::Integer t_order = 16;
    static constexpr sctl::Integer r_order = 5; 

  public:
    XsectionVis() = default;
    XsectionVis(const sctl::SlenderElemList<Real>& elem_lst, const sctl::Comm& comm = sctl::Comm::Self());
    const sctl::Vector<Real>& GetCoord() const;
    void SetCoord(const sctl::Vector<Real> new_coord);
    void WriteVTK(const std::string& fname, const sctl::Vector<Real>& F) const;
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


#include <utils_vis.cpp>

#endif 