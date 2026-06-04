// =============================================================================
// utils_vis.cpp
//
// Template implementation of the visualization helpers (VolumeVis, XsectionVis,
// CubeVolumeVisShifted) declared in utils_vis.hpp. Not a standalone translation
// unit: it is included at the bottom of utils_vis.hpp.
// =============================================================================

// Build interior target points for each element by sweeping its surface samples
// radially inward toward the per-cross-section centroid.
template <class Real> VolumeVis<Real>::VolumeVis(const sctl::SlenderElemList<Real>& elem_lst, const sctl::Comm& comm) : comm_(comm) {
  Nelem = elem_lst.Size();
  sctl::Vector<Real> s_param, sin_theta, cos_theta;
  for (sctl::Long i = 0; i < s_order; i++) {
    const Real t = i/(Real)(s_order-1);
    s_param.PushBack(t);
  }
  for (sctl::Long i = 0; i < t_order; i++) {
    const Real t = i/(Real)t_order;
    sin_theta.PushBack(sctl::sin<Real>(2*sctl::const_pi<Real>()*t));
    cos_theta.PushBack(sctl::cos<Real>(2*sctl::const_pi<Real>()*t));
  }
  for (sctl::Long elem_idx = 0; elem_idx < Nelem; elem_idx++) {
    const Real t_order_inv = 1/(Real)t_order;
    const Real r_order_inv = (1-1e-3)/(Real)(r_order-1);
    sctl::Vector<Real> X_, Xc(COORD_DIM);
    elem_lst.GetGeom(&X_, nullptr, nullptr, nullptr, nullptr, s_param, sin_theta, cos_theta, elem_idx);
    for (sctl::Long i = 0; i < s_order; i++) {
      Xc = 0;
      for (sctl::Long j = 0; j < t_order; j++) {
        for (sctl::Long l = 0; l < COORD_DIM; l++) {
          Xc[l] += X_[(i*t_order+j)*COORD_DIM+l] * t_order_inv;
        }
      }
      for (sctl::Long j = 0; j < t_order; j++) {
        for (sctl::Long k = 0; k < r_order; k++) {
          for (sctl::Long l = 0; l < COORD_DIM; l++) {
            coord.PushBack((X_[(i*t_order+j)*COORD_DIM+l]-Xc[l])*k*r_order_inv + Xc[l]);
          }
        }
      }
    }
  }
}

template <class Real> const sctl::Vector<Real>& VolumeVis<Real>::GetCoord() const {
  return coord;
}

template <class Real> void VolumeVis<Real>::WriteVTK(const std::string& fname, const sctl::Vector<Real>& F) const {
  sctl::VTUData vtu_data;
  GetVTUData(vtu_data, F);
  vtu_data.WriteVTK(fname, comm_);
}

// Emit hexahedral VTK cells over the (s, theta, r) sample lattice of each element.
template <class Real> void VolumeVis<Real>::GetVTUData(sctl::VTUData& vtu_data, const sctl::Vector<Real>& F) const {
  for (const auto& x : coord) vtu_data.coord.PushBack((float)x);
  for (const auto& x :     F) vtu_data.value.PushBack((float)x);
  for (sctl::Long l = 0; l < Nelem; l++) {
    const sctl::Long offset = l * s_order*t_order*r_order;
    for (sctl::Long i = 0; i < s_order-1; i++) {
      for (sctl::Long j = 0; j < t_order; j++) {
        for (sctl::Long k = 0; k < r_order-1; k++) {
          auto idx = [this,&offset](sctl::Long i, sctl::Long j, sctl::Long k) {
            return offset+(i*t_order+(j%t_order))*r_order+k;
          };
          vtu_data.connect.PushBack(idx(i+0,j+0,k+0));
          vtu_data.connect.PushBack(idx(i+0,j+0,k+1));
          vtu_data.connect.PushBack(idx(i+0,j+1,k+1));
          vtu_data.connect.PushBack(idx(i+0,j+1,k+0));
          vtu_data.connect.PushBack(idx(i+1,j+0,k+0));
          vtu_data.connect.PushBack(idx(i+1,j+0,k+1));
          vtu_data.connect.PushBack(idx(i+1,j+1,k+1));
          vtu_data.connect.PushBack(idx(i+1,j+1,k+0));
          vtu_data.offset.PushBack(vtu_data.connect.Dim());;
          vtu_data.types.PushBack(12);
        }
      }
    }
  }
}

// Build an N^3 uniform grid on a cube of edge L centered at (0.5, 0.5, 0.5),
// partitioned across MPI ranks along the slowest index.
template <class Real> CubeVolumeVisShifted<Real>::CubeVolumeVisShifted(const sctl::Long N_, Real L, const sctl::Comm& comm_) : N(N_), comm(comm_) {
  const sctl::Long pid = comm.Rank();
  const sctl::Long Np = comm.Size();

  const sctl::Long NN = sctl::pow<COORD_DIM-1,sctl::Long>(N);
  const sctl::Long a = (N-1)*(pid+0)/Np;
  const sctl::Long b = (N-1)*(pid+1)/Np;
  N0 = b-a+1;
  if (N0<2) return;

  coord.ReInit(N0 * NN * COORD_DIM);
  for (sctl::Long i = 0; i < N0; i++) {
    for (sctl::Long j = 0; j < NN; j++) {
      for (sctl::Long k = 0; k < COORD_DIM; k++) {
        sctl::Long idx = ((i+a)*NN+j);
        coord[(i*NN+j)*COORD_DIM+k] = (((idx/sctl::pow<sctl::Long>(N,k)) % N)/(Real)(N-1) -0.5) * L + 0.5; // Assuming center at (0.5,0.5,0.5).
      }
    }
  }
}

template <class Real> const sctl::Vector<Real>& CubeVolumeVisShifted<Real>::GetCoord() const {
  return coord;
}

// Emit hexahedral VTK cells over this rank's slab of the uniform grid.
template <class Real> void CubeVolumeVisShifted<Real>::GetVTUData(sctl::VTUData& vtu_data, const sctl::Vector<Real>& F) const {
  for (const auto& x : coord) vtu_data.coord.PushBack((float)x);
  for (const auto& x :     F) vtu_data.value.PushBack((float)x);
  for (sctl::Long i = 0; i < N0-1; i++) {
    for (sctl::Long j = 0; j < N-1; j++) {
      for (sctl::Long k = 0; k < N-1; k++) {
        auto idx = [this](sctl::Long i, sctl::Long j, sctl::Long k) {
          return (i*N+j)*N+k;
        };
        vtu_data.connect.PushBack(idx(i+0,j+0,k+0));
        vtu_data.connect.PushBack(idx(i+0,j+0,k+1));
        vtu_data.connect.PushBack(idx(i+0,j+1,k+1));
        vtu_data.connect.PushBack(idx(i+0,j+1,k+0));
        vtu_data.connect.PushBack(idx(i+1,j+0,k+0));
        vtu_data.connect.PushBack(idx(i+1,j+0,k+1));
        vtu_data.connect.PushBack(idx(i+1,j+1,k+1));
        vtu_data.connect.PushBack(idx(i+1,j+1,k+0));
        vtu_data.offset.PushBack(vtu_data.connect.Dim());;
        vtu_data.types.PushBack(12);
      }
    }
  }
}
template <class Real> void CubeVolumeVisShifted<Real>::WriteVTK(const std::string& fname, const sctl::Vector<Real>& F) const {
  sctl::VTUData vtu_data;
  GetVTUData(vtu_data, F);
  vtu_data.WriteVTK(fname, comm);
}

// Like VolumeVis, but offsets the innermost samples slightly off the surface to
// avoid stagnation points in cross-section mixing visualizations.
template <class Real> XsectionVis<Real>::XsectionVis(const sctl::SlenderElemList<Real>& elem_lst, const sctl::Comm& comm) : comm_(comm) {
  Nelem = elem_lst.Size();
  sctl::Vector<Real> s_param, sin_theta, cos_theta;
  for (sctl::Long i = 0; i < s_order; i++) {
    const Real t = i/(Real)(s_order-1);
    s_param.PushBack(t);
  }
  for (sctl::Long i = 0; i < t_order; i++) {
    const Real t = i/(Real)t_order;
    sin_theta.PushBack(sctl::sin<Real>(2*sctl::const_pi<Real>()*t));
    cos_theta.PushBack(sctl::cos<Real>(2*sctl::const_pi<Real>()*t));
  }
  for (sctl::Long elem_idx = 0; elem_idx < Nelem; elem_idx++) {
    const Real t_order_inv = 1/(Real)t_order;
    const Real r_order_inv = (1-1e-2)/(Real)(r_order-1); // offset points off the surface to avoid stagnation points in mixing visualizations.
    sctl::Vector<Real> X_, Xc(COORD_DIM);
    elem_lst.GetGeom(&X_, nullptr, nullptr, nullptr, nullptr, s_param, sin_theta, cos_theta, elem_idx);
    for (sctl::Long i = 0; i < s_order; i++) {
      Xc = 0;
      for (sctl::Long j = 0; j < t_order; j++) {
        for (sctl::Long l = 0; l < COORD_DIM; l++) {
          Xc[l] += X_[(i*t_order+j)*COORD_DIM+l] * t_order_inv;
        }
      }
      for (sctl::Long j = 0; j < t_order; j++) {
        for (sctl::Long k = 0; k < r_order; k++) {
          for (sctl::Long l = 0; l < COORD_DIM; l++) {
            coord.PushBack((X_[(i*t_order+j)*COORD_DIM+l]-Xc[l])*k*r_order_inv + Xc[l]);
          }
        }
      }
    }
  }
}

template <class Real> const sctl::Vector<Real>& XsectionVis<Real>::GetCoord() const {
  return coord;
}

template <class Real> void XsectionVis<Real>::SetCoord(const sctl::Vector<Real> new_coord) {
  coord = new_coord;
}

template <class Real> void XsectionVis<Real>::WriteVTK(const std::string& fname, const sctl::Vector<Real>& F) const {
  sctl::VTUData vtu_data;
  GetVTUData(vtu_data, F);
  vtu_data.WriteVTK(fname, comm_);
}

// Emit hexahedral VTK cells over the (s, theta, r) sample lattice of each element.
template <class Real> void XsectionVis<Real>::GetVTUData(sctl::VTUData& vtu_data, const sctl::Vector<Real>& F) const {
  for (const auto& x : coord) vtu_data.coord.PushBack((float)x);
  for (const auto& x :     F) vtu_data.value.PushBack((float)x);
  for (sctl::Long l = 0; l < Nelem; l++) {
    const sctl::Long offset = l * s_order*t_order*r_order;
    for (sctl::Long i = 0; i < s_order-1; i++) {
      for (sctl::Long j = 0; j < t_order; j++) {
        for (sctl::Long k = 0; k < r_order-1; k++) {
          auto idx = [this,&offset](sctl::Long i, sctl::Long j, sctl::Long k) {
            return offset+(i*t_order+(j%t_order))*r_order+k;
          };
          vtu_data.connect.PushBack(idx(i+0,j+0,k+0));
          vtu_data.connect.PushBack(idx(i+0,j+0,k+1));
          vtu_data.connect.PushBack(idx(i+0,j+1,k+1));
          vtu_data.connect.PushBack(idx(i+0,j+1,k+0));
          vtu_data.connect.PushBack(idx(i+1,j+0,k+0));
          vtu_data.connect.PushBack(idx(i+1,j+0,k+1));
          vtu_data.connect.PushBack(idx(i+1,j+1,k+1));
          vtu_data.connect.PushBack(idx(i+1,j+1,k+0));
          vtu_data.offset.PushBack(vtu_data.connect.Dim());;
          vtu_data.types.PushBack(12);
        }
      }
    }
  }
}