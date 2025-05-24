
template <class Real> VolumeVis<Real>::VolumeVis(const sctl::SlenderElemList<Real>& elem_lst, const sctl::Comm& comm, const bool shortened) : comm_(comm) {
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
  if (shortened) {
    for (sctl::Long elem_idx = 1; elem_idx < Nelem-1; elem_idx++) {
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
  } else {
    for (sctl::Long elem_idx = 0; elem_idx < Nelem; elem_idx++) {
      const Real t_order_inv = 1/(Real)t_order;
      const Real r_order_inv = (1-1e-6)/(Real)(r_order-1);
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
}

template <class Real> const sctl::Vector<Real>& VolumeVis<Real>::GetCoord() const {
  return coord;
}

template <class Real> void VolumeVis<Real>::WriteVTK(const std::string& fname, const sctl::Vector<Real>& F) const {
  sctl::VTUData vtu_data;
  GetVTUData(vtu_data, F);
  vtu_data.WriteVTK(fname, comm_);
}

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

template <class Real> StokesBIO<Real>::StokesBIO(const Real SL_scal, const Real DL_scal, const sctl::Comm comm)
  : comm_(comm), SL_scal_(SL_scal), DL_scal_(DL_scal), LayerPotenSL(ker_FxU, false, comm), LayerPotenDL(ker_DxU, false, comm) {
  LayerPotenSL.SetAccuracy(1e-14);
  LayerPotenDL.SetAccuracy(1e-14);
  LayerPotenSL.SetFMMKer(ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU, ker_FxU);
  LayerPotenDL.SetFMMKer(ker_DxU, ker_DxU, ker_DxU, ker_FSxU, ker_FSxU, ker_FSxU, ker_FxU, ker_FxU);
};

template <class Real> void StokesBIO<Real>::SetAccuracy(Real tol) {
  LayerPotenSL.SetAccuracy(tol);
  LayerPotenDL.SetAccuracy(tol);
}

template <class Real> template <class ElemLstType> void StokesBIO<Real>::AddElemList(const ElemLstType& elem_lst, const std::string& name) {
  LayerPotenSL.AddElemList(elem_lst, name);
  LayerPotenDL.AddElemList(elem_lst, name);
}

template <class Real> template <class ElemLstType> const ElemLstType& StokesBIO<Real>::GetElemList(const std::string& name) const {
  return LayerPotenSL.template GetElemList<ElemLstType>(name);
}

template <class Real> void StokesBIO<Real>::DeleteElemList(const std::string& name) {
  LayerPotenSL.DeleteElemList(name);
  LayerPotenDL.DeleteElemList(name);
}

template <class Real> template <class ElemLstType> void StokesBIO<Real>::DeleteElemList() {
  LayerPotenSL.template DeleteElemList<ElemLstType>();
  LayerPotenDL.template DeleteElemList<ElemLstType>();
}

template <class Real> void StokesBIO<Real>::SetTargetCoord(const sctl::Vector<Real>& Xtrg) {
  LayerPotenSL.SetTargetCoord(Xtrg);
  LayerPotenDL.SetTargetCoord(Xtrg);
}

template <class Real> void StokesBIO<Real>::SetTargetNormal(const sctl::Vector<Real>& Xn_trg) {
  LayerPotenSL.SetTargetNormal(Xn_trg);
  LayerPotenDL.SetTargetNormal(Xn_trg);
}

template <class Real> sctl::Long StokesBIO<Real>::Dim(sctl::Integer k) const {
  return LayerPotenSL.Dim(k);
}

template <class Real> void StokesBIO<Real>::Setup() const {
  if (SL_scal_) LayerPotenSL.Setup();
  if (DL_scal_) LayerPotenDL.Setup();
}

template <class Real> void StokesBIO<Real>::ClearSetup() const {
  LayerPotenSL.ClearSetup();
  LayerPotenDL.ClearSetup();
}

template <class Real> void StokesBIO<Real>::ComputePotential(sctl::Vector<Real>& U, const sctl::Vector<Real>& F) const {
  sctl::Vector<Real> Us, Ud;
  if (SL_scal_) LayerPotenSL.ComputePotential(Us, F);
  if (DL_scal_) LayerPotenDL.ComputePotential(Ud, F);

  if (SL_scal_ && DL_scal_) U = Us * SL_scal_ + Ud * DL_scal_;
  else if (SL_scal_) U = Us * SL_scal_;
  else if (DL_scal_) U = Ud * DL_scal_;
  else U.SetZero();
}

template <class Real> void StokesBIO<Real>::SqrtScaling(sctl::Vector<Real>& U) const {
  LayerPotenSL.SqrtScaling(U);
}

template <class Real> void StokesBIO<Real>::InvSqrtScaling(sctl::Vector<Real>& U) const {
  LayerPotenSL.InvSqrtScaling(U);
}

template <class Real> sctl::SlenderElemList<Real> PeriodicGeom<Real>::build_straight(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const Real r, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  for (sctl::Long k0 = -nbr_range; k0 <= nbr_range; k0++) { // 1D periodic in x direction.
    for (sctl::Long i = 0; i < Nelem; i++) {
      ElemOrderVec.PushBack(ElemOrder);
      FourierOrderVec.PushBack(FourierOrder);
      const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
      for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
        const Real x = (i+nodes[j])/Nelem;
        Xc.PushBack(k0+x);
        Xc.PushBack(0.5);
        Xc.PushBack(0.5);
        eps.PushBack(r);

        orient.PushBack(0);
        orient.PushBack(0);
        orient.PushBack(1);
      }
    }
    if (ptcls.Dim()>0) {
      sctl::Long Nptcl = ptcls.Dim();
      // many_sphs(ptcls_Xcs, ptcls_rs, Xc, r, Nptcl); // update ptcls_Xcs and ptcls_rs for ptcl info and filter_target later.
      if (nbr_range == 0) {
        many_sphs(ptcls_Xcs, ptcls_rs, Xc, r, Nptcl);
      }
      // TODO: if ptcl_Xcs and ptcl_rs NOT initialized, create them using many_sphs. otherwise, use what's there.
      for (sctl::Long p=0; p<Nptcl; p++) {
        const sctl::Long Nelem_sphere = ptcls[p];
        for (sctl::Long i = 0; i < Nelem_sphere; i++) { // add a sphere
          ElemOrderVec.PushBack(ElemOrder);
          FourierOrderVec.PushBack(FourierOrder);
          const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
          for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
            const Real ptcl_r = ptcls_rs[p];
            const Real theta = sctl::const_pi<Real>() * (i+nodes[j])/Nelem_sphere;
            Real x, y, z, ex, ey, ez, eps_j;
            // Real s = s_dsp + nds[k]*panel_len[dsp[i]+j];
            if (geom_mode == 1) {
              // sphere
              sphere_geom(x, y, z, ex, ey, ez, eps_j, theta, ptcl_r);
            } else if (geom_mode == 2) {
              // spheroid
              spheroid_geom(x, y, z, ex, ey, ez, eps_j, theta, ptcl_r);

            } else if (geom_mode == 3) {
              // bacteria
              bacteria_geom(x, y, z, ex, ey, ez, eps_j, 2*theta, ptcl_r); // 2pi for circular particles.
            } else if (geom_mode == 4) {
              // loop
              loop_geom(x, y, z, ex, ey, ez, eps_j, 2*theta, ptcl_r);
            } else {
              SCTL_ASSERT(false); // not implemented
            }
            Xc.PushBack(k0 + ptcls_Xcs[p*3]+x);
            Xc.PushBack(ptcls_Xcs[p*3+1]+y);
            Xc.PushBack(ptcls_Xcs[p*3+2]+z);
            eps.PushBack(eps_j);

            orient.PushBack(ex);
            orient.PushBack(ey);
            orient.PushBack(ez);
          }
        }
      }
    }
  }
  sctl::SlenderElemList<Real> elem_lst(ElemOrderVec, FourierOrderVec, Xc, eps, orient);
  return elem_lst;
};

template <class Real> sctl::SlenderElemList<Real> PeriodicGeom<Real>::build_sinusoidal(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const Real r, const Real mag, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  for (sctl::Long k0 = -nbr_range; k0 <= nbr_range; k0++) { // 1D periodic in x direction.
    for (sctl::Long i = 0; i < Nelem; i++) {
      ElemOrderVec.PushBack(ElemOrder);
      FourierOrderVec.PushBack(FourierOrder);
      const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
      for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
        const Real x = (i+nodes[j])/Nelem;
        Xc.PushBack(k0 + x);
        Xc.PushBack(mag * sctl::cos<Real>(2*sctl::const_pi<Real>()*x) + 0.5);
        Xc.PushBack(0.5);
        eps.PushBack(r);

        orient.PushBack(0);
        orient.PushBack(0);
        orient.PushBack(1);
      }
    }
    if (ptcls.Dim()>0) {
      sctl::Long Nptcl = ptcls.Dim();
      if (nbr_range == 0) {
        many_sphs(ptcls_Xcs, ptcls_rs, Xc, r, Nptcl);
      }
      for (sctl::Long p=0; p<Nptcl; p++) {
        const sctl::Long Nelem_sphere = ptcls[p];
        for (sctl::Long i = 0; i < Nelem_sphere; i++) { // add a sphere
          ElemOrderVec.PushBack(ElemOrder);
          FourierOrderVec.PushBack(FourierOrder);
          const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
          for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
            Real ptcl_r = ptcls_rs[p];
            const Real theta = sctl::const_pi<Real>() * (i+nodes[j])/Nelem_sphere;

            Real x, y, z, ex, ey, ez, eps_j;
            // Real s = s_dsp + nds[k]*panel_len[dsp[i]+j];
            if (geom_mode == 1) {
              // sphere
              sphere_geom(x, y, z, ex, ey, ez, eps_j, theta, ptcl_r);
            } else if (geom_mode == 2) {
              // spheroid
              spheroid_geom(x, y, z, ex, ey, ez, eps_j, theta, ptcl_r);

            } else if (geom_mode == 3) {
              // bacteria
              bacteria_geom(x, y, z, ex, ey, ez, eps_j, 2*theta, ptcl_r); // 2pi for circular particles.
            } else if (geom_mode == 4) {
              // loop
              loop_geom(x, y, z, ex, ey, ez, eps_j, 2*theta, ptcl_r);
            } else {
              SCTL_ASSERT(false); // not implemented
            }

            Xc.PushBack(k0 + ptcls_Xcs[p*3]+x);
            Xc.PushBack(ptcls_Xcs[p*3+1]+y);
            Xc.PushBack(ptcls_Xcs[p*3+2]+z);
            eps.PushBack(eps_j);

            orient.PushBack(ex);
            orient.PushBack(ey);
            orient.PushBack(ez);
          }
        }
      }
    }
  }
  sctl::SlenderElemList<Real> elem_lst(ElemOrderVec, FourierOrderVec, Xc, eps, orient);
  return elem_lst;
};

template <class Real> sctl::SlenderElemList<Real> PeriodicGeom<Real>::build_sinusoidal_mpi(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const Real r, const Real mag, const sctl::Comm& comm){
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  for (sctl::Long k0 = -nbr_range; k0 <= nbr_range; k0++) { // 1D periodic in x direction.
    for (sctl::Long i = 0; i < Nelem; i++) {
      ElemOrderVec.PushBack(ElemOrder);
      FourierOrderVec.PushBack(FourierOrder);
      const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
      for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
        const Real x = (i+nodes[j])/Nelem;
        Xc.PushBack(k0 + x);
        Xc.PushBack(mag * sctl::cos<Real>(2*sctl::const_pi<Real>()*x) + 0.5);
        Xc.PushBack(0.5);
        eps.PushBack(r);

        orient.PushBack(0);
        orient.PushBack(0);
        orient.PushBack(1);
      }
    }
  }
  // distributed memory
  const sctl::Long Ncpu = comm.Size();
  sctl::Vector<sctl::Long> obj_elem_cnt(Ncpu);
  sctl::Vector<sctl::Long> obj_elem_dsp(Ncpu);
  Real panel_len = 1/(Real)Nelem; // TODO: objects here have total length 1?
  obj_elem_cnt = ceil(Xc.Dim() / 3 / Ncpu); // TODO: ceiling function?
  obj_elem_cnt[Ncpu-1] = Xc.Dim() / 3 - (Ncpu-1) * obj_elem_cnt[0]; //  last cpu catches all remaining points.
  // Q: do nodes on single element have to stay on the same process? -- divide by panels not by nodes?
  obj_elem_dsp = 0.;
  sctl::omp_par::scan(obj_elem_cnt.begin(),obj_elem_dsp.begin(),Ncpu);
  // TODO: simpler to just copy it over?
  sctl::RigidBodyList<Real> rblist(comm, Ncpu, eps[0], "bacteria");
  sctl::SlenderElemList<Real> elem_lst;
  rblist.InitElemList(obj_elem_cnt, obj_elem_dsp, elem_lst, ElemOrder, FourierOrder, Xc, eps, orient, comm);
  // sctl::SlenderElemList<Real> elem_lst(ElemOrderVec, FourierOrderVec, Xc, eps, orient);
  return elem_lst;
};

template <class Real> std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> PeriodicGeom<Real>::filter_target(const sctl::Vector<Real> X, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real> ptcls_rs, const sctl::Vector<Real> ptcls_Xcs, const int geom_mode) {
  const sctl::Long N = X.Dim()/3; // number of targets
  const sctl::Long Nptcl = ptcls.Dim(); // number of particles
  sctl::Vector<sctl::Long> filtered_inds(N); // indicate whether the target was inside particle or not.
  filtered_inds.SetZero();
  sctl::Vector<Real> Xout; // collection of targets outside all ptcls.
  for (sctl::Long i = 0; i < N; i++) {
    const auto x = X.begin() + i*3;
    bool outside = true;
    for (sctl::Long j = 0; j < Nptcl; j++) {
      const auto pXc = ptcls_Xcs.begin() + j*3;
      const auto pr = ptcls_rs.begin() + j;
      outside = (outside && outside_ptcl(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0],geom_mode));
      if (!outside) {
        filtered_inds[i] = 1; // =1 if inside.
        break;  
      }
    }
    if (outside) {
      Xout.PushBack(x[0]);
      Xout.PushBack(x[1]);
      Xout.PushBack(x[2]);
    }
  }
  // std::cout << "Xout dim = " << Xout.Dim() << std::endl;
  return std::make_tuple(Xout,filtered_inds);
}

template <class Real> bool PeriodicGeom<Real>::outside_ptcl(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real pr, const int geom_mode) {
  if (geom_mode==1) {
    const Real d = (x1-pXc1)*(x1-pXc1) + (x2-pXc2)*(x2-pXc2) + (x3-pXc3)*(x3-pXc3);
    if (d < pr*pr*1.05) { // with buffer layer.
      return false;
    } else {
      return true;
    }
  } else if (geom_mode==2) {
    // TODO
    return true;
  } else if (geom_mode == 3) {
    // TODO
    return true;
  } else if (geom_mode == 4) {
    // TODO
    return true;
  } else {
    SCTL_ASSERT(false); // not implemented
  }
  
}

// could be private function?
template <class Real> void PeriodicGeom<Real>::many_sphs(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& Channel_Xc, const Real Channel_r, const sctl::Long Nobj) {
  srand48(2);

  ptcls_Xcs.ReInit(0);
  ptcls_rs.ReInit(0);
  sctl::Long obj_panels = Channel_Xc.Dim() / 3 / (Nobj+2); // floor of Npanels / Nobjs for approximately how many channel Xc nodes per obj.
  for (sctl::Long i = 1; i < Nobj+1; i++) {
    // base X,Y,Z on Xc of channel
    Real X0 = Channel_Xc[i*obj_panels*3];
    Real Y0 = Channel_Xc[i*obj_panels*3 + 1];
    Real Z0 = Channel_Xc[i*obj_panels*3 + 2];
    { // Set random offsets 
      // Real X = X0 + drand48()*Channel_r*0.03; // random offset that won't deviate too much away from center line.
      Real X = X0;
      Real Y = Y0 + drand48()*Channel_r*0.03;
      Real Z = Z0 + drand48()*Channel_r*0.03;
      // std::cout << "panel length /2 is " << (Channel_Xc[obj_panels*3] - Channel_Xc[0])/2. << std::endl;
      Real max_r = std::min(Channel_r - std::max({X-X0, Y-Y0, Z-Z0}), (Channel_Xc[obj_panels*3] - Channel_Xc[0])/2.);
      Real r = drand48() * max_r * 0.8; // random scaled down of max_r, at most 0.8 to avoid close-to-touching.
      // std::cout <<"r is " << r <<std::endl;

      ptcls_Xcs.PushBack(X);
      ptcls_Xcs.PushBack(Y);
      ptcls_Xcs.PushBack(Z);
      ptcls_rs.PushBack(r);
    }
  }
  // std::cout << "in make spheres; ptcls Xc dim = " << ptcls_Xcs.Dim() <<"; r dim = " << ptcls_rs.Dim() << std::endl;

}


template <class Real> void PeriodicGeom<Real>::bacteria_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad){
  Real t = theta/sctl::const_pi<Real>()-1; // -1:1
  Real aspect = sctl::const_pi<Real>()*3/2+1;

  Real L = aspect-1+sctl::const_pi<Real>()/2;
  Real scal = loop_rad/(1+L-sctl::const_pi<Real>()/2);
  // Real scal = loop_rad/(1+L-sctl::const_pi<Real>()/2) * 0.7;
  if (L*(1+t) < sctl::const_pi<Real>()/2) z = scal * (-sctl::cos<Real>(L*(1+t)) - L+sctl::const_pi<Real>()/2);
  else if (L*(1-t) < sctl::const_pi<Real>()/2) z = scal * (sctl::cos<Real>(L*(1-t)) + L-sctl::const_pi<Real>()/2);
  else z = scal * L * t;

  y = 0;
  x = 0;
  ex = 1/sctl::sqrt<Real>(3.);
  ey = 1/sctl::sqrt<Real>(3.);
  ez = 1/sctl::sqrt<Real>(3.);

  if (L*(1+t) < sctl::const_pi<Real>()/2) r = scal * sctl::sin<Real>(L*(1+t));
  else if (L*(1-t) < sctl::const_pi<Real>()/2) r = scal * sctl::sin<Real>(L*(1-t));
  else r = scal;
};

template <class Real> void PeriodicGeom<Real>::loop_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad){
  x = loop_rad * sctl::cos<Real>(theta);
  y = loop_rad * sctl::sin<Real>(theta);
  z = 0;
  ex = 0;
  ey = 0;
  ez = 1;
  r = 0.025;
};

template <class Real> void PeriodicGeom<Real>::sphere_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad){
  x = loop_rad*sctl::cos<Real>(theta);
  y = 0;
  z = 0;
  r = loop_rad*sctl::sin<Real>(theta);
  ex = 0;
  ey = 0;
  ez = 1;
}

template <class Real> void PeriodicGeom<Real>::spheroid_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad){
  Real A = loop_rad;
  Real C = A / 2.;

  // TODO
}