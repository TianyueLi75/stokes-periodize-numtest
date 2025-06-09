
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

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::build_straight(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real r, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
  comm_ = comm;
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  for (sctl::Long i = 0; i < Nelem; i++) {
    ElemOrderVec.PushBack(ElemOrder);
    FourierOrderVec.PushBack(FourierOrder);
    const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
    for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
      const Real x = (i+nodes[j])/Nelem;
      Xc.PushBack(x);
      Xc.PushBack(0.5);
      Xc.PushBack(0.5);
      eps.PushBack(r);

      orient.PushBack(0);
      orient.PushBack(0);
      orient.PushBack(1);
    }
  }
  sctl::Long Nelem_ptcl_tot = 0;
  if (ptcls.Dim()>0) {
    sctl::Long Nptcl = ptcls.Dim();
    if (nbr_range == 0) {
      many_sphs(ptcls_Xcs, ptcls_rs, Xc, 0, eps, Nptcl);
    }
    add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
      Nelem_ptcl_tot += ptcls[ptcl_i];
    }
  }
  sctl::Vector<Real> NormalOrient;
  { // set NormalOrient 
    for (sctl::Long i = 0; i < Nelem + Nelem_ptcl_tot; i++) {
      for (sctl::Long j = 0; j < ElemOrder*FourierOrder*COORD_DIM; j++) {
        NormalOrient.PushBack(i < Nelem ? 1 : -1);
      }
    }
  }
  sctl::SlenderElemList<Real> elem_lst;
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient, nbr_range, peri_mode);
  return std::make_tuple(elem_lst,NormalOrient_);
};

// not used...
template <class Real> sctl::SlenderElemList<Real> PeriodicGeom<Real>::build_sinusoidal_serial(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real r, const Real mag, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  if (peri_mode == 1) {
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
          many_sphs(ptcls_Xcs, ptcls_rs, Xc, (mag==0.1? 1 : 2), eps, Nptcl);
        }
        add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
      }
    }
  } else if (peri_mode == 3) {
    for (sctl::Long k2 = -nbr_range; k2<= nbr_range; k2++) {
      for (sctl::Long k1 = -nbr_range; k1<= nbr_range; k1++) {
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
              many_sphs(ptcls_Xcs, ptcls_rs, Xc, (mag==0.1? 1 : 2), eps, Nptcl);
            }
            add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
          }
        }
      }
    }
  }
  sctl::SlenderElemList<Real> elem_lst(ElemOrderVec, FourierOrderVec, Xc, eps, orient);
  // sctl::Vector<Real> NormalOrient;
  // { // set NormalOrient
  //   constexpr sctl::Integer COORD_DIM = 3;
  //   sctl::Vector<sctl::Long> elem_wise_node_cnt;
  //   elem_lst.GetNodeCoord(nullptr, nullptr, &elem_wise_node_cnt);
  //   for (sctl::Long i = 0; i < elem_wise_node_cnt.Dim(); i++) {
  //     for (sctl::Long j = 0; j < elem_wise_node_cnt[i]*COORD_DIM; j++) {
  //       NormalOrient.PushBack(i < Nelem ? 1 : -1);
  //     }
  //   }
  // }
  return elem_lst;
};

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::build_sinusoidal(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real r, const Real mag, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
  comm_ = comm;
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  for (sctl::Long i = 0; i < Nelem; i++) {
    ElemOrderVec.PushBack(ElemOrder);
    FourierOrderVec.PushBack(FourierOrder);
    const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
    for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
      const Real x = (i+nodes[j])/Nelem;
      Xc.PushBack(x);
      Xc.PushBack(mag * sctl::cos<Real>(2*sctl::const_pi<Real>()*x) + 0.5);
      Xc.PushBack(0.5);
      eps.PushBack(r);

      orient.PushBack(0);
      orient.PushBack(0);
      orient.PushBack(1);
    }
  }
  sctl::Long Nelem_ptcl_tot = 0;
  if (ptcls.Dim()>0) {
    sctl::Long Nptcl = ptcls.Dim();
    if (nbr_range == 0) {
      many_sphs(ptcls_Xcs, ptcls_rs, Xc, (mag==0.1? 1 : 2), eps, Nptcl);
    }
    add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
      Nelem_ptcl_tot += ptcls[ptcl_i];
    }
  }
  // std::cout << "size of Xc with ptcl: " << Xc.Dim() << std::endl;
  sctl::Vector<Real> NormalOrient;
  { // set NormalOrient 
    for (sctl::Long i = 0; i < Nelem + Nelem_ptcl_tot; i++) {
      for (sctl::Long j = 0; j < ElemOrder*FourierOrder*COORD_DIM; j++) {
        NormalOrient.PushBack(i < Nelem ? 1 : -1);
      }
    }
  }
  sctl::SlenderElemList<Real> elem_lst;
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient, nbr_range, peri_mode);
  return std::make_tuple(elem_lst,NormalOrient_);
};

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::build_conv_div(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real r1, const Real r2, const sctl::Comm& comm, sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
  comm_ = comm;
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;

  auto get_r = [&r1, &r2](const Real& x) {
    if (x < 0.1) {
      return 2*r1 + r2;
    } else if (x < 0.3) {
      return r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-0.1) / 0.2)+r1+r2;
    } else if (x < 0.7) {
      return r2;
    } else if (x < 0.9) {
      return r1+r2 - r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-0.7) / 0.2);
    } else {
      return 2*r1 + r2;
    }
  };

  for (sctl::Long i = 0; i < Nelem; i++) {
    ElemOrderVec.PushBack(ElemOrder);
    FourierOrderVec.PushBack(FourierOrder);
    const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
    for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
      const Real x = (i+nodes[j])/Nelem;
      Xc.PushBack(x);
      Xc.PushBack(0.5);
      Xc.PushBack(0.5);
      eps.PushBack(get_r(x));

      orient.PushBack(0);
      orient.PushBack(0);
      orient.PushBack(1);
    }
  }
  sctl::Long Nelem_ptcl_tot = 0;
  if (ptcls.Dim()>0) {
    // sctl::Long Nptcl = ptcls.Dim();
    if (nbr_range == 0) {
      // many_sphs(ptcls_Xcs, ptcls_rs, Xc, 3, eps, Nptcl);
      packed_sphs_conv_div(ptcls_Xcs, ptcls_rs, r1, r2);
    }
    // std::cout << ptcls_Xcs.Dim() << "; " << ptcls_rs.Dim() << std::endl;
    sctl::Vector<sctl::Long> ptcls_(ptcls_rs.Dim());
    ptcls_ = 1;
    ptcls.Swap(ptcls_);
    add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
      Nelem_ptcl_tot += ptcls[ptcl_i];
    }
  }
  // std::cout << "size of Xc with ptcl: " << Xc.Dim() << std::endl;
  sctl::Vector<Real> NormalOrient;
  { // set NormalOrient 
    for (sctl::Long i = 0; i < Nelem + Nelem_ptcl_tot; i++) {
      for (sctl::Long j = 0; j < ElemOrder*FourierOrder*COORD_DIM; j++) {
        NormalOrient.PushBack(i < Nelem ? 1 : -1);
      }
    }
  }
  sctl::SlenderElemList<Real> elem_lst;
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient, nbr_range, peri_mode);
  return std::make_tuple(elem_lst,NormalOrient_);
}

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::build_spiral(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real r_spiral, const Real r_channel, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
  comm_ = comm;
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;

  // sctl::Long n = 3;
  sctl::Long n = 1;

  // varying radius
  auto get_r = [&r_channel](const Real& x, const Real& angle) {
    if (x < 0.1 || x > 0.9) {
      return r_channel;
    } else {
      return r_channel + 0.3*r_channel*sctl::cos<Real>(angle);
    }
  };

  // // constant radius
  // auto get_r = [&r_channel](const Real& x, const Real& angle) {
  //   return r_channel;
  // };

  auto get_angle = [&n](const Real& x) {
    // no flat side:
    // return 2.*sctl::const_pi<Real>() * n * x;
    // flat at x<0.1 and x>0.9
    if (x<0.1) {
      return 0;
    } else if (x<0.9) {
      return 2.*sctl::const_pi<Real>() * n * (x-0.1) / 0.8; 
    } else {
      return 2.*sctl;:const_pi<Real>() * n;
    }
  }

  for (sctl::Long i = 0; i < Nelem; i++) {
    ElemOrderVec.PushBack(ElemOrder);
    FourierOrderVec.PushBack(FourierOrder);
    const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
    for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
      const Real x = (i+nodes[j])/Nelem;
      Xc.PushBack(x);
      // const Real angle = 2. * sctl::const_pi<Real>() * n * x; 
      const Real angle = get_angle(x);
      Xc.PushBack(0.5 + r_spiral * sctl::cos<Real>(angle));
      Xc.PushBack(0.5 + r_spiral * sctl::sin<Real>(angle));
      eps.PushBack(get_r(x, angle));

      orient.PushBack(0);
      orient.PushBack(0);
      orient.PushBack(1);
    }
  }
  sctl::Long Nelem_ptcl_tot = 0;
  if (ptcls.Dim()>0) {
    sctl::Long Nptcl = ptcls.Dim();
    if (nbr_range == 0) {
      many_sphs(ptcls_Xcs, ptcls_rs, Xc, 4, eps, Nptcl);
    }
    add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
      Nelem_ptcl_tot += ptcls[ptcl_i];
    }
  }
  // std::cout << "size of Xc with ptcl: " << Xc.Dim() << std::endl;
  sctl::Vector<Real> NormalOrient;
  { // set NormalOrient -- TODO: there should be a more direct way to create this for specific cpu without creating the whole list.
    for (sctl::Long i = 0; i < Nelem + Nelem_ptcl_tot; i++) {
      for (sctl::Long j = 0; j < ElemOrder*FourierOrder*COORD_DIM; j++) {
        NormalOrient.PushBack(i < Nelem ? 1 : -1);
      }
    }
  }
  sctl::SlenderElemList<Real> elem_lst;
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient, nbr_range, peri_mode);
  return std::make_tuple(elem_lst,NormalOrient_);
}

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::build_only_ptcls(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range, const sctl::Integer peri_mode, const Real box_sidelen, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
  comm_ = comm;
  sctl::Vector<Real> Xc_box, eps_box, orient_box;
  sctl::Vector<sctl::Long> ElemOrderVec_box, FourierOrderVec_box;
  for (sctl::Long i = 0; i < Nelem; i++) {
    ElemOrderVec_box.PushBack(ElemOrder);
    FourierOrderVec_box.PushBack(FourierOrder);
    const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec_box[i]);
    for (sctl::Long j = 0; j < ElemOrderVec_box[i]; j++) {
      const Real x = (i+nodes[j])/Nelem;
      Xc_box.PushBack(x);
      Xc_box.PushBack(0.5);
      Xc_box.PushBack(0.5);
      eps_box.PushBack(box_sidelen / 2.);

      orient_box.PushBack(0);
      orient_box.PushBack(0);
      orient_box.PushBack(1);
    }
  }
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  sctl::Long Nelem_ptcl_tot = 0;
  if (ptcls.Dim()>0) {
    sctl::Long Nptcl = ptcls.Dim();
    if (nbr_range == 0) {
      many_sphs(ptcls_Xcs, ptcls_rs, Xc_box, 0, eps_box, Nptcl);
    }
    add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
      Nelem_ptcl_tot += ptcls[ptcl_i];
    }
  }
  sctl::Vector<Real> NormalOrient;
  { // set NormalOrient 
    for (sctl::Long i = 0; i < Nelem_ptcl_tot; i++) {
      for (sctl::Long j = 0; j < ElemOrder*FourierOrder*COORD_DIM; j++) {
        NormalOrient.PushBack(-1);
      }
    }
  }
  sctl::SlenderElemList<Real> elem_lst;
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient, nbr_range, peri_mode);
  return std::make_tuple(elem_lst,NormalOrient_);
}

template <class Real> void PeriodicGeom<Real>::add_particles(sctl::Vector<sctl::Long>& ElemOrderVec, sctl::Vector<sctl::Long>& FourierOrderVec, sctl::Vector<Real>& Xc, sctl::Vector<Real>& eps, sctl::Vector<Real>& orient, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& ptcls_Xcs, const int geom_mode) {
  sctl::Long Nptcl = ptcls.Dim();
  // std::cout << "in add particles, Nptcl is " << Nptcl << std::endl;
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
        if (geom_mode == 0) {
          // sphere
          sphere_geom(x, y, z, ex, ey, ez, eps_j, theta, ptcl_r);
        } else if (geom_mode == 1) {
          // spheroid
          spheroid_geom(x, y, z, ex, ey, ez, eps_j, theta, ptcl_r);

        } else if (geom_mode == 2) {
          // bacteria
          bacteria_geom(x, y, z, ex, ey, ez, eps_j, 2*theta, ptcl_r); // 2pi for circular particles.
        } else if (geom_mode == 3) {
          // loop
          loop_geom(x, y, z, ex, ey, ez, eps_j, 2*theta, ptcl_r);
        } else {
          SCTL_ASSERT(false); // not implemented
        }

        Xc.PushBack(ptcls_Xcs[p*3]+x);
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


template <class Real> sctl::Vector<Real> PeriodicGeom<Real>::InitElemList(sctl::SlenderElemList<Real>& elem_lst, const sctl::Vector<sctl::Long>& ElemOrder, const sctl::Vector<sctl::Long>& FourierOrder, const sctl::Vector<Real>& X, const sctl::Vector<Real>& R, const sctl::Vector<Real>& OrientVec, const sctl::Vector<Real>& NormalOrient, const sctl::Integer nbr_range, const sctl::Integer peri_mode) {
  // peri_mode = j for j-periodic

  const sctl::Long Nelem = ElemOrder.Dim();
  // std::cout << "Nelem is " << Nelem << std::endl;
  sctl::Long loc_elem_cnt, loc_elem_dsp;
  if (Nelem) { // Set loc_elem_cnt, loc_elem_dsp 
    // node_dsp.ReInit(Nelem);
    // sctl::Vector<sctl::Long> node_cnt(Nelem);
    sctl::Vector<sctl::Long> node_cnt(Nelem), node_dsp(Nelem); node_dsp = 0;
    for (sctl::Long i = 0; i < Nelem; i++) {
      node_cnt[i] = ElemOrder[i] * FourierOrder[i] * FourierOrder[i];
    }
    sctl::omp_par::scan(node_cnt.begin(), node_dsp.begin(), Nelem);
    const sctl::Long Nnodes = node_cnt[Nelem-1] + node_dsp[Nelem-1];

    const sctl::Long Np = comm_.Size();
    const sctl::Long rank = comm_.Rank();
    sctl::Long a = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+0)/Np) - node_dsp.begin();
    sctl::Long b = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+1)/Np) - node_dsp.begin();
    if (rank == Np - 1) b = Nelem;
    if (rank == 0) a = 0;
    loc_elem_cnt = b-a;
    loc_elem_dsp = a;

    if (0 && !comm_.Rank()) { // Print partitioning
      std::cout<<"Partitioning: ";
      for (sctl::Long i = 0; i < comm_.Size(); i++) {
        sctl::Long a = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(i+0)/Np) - node_dsp.begin();
        sctl::Long b = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(i+1)/Np) - node_dsp.begin();
        if (i == Np - 1) b = Nelem;
        if (i == 0) a = 0;
        std::cout<<b-a<<' ';
      }
      std::cout<<'\n';

      std::cout<<"Weight: ";
      for (sctl::Long i = 0; i < comm_.Size(); i++) {
        sctl::Long a = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(i+0)/Np) - node_dsp.begin();
        sctl::Long b = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(i+1)/Np) - node_dsp.begin();
        if (i == Np - 1) b = Nelem;
        if (i == 0) a = 0;
        std::cout<<node_dsp[b-1]+node_cnt[b-1]-node_dsp[a]<<' ';
      }
      std::cout<<'\n';
    }
  } else {
    loc_elem_cnt = 0;
    loc_elem_dsp = 0;
  }

  // { // Set elem_lst
  const sctl::Vector<sctl::Long> LocElemOrder(loc_elem_cnt, (sctl::Iterator<sctl::Long>)ElemOrder.begin() + loc_elem_dsp, false);
  const sctl::Vector<sctl::Long> LocFourierOrder(loc_elem_cnt, (sctl::Iterator<sctl::Long>)FourierOrder.begin() + loc_elem_dsp, false);

  sctl::Long dsp = 0, cnt = 0;
  for (sctl::Long i = 0; i < loc_elem_dsp; i++) dsp += ElemOrder[i];
  for (sctl::Long i = 0; i < loc_elem_cnt; i++) cnt += ElemOrder[loc_elem_dsp+i];
  const sctl::Vector<Real> X_(cnt*COORD_DIM, (sctl::Iterator<Real>)X.begin() + dsp*COORD_DIM, false);
  const sctl::Vector<Real> R_(cnt, (sctl::Iterator<Real>)R.begin() + dsp, false);
  const sctl::Vector<Real> OrientVec_(cnt*COORD_DIM, (sctl::Iterator<Real>)OrientVec.begin() + dsp*COORD_DIM, false);
  const sctl::Vector<Real> NormalOrient_(cnt*FourierOrder[0]*COORD_DIM, (sctl::Iterator<Real>)NormalOrient.begin() + dsp*FourierOrder[0]*COORD_DIM, false); // TODO: this assumes all Fourier orders are the same.
  // std::cout << "cnt = " << cnt <<"; dsp = " << dsp << "; size NormalOrient_ before nbr: " << NormalOrient_.Dim() <<std::endl;
  // std::cout << "rank " << comm_.Rank() <<", size of Normal Orient all size " << NormalOrient.Dim() << ", looking at between " << dsp*FourierOrder[0]*COORD_DIM << " and " << dsp*FourierOrder[0]*COORD_DIM + cnt*FourierOrder[0]*COORD_DIM << std::endl;

  sctl::Vector<sctl::Long> LocElemOrder_nbr = vec_nbr_copy(LocElemOrder,nbr_range, peri_mode);
  sctl::Vector<sctl::Long> LocFourierOrder_nbr = vec_nbr_copy(LocFourierOrder,nbr_range, peri_mode);
  sctl::Vector<Real> X_nbr = X_nbr_copy(X_,nbr_range, peri_mode);
  sctl::Vector<Real> R_nbr = vec_nbr_copy(R_,nbr_range, peri_mode);
  sctl::Vector<Real> OrientVec_nbr = vec_nbr_copy(OrientVec_,nbr_range, peri_mode);
  // std::cout << "locElemord dim, fourierorder dim, X dim, R dim, Orient dim:  " << LocElemOrder_nbr.Dim() << ", " << LocFourierOrder_nbr.Dim() << ", " << X_nbr.Dim() <<", " << R_nbr.Dim() << ", " << OrientVec_nbr.Dim() << std::endl;
  elem_lst.template Init<Real>(LocElemOrder_nbr, LocFourierOrder_nbr, X_nbr, R_nbr, OrientVec_nbr);  
  // }
  return NormalOrient_;
}

template <class Real> sctl::Vector<Real> PeriodicGeom<Real>::vec_nbr_copy(const sctl::Vector<Real> X, const sctl::Integer nbr_range, const sctl::Integer peri_mode) {
  // sctl::Long Nrepeat = 2*nbr_range*peri_mode + 1;
  sctl::Long Nrepeat = 2*nbr_range + 1;
  if (peri_mode==3) {
    Nrepeat = Nrepeat * Nrepeat * Nrepeat;
  }
  sctl::Long N = X.Dim();
  sctl::Vector<Real> X_nbr(Nrepeat*N); // repeat X Nrepeat times
  for (sctl::Long k = 0; k < Nrepeat; k++) {
    for (sctl::Long i = 0; i < N; i++) {
      X_nbr[k*N+i] = X[i];
    }
  }
  return X_nbr;
}

template <class Real> sctl::Vector<sctl::Long> PeriodicGeom<Real>::vec_nbr_copy(const sctl::Vector<sctl::Long> X, const sctl::Integer nbr_range, const sctl::Integer peri_mode) {
  sctl::Long Nrepeat = 2*nbr_range + 1;
  if (peri_mode==3) {
    Nrepeat = Nrepeat * Nrepeat * Nrepeat;
  }
  sctl::Long N = X.Dim();
  sctl::Vector<sctl::Long> X_nbr(Nrepeat*N); // repeat X Nrepeat times
  for (sctl::Long k = 0; k < Nrepeat; k++) {
    for (sctl::Long i = 0; i < N; i++) {
      X_nbr[k*N+i] = X[i];
    }
  }
  return X_nbr;
}

template <class Real> sctl::Vector<Real> PeriodicGeom<Real>::X_nbr_copy(const sctl::Vector<Real> X, const sctl::Integer nbr_range, const sctl::Integer peri_mode) {
  if (nbr_range>0) { // duplicate geomtry to add images
    sctl::Vector<Real> Xc_;
    if (peri_mode == 1) {
      for (sctl::Long k0 = -nbr_range; k0 <= nbr_range; k0++) {
        for (sctl::Long i = 0; i < X.Dim()/3; i++) { // shift in x
          Xc_.PushBack(X[i*3+0] + k0);
          Xc_.PushBack(X[i*3+1]);
          Xc_.PushBack(X[i*3+2]);
        }
      }
    } else if (peri_mode == 2) {
      for (sctl::Long k1 = -nbr_range; k1 <= nbr_range; k1++) {
        for (sctl::Long k0 = -nbr_range; k0 <= nbr_range; k0++) {
          for (sctl::Long i = 0; i < X.Dim()/3; i++) { // shift in x
            Xc_.PushBack(X[i*3+0] + k0);
            Xc_.PushBack(X[i*3+1] + k1);
            Xc_.PushBack(X[i*3+2]);
          }
        }
      }
    } else if (peri_mode == 3) {
      for (sctl::Long k2 = -nbr_range; k2 <= nbr_range; k2++) {
        for (sctl::Long k1 = -nbr_range; k1 <= nbr_range; k1++) {
          for (sctl::Long k0 = -nbr_range; k0 <= nbr_range; k0++) {
            for (sctl::Long i = 0; i < X.Dim()/3; i++) { // shift in x
              Xc_.PushBack(X[i*3+0] + k0);
              Xc_.PushBack(X[i*3+1] + k1);
              Xc_.PushBack(X[i*3+2] + k2);
            }
          }
        }
      }
    } else {
      SCTL_ASSERT(false);
    }
    return Xc_;
  } else {
    return X;
  }
}

template <class Real> std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> PeriodicGeom<Real>::filter_target(const sctl::Vector<Real> X, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real> ptcls_rs, const sctl::Vector<Real> ptcls_Xcs, const int geom_mode) {
  const sctl::Long N = X.Dim()/3; // number of targets
  const sctl::Long Nptcl = ptcls.Dim(); // number of particles
  // std::cout << Nptcl <<std::endl;
  sctl::Vector<sctl::Long> filtered_inds(N); // indicate whether the target was inside particle or not.
  filtered_inds.SetZero();
  sctl::Vector<Real> Xout; // collection of targets outside all ptcls.
  for (sctl::Long i = 0; i < N; i++) {
    const auto x = X.begin() + i*3;
    bool outside = true;
    for (sctl::Long j = 0; j < Nptcl; j++) {
      const auto pXc = ptcls_Xcs.begin() + j*3;
      const auto pr = ptcls_rs.begin() + j;
      if (geom_mode==0) {
        outside = (outside && outside_sphere(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0]));
      } else if (geom_mode==1) {
        outside = (outside && outside_spheroid(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0],1.1));
      } else if (geom_mode == 2) {
        outside = (outside && outside_bacteria(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0]));
      } else if (geom_mode == 3) {
        outside = (outside && outside_loop(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0]));
      } else {
        SCTL_ASSERT(false); // not implemented
      }
      
      if (!outside) {
        // std::cout << "rank " << comm_.Rank() << ", inside particle" << std::endl;
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
  return std::make_tuple(Xout,filtered_inds);
}

template <class Real> bool PeriodicGeom<Real>::outside_sphere(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real pr) {
  const Real d = (x1-pXc1)*(x1-pXc1) + (x2-pXc2)*(x2-pXc2) + (x3-pXc3)*(x3-pXc3);
  if (d < pr*pr*1.05) { // with buffer layer.
    return false;
  } else {
    return true;
  }
}

template <class Real> bool PeriodicGeom<Real>::outside_spheroid(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real a, const Real u0) {
  const Real A = a * sctl::sqrt(u0*u0-1);
  const Real C = a * u0;
  const Real A2inv = 1. / (A*A); // TODO: numerical instability?
  const Real C2inv = 1. / (C*C); 
  const Real x1sq = (x1-pXc1)*(x1-pXc1);
  const Real x2sq = (x2-pXc2)*(x2-pXc2);
  const Real x3sq = (x3-pXc3)*(x3-pXc3);
  if (x2sq * C2inv + (x1sq+x3sq) * A2inv < 1.05) { // with buffer layer.
    return false;
  } else {
    return true;
  }
}

template <class Real> bool PeriodicGeom<Real>::outside_bacteria(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real pr) {
  //TODO
  return true;
}

template <class Real> bool PeriodicGeom<Real>::outside_loop(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real pr) {
  //TODO
  return true;
}


template <class Real> void PeriodicGeom<Real>::many_sphs(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& Channel_Xc, const sctl::Long Channel_mode, const sctl::Vector<Real>& Channel_eps, const sctl::Long Nobj) {
  srand48(2);

  ptcls_Xcs.ReInit(0);
  ptcls_rs.ReInit(0);
  sctl::Long N_segments, Nnodes_per_segment;
  sctl::Vector<sctl::Long> Nptcl_per_segment;
  if (Nobj >= 18) {
    N_segments = 20;
    Nnodes_per_segment = Channel_Xc.Dim() / 3 / N_segments;
    Nptcl_per_segment.ReInit(N_segments); // includes first and last.
    Nptcl_per_segment = Channel_Xc.Dim() / 3 / Nnodes_per_segment; // rounds down
    for (sctl::Long i = 0; i < (Channel_Xc.Dim()/3) % Nnodes_per_segment; i++) { // even out ptcls
      Nptcl_per_segment[i+1] += 1;
    }
  } else {
    N_segments = Nobj+2;
    Nptcl_per_segment.ReInit(N_segments);
    Nnodes_per_segment = Channel_Xc.Dim() / 3 / N_segments; // floor of Npanels / Nobjs for approximately how many channel Xc nodes per obj.
    Nptcl_per_segment = 1;
  }

  for (sctl::Long i = 1; i < N_segments-1; i++) {
    Real X0 = Channel_Xc[i*Nnodes_per_segment*3];
    Real Y0 = Channel_Xc[i*Nnodes_per_segment*3 + 1];
    Real Z0 = Channel_Xc[i*Nnodes_per_segment*3 + 2];
    Real eps0 = Channel_eps[i*Nnodes_per_segment];
    // Real X = X0 + drand48()*Channel_r*0.03; // random offset that won't deviate too much away from center line.
    Real X = X0;
    for (sctl::Long j=0; j < Nptcl_per_segment[i]; j++) {
      const Real angle = (2*sctl::const_pi<Real>() * j) / Nptcl_per_segment[i];
      Real Y = Y0 + sctl::cos<Real>(angle) * drand48() * eps0 * 0.7;
      Real Z = Z0 + sctl::sin<Real>(angle) * drand48() * eps0 * 0.7; 
      // Real max_r = std::min(Channel_r - std::max({X-X0, Y-Y0, Z-Z0}), (Channel_Xc[obj_panels*3] - Channel_Xc[0])/2.);
      // Real r = drand48() * max_r * 0.8; // random scaled down of max_r, at most 0.8 to avoid close-to-touching.
      Real r = eps0 * 0.15;
      ptcls_Xcs.PushBack(X);
      ptcls_Xcs.PushBack(Y);
      ptcls_Xcs.PushBack(Z);
      ptcls_rs.PushBack(r);
      if (ptcls_rs.Dim() == Nobj) {
        return;
      }
    }
  }
}

// template <class Real> void PeriodicGeom<Real>::many_sphs(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, const Real box_length, const sctl::Long geom_mode) {
//    // TODO: packed suspension of spheres to be repeated in 2- or 3-periodic geom.
//   return;
// }

template <class Real> void PeriodicGeom<Real>::packed_sphs_conv_div(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, const Real r1, const Real r2) {
  srand48(2);

  ptcls_Xcs.ReInit(0);
  ptcls_rs.ReInit(0);

  auto get_r = [&r1, &r2](const Real& x) {
    if (x < 0.1) {
      return 2*r1 + r2;
    } else if (x < 0.3) {
      return r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-0.1) / 0.2)+r1+r2;
    } else if (x < 0.7) {
      return r2;
    } else if (x < 0.9) {
      return r1+r2 - r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-0.7) / 0.2);
    } else {
      return 2*r1 + r2;
    }
  };

  const Real ptcl_r = 0.15 * r2; // radius/size of each particle.
  // hexagonal close packing
  const Real dx = 2. * (2. * ptcl_r); // buffer room between spheres = 2*diameter
  const Real dy = sctl::sqrt<Real>(3.) * (2. * ptcl_r);
  const Real dz = sctl::sqrt<Real>(6.) / 3. * 2. * (2. * ptcl_r);

  for (Real x=dx; x < 1.-dx; x+=dx) {
    Real X = x;
    Real r_channel = get_r(x);
    sctl::Long k = 0;
    for (Real z=-r_channel+(2. * ptcl_r); z<r_channel - (2. * ptcl_r); z += dz) {
      sctl::Long row = 0;
      Real y_bounds = sctl::sqrt<Real>(r_channel*r_channel - z*z);
      for (Real y=-y_bounds+(2. * ptcl_r); y<y_bounds - (2. * ptcl_r); y += dy) {
        Real y_shift = (row % 2)==0? 0 : (2. * ptcl_r);
        Real z_shift = (k % 2)==0? 0 : (2. * ptcl_r);
        Real Y = y + y_shift;
        Real Z = z + z_shift;

        sctl::Long iter_cnt = 0;
        Real ptcl_r_loc = drand48() * ptcl_r;
        while (iter_cnt < 30 && ptcl_r_loc >= 1e-5 && (Y*Y + Z*Z >= (r_channel-(2. * ptcl_r_loc))*(r_channel-(2. * ptcl_r_loc)))) {
          // reduce radius until particle fits inside channel.
          ptcl_r_loc *= 0.9;
          iter_cnt += 1;
        } 
        if (iter_cnt < 30 && ptcl_r_loc >= 1e-5) {
          ptcls_Xcs.PushBack(X);
          ptcls_Xcs.PushBack(Y + 0.5);
          ptcls_Xcs.PushBack(Z + 0.5);
          ptcls_rs.PushBack(ptcl_r_loc);
        }
        
      }
    } 

  }
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
  Real u0 = 1.1;
  x = 0;
  y = loop_rad * u0 * sctl::cos<Real>(theta);
  z = 0;
  r = loop_rad * sctl::sqrt<Real>(u0*u0 - 1) * sctl::sin<Real>(theta);
  ex = 0;
  ey = 1;
  ez = 0;
}