template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::build_straight(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const Real r, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
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
    if (Nptcl == 1) {
      // Special test case of one particle on centerline
      ptcls_Xcs.ReInit(3);
      ptcls_rs.ReInit(1);
      ptcls_Xcs = 0.5;
      ptcls_rs = r / 2.;
    } else {
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
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient);
  return std::make_tuple(elem_lst,NormalOrient_);
};

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::build_sinusoidal(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const Real r, const Real mag, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
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
    many_sphs(ptcls_Xcs, ptcls_rs, Xc, (mag==0.1? 1 : 2), eps, Nptcl);
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
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient);
  return std::make_tuple(elem_lst,NormalOrient_);
};

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::build_conv_div(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const Real r1, const Real r2, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_u0s, sctl::Vector<sctl::Long>& ptcls_ifprolate, const sctl::Long ptcl_ord, const sctl::Long N){
  comm_ = comm;
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;

  auto get_r = [&r1, &r2](const Real& x) {
    Real cutoff1 = 0.25; // end of large radius left
    Real cutoff2 = 0.35; // start of small radius left
    Real cutoff3 = 0.65; // end of small radius right
    Real cutoff4 = 0.75; // start of large radius right
    if (x < cutoff1) {
      return 2*r1 + r2;
    } else if (x < cutoff2) {
      return r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-cutoff1) / (cutoff2-cutoff1))+r1+r2;
    } else if (x < cutoff3) {
      return r2;
    } else if (x < cutoff4) {
      return r1+r2 - r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-cutoff3) / (cutoff4-cutoff3));
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
      orient.PushBack(1);
      orient.PushBack(0);
    }
  }
  sctl::Long Nelem_ptcl_tot = 0;
  sctl::Vector<Real> ptcls_thetas;
  sctl::Vector<Real> ptcls_phis;
  if (ptcls.Dim()>0) {
    // packed_sphs_conv_div(ptcls_Xcs, ptcls_rs, r1, r2); // Fit grid of small spheres inside the channel.
    packed_spheroids_conv_div(ptcls_Xcs, ptcls_rs, ptcls_u0s, ptcls_thetas, ptcls_phis, ptcls_ifprolate, r1, r2, N);
    // packed_spheres_conv_div(ptcls_Xcs, ptcls_rs, ptcls_thetas, ptcls_phis, r1, r2); // DEBUGGING VSLIP: use all spheres and no rotation

    sctl::Vector<sctl::Long> ptcls_(ptcls_rs.Dim());
    ptcls_ = ptcl_ord;
    ptcls.Swap(ptcls_);
    
    // // DEBUGGING: remove rotation
    // ptcls_thetas.ReInit(ptcls_rs.Dim());
    // ptcls_phis.ReInit(ptcls_rs.Dim());
    // ptcls_thetas.SetZero();
    // ptcls_phis.SetZero();
    // ///////////////////////////
    
    add_spheroids_rotated(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, ptcls_u0s, ptcls_thetas, ptcls_phis, ptcls_ifprolate);
    // add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, 0); // CHOCO FEB 2026: DEBUG VSLIP with spheres
    // add_particles_rotated(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, 0, ptcls_thetas, ptcls_phis); 
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
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient, true); 
  return std::make_tuple(elem_lst,NormalOrient_,ptcls_thetas,ptcls_phis);
}

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::build_conv_div_sph(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const Real r1, const Real r2, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const sctl::Long ptcl_ord){
  comm_ = comm;
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;

  auto get_r = [&r1, &r2](const Real& x) {
    Real cutoff1 = 0.25; // end of large radius left
    Real cutoff2 = 0.35; // start of small radius left
    Real cutoff3 = 0.65; // end of small radius right
    Real cutoff4 = 0.75; // start of large radius right
    if (x < cutoff1) {
      return 2*r1 + r2;
    } else if (x < cutoff2) {
      return r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-cutoff1) / (cutoff2-cutoff1))+r1+r2;
    } else if (x < cutoff3) {
      return r2;
    } else if (x < cutoff4) {
      return r1+r2 - r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-cutoff3) / (cutoff4-cutoff3));
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
      orient.PushBack(1);
      orient.PushBack(0);
    }
  }
  sctl::Long Nelem_ptcl_tot = 0;
  sctl::Vector<Real> ptcls_thetas;
  sctl::Vector<Real> ptcls_phis;
  if (ptcls.Dim()>0) {
    packed_sphs_conv_div(ptcls_Xcs, ptcls_rs, r1, r2); // Fit grid of small spheres inside the channel.

    sctl::Vector<sctl::Long> ptcls_(ptcls_rs.Dim());
    ptcls_ = ptcl_ord;
    ptcls.Swap(ptcls_);
    
    add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, 0); 
    for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
      Nelem_ptcl_tot += ptcls[ptcl_i];
    }
  }
  sctl::Vector<Real> NormalOrient;
  { 
    for (sctl::Long i = 0; i < Nelem + Nelem_ptcl_tot; i++) {
      for (sctl::Long j = 0; j < ElemOrder*FourierOrder*COORD_DIM; j++) {
        NormalOrient.PushBack(i < Nelem ? 1 : -1);
      }
    }
  }
  sctl::SlenderElemList<Real> elem_lst;
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient, true); 
  return std::make_tuple(elem_lst,NormalOrient_,ptcls_thetas,ptcls_phis);
}


template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::build_trefoil(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const sctl::Long ptcl_ord, const int geom_mode){
  comm_ = comm;
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  Real r_min = 0.01;
  Real r_max = 0.04;

  // varying helix radius
  auto get_r = [&r_min,&r_max](const Real& x) {
    Real angle = sctl::const_pi<Real>() * (16.*x - 28./3.); // =8*(t-pi/6), t = (x-0.5)*2pi
    return r_min + (r_max - r_min) * (0.5 * sctl::sin<Real>(angle) + 0.5);
  };

  auto get_xyz = [](const Real& x) {
    const Real xminus = x-0.5;
    const Real x4pi = 4.*sctl::const_pi<Real>()*xminus;
    const Real x8pi = 2.*x4pi;
    const Real xminus2 = xminus * xminus;
    const Real xminus5 = xminus2 * xminus2 * xminus;
    Real xcoeff = xminus2 * 4. - 1.;
    xcoeff = xcoeff / 5.;
    Real x_ = 0.5 * xminus * sctl::cos<Real>(x4pi) + 8. * xminus5 + 0.5;
    Real y_ = sctl::sin<Real>(x4pi) * xcoeff + 0.5;
    Real z_ = sctl::sin<Real>(x8pi) * xcoeff + 0.5;
    // std::cout << "inside getxyz, x = " << x_ << ", y = " << y_ << ", z = " << z_ << std::endl;

    return std::make_tuple(x_,y_,z_);
  };

  for (sctl::Long i = 0; i < Nelem; i++) {
    ElemOrderVec.PushBack(ElemOrder);
    FourierOrderVec.PushBack(FourierOrder);
    const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
    for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
      const Real x = (i+nodes[j])/Nelem;
      std::tuple<Real,Real,Real> xyz_j = get_xyz(x);
      Xc.PushBack(std::get<0>(xyz_j));
      Xc.PushBack(std::get<1>(xyz_j));
      Xc.PushBack(std::get<2>(xyz_j));
      eps.PushBack(get_r(x)); 
      // if (i < 10) {
      //   std::cout << "Xc = " << std::get<0>(xyz_j) << ", " << std::get<1>(xyz_j) << ", " << std::get<2>(xyz_j) << "; radius is " << get_r(x) << std::endl;
      // }

      orient.PushBack(0);
      orient.PushBack(0);
      orient.PushBack(1);
    }

  }

  // TODO: particles in channel in lattice.
  sctl::Long Nelem_ptcl_tot = 0;
  if (ptcls.Dim()>0) {
    // sctl::Long Nptcl = ptcls.Dim();
    packed_sphs_trefoil(ptcls_Xcs, ptcls_rs, r_min, r_max);
    // std::cout << ptcls_Xcs.Dim() << "; " << ptcls_rs.Dim() << std::endl;
    sctl::Vector<sctl::Long> ptcls_(ptcls_rs.Dim());
    ptcls_ = ptcl_ord;
    ptcls.Swap(ptcls_);
    add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
      Nelem_ptcl_tot += ptcls[ptcl_i];
    }
  }
  // std::cout << "size of Xc with ptcl: " << Xc.Dim() << std::endl;
  sctl::Vector<Real> NormalOrient;
  { // set normal to surface
    for (sctl::Long i = 0; i < Nelem + Nelem_ptcl_tot; i++) {
      for (sctl::Long j = 0; j < ElemOrder*FourierOrder*COORD_DIM; j++) {
        NormalOrient.PushBack(i < Nelem ? 1 : -1);
      }
    }
  }
  // std::cout << "right before init elem list" << std::endl;
  sctl::SlenderElemList<Real> elem_lst;
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient);
  return std::make_tuple(elem_lst,NormalOrient_);
}

template <class Real> std::tuple<bool, Real, Real, Real> PeriodicGeom<Real>::in_trefoil(Real a, Real b, Real c) {
    Real r_min = 0.01;
    Real r_max = 0.04;

    if (a>1+1e-5 || a < -1e-5) { // shift x to within [0,1].
        a = a - std::floor(a);
    }

    auto get_r = [&r_min,&r_max](const Real& x) {
        Real angle = sctl::const_pi<Real>() * (16.*x - 28./3.); // =8*(t-pi/6), t = (x-0.5)*2pi
        return r_min + (r_max - r_min) * (0.5 * sctl::sin<Real>(angle) + 0.5);
    };

    auto get_xyz = [](const Real& x) {
        const Real xminus = x-0.5;
        const Real x4pi = 4.*sctl::const_pi<Real>()*xminus;
        const Real x8pi = 2.*x4pi;
        const Real xminus2 = xminus * xminus;
        const Real xminus5 = xminus2 * xminus2 * xminus;
        Real xcoeff = xminus2 * 4. - 1.;
        xcoeff = xcoeff / 5.;
        Real x_ = 0.5 * xminus * sctl::cos<Real>(x4pi) + 8. * xminus5 + 0.5;
        Real y_ = sctl::sin<Real>(x4pi) * xcoeff + 0.5;
        Real z_ = sctl::sin<Real>(x8pi) * xcoeff + 0.5;

        return std::make_tuple(x_,y_,z_);
    };

    Real min_dist2 = 10.;
    Real closest_x = 0.;
    const int N = 4000; // resolution of the sampling
    for (int i = 0; i <= N; i++) {
        Real x = (Real)i / N; // TODO: account for distributed memory for x \in (a,b) instead of (0,1).
        std::tuple<Real,Real,Real> xchere = get_xyz(x);
        Real cx = std::get<0>(xchere);
        Real cy = std::get<1>(xchere);
        Real cz = std::get<2>(xchere);

        Real dx = cx - a;
        Real dy = cy - b;
        Real dz = cz - c;

        Real dist2 = dx*dx + dy*dy + dz*dz;

        if (dist2 < min_dist2) {
            min_dist2 = dist2;
            closest_x = x;
        }
    }

    Real r = get_r(closest_x);
    bool is_in_trefoil = (min_dist2 <= r*r);
    std::tuple<Real,Real,Real> closest_xyz = get_xyz(closest_x);
    Real xc = std::get<0>(closest_xyz);
    Real yc = std::get<1>(closest_xyz);
    Real zc = std::get<2>(closest_xyz);
    return std::make_tuple(is_in_trefoil, xc, yc, zc);
}

template <class Real> sctl::SlenderElemList<Real> PeriodicGeom<Real>::free_ptcls(const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs){
  comm_ = comm;
  sctl::Long Nelem_ptcl_tot = 0;
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs,0);
  std::cout << "Comm.Rank = " << comm.Rank() << " size of Xc list after adding particles: " << Xc.Dim() << "; size of eps: " << eps.Dim() << std::endl;
  for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
    Nelem_ptcl_tot += ptcls[ptcl_i];
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
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient);
  std::cout << "Comm.Rank = " << comm.Rank() << " size of elem here is " << elem_lst.Size() << std::endl;
  return elem_lst;
}


template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::many_ptcls1(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
  comm_ = comm;
  sctl::Long Nelem_ptcl_tot = 0;
  SCTL_ASSERT(ptcls.Dim()==0);
  SCTL_ASSERT(ptcls_Xcs.Dim()==0);
  SCTL_ASSERT(ptcls_rs.Dim()==0);
  
  // One sphere centered in space.
  ptcls_Xcs.PushBack(0.5);
  ptcls_Xcs.PushBack(0.5);
  ptcls_Xcs.PushBack(0.5);
  // ptcls_rs.PushBack(0.3);
  ptcls_rs.PushBack(0.1);

  ptcls.PushBack(Nelem);
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
  for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
    Nelem_ptcl_tot += ptcls[ptcl_i];
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
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient);
  return std::make_tuple(elem_lst,NormalOrient_);
}

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::many_ptcls3(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
  comm_ = comm;
  sctl::Long Nelem_ptcl_tot = 0;
  SCTL_ASSERT(ptcls.Dim()==0);
  SCTL_ASSERT(ptcls_Xcs.Dim()==0);
  SCTL_ASSERT(ptcls_rs.Dim()==0);
  // 3 spheres in space
  ptcls_Xcs.PushBack(0.2);
  ptcls_Xcs.PushBack(0.2);
  ptcls_Xcs.PushBack(0.2);
  ptcls_rs.PushBack(0.06);

  ptcls_Xcs.PushBack(0.7);
  ptcls_Xcs.PushBack(0.25);
  ptcls_Xcs.PushBack(0.65);
  ptcls_rs.PushBack(0.2);

  ptcls_Xcs.PushBack(0.5);
  ptcls_Xcs.PushBack(0.25);
  ptcls_Xcs.PushBack(0.3);
  ptcls_rs.PushBack(0.1);
  
  ptcls.PushBack(Nelem);
  ptcls.PushBack(Nelem);
  ptcls.PushBack(Nelem);
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
  for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
    Nelem_ptcl_tot += ptcls[ptcl_i];
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
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient);
  return std::make_tuple(elem_lst,NormalOrient_);
}

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::many_spheroids3(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs){
  const sctl::Long geom_mode = 1; // spheroids
  
  comm_ = comm;
  sctl::Long Nelem_ptcl_tot = 0;
  SCTL_ASSERT(ptcls.Dim()==0);
  SCTL_ASSERT(ptcls_Xcs.Dim()==0);
  SCTL_ASSERT(ptcls_rs.Dim()==0);
  // 3 spheroids in space
  ptcls_Xcs.PushBack(0.2);
  ptcls_Xcs.PushBack(0.2);
  ptcls_Xcs.PushBack(0.2);
  ptcls_rs.PushBack(0.06);

  ptcls_Xcs.PushBack(0.7);
  ptcls_Xcs.PushBack(0.25);
  ptcls_Xcs.PushBack(0.65);
  ptcls_rs.PushBack(0.2);

  ptcls_Xcs.PushBack(0.5);
  ptcls_Xcs.PushBack(0.65);
  ptcls_Xcs.PushBack(0.35);
  ptcls_rs.PushBack(0.2);

  ptcls.PushBack(Nelem);
  ptcls.PushBack(Nelem);
  ptcls.PushBack(Nelem);

  // centerline orientation of the particle
  sctl::Vector<Real> ptcls_thetas(3);
  sctl::Vector<Real> ptcls_phis(3);
  ptcls_thetas[0] = 0;
  ptcls_thetas[1] = sctl::const_pi<Real>() / 10.;
  ptcls_thetas[2] = 5. * sctl::const_pi<Real>() / 3.;
  ptcls_phis[0] = 0;
  ptcls_phis[1] = 0;
  ptcls_phis[2] = sctl::const_pi<Real>() / 5.;
  
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  add_particles_rotated(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode, ptcls_thetas, ptcls_phis);
  for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
    Nelem_ptcl_tot += ptcls[ptcl_i];
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
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient);
  return std::make_tuple(elem_lst,NormalOrient_,ptcls_thetas,ptcls_phis);
}

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::many_loops3(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs){
  const sctl::Long geom_mode = 3; // loops
  
  comm_ = comm;
  sctl::Long Nelem_ptcl_tot = 0;
  SCTL_ASSERT(ptcls.Dim()==0);
  SCTL_ASSERT(ptcls_Xcs.Dim()==0);
  SCTL_ASSERT(ptcls_rs.Dim()==0);
  // 3 spheroids in space
  ptcls_Xcs.PushBack(0.2);
  ptcls_Xcs.PushBack(0.2);
  ptcls_Xcs.PushBack(0.2);
  ptcls_rs.PushBack(0.06);

  ptcls_Xcs.PushBack(0.7);
  ptcls_Xcs.PushBack(0.25);
  ptcls_Xcs.PushBack(0.65);
  ptcls_rs.PushBack(0.2);

  ptcls_Xcs.PushBack(0.5);
  ptcls_Xcs.PushBack(0.65);
  ptcls_Xcs.PushBack(0.35);
  ptcls_rs.PushBack(0.2);

  ptcls.PushBack(Nelem);
  ptcls.PushBack(Nelem);
  ptcls.PushBack(Nelem);

  // centerline orientation of the particle
  sctl::Vector<Real> ptcls_thetas(3);
  sctl::Vector<Real> ptcls_phis(3);
  ptcls_thetas[0] = 0;
  ptcls_thetas[1] = sctl::const_pi<Real>() / 10.;
  ptcls_thetas[2] = 5. * sctl::const_pi<Real>() / 3.;
  ptcls_phis[0] = 0;
  ptcls_phis[1] = 0;
  ptcls_phis[2] = sctl::const_pi<Real>() / 5.;
  
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  add_particles_rotated(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode, ptcls_thetas, ptcls_phis);
  for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
    Nelem_ptcl_tot += ptcls[ptcl_i];
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
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient);
  return std::make_tuple(elem_lst,NormalOrient_,ptcls_thetas,ptcls_phis);
}

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::many_ptcls2(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Comm& comm, const sctl::Long Nptcl, sctl::Vector<sctl::Long>& ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
  comm_ = comm;
  sctl::Long Nelem_ptcl_tot = 0;
  // std::cout << "create matrix" << std::endl;
  std::string data_filename;
  sctl::Matrix<Real> Xc_from_file(Nptcl,4);
  if (Nptcl==27 || Nptcl==64||Nptcl==343|| Nptcl==512 ||Nptcl==729 ||Nptcl==1331||Nptcl==1728||Nptcl==2197) { // Nptcl = 125 and Nptcl = 1000 are ignored because they overlap with _larger files, which are used for actual scaling data.
    data_filename = "data/sphere_data_"+std::to_string(Nptcl)+"_grid.txt";
  } else {
    data_filename = "data/sphere_data_"+std::to_string(Nptcl)+"_larger.txt";
  }
  std::ifstream infile(data_filename);
  if (!infile) {
      std::cerr << "Error opening file " << data_filename << std::endl;
      SCTL_ASSERT(false);
  }
  for (sctl::Long row=0; row < Nptcl; row++) {
    for (sctl::Long col=0; col < 4; col++) {
      if (!(infile >> Xc_from_file(row,col))) {
        std::cerr << "not enough entries in data file" << std::endl;
      }
    }
  }
  for (sctl::Long i=0; i<Nptcl; i++) {
    ptcls_Xcs.PushBack(Xc_from_file(i,0));
    ptcls_Xcs.PushBack(Xc_from_file(i,1));
    ptcls_Xcs.PushBack(Xc_from_file(i,2));
    ptcls_rs.PushBack(Xc_from_file(i,3));
    ptcls.PushBack(Nelem);
  }
  sctl::Vector<Real> Xc, eps, orient;
  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  add_particles(ElemOrderVec, FourierOrderVec, Xc, eps, orient, ElemOrder, FourierOrder, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
  for (sctl::Long ptcl_i = 0; ptcl_i < ptcls.Dim(); ptcl_i++) {
    Nelem_ptcl_tot += ptcls[ptcl_i];
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
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient);
  return std::make_tuple(elem_lst,NormalOrient_);
}

/**
    Given arrays of properties for spheroids, return the SlenderELem List representing this setup.
*/
template <class Real> sctl::SlenderElemList<Real> PeriodicGeom<Real>::spheroid_system(                                                           
                                                                const sctl::Long Nelem_ptcl, 
                                                                const sctl::Long ElemOrder, 
                                                                const sctl::Long FourierOrder, 
                                                                const sctl::Vector<Real> Xcenter_lst, 
                                                                const sctl::Vector<sctl::Long> if_prolate_lst, // true if prolate
                                                                const sctl::Vector<Real> u0_lst, 
                                                                const sctl::Vector<Real> r_lst, 
                                                                const sctl::Vector<Real> theta_lst, 
                                                                const sctl::Vector<Real> phi_lst, 
                                                                const sctl::Comm comm) 
{
    const sctl::Long Nptcls = u0_lst.Dim();
    sctl::Vector<sctl::Long> ElemOrderVec(Nptcls * Nelem_ptcl);
    sctl::Vector<sctl::Long> FourierOrderVec(Nptcls * Nelem_ptcl);
    ElemOrderVec = ElemOrder;
    FourierOrderVec = FourierOrder;
    
    const auto prolate_geom = [](Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& cylindrical_radius, const Real polar_angle, const Real size, const Real u0) {
        SCTL_ASSERT(u0 > 1.); 
        x = size * u0 * sctl::cos<Real>(polar_angle);
        y = 0.;
        z = 0.;
        cylindrical_radius = size * sctl::sqrt<Real>(u0*u0 - 1) * sctl::sin<Real>(polar_angle);
        ex = 0.;
        ey = 1.;
        ez = 0.;
    };

    const auto oblate_geom = [](Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& cylindrical_radius, const Real polar_angle, const Real size, const Real u0) {
        x = size * u0 * sctl::cos<Real>(polar_angle); 
        y = 0.;
        z = 0.;
        cylindrical_radius = size * sctl::sqrt<Real>(u0*u0 + 1) * sctl::sin<Real>(polar_angle);
        ex = 0.;
        ey = 1.;
        ez = 0.;
    };

    sctl::Vector<Real> Xc, eps, orient;
    for (sctl::Long pid = 0; pid < Nptcls; pid ++) {
        const Real ptcl_size = r_lst[pid];
        const Real ptcl_u0 = u0_lst[pid];
        const Real theta_rotate = theta_lst[pid];
        const Real phi_rotate = phi_lst[pid];
        const sctl::Long if_prolate_here = if_prolate_lst[pid];
        const Real cos_theta_rotate = sctl::cos<Real>(theta_rotate);
        const Real sin_theta_rotate = sctl::sin<Real>(theta_rotate);
        const Real cos_phi_rotate = sctl::cos<Real>(phi_rotate);
        const Real sin_phi_rotate = sctl::sin<Real>(phi_rotate);
        for (sctl::Long i=0; i < Nelem_ptcl; i++) {
            const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
            for (sctl::Long j=0; j<ElemOrderVec[i]; j++) {
                const Real theta = sctl::const_pi<Real>() * (i+nodes[j])/Nelem_ptcl;
                Real x, y, z, ex, ey, ez, eps_j;
                if (if_prolate_here) {
                    prolate_geom(x,y,z,ex,ey,ez,eps_j,theta,ptcl_size,ptcl_u0);
                } else {
                    oblate_geom(x,y,z,ex,ey,ez,eps_j,theta,ptcl_size,ptcl_u0);
                }
                Real x_rotated = cos_theta_rotate * cos_phi_rotate * x - sin_phi_rotate * y + sin_theta_rotate * cos_phi_rotate * z;
                Real y_rotated = cos_theta_rotate * sin_phi_rotate * x + cos_phi_rotate * y + sin_theta_rotate * sin_phi_rotate * z;
                Real z_rotated = -sin_theta_rotate * x + cos_theta_rotate * z;
                Real ex_rotated = cos_theta_rotate * cos_phi_rotate * ex - sin_phi_rotate * ey + sin_theta_rotate * cos_phi_rotate * ez;
                Real ey_rotated = cos_theta_rotate * sin_phi_rotate * ex + cos_phi_rotate * ey + sin_theta_rotate * sin_phi_rotate * ez;
                Real ez_rotated = -sin_theta_rotate * ex + cos_theta_rotate * ez;

                Xc.PushBack(Xcenter_lst[pid*3+0]+x_rotated);
                Xc.PushBack(Xcenter_lst[pid*3+1]+y_rotated);
                Xc.PushBack(Xcenter_lst[pid*3+2]+z_rotated);
                eps.PushBack(eps_j);
                orient.PushBack(ex_rotated);
                orient.PushBack(ey_rotated);
                orient.PushBack(ez_rotated);
            }
        }
    }

    const auto init_elem_lst = [](sctl::SlenderElemList<Real>& elem_lst, const sctl::Vector<sctl::Long>& ElemOrderVec, const sctl::Vector<sctl::Long>& FourierOrderVec, const sctl::Vector<Real>& X, const sctl::Vector<Real>& R, const sctl::Vector<Real>& OrientVec, const sctl::Comm comm) {
        const sctl::Long Nelem_tot = ElemOrderVec.Dim();
        sctl::Long loc_elem_cnt, loc_elem_dsp;
        if (Nelem_tot) { // Set loc_elem_cnt, loc_elem_dsp 
            sctl::Vector<sctl::Long> node_cnt(Nelem_tot), node_dsp(Nelem_tot); node_dsp = 0;
            for (sctl::Long i = 0; i < Nelem_tot; i++) {
            node_cnt[i] = ElemOrderVec[i] * FourierOrderVec[i] * FourierOrderVec[i];
            }
            sctl::omp_par::scan(node_cnt.begin(), node_dsp.begin(), Nelem_tot);
            const sctl::Long Nnodes = node_cnt[Nelem_tot-1] + node_dsp[Nelem_tot-1];

            const sctl::Long Np = comm.Size();
            const sctl::Long rank = comm.Rank();
            sctl::Long a = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+0)/Np) - node_dsp.begin();
            sctl::Long b = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+1)/Np) - node_dsp.begin();
            if (rank == Np - 1) b = Nelem_tot;
            if (rank == 0) a = 0;
            loc_elem_cnt = b-a;
            loc_elem_dsp = a;
        } else {
            loc_elem_cnt = 0;
            loc_elem_dsp = 0;
        }

        const sctl::Vector<sctl::Long> LocElemOrder(loc_elem_cnt, (sctl::Iterator<sctl::Long>)ElemOrderVec.begin() + loc_elem_dsp, false);
        const sctl::Vector<sctl::Long> LocFourierOrder(loc_elem_cnt, (sctl::Iterator<sctl::Long>)FourierOrderVec.begin() + loc_elem_dsp, false);

        sctl::Long dsp = 0, cnt = 0;
        for (sctl::Long i = 0; i < loc_elem_dsp; i++) dsp += ElemOrderVec[i];
        for (sctl::Long i = 0; i < loc_elem_cnt; i++) cnt += ElemOrderVec[loc_elem_dsp+i];
        const sctl::Vector<Real> X_(cnt*3, (sctl::Iterator<Real>)X.begin() + dsp*3, false);
        const sctl::Vector<Real> R_(cnt, (sctl::Iterator<Real>)R.begin() + dsp, false);
        const sctl::Vector<Real> OrientVec_(cnt*3, (sctl::Iterator<Real>)OrientVec.begin() + dsp*3, false);

        elem_lst.template Init<Real>(LocElemOrder, LocFourierOrder, X_, R_, OrientVec_);  
    };

    sctl::SlenderElemList<Real> elem_lst;
    init_elem_lst(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, comm);
  
    return elem_lst;

}

template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::build_only_ptcls(const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const Real box_sidelen, const sctl::Comm& comm, const sctl::Vector<sctl::Long> ptcls, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_Xcs, const int geom_mode){
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
    many_sphs(ptcls_Xcs, ptcls_rs, Xc_box, 0, eps_box, Nptcl);
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
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient);
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

        } else if (geom_mode == 3) {
          // loop
          loop_geom(x, y, z, ex, ey, ez, eps_j, 2*theta, ptcl_r, 0.05);
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

template <class Real> void PeriodicGeom<Real>::add_spheroids(sctl::Vector<sctl::Long>& ElemOrderVec, sctl::Vector<sctl::Long>& FourierOrderVec, sctl::Vector<Real>& Xc, sctl::Vector<Real>& eps, sctl::Vector<Real>& orient, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& ptcls_Xcs, const sctl::Vector<Real>& ptcls_u0s, const sctl::Vector<sctl::Long>& ptcls_ifprolate) {
  const auto prolate_geom = [](Real& x, Real& y, Real& z, Real& cylindrical_radius, const Real polar_angle, const Real size, const Real u0) {
        SCTL_ASSERT(u0 > 1.); 
        x = size * u0 * sctl::cos<Real>(polar_angle);
        y = 0.;
        z = 0.;
        cylindrical_radius = size * sctl::sqrt<Real>(u0*u0 - 1) * sctl::sin<Real>(polar_angle);
    };

    const auto oblate_geom = [](Real& x, Real& y, Real& z, Real& cylindrical_radius, const Real polar_angle, const Real size, const Real u0) {
        x = size * u0 * sctl::cos<Real>(polar_angle); 
        y = 0.;
        z = 0.;
        cylindrical_radius = size * sctl::sqrt<Real>(u0*u0 + 1) * sctl::sin<Real>(polar_angle);
    };
  
  sctl::Long Nptcl = ptcls.Dim();
  std::cout << "in add particles, Nptcl is " << Nptcl << std::endl;
  for (sctl::Long p=0; p<Nptcl; p++) {
    const sctl::Long Nelem_sphere = ptcls[p];
    Real ptcl_size = ptcls_rs[p];
    Real ptcl_u0 = ptcls_u0s[p];
    sctl::Long if_prolate_here = ptcls_ifprolate[p];
    for (sctl::Long i = 0; i < Nelem_sphere; i++) { // add a sphere
      ElemOrderVec.PushBack(ElemOrder);
      FourierOrderVec.PushBack(FourierOrder);
      const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
      for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
        const Real theta = sctl::const_pi<Real>() * (i+nodes[j])/Nelem_sphere;
        Real x,y,z,eps_j;
        if (if_prolate_here) {
            prolate_geom(x,y,z,eps_j,theta,ptcl_size,ptcl_u0);
        } else {
            oblate_geom(x,y,z,eps_j,theta,ptcl_size,ptcl_u0);
        }
        // std::cout << "x,y,z: " << x << ", " << y  << ", " << z << ", eps = " << eps_j << std::endl;

        Xc.PushBack(ptcls_Xcs[p*3]+x);
        Xc.PushBack(ptcls_Xcs[p*3+1]+y);
        Xc.PushBack(ptcls_Xcs[p*3+2]+z);
        eps.PushBack(eps_j);

        orient.PushBack(0); // placeholders
        orient.PushBack(0);
        orient.PushBack(1);
      }
    }
  }
}

template <class Real> void PeriodicGeom<Real>::add_spheroids_rotated(sctl::Vector<sctl::Long>& ElemOrderVec, sctl::Vector<sctl::Long>& FourierOrderVec, sctl::Vector<Real>& Xc, sctl::Vector<Real>& eps, sctl::Vector<Real>& orient, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& ptcls_Xcs, const sctl::Vector<Real>& ptcls_u0s, const sctl::Vector<Real>& ptcls_thetas, const sctl::Vector<Real>& ptcls_phis, const sctl::Vector<sctl::Long>& ptcls_ifprolate) {
  const auto prolate_geom = [](Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& cylindrical_radius, const Real polar_angle, const Real size, const Real u0) {
        SCTL_ASSERT(u0 > 1.); 
        x = size * u0 * sctl::cos<Real>(polar_angle);
        y = 0.;
        z = 0.;
        cylindrical_radius = size * sctl::sqrt<Real>(u0*u0 - 1) * sctl::sin<Real>(polar_angle);
        ex = 0.;
        ey = 1.;
        ez = 0.;
    };

    const auto oblate_geom = [](Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& cylindrical_radius, const Real polar_angle, const Real size, const Real u0) {
        x = size * u0 * sctl::cos<Real>(polar_angle); 
        y = 0.;
        z = 0.;
        cylindrical_radius = size * sctl::sqrt<Real>(u0*u0 + 1) * sctl::sin<Real>(polar_angle);
        ex = 0.;
        ey = 1.;
        ez = 0.;
    };
  
  sctl::Long Nptcl = ptcls.Dim();
  // std::cout << "in add particles, Nptcl is " << Nptcl << std::endl;
  for (sctl::Long p=0; p<Nptcl; p++) {
    const sctl::Long Nelem_sphere = ptcls[p];
    Real ptcl_size = ptcls_rs[p];
    Real ptcl_u0 = ptcls_u0s[p];
    const Real theta_rotate = ptcls_thetas[p];
    const Real phi_rotate = ptcls_phis[p];
    const Real cos_theta_rotate = sctl::cos<Real>(theta_rotate);
    const Real sin_theta_rotate = sctl::sin<Real>(theta_rotate);
    const Real cos_phi_rotate = sctl::cos<Real>(phi_rotate);
    const Real sin_phi_rotate = sctl::sin<Real>(phi_rotate);
    sctl::Long if_prolate_here = ptcls_ifprolate[p];
    for (sctl::Long i = 0; i < Nelem_sphere; i++) { // add a sphere
      ElemOrderVec.PushBack(ElemOrder);
      FourierOrderVec.PushBack(FourierOrder);
      const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
      for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
        const Real theta = sctl::const_pi<Real>() * (i+nodes[j])/Nelem_sphere;
        Real x,y,z,ex,ey,ez,eps_j;
        if (if_prolate_here) {
            prolate_geom(x,y,z,ex,ey,ez,eps_j,theta,ptcl_size,ptcl_u0);
        } else {
            oblate_geom(x,y,z,ex,ey,ez,eps_j,theta,ptcl_size,ptcl_u0);
        }
        Real x_rotated = cos_theta_rotate * cos_phi_rotate * x - sin_phi_rotate * y + sin_theta_rotate * cos_phi_rotate * z;
        Real y_rotated = cos_theta_rotate * sin_phi_rotate * x + cos_phi_rotate * y + sin_theta_rotate * sin_phi_rotate * z;
        Real z_rotated = -sin_theta_rotate * x + cos_theta_rotate * z;
                
        Real ex_rotated = cos_theta_rotate * cos_phi_rotate * ex - sin_phi_rotate * ey + sin_theta_rotate * cos_phi_rotate * ez;
        Real ey_rotated = cos_theta_rotate * sin_phi_rotate * ex + cos_phi_rotate * ey + sin_theta_rotate * sin_phi_rotate * ez;
        Real ez_rotated = -sin_theta_rotate * ex + cos_theta_rotate * ez;

        Xc.PushBack(ptcls_Xcs[p*3]+x_rotated);
        Xc.PushBack(ptcls_Xcs[p*3+1]+y_rotated);
        Xc.PushBack(ptcls_Xcs[p*3+2]+z_rotated);
        eps.PushBack(eps_j);

        orient.PushBack(ex_rotated);
        orient.PushBack(ey_rotated);
        orient.PushBack(ez_rotated);
      }
    }
  }
}

// 
template <class Real> std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> PeriodicGeom<Real>::loops_system(const sctl::Vector<sctl::Long> ptcls, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Vector<Real>& ptcls_Xcs, const sctl::Vector<Real>& ptcls_major_rs, const sctl::Vector<Real>& ptcls_minor_rs, const sctl::Vector<Real>& ptcls_thetas, const sctl::Vector<Real>& ptcls_phis, sctl::Comm& comm) {
  comm_ = comm;
  sctl::Long Nptcl = ptcls.Dim();
  // std::cout << "in add particles, Nptcl is " << Nptcl << std::endl;

  sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
  sctl::Vector<Real> Xc, eps, orient;
  for (sctl::Long p=0; p<Nptcl; p++) {
    const sctl::Long Nelem_sphere = ptcls[p];
    const Real theta_rotate = ptcls_thetas[p];
    const Real phi_rotate = ptcls_phis[p];
    const Real major_r = ptcls_major_rs[p];
    const Real minor_r = ptcls_minor_rs[p];

    const Real cos_theta_rotate = sctl::cos<Real>(theta_rotate);
    const Real sin_theta_rotate = sctl::sin<Real>(theta_rotate);
    const Real cos_phi_rotate = sctl::cos<Real>(phi_rotate);
    const Real sin_phi_rotate = sctl::sin<Real>(phi_rotate);

    for (sctl::Long i = 0; i < Nelem_sphere; i++) { // add a sphere
      ElemOrderVec.PushBack(ElemOrder);
      FourierOrderVec.PushBack(FourierOrder);
      const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrder);
      for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
        const Real theta = sctl::const_pi<Real>() * (i+nodes[j])/Nelem_sphere;
        Real x, y, z, ex, ey, ez, eps_j;
        loop_geom(x, y, z, ex, ey, ez, eps_j, 2*theta, major_r, minor_r);
        // std::cout << "loop geom with major r " << major_r << ", minor r " << minor_r << ", eps_j = " << eps_j << std::endl;

        Real x_rotated = cos_theta_rotate * cos_phi_rotate * x - sin_phi_rotate * y + sin_theta_rotate * cos_phi_rotate * z;
        Real y_rotated = cos_theta_rotate * sin_phi_rotate * x + cos_phi_rotate * y + sin_theta_rotate * sin_phi_rotate * z;
        Real z_rotated = -sin_theta_rotate * x + cos_theta_rotate * z;
        Real ex_rotated = cos_theta_rotate * cos_phi_rotate * ex - sin_phi_rotate * ey + sin_theta_rotate * cos_phi_rotate * ez;
        Real ey_rotated = cos_theta_rotate * sin_phi_rotate * ex + cos_phi_rotate * ey + sin_theta_rotate * sin_phi_rotate * ez;
        Real ez_rotated = -sin_theta_rotate * ex + cos_theta_rotate * ez;

        Xc.PushBack(ptcls_Xcs[p*3]+x_rotated);
        Xc.PushBack(ptcls_Xcs[p*3+1]+y_rotated);
        Xc.PushBack(ptcls_Xcs[p*3+2]+z_rotated);
        eps.PushBack(eps_j);

        orient.PushBack(ex_rotated); 
        orient.PushBack(ey_rotated);
        orient.PushBack(ez_rotated);
      }
    }
  }
  // std::cout << "size of Xc: " << Xc.Dim() << ", size of eps" << eps.Dim() << std::endl;

  sctl::Vector<Real> NormalOrient(Xc.Dim() * FourierOrder);
  NormalOrient = -1.;
  sctl::SlenderElemList<Real> elem_lst;
  sctl::Vector<Real> NormalOrient_ = InitElemList(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, NormalOrient);
  return std::make_tuple(elem_lst,NormalOrient_);
}

template <class Real> void PeriodicGeom<Real>::add_particles_rotated(sctl::Vector<sctl::Long>& ElemOrderVec, sctl::Vector<sctl::Long>& FourierOrderVec, sctl::Vector<Real>& Xc, sctl::Vector<Real>& eps, sctl::Vector<Real>& orient, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real>& ptcls_rs, const sctl::Vector<Real>& ptcls_Xcs, const int geom_mode, const sctl::Vector<Real> ptcls_thetas, const sctl::Vector<Real> ptcls_phis) {
  sctl::Long Nptcl = ptcls.Dim();
  for (sctl::Long p=0; p<Nptcl; p++) {
    const sctl::Long Nelem_sphere = ptcls[p];
    const Real theta_rotate = ptcls_thetas[p];
    const Real phi_rotate = ptcls_phis[p];
    const Real cos_theta_rotate = sctl::cos<Real>(theta_rotate);
    const Real sin_theta_rotate = sctl::sin<Real>(theta_rotate);
    const Real cos_phi_rotate = sctl::cos<Real>(phi_rotate);
    const Real sin_phi_rotate = sctl::sin<Real>(phi_rotate);
    // std::cout << "DEBUG: theta is " << theta_rotate << ", phi is " << phi_rotate <<"; cos theta, e.g. = " << cos_theta_rotate << std::endl;
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
        } else if (geom_mode == 3) {
          // loop
          loop_geom(x, y, z, ex, ey, ez, eps_j, 2*theta, ptcl_r, 0.05);
        } else {
          SCTL_ASSERT(false); // not implemented
        }
        Real x_rotated = cos_theta_rotate * cos_phi_rotate * x - sin_phi_rotate * y + sin_theta_rotate * cos_phi_rotate * z;
        Real y_rotated = cos_theta_rotate * sin_phi_rotate * x + cos_phi_rotate * y + sin_theta_rotate * sin_phi_rotate * z;
        Real z_rotated = -sin_theta_rotate * x + cos_theta_rotate * z;
        Real ex_rotated = cos_theta_rotate * cos_phi_rotate * ex - sin_phi_rotate * ey + sin_theta_rotate * cos_phi_rotate * ez;
        Real ey_rotated = cos_theta_rotate * sin_phi_rotate * ex + cos_phi_rotate * ey + sin_theta_rotate * sin_phi_rotate * ez;
        Real ez_rotated = -sin_theta_rotate * ex + cos_theta_rotate * ez;
        // std::cout << "DEBUG, centerline node was " << x << ", " << y <<", " << z << ", rotated to be "<<x_rotated <<", "<<y_rotated<<", "<<z_rotated<<std::endl;

        Xc.PushBack(ptcls_Xcs[p*3]+x_rotated);
        Xc.PushBack(ptcls_Xcs[p*3+1]+y_rotated);
        Xc.PushBack(ptcls_Xcs[p*3+2]+z_rotated);
        eps.PushBack(eps_j);

        orient.PushBack(ex_rotated);
        orient.PushBack(ey_rotated);
        orient.PushBack(ez_rotated);
      }
    }
  }
}

template <class Real> sctl::Vector<Real> PeriodicGeom<Real>::InitElemList(sctl::SlenderElemList<Real>& elem_lst, const sctl::Vector<sctl::Long>& ElemOrder, const sctl::Vector<sctl::Long>& FourierOrder, const sctl::Vector<Real>& X, const sctl::Vector<Real>& R, const sctl::Vector<Real>& OrientVec, const sctl::Vector<Real>& NormalOrient, bool use_orient) {
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
  } else {
    loc_elem_cnt = 0;
    loc_elem_dsp = 0;
  }

  const sctl::Vector<sctl::Long> LocElemOrder(loc_elem_cnt, (sctl::Iterator<sctl::Long>)ElemOrder.begin() + loc_elem_dsp, false);
  const sctl::Vector<sctl::Long> LocFourierOrder(loc_elem_cnt, (sctl::Iterator<sctl::Long>)FourierOrder.begin() + loc_elem_dsp, false);

  sctl::Long dsp = 0, cnt = 0;
  for (sctl::Long i = 0; i < loc_elem_dsp; i++) dsp += ElemOrder[i];
  for (sctl::Long i = 0; i < loc_elem_cnt; i++) cnt += ElemOrder[loc_elem_dsp+i];
  const sctl::Vector<Real> X_(cnt*COORD_DIM, (sctl::Iterator<Real>)X.begin() + dsp*COORD_DIM, false);
  const sctl::Vector<Real> R_(cnt, (sctl::Iterator<Real>)R.begin() + dsp, false);
  const sctl::Vector<Real> OrientVec_(cnt*COORD_DIM, (sctl::Iterator<Real>)OrientVec.begin() + dsp*COORD_DIM, false);
  const sctl::Vector<Real> NormalOrient_(cnt*FourierOrder[0]*COORD_DIM, (sctl::Iterator<Real>)NormalOrient.begin() + dsp*FourierOrder[0]*COORD_DIM, false); // TODO: this assumes all Fourier orders are the same.

  if (use_orient) {
    elem_lst.template Init<Real>(LocElemOrder, LocFourierOrder, X_, R_, OrientVec_);  
  } else {
    elem_lst.template Init<Real>(LocElemOrder, LocFourierOrder, X_, R_);  
  }
  

  return NormalOrient_;
}

template <class Real> std::tuple<sctl::Long,sctl::Long> PeriodicGeom<Real>::GetGlobalIdx(const sctl::Vector<sctl::Long>& ElemOrder, const sctl::Vector<sctl::Long>& FourierOrder, const sctl::Comm& comm) {
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

    const sctl::Long Np = comm.Size();
    const sctl::Long rank = comm.Rank();
    sctl::Long a = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+0)/Np) - node_dsp.begin();
    sctl::Long b = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+1)/Np) - node_dsp.begin();
    if (rank == Np - 1) b = Nelem;
    if (rank == 0) a = 0;
    loc_elem_cnt = b-a;
    loc_elem_dsp = a;

    if (0 && !comm.Rank()) { // Print partitioning
      std::cout<<"Partitioning: ";
      for (sctl::Long i = 0; i < comm.Size(); i++) {
        sctl::Long a = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(i+0)/Np) - node_dsp.begin();
        sctl::Long b = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(i+1)/Np) - node_dsp.begin();
        if (i == Np - 1) b = Nelem;
        if (i == 0) a = 0;
        std::cout<<b-a<<' ';
      }
      std::cout<<'\n';

      std::cout<<"Weight: ";
      for (sctl::Long i = 0; i < comm.Size(); i++) {
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
  return std::make_tuple(loc_elem_cnt,loc_elem_dsp);
}

template <class Real> sctl::Vector<Real> PeriodicGeom<Real>::vec_nbr_copy(const sctl::Vector<Real> X, const sctl::Integer nbr_range, const sctl::Integer peri_mode) {
  // sctl::Long Nrepeat = 2*nbr_range*peri_mode + 1;
  sctl::Long Nrepeat = 2*nbr_range + 1;
  if (peri_mode==3) {
    Nrepeat = Nrepeat * Nrepeat * Nrepeat;
  } else if (peri_mode==2) {
    Nrepeat = Nrepeat * Nrepeat;
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
  sctl::Long N = X.Dim();
  if (nbr_range>0) { // duplicate geomtry to add images
    sctl::Long Nrepeat = 2*nbr_range + 1;
    if (peri_mode==3) {
      Nrepeat = Nrepeat * Nrepeat * Nrepeat;
    } else if (peri_mode==2) {
      Nrepeat = Nrepeat * Nrepeat;
    }
    sctl::Vector<Real> Xc_(N*Nrepeat);
    sctl::Long Xcptr = 0;
    if (peri_mode == 1) {
      for (sctl::Long k0 = -nbr_range; k0 <= nbr_range; k0++) {
        Xcptr = (k0+nbr_range)*N;
        for (sctl::Long i = 0; i < N/3; i++) { // shift in x
          Xc_[Xcptr + i*3+0] = X[i*3+0] + k0;
          Xc_[Xcptr + i*3+1] = X[i*3+1];
          Xc_[Xcptr + i*3+2] = X[i*3+2];
        }
      }
    } else if (peri_mode == 2) {
      for (sctl::Long k1 = -nbr_range; k1 <= nbr_range; k1++) {
        for (sctl::Long k0 = -nbr_range; k0 <= nbr_range; k0++) {
          Xcptr = (k0+nbr_range)*N + (k1+nbr_range)*(2*nbr_range+1)*N;
          for (sctl::Long i = 0; i < X.Dim()/3; i++) { // shift in x
            Xc_[Xcptr + i*3+0] = X[i*3+0] + k0;
            Xc_[Xcptr + i*3+1] = X[i*3+1] + k1;
            Xc_[Xcptr + i*3+2] = X[i*3+2];
          }
        }
      }
    } else if (peri_mode == 3) {
      for (sctl::Long k2 = -nbr_range; k2 <= nbr_range; k2++) {
        for (sctl::Long k1 = -nbr_range; k1 <= nbr_range; k1++) {
          for (sctl::Long k0 = -nbr_range; k0 <= nbr_range; k0++) {
            Xcptr = (k0+nbr_range)*N + (k1+nbr_range)*(2*nbr_range+1)*N + (k2+nbr_range)*(2*nbr_range+1)*(2*nbr_range+1)*N;
            for (sctl::Long i = 0; i < X.Dim()/3; i++) {
              Xc_[Xcptr + i*3+0] = X[i*3+0] + k0;
              Xc_[Xcptr + i*3+1] = X[i*3+1] + k1;
              Xc_[Xcptr + i*3+2] = X[i*3+2] + k2;
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

// Replaced by XsectionVis
template <class Real> sctl::Vector<Real> PeriodicGeom<Real>::form_targets(const sctl::Long r_ord, const sctl::Long azi_ord, const sctl::SlenderElemList<Real>& elem_lst, const sctl::Comm& comm) {
  // Create streakline targets at Ngroups cross sections, divided evenly among processes in comm, and return the local list.
  sctl::Vector<Real> Xtrgs;

  // elems in elem_lst stored distributively, so use Ngroups_per_process.
  sctl::Vector<Real> s_param, sin_theta, cos_theta;
  s_param.PushBack(0.); // only take starting value of panel
  // s_param.PushBack(1);
  for (sctl::Long i = 0; i < azi_ord; i++) {
    const Real t = i/(Real)azi_ord;
    sin_theta.PushBack(sctl::sin<Real>(2*sctl::const_pi<Real>()*t));
    cos_theta.PushBack(sctl::cos<Real>(2*sctl::const_pi<Real>()*t));
  }
  if (!comm.Rank()) {
    std::cout << "length of sin theta is " << sin_theta.Dim() << std::endl;
  }

  sctl::Long Nelem = elem_lst.Size();
  for (sctl::Long elem_idx = 0; elem_idx < Nelem; elem_idx++) {
    const Real t_order_inv = 1/(Real)azi_ord;
    const Real r_order_inv = (1-1e-3)/(Real)(r_ord-1);
    sctl::Vector<Real> X_, Xc(3);
    elem_lst.GetGeom(&X_, nullptr, nullptr, nullptr, nullptr, s_param, sin_theta, cos_theta, elem_idx);
    // std::cout << "on Rank " << comm.Rank() << ", size of elem lst target points is " << X_.Dim() << std::endl;
    Xc = 0;
    for (sctl::Long j = 0; j < azi_ord; j++) {
      for (sctl::Long l = 0; l < 3; l++) {
        Xc[l] += X_[j*3+l] * t_order_inv;
      }
    }
    for (sctl::Long j = 0; j < azi_ord; j++) {
      for (sctl::Long k = 0; k < azi_ord; k++) {
        for (sctl::Long l = 0; l < 3; l++) {
          Xtrgs.PushBack((X_[j*3+l]-Xc[l])*k*r_order_inv + Xc[l]);
        }
      }
    }
  }
  
  return Xtrgs;
  
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
        Real u0 = 1.1; // todo: take as input..
        int if_prolate = 1;
        outside = (outside && outside_spheroid(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0],u0,if_prolate));
      } else if (geom_mode == 3) {
        outside = (outside && outside_loop(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0],0.05)); // ptcls_rs for loops store the size of the loop, so loop_rad, rather than the "thickness", which is hardcoded to be 0.025.
      } else {
        SCTL_ASSERT(false); // not implemented
      }
      
      if (!outside) {
        // std::cout << "rank " << comm_.Rank() << ", inside particle" << std::endl;
        // std::cout << "inside particle, remove." << std::endl;
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

template <class Real> std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> PeriodicGeom<Real>::filter_target_rotated(const sctl::Vector<Real> X, const sctl::Vector<sctl::Long> ptcls, const sctl::Vector<Real> ptcls_rs, const sctl::Vector<Real> ptcls_Xcs, const int geom_mode, const sctl::Vector<Real> ptcls_thetas, const sctl::Vector<Real> ptcls_phis) {
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
      const auto ptheta = ptcls_thetas.begin() + j;
      const auto pphi = ptcls_phis.begin() + j;
      if (geom_mode==0) {
        outside = (outside && outside_sphere(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0]));
      } else if (geom_mode==1) {
        Real u0 = 1.1;
        int if_prolate = 1;
        outside = (outside && outside_spheroid_rotated(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0],u0,if_prolate,ptheta[0],pphi[0]));
      } else if (geom_mode == 3) {
        outside = (outside && outside_loop_rotated(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0],0.05,ptheta[0],pphi[0])); // hardcoded geometry r.n.
      } else {
        SCTL_ASSERT(false); // not implemented
      }
      
      if (!outside) {
        // std::cout << "rank " << comm_.Rank() << ", inside particle" << std::endl;
        // std::cout << "target "<< i << " inside particle, remove." << std::endl;
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

template <class Real> std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> PeriodicGeom<Real>::filter_spheroids(
                                                                                                const sctl::Vector<Real> X, 
                                                                                                const sctl::Vector<Real> r_all, 
                                                                                                const sctl::Vector<Real> u0_all, 
                                                                                                const sctl::Vector<Real> Xcenter_all, 
                                                                                                const sctl::Vector<sctl::Long> if_prolate_all) 
{
    const sctl::Long N = X.Dim()/3; // number of targets
    const sctl::Long Nptcl = r_all.Dim(); // number of particles
    sctl::Vector<sctl::Long> filtered_inds(N); // indicate whether the target was inside particle or not.
    filtered_inds.SetZero();
    sctl::Vector<Real> Xout; // collection of targets outside all ptcls.
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        bool outside = true;
        for (sctl::Long j = 0; j < Nptcl; j++) {
            const auto pXc = Xcenter_all.begin() + j*3;
            const auto pr = r_all.begin() + j;
            const auto pu0 = u0_all.begin() + j;
            const auto p_if_prolate = if_prolate_all.begin()+j;
            outside = (outside && outside_spheroid(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0],pu0[0],p_if_prolate[0]));
            
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
    return std::make_tuple(Xout,filtered_inds);
}

template <class Real> std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> PeriodicGeom<Real>::filter_spheroids_rotated(
                                                                                                const sctl::Vector<Real> X, 
                                                                                                const sctl::Vector<Real> r_all, 
                                                                                                const sctl::Vector<Real> u0_all, 
                                                                                                const sctl::Vector<Real> Xcenter_all, 
                                                                                                const sctl::Vector<sctl::Long> if_prolate_all,
                                                                                                const sctl::Vector<Real> theta_all, 
                                                                                                const sctl::Vector<Real> phi_all) 
{
    const sctl::Long N = X.Dim()/3; // number of targets
    const sctl::Long Nptcl = r_all.Dim(); // number of particles
    sctl::Vector<sctl::Long> filtered_inds(N); // indicate whether the target was inside particle or not.
    filtered_inds.SetZero();
    sctl::Vector<Real> Xout; // collection of targets outside all ptcls.
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        bool outside = true;
        for (sctl::Long j = 0; j < Nptcl; j++) {
            const auto pXc = Xcenter_all.begin() + j*3;
            const auto pr = r_all.begin() + j;
            const auto pu0 = u0_all.begin() + j;
            const auto ptheta = theta_all.begin() + j;
            const auto pphi = phi_all.begin() + j;
            const auto p_if_prolate = if_prolate_all.begin()+j;
            outside = (outside && outside_spheroid_rotated(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],pr[0],pu0[0],p_if_prolate[0],ptheta[0],pphi[0]));
            
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
    return std::make_tuple(Xout,filtered_inds);
}

template <class Real> std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> PeriodicGeom<Real>::filter_loops_rotated(
                                                                                                const sctl::Vector<Real> X, 
                                                                                                const sctl::Vector<Real> major_r_all, 
                                                                                                const sctl::Vector<Real> minor_r_all, 
                                                                                                const sctl::Vector<Real> Xcenter_all, 
                                                                                                const sctl::Vector<Real> theta_all, 
                                                                                                const sctl::Vector<Real> phi_all) 
{
    const sctl::Long N = X.Dim()/3; // number of targets
    const sctl::Long Nptcl = major_r_all.Dim(); // number of particles
    sctl::Vector<sctl::Long> filtered_inds(N); // indicate whether the target was inside particle or not.
    filtered_inds.SetZero();
    sctl::Vector<Real> Xout; // collection of targets outside all ptcls.
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        bool outside = true;
        for (sctl::Long j = 0; j < Nptcl; j++) {
            const auto pXc = Xcenter_all.begin() + j*3;
            const auto major_r = major_r_all.begin() + j;
            const auto minor_r = minor_r_all.begin() + j;
            const auto ptheta = theta_all.begin() + j;
            const auto pphi = phi_all.begin() + j;
            outside = (outside && outside_loop_rotated(x[0],x[1],x[2],pXc[0],pXc[1],pXc[2],major_r[0],minor_r[0],ptheta[0],pphi[0]));
            
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

template <class Real> bool PeriodicGeom<Real>::outside_spheroid(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real a, const Real u0, const int if_prolate) {
  const Real v1 = x1-pXc1;
  const Real v2 = x2-pXc2;
  const Real v3 = x3-pXc3;
  const Real buffer = 1.25;

  Real A;
  if (if_prolate) {
    A = a * sctl::sqrt(u0*u0-1);
  } else {
    A = a * sctl::sqrt(u0*u0+1);
  }
  const Real C = a * u0;
  const Real A2inv = 1. / (A*A); // TODO: numerical instability?
  const Real C2inv = 1. / (C*C); 
  const Real x1sq = v1*v1;
  const Real x2sq = v2*v2;
  const Real x3sq = v3*v3;
  if (x1sq * C2inv + (x2sq+x3sq) * A2inv < buffer) { 
    return false;
  } else {
    return true;
  }
}

template <class Real> bool PeriodicGeom<Real>::outside_spheroid_rotated(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real a, const Real u0, const int if_prolate, const Real ptheta, const Real pphi) {
  const Real v1 = x1-pXc1;
  const Real v2 = x2-pXc2;
  const Real v3 = x3-pXc3;
  const Real buffer = 1.05;
  
  // Counter rotate target -- R^{-1} = R^T
  const Real cos_theta_rotate = sctl::cos<Real>(ptheta);
  const Real sin_theta_rotate = sctl::sin<Real>(ptheta);
  const Real cos_phi_rotate = sctl::cos<Real>(pphi);
  const Real sin_phi_rotate = sctl::sin<Real>(pphi);
  Real v1_rotated = cos_theta_rotate * cos_phi_rotate * v1 + sin_phi_rotate * cos_theta_rotate * v2 - sin_theta_rotate * v3;
  Real v2_rotated = - sin_phi_rotate * v1 + cos_phi_rotate * v2;
  Real v3_rotated = cos_phi_rotate * sin_theta_rotate * v1 + sin_phi_rotate * sin_theta_rotate * v2 + cos_theta_rotate * v3; 
  // std::cout << "difference vector before rotation: " << v1 << ", " << v2 << ", " << v3 << "; " << std::endl;
  // std::cout << "Rotation by theta = " << ptheta << ", phi = " << pphi << "; " << std::endl;
  // std::cout << "New rotated vector is " << v1_rotated<<", " << v2_rotated << ", " << v3_rotated << ". " << std::endl;

  Real A;
  if (if_prolate) {
    A = a * sctl::sqrt(u0*u0-1);
  } else {
    A = a * sctl::sqrt(u0*u0+1);
  }
  const Real C = a * u0;
  const Real A2inv = 1. / (A*A); // TODO: numerical instability?
  const Real C2inv = 1. / (C*C); 
  const Real x1sq = v1_rotated*v1_rotated; // should be same as v1^2, etc
  const Real x2sq = v2_rotated*v2_rotated;
  const Real x3sq = v3_rotated*v3_rotated;
  if (x1sq * C2inv + (x2sq+x3sq) * A2inv < buffer) { // with buffer layer.
    // std::cout << "INSIDE spheroid centered at = " << pXc1 << ", " << pXc2 << ", " << pXc3 << "; major axis is " << x1sq * C2inv << ", minor: " << (x2sq+x3sq) * A2inv << std::endl;
    return false;
  } else {
    // std::cout << "OUTSIDE spheroid centered at = " << pXc1 << ", " << pXc2 << ", " << pXc3 << "; major axis is " << x1sq * C2inv << ", minor: " << (x2sq+x3sq) * A2inv << std::endl;
    return true;
  }
}

template <class Real> bool PeriodicGeom<Real>::outside_loop(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real major_r, const Real minor_r) {
  
  const Real v1 = x1-pXc1;
  const Real v2 = x2-pXc2;
  const Real v3 = x3-pXc3;
  const Real buffer = 1.25;

  bool A = v1*v1 + v2*v2  > (major_r-minor_r*buffer)*(major_r-minor_r*buffer);
  bool B = v1*v1 + v2*v2  < (major_r+minor_r*buffer)*(major_r+minor_r*buffer); // x-y direction within pr ring
  bool C = v3 < minor_r * buffer && v3 > -minor_r * buffer; // z direction between [-pr, pr]
  bool inside = A && B && C; // includes more points than necessary to be "inside"
  // std::cout << "target at " << x1 << ", " << x2 << ", " << x3 << ", particle center at " << pXc1 << ", " << pXc2 << ", " << pXc3 << ". Large radius of loop is " << loop_rad << ", radius of the ringlet is " << pr << std::endl;
  // std::cout << "x,y distance ^2 to center: " << v1*v1+v2*v2 << ", radius min ^2 = " << (loop_rad-pr)*(loop_rad-pr) << ", max ^2 = " << (loop_rad+pr)*(loop_rad+pr) << std::endl;

  return !inside;
}

template <class Real> bool PeriodicGeom<Real>::outside_loop_rotated(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real major_r, const Real minor_r, const Real ptheta, const Real pphi) {
  const Real v1 = x1-pXc1;
  const Real v2 = x2-pXc2;
  const Real v3 = x3-pXc3;
  const Real buffer = 1.25;

  // Counter rotate target -- R^{-1} = R^T
  const Real cos_theta_rotate = sctl::cos<Real>(ptheta);
  const Real sin_theta_rotate = sctl::sin<Real>(ptheta);
  const Real cos_phi_rotate = sctl::cos<Real>(pphi);
  const Real sin_phi_rotate = sctl::sin<Real>(pphi);
  Real v1_rotated = cos_theta_rotate * cos_phi_rotate * v1 + sin_phi_rotate * cos_theta_rotate * v2 - sin_theta_rotate * v3;
  Real v2_rotated = - sin_phi_rotate * v1 + cos_phi_rotate * v2;
  Real v3_rotated = cos_phi_rotate * sin_theta_rotate * v1 + sin_phi_rotate * sin_theta_rotate * v2 + cos_theta_rotate * v3; 

  bool A = v1_rotated*v1_rotated + v2_rotated*v2_rotated  > (major_r-minor_r*buffer)*(major_r-minor_r*buffer); // decrease inner radius to filter out more close eval points.
  bool B = v1_rotated*v1_rotated + v2_rotated*v2_rotated  < (major_r+minor_r*buffer)*(major_r+minor_r*buffer); // x-y direction within pr ring
  bool C = v3_rotated < minor_r * buffer && v3_rotated > -minor_r * buffer; // z direction between [-pr, pr]
  bool inside = A && B && C; // includes more points than necessary to be "inside"
  // std::cout << "target at " << x1 << ", " << x2 << ", " << x3 << ", particle center at " << pXc1 << ", " << pXc2 << ", " << pXc3 << ". Large radius of loop is " << loop_rad << ", radius of the ringlet is " << pr << std::endl;
  // std::cout << "x,y distance ^2 to center: " << v1_rotated*v1_rotated+v2_rotated*v2_rotated << ", radius min ^2 = " << (loop_rad-pr)*(loop_rad-pr) << ", max ^2 = " << (loop_rad+pr)*(loop_rad+pr) << std::endl;
  // std::cout << "z difference value: " << v3_rotated << ", one sided radius for z bound: " << pr << std::endl;

  return !inside;
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

template <class Real> void PeriodicGeom<Real>::packed_sphs_conv_div(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, const Real r1, const Real r2) {
  srand48(2);

  ptcls_Xcs.ReInit(0);
  ptcls_rs.ReInit(0);

  // auto get_r = [&r1, &r2](const Real& x) {
  //   if (x < 0.1) {
  //     return 2*r1 + r2;
  //   } else if (x < 0.3) {
  //     return r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-0.1) / 0.2)+r1+r2;
  //   } else if (x < 0.7) {
  //     return r2;
  //   } else if (x < 0.9) {
  //     return r1+r2 - r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-0.7) / 0.2);
  //   } else {
  //     return 2*r1 + r2;
  //   }
  // };

  auto get_r = [&r1, &r2](const Real& x) {
    Real cutoff1 = 0.25; // end of large radius left
    Real cutoff2 = 0.35; // start of small radius left
    Real cutoff3 = 0.65; // end of small radius right
    Real cutoff4 = 0.75; // start of large radius right
    if (x < cutoff1) {
      return 2*r1 + r2;
    } else if (x < cutoff2) {
      return r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-cutoff1) / (cutoff2-cutoff1))+r1+r2;
    } else if (x < cutoff3) {
      return r2;
    } else if (x < cutoff4) {
      return r1+r2 - r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-cutoff3) / (cutoff4-cutoff3));
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
        while (iter_cnt < 30 && ptcl_r_loc >= 0.005 && (Y*Y + Z*Z >= (r_channel-(2. * ptcl_r_loc))*(r_channel-(2. * ptcl_r_loc)))) {
          // reduce radius until particle fits inside channel.
          ptcl_r_loc *= 0.95;
          iter_cnt += 1;
        } 
        if (iter_cnt < 30 && ptcl_r_loc >= 0.005) {
          ptcls_Xcs.PushBack(X);
          ptcls_Xcs.PushBack(Y + 0.5);
          ptcls_Xcs.PushBack(Z + 0.5);
          ptcls_rs.PushBack(ptcl_r_loc);
        }
        
      }
    } 

  }
}

template <class Real> void PeriodicGeom<Real>::packed_spheroids_conv_div(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_u0s, sctl::Vector<Real>& ptcls_thetas, sctl::Vector<Real>& ptcls_phis, sctl::Vector<sctl::Long>& ptcls_ifprolate, const Real r1, const Real r2, const sctl::Long N) {
  // const int N = 2560; // use spheroid system in unit box for a densely packed system
  // const int N = 80;
  const Real buffer = 1.15;

  // Clear all inputting arrays for populating later
  ptcls_Xcs.ReInit(0);
  ptcls_rs.ReInit(0);
  ptcls_u0s.ReInit(0);
  ptcls_ifprolate.ReInit(0);

  // Get original system from file
  sctl::Vector<Real> Xcenter_all, u0_all, r_all, theta_all, phi_all;
  sctl::Vector<sctl::Long> if_prolate_all;
  std::string filename_Xcenter, filename_u0, filename_size, filename_ifprolate;
  filename_Xcenter = "data/bie_spheroids/Xcenter_"+std::to_string(N)+".txt";
  filename_u0 = "data/bie_spheroids/u0_"+std::to_string(N)+".txt";
  filename_size = "data/bie_spheroids/size_"+std::to_string(N)+".txt";
  filename_ifprolate = "data/bie_spheroids/ifprolate_"+std::to_string(N)+".txt";

  const auto read_file_double = [](const std::string filename, sctl::Vector<Real>& Out, const sctl::Long expected_length) {
    Out.ReInit(expected_length);
    Out.SetZero();
    
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error opening file " << filename << std::endl;
        SCTL_ASSERT(false);
    }

    for (sctl::Long i=0; i<expected_length; i++ ) {
        if (!(infile >> Out[i])) {
            std::cerr << "File cut short at " << i << ", before expected length " << expected_length << std::endl;
        }
    }
  };
  const auto read_file_int = [](const std::string filename, sctl::Vector<sctl::Long>& Out, const sctl::Long expected_length) {
    Out.ReInit(expected_length);
    Out.SetZero();
    
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error opening file " << filename << std::endl;
        SCTL_ASSERT(false);
    }

    for (sctl::Long i=0; i<expected_length; i++ ) {
        if (!(infile >> Out[i])) {
            std::cerr << "File cut short at " << i << ", before expected length " << expected_length << std::endl;
        }
    }
  };

  read_file_double(filename_Xcenter, Xcenter_all, N * 3);
  read_file_double(filename_u0, u0_all, N);
  read_file_double(filename_size, r_all, N);
  read_file_int(filename_ifprolate, if_prolate_all, N);

  // Randomly give rotation angles
  srand48(2);
  for (sctl::Long i=0; i<N; i++) {
      const Real theta_rotate = drand48() * sctl::const_pi<Real>() * 2.;
      const Real phi_rotate = drand48() * sctl::const_pi<Real>();
      // const Real theta_rotate = 0.;
      // const Real phi_rotate = 0.;
      theta_all.PushBack(theta_rotate);
      phi_all.PushBack(phi_rotate);
  }

  // Conv-div channel
  const auto get_r = [&r1, &r2](const Real& x) {
    const Real cutoff1 = 0.25; // end of large radius left
    const Real cutoff2 = 0.35; // start of small radius left
    const Real cutoff3 = 0.65; // end of small radius right
    const Real cutoff4 = 0.75; // start of large radius right
    if (x < cutoff1) {
      return 2*r1 + r2;
    } else if (x < cutoff2) {
      return r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-cutoff1) / (cutoff2-cutoff1))+r1+r2;
    } else if (x < cutoff3) {
      return r2;
    } else if (x < cutoff4) {
      return r1+r2 - r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-cutoff3) / (cutoff4-cutoff3));
    } else {
      return 2*r1 + r2;
    }
  };

  // // Conv-div channel
  // const auto get_r = [&r1, &r2](const Real& x) {
  //   return 2*r1 + r2;
  // };

  // Populate arrays with only spheroids that are in the channel
  for (sctl::Long n=0; n<N; n++) {
    // Spheroid data
    sctl::Vector<Real> Xc_here(3, (sctl::Iterator<Real>) Xcenter_all.begin() + n*3, false);
    sctl::Vector<Real> u0_here(1, (sctl::Iterator<Real>) u0_all.begin() + n, false);
    sctl::Vector<Real> r_here(1, (sctl::Iterator<Real>) r_all.begin() + n, false);
    sctl::Vector<Real> theta_here(1, (sctl::Iterator<Real>) theta_all.begin() + n, false);
    sctl::Vector<Real> phi_here(1, (sctl::Iterator<Real>) phi_all.begin() + n, false);
    sctl::Vector<sctl::Long> if_prolate_here(1, (sctl::Iterator<sctl::Long>) if_prolate_all.begin() + n, false);

    Real channel_r_here = get_r(Xc_here[0]);
    Real major_axis;
    if (if_prolate_here[0]) { // prolate, use C
      major_axis = r_here[0] * u0_here[0];
    } else { // oblate, use A
      major_axis = r_here[0] * sctl::sqrt(u0_here[0]*u0_here[0]+1.);
    }
    Real Xc_r2 = (Xc_here[1]-0.5)*(Xc_here[1]-0.5) + (Xc_here[2]-0.5)*(Xc_here[2]-0.5);

    if ((sctl::sqrt<Real>(Xc_r2) + major_axis)*buffer < channel_r_here) {
      ptcls_Xcs.PushBack(Xc_here[0]);
      ptcls_Xcs.PushBack(Xc_here[1]);
      ptcls_Xcs.PushBack(Xc_here[2]);
      ptcls_rs.PushBack(r_here[0]);
      ptcls_u0s.PushBack(u0_here[0]);
      ptcls_thetas.PushBack(theta_here[0]);
      ptcls_phis.PushBack(phi_here[0]);
      ptcls_ifprolate.PushBack(if_prolate_here[0]);
    }
  }
}

// CHOCO FEB 2026: DEBUG VSLIP using spheres instead of spheroids.
template <class Real> void PeriodicGeom<Real>::packed_spheres_conv_div(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, sctl::Vector<Real>& ptcls_thetas, sctl::Vector<Real>& ptcls_phis, const Real r1, const Real r2) {
  const int N = 2000; 
  const Real buffer = 1.15;

  // Clear all inputting arrays for populating later
  ptcls_Xcs.ReInit(0);
  ptcls_rs.ReInit(0);
  ptcls_thetas.ReInit(0);
  ptcls_phis.ReInit(0);

  // Get original system from file
  sctl::Matrix<Real> Xc_from_file(N,4);
  std::string data_filename = "data/sphere_data_"+std::to_string(N)+"_larger.txt";
  std::ifstream infile(data_filename);
  if (!infile) {
      std::cerr << "Error opening file " << data_filename << std::endl;
      SCTL_ASSERT(false);
  }
  for (sctl::Long row=0; row < N; row++) {
    for (sctl::Long col=0; col < 4; col++) {
      if (!(infile >> Xc_from_file(row,col))) {
        std::cerr << "not enough entries in data file" << std::endl;
      }
    }
  }

  // Conv-div channel
  const auto get_r = [&r1, &r2](const Real& x) {
    const Real cutoff1 = 0.25; // end of large radius left
    const Real cutoff2 = 0.35; // start of small radius left
    const Real cutoff3 = 0.65; // end of small radius right
    const Real cutoff4 = 0.75; // start of large radius right
    if (x < cutoff1) {
      return 2*r1 + r2;
    } else if (x < cutoff2) {
      return r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-cutoff1) / (cutoff2-cutoff1))+r1+r2;
    } else if (x < cutoff3) {
      return r2;
    } else if (x < cutoff4) {
      return r1+r2 - r1*sctl::cos<Real>(sctl::const_pi<Real>() * (x-cutoff3) / (cutoff4-cutoff3));
    } else {
      return 2*r1 + r2;
    }
  };

  srand48(2);
  for (sctl::Long i=0; i<N; i++) {
    Real x1 = Xc_from_file(i,0);
    Real x2 = Xc_from_file(i,1);
    Real x3 = Xc_from_file(i,2);
    Real r_here = Xc_from_file(i,3);
    Real channel_r_here = get_r(x1);
    Real Xc_r2 = (x2-0.5)*(x2-0.5) + (x3-0.5)*(x3-0.5);
    // std::cout <<" particle " << i << ", centered at " << x1 << ","<<x2<<","<<x3<<", radius " << r_here << "; channel radius here is " << channel_r_here << std::endl;
    if ((sctl::sqrt<Real>(Xc_r2) + r_here) * buffer < channel_r_here) {
      ptcls_Xcs.PushBack(x1);
      ptcls_Xcs.PushBack(x2);
      ptcls_Xcs.PushBack(x3);
      ptcls_rs.PushBack(r_here);

      const Real theta_rotate = drand48() * sctl::const_pi<Real>() * 2.;
      const Real phi_rotate = drand48() * sctl::const_pi<Real>();
      ptcls_thetas.PushBack(theta_rotate);
      ptcls_phis.PushBack(phi_rotate);
    }
  }
}

template <class Real> void PeriodicGeom<Real>::packed_sphs_trefoil(sctl::Vector<Real>& ptcls_Xcs, sctl::Vector<Real>& ptcls_rs, const Real r_min, const Real r_max) {
  srand48(2);

  ptcls_Xcs.ReInit(0);
  ptcls_rs.ReInit(0);

  auto get_r = [&r_min,&r_max](const Real& x) {
    Real angle = sctl::const_pi<Real>() * (16.*x - 28./3.); // =8*(t-pi/6), t = (x-0.5)*2pi
    return r_min + (r_max - r_min) * (0.5 * sctl::sin<Real>(angle) + 0.5);
  };

  auto get_xyz = [](const Real& x) {
    const Real xminus = x-0.5;
    const Real x4pi = 4.*sctl::const_pi<Real>()*xminus;
    const Real x8pi = 2.*x4pi;
    const Real xminus2 = xminus * xminus;
    const Real xminus5 = xminus2 * xminus2 * xminus;
    Real xcoeff = xminus2 * 4. - 1.;
    xcoeff = xcoeff / 5.;
    Real x_ = 0.5 * xminus * sctl::cos<Real>(x4pi) + 8. * xminus5 + 0.5;
    Real y_ = sctl::sin<Real>(x4pi) * xcoeff + 0.5;
    Real z_ = sctl::sin<Real>(x8pi) * xcoeff + 0.5;

    return std::make_tuple(x_,y_,z_);
  };

  const Real ptcl_r = 0.15 * r_max; // radius/size of each particle.
  // hexagonal close packing
  const Real dx = 0.5 * r_max; // buffer room between centers of spheres = 1.2*diameter
  const Real dr = 1.2 * (2. * ptcl_r); // buffer in r direction also 
  const Real dtheta = sctl::const_pi<Real>() / 3.; // 2pi/6 so 6 spheres in theta direction

  for (Real x=dx; x < 1.-dx; x+=dx) {
    Real r_channel = get_r(x);
    for (Real r = dr; r < r_channel-1.2*ptcl_r; r += dr) {
      for (Real theta = dtheta; theta < 2*sctl::const_pi<Real>()-dtheta; theta += dtheta) {
        Real Y = r * sctl::cos<Real>(theta);
        Real Z = r * sctl::sin<Real>(theta);
        Real ptcl_r_loc = (drand48() * 0.4 + 0.5) * ptcl_r; // shrink ratio between 0.5 and 0.9
        std::tuple<Real,Real,Real> xyz_j = get_xyz(x);
        if (ptcl_r_loc >= 1e-5) {
          ptcls_Xcs.PushBack(std::get<0>(xyz_j));
          ptcls_Xcs.PushBack(Y + std::get<1>(xyz_j));
          ptcls_Xcs.PushBack(Z + std::get<2>(xyz_j));
          ptcls_rs.PushBack(ptcl_r_loc);
        }
      }
    }
  }
}

template <class Real> void PeriodicGeom<Real>::loop_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real major_r, const Real minor_r){
  x = major_r * sctl::cos<Real>(theta);
  y = major_r * sctl::sin<Real>(theta);
  z = 0;
  ex = 0;
  ey = 0;
  ez = 1;
  r = minor_r;
};

template <class Real> void PeriodicGeom<Real>::sphere_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad){
  x = loop_rad*sctl::cos<Real>(theta);
  y = 0;
  z = 0;
  r = loop_rad*sctl::sin<Real>(theta);
  ex = 0;
  ey = 1;
  ez = 0;
}

// DEPRECATED
template <class Real> void PeriodicGeom<Real>::spheroid_geom(Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta, const Real loop_rad){
  Real u0 = 1.1;
  x = loop_rad * u0 * sctl::cos<Real>(theta);
  y = 0.;
  z = 0.;
  r = loop_rad * sctl::sqrt<Real>(u0*u0 - 1) * sctl::sin<Real>(theta);
  ex = 0;
  ey = 0;
  ez = 1;
}


// Use FMM to approximate field generated by infinite copies of Xsrc's in <peri_mode>-periodic geometry; approximation done by adding <Ncopy> many source (with strength <sigma>) contributions using FMM.
template <class Real> sctl::Vector<Real> PeriodicGeom<Real>::exact_field_fmm(const sctl::Vector<Real>& Xtrg, const sctl::Vector<Real>& Xsrc, const sctl::Vector<Real>& sigma, const sctl::Long Ncopy, const sctl::Integer peri_mode) {
    sctl::Stokes3D_FxU ker;
    Real tol_ = 1e-14; // arbitrary multipole accuracy for fmm
  
    const sctl::Long N = Xtrg.Dim()/3; // Number of targets.
    sctl::Vector<Real> U(N*3); // output vector.
    U.SetZero();
    // Setup FMM
    fmm.SetKernels(ker, ker, ker);
    fmm.AddSrc("Src", ker, ker);
    fmm.AddTrg("Trg", ker, ker);
    fmm.SetKernelS2T("Src", "Trg", ker);
    fmm.SetAccuracy((sctl::Integer)(sctl::log(tol_)/sctl::log(0.1))+1);
    fmm.SetTrgCoord("Trg", Xtrg);
    
    // Vector of all sources for approximation
    PeriodicGeom<Real> obj;
    sctl::Vector<Real> Xsrc_, sigma_;
    if (peri_mode==3) {
      // 3-Peri requires too large memory for all copies to be loaded at the same time.
      Xsrc_ = obj.X_nbr_copy(Xsrc,Ncopy,2);
      sigma_ = obj.vec_nbr_copy(sigma,Ncopy,2);
      sctl::Long Nsrc2D = Xsrc_.Dim()/3;
      sctl::Vector<Real> Zshift_(Xsrc.Dim());
      Zshift_.SetZero();
      for (int i=0; i<Xsrc.Dim()/3; i++) {
        Zshift_[i*3+2] = 1.;
      }
      // std::cout << "size of Zshift_ before copies: " << Zshift_.Dim() << std::endl;
      sctl::Vector<Real> Zshift = obj.vec_nbr_copy(Zshift_,Ncopy,2);
      sctl::Vector<Real> Uloc(Xtrg.Dim());
      sctl::Vector<Real> Xsrc_plane;
      fmm.SetSrcDensity("Src", sigma_);
      for (int k3=-Ncopy; k3<Ncopy; k3++) {
        std::cout << "K3 = " << k3 << std::endl;
        Uloc.SetZero(); // reset U
        // std::cout << "size of Xsrc: " << Xsrc_.Dim() << ", size of Zshift: " << Zshift.Dim() << std::endl;
        Xsrc_plane = Xsrc_ + k3 * Zshift; // Shift 2D grid by k3*(0,0,1) for each point.
        // std::cout << "Check shift, Xsrc_[0] was " << Xsrc_[0] << ", " << Xsrc_[1] << ", " << Xsrc_[2] << ", is (after shift) " << Xsrc_plane[0] << ", " << Xsrc_plane[1] << ", " << Xsrc_plane[2] << ". " << std::endl;
        fmm.SetSrcCoord("Src", Xsrc_plane);
        fmm.Eval(Uloc, "Trg");
        U += Uloc;
      }
    } else {
      std::cout << "peri mode < 3" << std::endl;
      Xsrc_ = obj.X_nbr_copy(Xsrc,Ncopy,peri_mode);
      sigma_ = obj.vec_nbr_copy(sigma,Ncopy,peri_mode);
      fmm.SetSrcDensity("Src", sigma_);
      fmm.SetSrcCoord("Src", Xsrc_); 
      fmm.Eval(U, "Trg");
    }
    return U;
}