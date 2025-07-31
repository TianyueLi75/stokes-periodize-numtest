template <class Real> const sctl::Vector<Real>& Periodize3D<Real>::GetProxySurf() {
  static const sctl::Vector<Real> proxy_surf = [](){
    sctl::Vector<Real> X;
    std::string data_file = "data/Mbc_ue2dc_1d_l"+std::to_string(level)+"_m"+std::to_string(m0)+".mat";
    X.template Read<PrecompReal>("data/dn_equiv_surf_l30_m20.mat");
    return X;
  }();
  return proxy_surf;
}

template <class Real> void Periodize3D<Real>::EvalFarField(sctl::Vector<Real>& U_far, const sctl::Vector<Real>& Xt, const sctl::Vector<Real>& U_proxy) {
  const auto& Mbc0 = GetMat_UC2DE0();
  const auto& Mbc1 = GetMat_UC2DE1();
  const sctl::Long N = Mbc0.Dim(0);
  SCTL_ASSERT(U_proxy.Dim() == N);

  // Compute the equivalent density at proxy points
  auto proxy_density = (sctl::Matrix<Real>(1,N,(sctl::Iterator<Real>)U_proxy.begin(),false) * Mbc0) * Mbc1;

  // Evaluate the potential from proxy points at the targets Xt
  U_far = 0;
  static const sctl::Stokes3D_FxU stokeslet;
  stokeslet.template Eval<Real,true>(U_far, Xt, GetProxySurf(), sctl::Vector<Real>(), sctl::Vector<Real>(N,proxy_density.begin(),false));
}

template <class Real> const sctl::Matrix<Real>& Periodize3D<Real>::GetMat_UC2DE0() {
  static sctl::Matrix<Real> Mbc = [](){
    std::string data1 = "data/M_uc2ue0_l"+std::to_string(level)+"_m"+std::to_string(m0)+".mat";
    std::string data2 = "data/M_uc2ue1_l"+std::to_string(level)+"_m"+std::to_string(m0)+".mat";
    std::string data3 = "data/Mbc_ue2dc_l"+std::to_string(level)+"_m"+std::to_string(m0)+".mat";
    std::string data = "data/M_dc2de0_l"+std::to_string(level)+"_m"+std::to_string(m0)+".mat";
    sctl::Matrix<Real> Mbc_ue2dc, M_dc2de0, M_uc2ue0, M_uc2ue1;
    M_uc2ue0.template Read<PrecompReal>("data/M_uc2ue0_l30_m20.mat");
    M_uc2ue1.template Read<PrecompReal>("data/M_uc2ue1_l30_m20.mat");
    Mbc_ue2dc.template Read<PrecompReal>("data/Mbc_ue2dc_l30_m20.mat");
    M_dc2de0.template Read<PrecompReal>("data/M_dc2de0_l30_m20.mat");
    return (M_uc2ue0 * (M_uc2ue1 * Mbc_ue2dc)) * M_dc2de0;
  }();
  return Mbc;
}

template <class Real> const sctl::Matrix<Real>& Periodize3D<Real>::GetMat_UC2DE1() {
  static sctl::Matrix<Real> Mbc = [](){
    sctl::Matrix<Real> M_dc2de1;
    std::string data_file = "data/M_dc2de1_l"+std::to_string(level)+"_m"+std::to_string(m0)+".mat";
    M_dc2de1.template Read<PrecompReal>(data_file.c_str());
    return M_dc2de1;
  }();
  return Mbc;
}