#include <tuple>

template <class Real> class PeriodizeOp {
  static constexpr Real tol = sctl::machine_eps<Real>()*64; // tolerance for pseudo-inverse
  static constexpr sctl::Integer COORD_DIM = 3;
  static constexpr sctl::Long m0 = 20; // multipole order
  static constexpr sctl::Long level = 30; // levels of tree expansion and evaluation

  using KerM2M = sctl::Stokes3D_FxU;
  using KerM2L = sctl::Stokes3D_FxU;
  using KerL2L = sctl::Stokes3D_FxU;

  public:

    static constexpr sctl::Long Nsurf() {
      return 6*(m0-1)*(m0-1)+2;
    }
    static sctl::Vector<Real> uc_surf(const Real box_length=1, const sctl::Vector<Real>& X0 = sctl::Vector<Real>{0,0,0}) {
      return proxy_surf(box_length*2.95, X0);
    }
    static sctl::Vector<Real> ue_surf(const Real box_length=1, const sctl::Vector<Real>& X0 = sctl::Vector<Real>{0,0,0}) {
      return proxy_surf(box_length*1.05, X0);
    }
    static sctl::Vector<Real> dc_surf(const Real box_length=1, const sctl::Vector<Real>& X0 = sctl::Vector<Real>{0,0,0}) {
      return proxy_surf(box_length*1.05, X0);
    }
    static sctl::Vector<Real> de_surf(const Real box_length=1, const sctl::Vector<Real>& X0 = sctl::Vector<Real>{0,0,0}) {
      return proxy_surf(box_length*2.95, X0);
    }

    static std::tuple<sctl::Matrix<Real>, sctl::Matrix<Real>> UC2UE(const Real box_length=1) {
      sctl::Profile::Scoped prof(__FUNCTION__);
      const auto Xc = uc_surf(box_length);
      const auto Xe = ue_surf(box_length);
      const KerM2M ker_m2m;

      sctl::Matrix<Real> Me2c, U,S,Vt, Mc2e0, Mc2e1;
      ker_m2m.KernelMatrix<Real,true>(Me2c, Xc, Xe, sctl::Vector<Real>());
      sctl::Matrix<Real>(Me2c).SVD(U,S,Vt);

      Real max_val = 0;
      for (sctl::Long i = 0; i < std::min(S.Dim(0),S.Dim(1)); i++) max_val = std::max<Real>(max_val, sctl::fabs(S[i][i]));
      for (sctl::Long i = 0; i < std::min(S.Dim(0),S.Dim(1)); i++) S[i][i] = (sctl::fabs(S[i][i]) < max_val*tol ? 0 : 1/S[i][i]);

      Mc2e0 = Vt.Transpose();
      Mc2e1 = S * U.Transpose();
      return std::make_tuple(Mc2e0, Mc2e1);
    }

    static std::tuple<sctl::Matrix<Real>, sctl::Matrix<Real>> DC2DE(const Real box_length=1) {
      sctl::Profile::Scoped prof(__FUNCTION__);
      const auto Xc = dc_surf(box_length);
      const auto Xe = de_surf(box_length);
      const KerL2L ker_l2l;

      sctl::Matrix<Real> Me2c, U,S,Vt, Mc2e0, Mc2e1;
      ker_l2l.KernelMatrix<Real,true>(Me2c, Xc, Xe, sctl::Vector<Real>());
      sctl::Matrix<Real>(Me2c).SVD(U,S,Vt);

      Real max_val = 0;
      for (sctl::Long i = 0; i < std::min(S.Dim(0),S.Dim(1)); i++) max_val = std::max<Real>(max_val, sctl::fabs(S[i][i]));
      for (sctl::Long i = 0; i < std::min(S.Dim(0),S.Dim(1)); i++) S[i][i] = (sctl::fabs(S[i][i]) < max_val*tol ? 0 : 1/S[i][i]);

      Mc2e0 = Vt.Transpose() * S;
      Mc2e1 = U.Transpose();
      return std::make_tuple(Mc2e0, Mc2e1);
    }

    static sctl::Matrix<Real> BC_UE2DC() {
      static const auto M = BC_UE2DC_helper();
      return M;
    }

  private:

    static sctl::Vector<Real> proxy_surf(const Real box_length, const sctl::Vector<Real>& Xc) {
      static const sctl::Vector<Real> X0 = []() {
        sctl::Vector<Real> X;
        for (sctl::Long i0 = 0; i0 < m0; i0++) {
          for (sctl::Long i1 = 0; i1 < m0; i1++) {
            for (sctl::Long i2 = 0; i2 < m0; i2++) {
              if (i0==0 || i0==m0-1 || i1==0 || i1==m0-1 || i2==0 || i2==m0-1) {
                const Real x = i0/(Real)(m0-1);
                const Real y = i1/(Real)(m0-1);
                const Real z = i2/(Real)(m0-1);
                X.PushBack(x-0.5);
                X.PushBack(y-0.5);
                X.PushBack(z-0.5);
              }
            }
          }
        }
        return X;
      }();
      SCTL_ASSERT(X0.Dim() == Nsurf()*COORD_DIM);

      sctl::Vector<Real> X(Nsurf()*COORD_DIM);
      for (sctl::Long i = 0; i < Nsurf(); i++) {
        for (sctl::Long k = 0; k < COORD_DIM; k++) {
          X[i*COORD_DIM+k] = X0[i*COORD_DIM+k]*box_length + Xc[k];
        }
      }
      return X;
    }

    static sctl::Matrix<Real> BC_UE2DC_helper() {
      sctl::Profile::Scoped prof(__FUNCTION__);

      std::string data_file = "data/Mbc_ue2dc_1d_l"+std::to_string(level)+"_m"+std::to_string(m0)+".mat";
      sctl::Matrix<Real> M;
      // M.template Read<sctl::QuadReal>("data/Mbc_ue2dc_1d_l30_m20.mat");
      M.template Read<sctl::QuadReal>(data_file.c_str());
      if (M.Dim(0) || M.Dim(1)) return M;

      const KerM2L ker_m2l;
      const sctl::Integer kdim[2] = {KerM2L::SrcDim(), KerM2L::TrgDim()};
      const auto X0 = dc_surf();

      sctl::Matrix<Real> M2M(Nsurf()*kdim[0], Nsurf()*kdim[0]);
      M.ReInit(Nsurf()*kdim[0], Nsurf()*kdim[1]);
      M2M.SetZero();
      M.SetZero();

      for (sctl::Long i = 0; i < M2M.Dim(0); i++) { // M2M <-- Identity
        M2M[i][i] = 1;
      }

      for (sctl::Long l = 0; l < level; l++) { // tree-code (hierarchical) summation
        std::cout<<"level = "<<l<<'\n';
        const sctl::Long box_length = ((sctl::Long)1) << l;

        sctl::Matrix<Real> M0, M1;
        ker_m2l.KernelMatrix<Real,true>(M0, X0, ue_surf(box_length, sctl::Vector<Real>{ (Real)1.5*box_length+0.5,0,0}), sctl::Vector<Real>());
        ker_m2l.KernelMatrix<Real,true>(M1, X0, ue_surf(box_length, sctl::Vector<Real>{-(Real)1.5*box_length-0.5,0,0}), sctl::Vector<Real>());
        M += M2M * (M0 + M1);

        ker_m2l.KernelMatrix<Real,true>(M0, uc_surf(2*box_length), ue_surf(box_length, sctl::Vector<Real>{ (Real)0.5*box_length,0,0}), sctl::Vector<Real>());
        ker_m2l.KernelMatrix<Real,true>(M1, uc_surf(2*box_length), ue_surf(box_length, sctl::Vector<Real>{-(Real)0.5*box_length,0,0}), sctl::Vector<Real>());
        const auto [M_uc2ue0, M_uc2ue1] = PeriodizeOp<Real>::UC2UE(2*box_length);
        M2M = M2M * (((M0+M1) * M_uc2ue0) * M_uc2ue1);
      }
      //for (sctl::Long i = -1000; i <= 1000; i++) { // direct summation
      //  if (abs(i) >= 2) {
      //    sctl::Matrix<Real> M_;
      //    ker_m2l.KernelMatrix<Real,true>(M_, X0, ue_surf(1, sctl::Vector<Real>{(Real)i,0,0}), sctl::Vector<Real>());
      //    M += M_;
      //  }
      //}

      // M.template Write<sctl::QuadReal>("data/Mbc_ue2dc_1d_l30_m20.mat");
      M.template Write<sctl::QuadReal>(data_file.c_str());
      return M;
    }
};

////////////// Periodize1D /////////////////////////////
template <class Real> const sctl::Vector<Real>& Periodize1D<Real>::GetProxySurf() {
  static const auto X = PeriodizeOp<Real>::uc_surf(1, sctl::Vector<Real>{0.5,0.5,0.5});
  return X;
}

template <class Real> sctl::Vector<Real> Periodize1D<Real>::GetProxySurf(const sctl::Long level_in, const sctl::Long m0_in) {
  // sctl::Vector<Real> proxy_surf = [&level_in,&m0_in](){
  std::string data_dn = "data/dn_equiv_surf_1d_l"+std::to_string(level_in)+"_m"+std::to_string(m0_in)+".mat";
  // std::string data_up = "data/up_check_surf_l"+std::to_string(level_in)+"_m"+std::to_string(m0_in)+".mat";
  sctl::Vector<Real> X;
  X.template Read<PrecompReal>(data_dn.c_str());
  if (!X.Dim()) {
    X = PeriodizeOp<Real>::uc_surf(1, sctl::Vector<Real>{0.5,0.5,0.5});
    X.template Write<PrecompReal>(data_dn.c_str());
    // X.template Write<PrecompReal>(data_up.c_str());
  }
  // std::cout << "size of proxy surf is " << X.Dim() << std::endl;
  return X;
  // }();
  // std::cout << "reading proxy surf or making it" << std::endl;
  // return proxy_surf;
}


template <class Real> void Periodize1D<Real>::EvalFarField(sctl::Vector<Real>& U_far, const sctl::Vector<Real>& Xt, const sctl::Vector<Real>& U_proxy) {
  // std::cout << "in eval far field, getting matrices. " << std::endl;
  const auto& Mbc0 = GetMat_UC2DE0();
  const auto& Mbc1 = GetMat_UC2DE1();
  const sctl::Long N = Mbc0.Dim(0);
  SCTL_ASSERT(U_proxy.Dim() == N);

  // Compute the equivalent density at proxy points
  auto proxy_density = (sctl::Matrix<Real>(1,N,(sctl::Iterator<Real>)U_proxy.begin(),false) * Mbc0) * Mbc1;

  // std::cout << "evaluate potential from proxy to Xt" << std::endl;
  // Evaluate the potential from proxy points at the targets Xt
  U_far = 0;
  static const sctl::Stokes3D_FxU stokeslet;
  stokeslet.template Eval<Real,true>(U_far, Xt, GetProxySurf(), sctl::Vector<Real>(), sctl::Vector<Real>(N,proxy_density.begin(),false));
}

template <class Real> void Periodize1D<Real>::EvalFarField(sctl::Vector<Real>& U_far, const sctl::Vector<Real>& Xt, const sctl::Vector<Real>& U_proxy, const sctl::Long level_in, const sctl::Long m0_in) {
  // std::cout << "in eval far field, getting matrices. " << std::endl;
  const auto Mbc0 = GetMat_UC2DE0(level_in, m0_in);
  const auto Mbc1 = GetMat_UC2DE1(level_in, m0_in);
  const sctl::Long N = Mbc0.Dim(0);
  // std::cout << N << ", " << U_proxy.Dim() << ", l=" << level_in << ", m0=" << m0_in << std::endl;
  SCTL_ASSERT(U_proxy.Dim() == N);

  // Compute the equivalent density at proxy points
  auto proxy_density = (sctl::Matrix<Real>(1,N,(sctl::Iterator<Real>)U_proxy.begin(),false) * Mbc0) * Mbc1;

  // std::cout << "evaluate potential from proxy to Xt" << std::endl;
  // Evaluate the potential from proxy points at the targets Xt
  U_far = 0;
  static const sctl::Stokes3D_FxU stokeslet;
  stokeslet.template Eval<Real,true>(U_far, Xt, GetProxySurf(level_in, m0_in), sctl::Vector<Real>(), sctl::Vector<Real>(N,proxy_density.begin(),false));
}

template <class Real> const sctl::Matrix<Real>& Periodize1D<Real>::GetMat_UC2DE0() {
  static sctl::Matrix<Real> Mbc = [](){
    const auto [M_uc2ue0, M_uc2ue1] = PeriodizeOp<Real>::UC2UE();
    const auto [M_dc2de0, M_dc2de1] = PeriodizeOp<Real>::DC2DE();
    const auto Mbc_ue2dc = PeriodizeOp<Real>::BC_UE2DC();
    return (M_uc2ue0 * (M_uc2ue1 * Mbc_ue2dc)) * M_dc2de0;
  }();
  return Mbc;
}

template <class Real> sctl::Matrix<Real> Periodize1D<Real>::GetMat_UC2DE0(const sctl::Long level_in, const sctl::Long m0_in) {
  // sctl::Matrix<Real> Mbc = [&level_in,&m0_in](){
  std::string data1 = "data/M_uc2ue0_1d_l"+std::to_string(level_in)+"_m"+std::to_string(m0_in)+".mat";
  std::string data2 = "data/M_uc2ue1_1d_l"+std::to_string(level_in)+"_m"+std::to_string(m0_in)+".mat";
  std::string data3 = "data/Mbc_ue2dc_1d_l"+std::to_string(level_in)+"_m"+std::to_string(m0_in)+".mat";
  std::string data4 = "data/M_dc2de0_1d_l"+std::to_string(level_in)+"_m"+std::to_string(m0_in)+".mat";
  sctl::Matrix<Real> Mbc_ue2dc, M_dc2de0, M_uc2ue0, M_uc2ue1;
  M_uc2ue0.template Read<PrecompReal>(data1.c_str());

  if (M_uc2ue0.Dim(0) || M_uc2ue0.Dim(1)) {
    // std::cout << " successfully read file." << std::endl;
    M_uc2ue1.template Read<PrecompReal>(data2.c_str());
    Mbc_ue2dc.template Read<PrecompReal>(data3.c_str());
    M_dc2de0.template Read<PrecompReal>(data4.c_str());
  } else {
    // std::cout << " Didn't read file, making now." << std::endl;
    std::tuple<sctl::Matrix<Real>, sctl::Matrix<Real>> tpl1 = PeriodizeOp<Real>::UC2UE();
    M_uc2ue0 = std::get<0>(tpl1);
    M_uc2ue1 = std::get<1>(tpl1);
    std::tuple<sctl::Matrix<Real>, sctl::Matrix<Real>> tpl2 = PeriodizeOp<Real>::DC2DE();
    M_dc2de0 = std::get<0>(tpl2);
    Mbc_ue2dc = PeriodizeOp<Real>::BC_UE2DC();
    M_uc2ue0.template Write<PrecompReal>(data1.c_str());
    M_uc2ue1.template Write<PrecompReal>(data2.c_str());
    Mbc_ue2dc.template Write<PrecompReal>(data3.c_str());
    M_dc2de0.template Write<PrecompReal>(data4.c_str());
  }
  // std::cout << level_in << ", " << m0_in << std::endl;
  // std::cout << M_uc2ue0.Dim(0)<< ", "  << M_uc2ue1.Dim(0)<< ", "  << Mbc_ue2dc.Dim(0)<< ", "  << M_dc2de0.Dim(0) << std::endl;
  return (M_uc2ue0 * (M_uc2ue1 * Mbc_ue2dc)) * M_dc2de0;
  // }();
  // return Mbc;
}

template <class Real> const sctl::Matrix<Real>& Periodize1D<Real>::GetMat_UC2DE1() {
  static sctl::Matrix<Real> Mbc = [](){
    const auto [M_dc2de0, M_dc2de1] = PeriodizeOp<Real>::DC2DE();
    return M_dc2de1;
  }();
  return Mbc;
}

template <class Real> sctl::Matrix<Real> Periodize1D<Real>::GetMat_UC2DE1(const sctl::Long level_in, const sctl::Long m0_in) {
  // sctl::Matrix<Real> Mbc = [&level_in,&m0_in](){
  std::string data_file = "data/M_dc2de1_1d_l"+std::to_string(level_in)+"_m"+std::to_string(m0_in)+".mat";
  sctl::Matrix<Real> M_dc2de1;
  M_dc2de1.template Read<PrecompReal>(data_file.c_str());

  if (M_dc2de1.Dim(0)==0 && M_dc2de1.Dim(1)==0) {
    std::cout << "Writing data file" << std::endl;
    std::tuple<sctl::Matrix<Real>, sctl::Matrix<Real>> tpl2 = PeriodizeOp<Real>::DC2DE();
    M_dc2de1 = std::get<1>(tpl2);
    M_dc2de1.template Write<PrecompReal>(data_file.c_str());
  }
  return M_dc2de1;
  // }();
  // return Mbc;
}

////////////// Periodize3D /////////////////////////////
template <class Real> const sctl::Vector<Real>& Periodize3D<Real>::GetProxySurf() {
  static const sctl::Vector<Real> proxy_surf = [](){
    sctl::Vector<Real> X;
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
    M_dc2de1.template Read<PrecompReal>("data/M_dc2de1_l30_m20.mat");
    return M_dc2de1;
  }();
  return Mbc;
}