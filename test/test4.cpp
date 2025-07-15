#include "periodize.hpp"
#include "utils.hpp"

/**
 * Background flow with unit pressure gradient along X-axis.
 */
template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
  const sctl::Long N = X.Dim()/3;
  sctl::Vector<Real> U(N*3);
  for (sctl::Long i = 0; i < N; i++) {
    const auto x = X.begin() + i*3;
    U[i*3+0] = -((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
    U[i*3+1] = 0;
    U[i*3+2] = 0;
  }
  return U;
}

// // Background flow that is constant in only x direction.
// template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
//   const sctl::Long N = X.Dim()/3;
//   sctl::Vector<Real> U(3*N);
//   for (sctl::Long i = 0; i < N; i++) {
//     U[i*3+0] = 1.;
//     U[i*3+1] = 0.;
//     U[i*3+2] = 0.;
//   }
//   return U;
// }

template <class Real> void test(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Integer peri_mode, sctl::Comm comm) {

  // Combine single-layer and double-layer kernels in these proportions
  const Real SL_scal = 1.0;
  const Real DL_scal = 1.0;

  const Real tol = 1e-15;
  const Real gmres_tol = 1e-10;
  const sctl::Long ElemOrder = 10;

  PeriodicGeom<Real> obj;
  sctl::Vector<sctl::Long> ptcls;
  sctl::Vector<Real> ptcls_Xcs;
  sctl::Vector<Real> ptcls_rs;
  sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
  sctl::Vector<Real> NormalOrient;
  std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 0, 1, comm, 25, ptcls, ptcls_rs, ptcls_Xcs, 0);
  elem_lst0 = std::get<0>(build0);
  std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 1, peri_mode, 25, comm, ptcls, ptcls_rs, ptcls_Xcs, 0);
  elem_lst_nbr = std::get<0>(build_nbr);
  NormalOrient = std::get<1>(build_nbr);
  const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); 
  std::cout << "periodic mode is " << peri_mode << ", Nrepeat is " << Nrepeat << std::endl;

  sctl::Vector<Real> X0; // target coordinates
  elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
//   std::cout << "size of X0: " <<  X0.Dim() << std::endl;
  sctl::Vector<Real> X_proxy;
  if (peri_mode == 1) {
    X_proxy = Periodize1D<Real>::GetProxySurf(); // proxy points coordinates
  } else if (peri_mode == 3) {
    X_proxy = Periodize3D<Real>::GetProxySurf(); // proxy points coordinates
  } else {
    SCTL_ASSERT(false);
  }
  
  std::string nbr_vis = "vis/Snbr_dense_"+std::to_string(peri_mode)+"_periodic";
  sctl::Vector<Real> Xnbr,Xnbrn;
  elem_lst_nbr.GetNodeCoord(&Xnbr, &Xnbrn, nullptr);
  elem_lst_nbr.WriteVTK(nbr_vis,Xnbr,comm); // visualization with particle inside.
//   std::cout << "after Sptcl" << std::endl;

  StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
  LayerPotenOp0.AddElemList(elem_lst_nbr);
  LayerPotenOp0.SetTargetCoord(X0);
  LayerPotenOp0.SetAccuracy(tol);

  StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal, comm); // potential from elem_lst0 to proxy points
  LayerPotenOp_proxy.AddElemList(elem_lst0);
  LayerPotenOp_proxy.SetTargetCoord(X_proxy);
  LayerPotenOp_proxy.SetAccuracy(tol);

//   std::cout << "after making LPO" << std::endl;

  // periodized layer potential operator
  const auto BIO = [&DL_scal,&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat,NormalOrient,&peri_mode, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    const sctl::Long N = sigma.Dim();
    // std::cout << "in BIO, dim of sigma is " << N << std::endl;

    sctl::Vector<Real> sigma_nbr(Nrepeat*N); // repeat sigma Nrepeat times
    for (sctl::Long k = 0; k < Nrepeat; k++) {
      for (sctl::Long i = 0; i < N; i++) {
        sigma_nbr[k*N+i] = sigma[i];
      }
    }

    U->SetZero();
    LayerPotenOp0.ComputePotential(*U, sigma_nbr);
    // if (DL_scal && U->Dim() == N) (*U) -= sigma*0.5*NormalOrient * DL_scal; // for double-layer
    if (DL_scal && U->Dim() == N) {
    //   std::cout << "self to self, dim of NormalOrient is " << NormalOrient.Dim() << std::endl;
      (*U) -= sigma*0.5*NormalOrient * DL_scal;
    }

    { // Add far-field
      sctl::Vector<Real> U_proxy, U_far;
      LayerPotenOp_proxy.ComputePotential(U_proxy, sigma);
      // std::cout << "right before eval far" << std::endl;
      if (peri_mode==1) {
        // 1-periodic
        Periodize1D<Real>::EvalFarField(U_far, X0, U_proxy);
      } else if (peri_mode==3) {
        // 3-periodic
        Periodize3D<Real>::EvalFarField(U_far, X0, U_proxy);
      } else {
        std::cout << "2-periodic not yet implemented." << std::endl;
        SCTL_ASSERT(false);
      }
      
      (*U) += U_far;
    } 
    // comm.Barrier();
  };

  // Solve for sigma to satisfy no-slip boundary conditions: BIO(sigma) + bg_flow = 0
  sctl::Vector<Real> sigma;
  sctl::GMRES<Real> solver(comm);
  sctl::Vector<Real> U0 = -bg_flow(X0);
  // std::cout << "U0:" << std::endl;
  // for (int i=0; i<U0.Dim(); i++) {
  //   std::cout << U0[i] << std::endl;
  // }
  solver(&sigma, BIO, U0, gmres_tol);

  { // Evaluate in interior, and write visualization
    std::cout << "Rank " << comm.Rank()<< " calculating target points." << std::endl;
    PeriodicGeom<Real> trg;
    // const sctl::Long Nelem_trg = 4;
    // const sctl::Long FourierOrder_trg = 16;
    // sctl::SlenderElemList<Real> elem_lst_trg;
    // sctl::Vector<sctl::Long> ptcls_trg;
    // sctl::Vector<Real> ptcls_Xcs_trg;
    // sctl::Vector<Real> ptcls_rs_trg;
    // std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_straight(Nelem_trg, ElemOrder, FourierOrder_trg, 0, peri_mode, 0.9/2., comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, 0);
    // elem_lst_trg = std::get<0>(build_trg);
    
    CubeVolumeVisShifted<Real> vol_vis(50, 1.0, comm);
    // VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
    // X0 = vol_vis.GetCoord();
    sctl::Vector<Real> X0_all = vol_vis.GetCoord();
    // std::cout << "size of X0 all is " << X0_all.Dim();
    sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
    X0 = std::get<0>(trg_tuple);
    filtered_inds = std::get<1>(trg_tuple);
    // std::cout << "size of X0 is " << X0.Dim() << std::endl;

    LayerPotenOp0.SetTargetCoord(X0);
    sctl::Vector<Real> U;
    BIO(&U, sigma);
    U += bg_flow(X0);
    // sctl::Vector<Real> U_vis = U;
    sctl::Vector<Real> U_vis(X0_all.Dim());
    U_vis = 0.;
    sctl::Long X1_ptr = 0;
    for (sctl::Long i=0; i<X0_all.Dim()/3; i++) {
      if (filtered_inds[i] == 0) {
        U_vis[i*3] = U[X1_ptr*3];
        U_vis[i*3+1] = U[X1_ptr*3+1];
        U_vis[i*3+2] = U[X1_ptr*3+2];
        X1_ptr += 1;
      }
    }
    // sctl::Vector<Real> bgU = bg_flow(X0);
    // vol_vis.WriteVTK("vis/backgroundU", bgU);

    std::string filename_vis = "vis/U_dense_"+std::to_string(peri_mode)+"_periodic";
    vol_vis.WriteVTK(filename_vis, U_vis);

  }
}


int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    //sctl::Profile::Enable(true);
    sctl::Comm comm = sctl::Comm::World();
    long Nelem_ptcl = std::stol(argv[1]); // number of elements
    long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
    int peri_mode = std::stoi(argv[3]); // what kind of periodicity does the system have; peri_mode = j for j-periodic.
    test<Real>(Nelem_ptcl, FourierOrder, peri_mode, comm);
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}
