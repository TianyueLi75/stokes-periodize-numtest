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

template <class Real> void test(sctl::Long Nelem_ptcl, sctl::Long FourierOrder, sctl::Integer peri_mode, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode) {

  // Combine single-layer and double-layer kernels in these proportions
  const Real SL_scal = 1.0;
  const Real DL_scal = 1.0;

  const Real tol = 1e-15;
  const Real gmres_tol = 1e-10;
  const sctl::Long ElemOrder = 10;

  PeriodicGeom<Real> obj;
  sctl::Vector<sctl::Long> ptcls;
  if (Nptcl>0) {
    ptcls.ReInit(Nptcl);
    ptcls = 1;
  }
  sctl::Vector<Real> ptcls_Xcs;
  sctl::Vector<Real> ptcls_rs;
  sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
  sctl::Vector<Real> NormalOrient;
  Real box_sidelen = 1.;
  std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_only_ptcls(Nelem_ptcl, ElemOrder, FourierOrder, 0, 1, box_sidelen, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
  elem_lst0 = std::get<0>(build0);
  std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_only_ptcls(Nelem_ptcl, ElemOrder, FourierOrder, 1, peri_mode, box_sidelen, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
  elem_lst_nbr = std::get<0>(build_nbr);
  NormalOrient = std::get<1>(build_nbr);
  const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); // should be 3
  Nptcl = ptcls_rs.Dim(); 
  std::cout << "periodic mode is " << peri_mode << ", Nrepeat is " << Nrepeat << std::endl;

  sctl::Vector<Real> X0; // target coordinates
  elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
  const auto X_proxy = Periodize<Real>::GetProxySurf(); // proxy points coordinates

  sctl::Vector<Real> Xnbr,Xnbrn;
  elem_lst_nbr.GetNodeCoord(&Xnbr, &Xnbrn, nullptr);
  // elem_lst_nbr.WriteVTK("vis/S-nbr-normal",Xnbrn,comm);
  elem_lst_nbr.WriteVTK("vis/S-ptcl-only",Xnbr,comm); // visualization with particle inside.

  StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
  LayerPotenOp0.AddElemList(elem_lst_nbr);
  LayerPotenOp0.SetTargetCoord(X0);
  LayerPotenOp0.SetAccuracy(tol);

  // TODO: if 3-periodic, initiate using periodize-example code.
  StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal, comm); // potential from elem_lst0 to proxy points
  LayerPotenOp_proxy.AddElemList(elem_lst0);
  LayerPotenOp_proxy.SetTargetCoord(X_proxy);
  LayerPotenOp_proxy.SetAccuracy(tol);

  // periodized layer potential operator
  const auto BIO = [&DL_scal,&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat,NormalOrient, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    const sctl::Long N = sigma.Dim();
    // std::cout << "in BIO, dim of sigma is " << N << std::endl;

    sctl::Vector<Real> sigma_nbr(Nrepeat*N); // repeat sigma Nrepeat times
    for (sctl::Long k = 0; k < Nrepeat; k++) {
      for (sctl::Long i = 0; i < N; i++) {
        sigma_nbr[k*N+i] = sigma[i];
      }
    }

    // std::cout << " before U set zero" << std::endl;
    U->SetZero();
    // std::cout << "after U set zero" << std::endl;
    // std::cout << "wavy 2 debug: size of LayerPotenOp is "<< LayerPotenOp0.Dim(1) << ", " << LayerPotenOp0.Dim(0) << std::endl;
    LayerPotenOp0.ComputePotential(*U, sigma_nbr);
    // if (DL_scal && U->Dim() == N) (*U) -= sigma*0.5*NormalOrient * DL_scal; // for double-layer
    if (DL_scal && U->Dim() == N) {
      // std::cout << "self to self, dim of NormalOrient is " << NormalOrient.Dim() << std::endl;
      (*U) -= sigma*0.5*NormalOrient * DL_scal;
    }

    { // Add far-field
      sctl::Vector<Real> U_proxy, U_far;
      LayerPotenOp_proxy.ComputePotential(U_proxy, sigma);
      Periodize<Real>::EvalFarField(U_far, X0, U_proxy);
      (*U) += U_far;
    } 
    // comm.Barrier();
  };

  // Solve for sigma to satisfy no-slip boundary conditions: BIO(sigma) + bg_flow = 0
  sctl::Vector<Real> sigma;
  sctl::GMRES<Real> solver(comm);
  solver(&sigma, BIO, -bg_flow(X0), gmres_tol);

  { // Evaluate in interior, and write visualization
    std::cout << "Rank " << comm.Rank()<< " calculating target points." << std::endl;
    
    VolumeVis<Real> vol_vis(elem_lst0, comm); 
    sctl::Vector<Real> X0 = vol_vis.GetCoord();
    std::cout << "size of X0 is " << X0.Dim();
    LayerPotenOp0.SetTargetCoord(X0);
    sctl::Vector<Real> U,U2;
    BIO(&U, sigma);
    std::cout << "size of U before bgflow: " << U.Dim() << std::endl;
    U2 = bg_flow(X0);
    std::cout << "size of bgflow(X0): " << U2.Dim() << std::endl;
    U += bg_flow(X0);

    std::string filename_vis = "vis/U_ptcl_only_"+std::to_string(peri_mode)+"_periodic";
    vol_vis.WriteVTK(filename_vis, U);

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
    long Nptcl = std::stol(argv[4]); // number of particles inside
    long geom_mode = std::stol(argv[5]); // =0: spheres; =1: spheroids; =3: bacteria; =4: loop.
    test<Real>(Nelem_ptcl, FourierOrder, peri_mode, comm, Nptcl, geom_mode);
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}
