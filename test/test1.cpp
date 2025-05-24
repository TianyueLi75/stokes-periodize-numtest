// export OMP_NUM_THREADS=16; time make DEBUG=0 -B bin/test1 && time mpirun -n 1 --map-by slot:pe=$OMP_NUM_THREADS ./bin/test1

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

template <class Real> void test(sctl::Comm comm) {
  // std::cout << "size of overall comm is " << comm.Size() <<"; rank in overall comm is "<< comm.Rank() << std::endl;

  // Combine single-layer and double-layer kernels in these proportions
  const Real SL_scal = 1.0;
  const Real DL_scal = 1.0;

  const Real tol = 1e-14;
  const Real gmres_tol = 1e-8;
  // const sctl::Long gmres_max_iter = 50;
  const sctl::Long Nelem_channel = 4;
  const sctl::Long ElemOrder = 10;
  const sctl::Long FourierOrder = 16;

  PeriodicGeom<Real> obj;
  sctl::Vector<sctl::Long> ptcls; // no length initialization -- no particles.
  // sctl::Vector<sctl::Long> ptcls(1);
  // ptcls = 4;
  int geom_mode = 3;
  sctl::Vector<Real> ptcls_Xcs;
  sctl::Vector<Real> ptcls_rs;
  const auto elem_lst0 = obj.build_sinusoidal(Nelem_channel, ElemOrder, FourierOrder, 0, 0.2, 0.1, ptcls, ptcls_rs, ptcls_Xcs,geom_mode);
  const auto elem_lst_nbr = obj.build_sinusoidal(Nelem_channel, ElemOrder, FourierOrder, 1, 0.2, 0.1, ptcls, ptcls_rs, ptcls_Xcs,geom_mode);
  const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); // should be 3
  // std::cout << "just after build sinusoidal; ptcl Xc dim = " << ptcls_Xcs.Dim() << "; ptcl rs dim = " << ptcls_rs.Dim() << std::endl;

  sctl::Vector<Real> X0; // target coordinates
  elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
  const auto X_proxy = Periodize<Real>::GetProxySurf(); // proxy points coordinates

  sctl::Vector<Real> NormalOrient; // normal orientation (interior=1, exterior=-1)
  { // set NormalOrient
    constexpr sctl::Integer COORD_DIM = 3;
    sctl::Vector<sctl::Long> elem_wise_node_cnt;
    elem_lst0.GetNodeCoord(nullptr, nullptr, &elem_wise_node_cnt);
    for (sctl::Long i = 0; i < elem_wise_node_cnt.Dim(); i++) {
      for (sctl::Long j = 0; j < elem_wise_node_cnt[i]*COORD_DIM; j++) {
        NormalOrient.PushBack(i < Nelem_channel ? 1 : -1);
      }
    }
  }

  // elem_lst0.WriteVTK("vis/X_utils", X0, comm);

  StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
  LayerPotenOp0.AddElemList(elem_lst_nbr);
  LayerPotenOp0.SetTargetCoord(X0);
  LayerPotenOp0.SetAccuracy(tol);

  StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal, comm); // potential from elem_lst0 to proxy points
  LayerPotenOp_proxy.AddElemList(elem_lst0);
  LayerPotenOp_proxy.SetTargetCoord(X_proxy);
  LayerPotenOp_proxy.SetAccuracy(tol);

  // periodized layer potential operator
  const auto BIO = [&DL_scal,&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat,NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    const sctl::Long N = sigma.Dim();

    sctl::Vector<Real> sigma_nbr(Nrepeat*N); // repeat sigma Nrepeat times
    for (sctl::Long k = 0; k < Nrepeat; k++) {
      for (sctl::Long i = 0; i < N; i++) {
        sigma_nbr[k*N+i] = sigma[i];
      }
    }

    U->SetZero();
    LayerPotenOp0.ComputePotential(*U, sigma_nbr);
    if (DL_scal && U->Dim() == N) (*U) -= sigma*0.5*NormalOrient * DL_scal; // for double-layer

    { // Add far-field
      sctl::Vector<Real> U_proxy, U_far;
      LayerPotenOp_proxy.ComputePotential(U_proxy, sigma);
      Periodize<Real>::EvalFarField(U_far, X0, U_proxy);
      (*U) += U_far;
    }
  };

  // Solve for sigma to satisfy no-slip boundary conditions: BIO(sigma) + bg_flow = 0
  sctl::Vector<Real> sigma;
  sctl::GMRES<Real> solver(comm);
  // solver(&sigma, BIO, -bg_flow(X0), gmres_tol, gmres_max_iter);
  solver(&sigma, BIO, -bg_flow(X0), gmres_tol);
  // elem_lst0.WriteVTK("vis/sigma", sigma, comm);

  { // Evaluate in interior, and write visualization
    // Note: BIO defined to take in X0 as target.
    VolumeVis<Real> vol_vis(elem_lst0, comm, false);
    sctl::Vector<Real> X0_all = vol_vis.GetCoord(); // set new target coordinates
    sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
    if (!ptcls.Dim()) {
      std::cout << "no particles." << std::endl;
      X0 = X0_all;
      filtered_inds = 1;
    } else {
      std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = obj.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs,geom_mode);
      X0 = std::get<0>(trg_tuple);
      filtered_inds = std::get<1>(trg_tuple);
    }
    LayerPotenOp0.SetTargetCoord(X0);
    sctl::Vector<Real> U;
    BIO(&U, sigma);
    U += bg_flow(X0);
    sctl::Vector<Real> U_vis(X0_all.Dim());
    if (!ptcls.Dim()) {
      std::cout << "no particles." << std::endl;
      U_vis = U;
    } else {
      U_vis = 0.;
      sctl::Long X1_ptr = 0;
      for (sctl::Long i=0; i<X0_all.Dim()/3; i++) {
        // std::cout << "target id " << i << std::endl; 
        if (filtered_inds[i] == 0) {
          // std::cout << "not inside, inputting " << X1_ptr << "th target" << std::endl;
          U_vis[i*3] = U[X1_ptr*3];
          U_vis[i*3+1] = U[X1_ptr*3+1];
          U_vis[i*3+2] = U[X1_ptr*3+2];
          X1_ptr += 1;
        }
      }
    }
    vol_vis.WriteVTK("vis/U_utils", U_vis);
  }
}

int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    //sctl::Profile::Enable(true);
    sctl::Comm comm = sctl::Comm::World();
    // std::cout << "size of comm is "<< comm.Size() << std::endl;
    test<Real>(comm);
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}

