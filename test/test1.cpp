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

template <class Real> void test(sctl::Long Nelem_channel, sctl::Long FourierOrder, bool write_ref, sctl::Comm comm, sctl::Long Nptcl, sctl::Long channel_mode, sctl::Long geom_mode) {

  // Combine single-layer and double-layer kernels in these proportions
  const Real SL_scal = 1.0;
  const Real DL_scal = 1.0;

  const Real tol = 1e-15;
  const Real gmres_tol = 1e-10;
  const sctl::Long ElemOrder = 10;

  // std::cout << "Number of threads before is "<< omp_get_num_threads() << std::endl;
  // omp_set_num_threads(16);
  // std::cout << "Number of threads after is "<< omp_get_num_threads() << std::endl;

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
  // peri_mode = 1 for 1-periodic.
  // TODO: remove sinusoidal with mag 0.3.
  if (channel_mode == 0) {// straight channel
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_straight(Nelem_channel, ElemOrder, FourierOrder, 0, 3, 0.25, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst0 = std::get<0>(build0);
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_straight(Nelem_channel, ElemOrder, FourierOrder, 1, 3, 0.25, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst_nbr = std::get<0>(build_nbr);
    NormalOrient = std::get<1>(build_nbr); // only needed for self eval which uses nbr object.
  } else if (channel_mode == 1) {// sinusoidal with mag = 0.1
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_sinusoidal(Nelem_channel, ElemOrder, FourierOrder, 0, 1, 0.2, 0.1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst0 = std::get<0>(build0);
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_sinusoidal(Nelem_channel, ElemOrder, FourierOrder, 1, 1, 0.2, 0.1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst_nbr = std::get<0>(build_nbr);
    NormalOrient = std::get<1>(build_nbr);
  } else if (channel_mode == 2) {// sinusoidal with mag = 0.3
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_sinusoidal(Nelem_channel, ElemOrder, FourierOrder, 0, 1, 0.1, 0.25, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst0 = std::get<0>(build0);
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_sinusoidal(Nelem_channel, ElemOrder, FourierOrder, 1, 1, 0.1, 0.25, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);  
    elem_lst_nbr = std::get<0>(build_nbr);
    NormalOrient = std::get<1>(build_nbr);
  } else if (channel_mode == 3) {
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_conv_div(Nelem_channel, ElemOrder, FourierOrder, 0, 1, 0.1, 0.1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst0 = std::get<0>(build0);
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_conv_div(Nelem_channel, ElemOrder, FourierOrder, 1, 1, 0.1, 0.1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);  
    elem_lst_nbr = std::get<0>(build_nbr);
    NormalOrient = std::get<1>(build_nbr);
  } else if (channel_mode == 4) {
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_spiral(Nelem_channel, ElemOrder, FourierOrder, 0, 1, 0.2, 0.05, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst0 = std::get<0>(build0);
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_spiral(Nelem_channel, ElemOrder, FourierOrder, 1, 1, 0.2, 0.05, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);  
    elem_lst_nbr = std::get<0>(build_nbr);
    NormalOrient = std::get<1>(build_nbr);
  } else {
    SCTL_ASSERT(false); // not implemented
  }
  // std::cout << "Size of elem_lst_nbr is " << elem_lst_nbr.Size() << ", Size of elem_lst0 is " << elem_lst0.Size() <<std::endl;
  const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); // should be 3
  Nptcl = ptcls_rs.Dim(); // Number of particles could have changed after initializing.

  sctl::Vector<Real> X0; // target coordinates
  elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
  const auto X_proxy = Periodize1D<Real>::GetProxySurf(); // proxy points coordinates

  // for (long i=0; i < ptcls_Xcs.Dim()/3; i++) {
  //   std::cout << "rank " << comm.Rank() << " stores ptcl at " << ptcls_Xcs[i*3] << ", " << ptcls_Xcs[i*3+1] << ", " << ptcls_Xcs[i*3+2] <<std::endl;
  // }

  // sctl::Vector<Real> Xnbr,Xnbrn;
  // elem_lst_nbr.GetNodeCoord(&Xnbr, Xnbrn, nullptr);
  // // elem_lst_nbr.WriteVTK("vis/S-nbr-normal",Xnbrn,comm);
  elem_lst0.WriteVTK("vis/S-ptcl",X0,comm); // visualization with particle inside.

  StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
  LayerPotenOp0.AddElemList(elem_lst_nbr);
  LayerPotenOp0.SetTargetCoord(X0);
  LayerPotenOp0.SetAccuracy(tol);

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
      // std::cout << "in far eval" << std::endl;
      sctl::Vector<Real> U_proxy, U_far;
      LayerPotenOp_proxy.ComputePotential(U_proxy, sigma);
      Periodize1D<Real>::EvalFarField(U_far, X0, U_proxy);
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
    //TODO: should split target up by comm.Rank() not element, since Vol_vis should be shortened..
    PeriodicGeom<Real> trg;
    const sctl::Long Nelem_trg = 4;
    const sctl::Long FourierOrder_trg = 16;
    sctl::SlenderElemList<Real> elem_lst_trg;
    sctl::Vector<sctl::Long> ptcls_trg;
    sctl::Vector<Real> ptcls_Xcs_trg;
    sctl::Vector<Real> ptcls_rs_trg;
    // peri_mode = 1 for 1-periodic
    if (channel_mode == 0) {
      std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_straight(Nelem_trg, ElemOrder, FourierOrder_trg, 0, 1, 0.25, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, geom_mode);
      elem_lst_trg = std::get<0>(build_trg);
    } else if (channel_mode == 1) {
      std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_sinusoidal(Nelem_trg, ElemOrder, FourierOrder_trg, 0, 1, 0.2, 0.1, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, geom_mode);
      elem_lst_trg = std::get<0>(build_trg);
    } else if (channel_mode == 2) {
      std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_sinusoidal(Nelem_trg, ElemOrder, FourierOrder_trg, 0, 1, 0.1, 0.25, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, geom_mode);
      elem_lst_trg = std::get<0>(build_trg);
    } else if (channel_mode == 3) {
      std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_conv_div(Nelem_channel, ElemOrder, FourierOrder, 0, 1, 0.1, 0.1, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, geom_mode);
      elem_lst_trg = std::get<0>(build_trg);
    } else if (channel_mode == 4) {
      std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_spiral(Nelem_channel, ElemOrder, FourierOrder, 0, 1, 0.2, 0.02, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, geom_mode);
      elem_lst_trg = std::get<0>(build_trg);
    } else {
      SCTL_ASSERT(false); // not implemented
    }
    // std::cout << "dim of elem_lst_trg is " << elem_lst_trg.Size() << std::endl;

    VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
    sctl::Vector<Real> X0_all = vol_vis.GetCoord(); // set new target coordinates
    // std::cout << "dim of X0_all is " << X0_all.Dim() << std::endl;
    sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
    if (!ptcls.Dim()) {
      X0 = X0_all;
      filtered_inds = 1;
    } else {
      std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
      X0 = std::get<0>(trg_tuple);
      filtered_inds = std::get<1>(trg_tuple);
    }
    // std::cout << "dim of X0 is " << X0.Dim() << std::endl;
    LayerPotenOp0.SetTargetCoord(X0);
    sctl::Vector<Real> U;
    BIO(&U, sigma);
    U += bg_flow(X0);
    sctl::Vector<Real> U_vis(X0_all.Dim());
    if (!ptcls.Dim()) {
      U_vis = U;
    } else {
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
    }
    sctl::Vector<sctl::Long> size_loc(1);
    size_loc[0] = X0_all.Dim();
    sctl::Vector<sctl::Long> size_all(1);
    comm.Allreduce((sctl::Iterator<sctl::Long>) size_loc.begin(), (sctl::Iterator<sctl::Long>) size_all.begin(), 1, sctl::CommOp::SUM);
    std::cout << "rank " << comm.Rank() << " size loc = " << size_loc[0] << ", size all is " << size_all[0] << std::endl;
    std::string filename = "U_8_16_"+std::to_string(channel_mode)+"_"+std::to_string(geom_mode)+"_"+std::to_string(Nptcl);
    std::string filename_out = "out/"+filename+".txt";
    std::string filename_vis = "vis/"+filename;
    if (write_ref) {
      vol_vis.WriteVTK(filename_vis, U_vis);
      // sctl::Vector<Real> U_vis_all(size_all[0]);
      // comm.Allgather((sctl::Iterator<Real>) U_vis.begin(), size_loc[0], (sctl::Iterator<Real>) U_vis_all.begin(), size_all[0]);
      // std::cout << U_vis_all.Dim() << std::endl;
      // if (!comm.Rank()) {
      //   U_vis_all.Write(filename_out.c_str());
      // }
    } else {
      sctl::Vector<Real> U_ref;
      if (!comm.Rank()) {
        U_ref.Read(filename_out.c_str());
      }
      comm.PartitionN(U_ref,size_loc[0]);
      // std::cout << "dim of U ref is " << U_ref.Dim() << ", dim of U vis is " << U_vis.Dim() << std::endl;
      const auto err = U_vis - U_ref;
      double max_err = 0;
      for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
      sctl::Vector<Real> err_loc(1);
      err_loc[0] = max_err;
      sctl::Vector<Real> err_all(1);
      err_all[0] = 0;
      comm.Allreduce((sctl::Iterator<sctl::Long>) err_loc.begin(), (sctl::Iterator<sctl::Long>) err_all.begin(), 1, sctl::CommOp::MAX);
      if (!comm.Rank()) {
        std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << std::endl;
      }
    }
    // if (!comm.Rank()) {
    //   if (write_ref) {
    //     vol_vis.WriteVTK(filename_vis, U_vis);
    //     U_vis.Write(filename_out.c_str());
    //   } else {
    //     sctl::Vector<Real> U_ref;
    //     U_ref.Read(filename_out.c_str());
    //     // 
    //     // std::cout << "dim of U ref is " << U_ref.Dim() << ", dim of U vis is " << U_vis.Dim() << std::endl;
    //     const auto err = U_vis_all - U_ref;
    //     // std::cout << "dim of U ref is " << U_ref.Dim() << ", dim of U vis is " << U_vis.Dim() << std::endl;
    //     // sctl::Vector<Real> U_ref_loc(U_vis.Dim(),(sctl::Iterator<Real>) U_ref.begin(),false);
    //     // const auto err = U_vis - U_ref_loc;
    //     double max_err = 0;
    //     for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
    //     std::cout<<"Max error = "<< max_err <<'\n';
    //   }
    // }

  }
}

int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    //sctl::Profile::Enable(true);
    sctl::Comm comm = sctl::Comm::World();
    long Nelem_channel = std::stol(argv[1]); // number of elements
    long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
    int write_ref = std::stoi(argv[3]); // whether the parameters are considered ''artifical true solution''.
    long Nptcl = std::stol(argv[4]); // number of particles inside
    long channel_mode = std::stol(argv[5]); // =0: straight; =1: sinusoidal mag=0.1; =2: sinusoidal mag=0.3
    long geom_mode = std::stol(argv[6]); // =0: spheres; =1: spheroids; =3: bacteria; =4: loop.
    test<Real>(Nelem_channel, FourierOrder, (write_ref==1), comm, Nptcl, channel_mode, geom_mode);
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}

