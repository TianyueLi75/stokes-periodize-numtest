#include "periodize.hpp"
#include "utils.hpp"
// Test script focusing on particle-only 1 and 3 periodic solvers via "manufactured solutions":
//     Randomly place Stokeslets inside each particle, with random forces, assumed to be periodic like the particles
//     Evaluate flow field from Stokeslet on surface of each particle as imposed BC
//     Solve BIE for periodic spheres given BC
//     Evaluate flow field at targets exterior to spheres through BIE solution and compare with exact flow field from Stokeslets.

// /**
//  * Background flow with unit pressure gradient along X-axis.
//  */
// template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
//     const Real pdrive = 1;
//     const sctl::Long N = X.Dim()/3;
//     sctl::Vector<Real> U(N*3);
//     for (sctl::Long i = 0; i < N; i++) {
//         const auto x = X.begin() + i*3;
//         U[i*3+0] = -pdrive * ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
//         U[i*3+1] = 0;
//         U[i*3+2] = 0;
//     }
//     return U;
// }

/**
 * Calculate flow field at target Xtrg give Stokeslets at Xsrc and strength sigma, repeated Ncopy times
 * Output in order [Ex1; Ey1; Ez1; Ex2; Ey2; Ez2;...]
 */
template <class Real> sctl::Vector<Real> exact_field(const sctl::Vector<Real>& Xtrg, const sctl::Vector<Real>& Xsrc, const sctl::Vector<Real>& sigma, const sctl::Long Ncopy, const sctl::Integer peri_mode) {
    sctl::Stokes3D_FxU ker;
  
    const sctl::Long N = Xtrg.Dim()/3;
    sctl::Vector<Real> U(N*3);
    U = 0.;
    PeriodicGeom<Real> obj;
    sctl::Vector<Real> Xsrc_ = obj.X_nbr_copy(Xsrc,Ncopy,peri_mode);
    sctl::Vector<Real> sigma_ = obj.vec_nbr_copy(sigma,Ncopy,peri_mode);
    ker.Eval(U,Xtrg,Xsrc_,Xsrc_,sigma_);
    return U;
}

template <class Real> void test(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Integer peri_mode, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, sctl::Long Ncopy) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real tol = 1e-15;
    const Real gmres_tol = 1e-13;
    const sctl::Long ElemOrder = 10;

    // TODO: currently specifically using one sphere, but should be able to look at other constructions.
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 1, peri_mode, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst_nbr = std::get<0>(build_nbr);
        NormalOrient = std::get<1>(build_nbr);
    } else { // TODO: once many_ptcl2() gets Nptcl argument, change that.
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 0, 1, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 1, peri_mode, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst_nbr = std::get<0>(build_nbr);
        NormalOrient = std::get<1>(build_nbr);
    }
    const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); 
    Nptcl = ptcls_rs.Dim(); 
    std::cout << "periodic mode is " << peri_mode << ", Nrepeat is " << Nrepeat << std::endl;

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Vector<Real> X_proxy;
    if (peri_mode == 1) {
        X_proxy = Periodize1D<Real>::GetProxySurf(); // proxy points coordinates
    } else if (peri_mode == 3) {
        X_proxy = Periodize3D<Real>::GetProxySurf(); // proxy points coordinates
    } else {
        SCTL_ASSERT(false);
    }

    // Create point charges at random locations close to particle center, by a distance of at most 0.2r.
    // Currently one charge per particle.
    sctl::Vector<Real> Xsrc(2*Nptcl*3);
    sctl::Vector<Real> Stokeslet_sigma(2*Nptcl*3);
    srand48(2);
    for (sctl::Long i=0; i<Nptcl; i++) {

        const Real disp = 0.2 * ptcls_rs[i];
        const Real disp_y = disp * drand48();
        const Real disp_z = disp * drand48();
        const Real rand_mag = drand48()-0.5;
        Xsrc[i*6+0] = ptcls_Xcs[i*3+0];
        Xsrc[i*6+1] = ptcls_Xcs[i*3+1] + disp_y;
        Xsrc[i*6+2] = ptcls_Xcs[i*3+2] + disp_z;
        Stokeslet_sigma[i*6+0] = 0.; 
        Stokeslet_sigma[i*6+1] = -rand_mag * disp_y; 
        Stokeslet_sigma[i*6+2] = -rand_mag * disp_z; 
        Xsrc[i*6+3] = ptcls_Xcs[i*3+0];
        Xsrc[i*6+4] = ptcls_Xcs[i*3+1] - disp_y;
        Xsrc[i*6+5] = ptcls_Xcs[i*3+2] - disp_z;
        Stokeslet_sigma[i*6+3] = 0.; 
        Stokeslet_sigma[i*6+4] = rand_mag * disp_y; 
        Stokeslet_sigma[i*6+5] = rand_mag * disp_z; 
        
        // DEBUG
        // std::cout << "location of " << i*2+0 << "-th charge is (" << Xsrc[i*6+0] << ", " << Xsrc[i*6+1] << ", " << Xsrc[i*6+2] << "); strength is " << Stokeslet_sigma[i*6+0] << ", " << Stokeslet_sigma[i*6+1] << ", " << Stokeslet_sigma[i*6+2] << std::endl;
        // std::cout << "location of " << i*2+1 << "-th charge is (" << Xsrc[i*6+3] << ", " << Xsrc[i*6+4] << ", " << Xsrc[i*6+5] << "); strength is " << Stokeslet_sigma[i*6+3] << ", " << Stokeslet_sigma[i*6+4] << ", " << Stokeslet_sigma[i*6+5] << std::endl;
    }
    sctl::Vector<Real> field_on_surf = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy, peri_mode);

    if (write_ref) {
        elem_lst0.WriteVTK("vis/ptcl_density",field_on_surf,comm);
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst_nbr);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);

    StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal, comm); // potential from elem_lst0 to proxy points
    LayerPotenOp_proxy.AddElemList(elem_lst0);
    LayerPotenOp_proxy.SetTargetCoord(X_proxy);
    LayerPotenOp_proxy.SetAccuracy(tol);

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
        if (DL_scal && U->Dim() == N) {
            (*U) -= sigma*0.5*NormalOrient * DL_scal;
        }

        { // Add far-field
            sctl::Vector<Real> U_proxy, U_far;
            LayerPotenOp_proxy.ComputePotential(U_proxy, sigma);
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

    sctl::Vector<Real> sigma;
    sctl::GMRES<Real> solver(comm);
    solver(&sigma, BIO, field_on_surf, gmres_tol);
    sctl::Profile::print();

    { // Evaluate in interior, and write visualization
        // std::cout << "Rank " << comm.Rank()<< " calculating target points." << std::endl;
        PeriodicGeom<Real> trg;    
        CubeVolumeVisShifted<Real> vol_vis(20, 1.0, comm);
        // VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        // X0 = vol_vis.GetCoord();
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        sctl::Profile::print();
        // Exact solution from Xsrc and Stokeslet_sigma
        sctl::Vector<Real> field_on_trg = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy, peri_mode);

        // get max abs error
        const auto err = U - field_on_trg;
        double max_err = 0;
        Real max_u = 0.;
        for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
        for (const auto e : field_on_trg) max_u = std::max<Real>(max_u, sctl::fabs(e));
        
        // std::cout<< "Rank " << comm.Rank() << " Max error = "<< std::setprecision(10) << max_err << ", max u = " << max_u << std::endl;

        sctl::Vector<Real> err_loc(1);
        err_loc[0] = max_err;
        sctl::Vector<Real> err_all(1);
        err_all[0] = 0;
        // comm.Allreduce((sctl::Iterator<sctl::Long>) err_loc.begin(), (sctl::Iterator<sctl::Long>) err_all.begin(), 1, sctl::CommOp::MAX);
        comm.Allreduce((sctl::Iterator<Real>) err_loc.begin(), (sctl::Iterator<Real>) err_all.begin(), 1, sctl::CommOp::MAX);
        
        sctl::Vector<Real> u_loc(1);
        u_loc[0] = max_u;
        sctl::Vector<Real> u_all(1);
        u_all[0] = 0.;
        comm.Allreduce((sctl::Iterator<Real>) u_loc.begin(), (sctl::Iterator<Real>) u_all.begin(), 1, sctl::CommOp::MAX);

        if (!comm.Rank()) {
            std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << std::endl;
            std::cout<<"Max relative error = "<< std::setprecision(10) << err_all[0] / u_all[0] << std::endl;
        }

        
        if (write_ref) {
            // For visualization
            sctl::Vector<Real> U_vis(X0_all.Dim());
            sctl::Vector<Real> U_vis_exact(X0_all.Dim());
            sctl::Vector<Real> err_vis(X0_all.Dim());
            U_vis = 0.;
            U_vis_exact = 0.;
            err_vis = 0.;
            sctl::Long X1_ptr = 0;
            for (sctl::Long i=0; i<X0_all.Dim()/3; i++) {
                if (filtered_inds[i] == 0) {
                    U_vis[i*3] = U[X1_ptr*3];
                    U_vis[i*3+1] = U[X1_ptr*3+1];
                    U_vis[i*3+2] = U[X1_ptr*3+2];

                    U_vis_exact[i*3] = field_on_trg[X1_ptr*3];
                    U_vis_exact[i*3+1] = field_on_trg[X1_ptr*3+1];
                    U_vis_exact[i*3+2] = field_on_trg[X1_ptr*3+2];

                    err_vis[i*3] = err[X1_ptr*3];
                    err_vis[i*3+1] = err[X1_ptr*3+1];
                    err_vis[i*3+2] = err[X1_ptr*3+2];
                    X1_ptr += 1;
                }
            }
            vol_vis.WriteVTK("vis/exact_soln", U_vis_exact);
            vol_vis.WriteVTK("vis/BIE_soln", U_vis);
            vol_vis.WriteVTK("vis/err", err_vis);
        }
    }
}


int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    sctl::Profile::Enable(true);
    sctl::Comm comm = sctl::Comm::World();
    long Nelem_ptcl = std::stol(argv[1]); // number of elements
    long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
    int write_ref = std::stol(argv[3]);
    int peri_mode = std::stoi(argv[4]); // what kind of periodicity does the system have; peri_mode = j for j-periodic.
    long Nptcl = std::stol(argv[5]); // number of particles inside
    long geom_mode = std::stol(argv[6]); // =0: spheres; =1: spheroids; =3: bacteria; =4: loop.
    long Ncopy = std::stol(argv[7]); // Number of copies on each side to add to sources 

    // sctl::Profile::Tic("test");
    test<Real>(Nelem_ptcl, FourierOrder, (write_ref==1), peri_mode, comm, Nptcl, geom_mode, Ncopy);
    // sctl::Profile::Toc();
    // sctl::Profile::print();
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}
