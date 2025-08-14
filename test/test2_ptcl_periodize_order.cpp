#include "periodize.hpp"
#include "utils.hpp"
// Test script focusing on particle-only 1 and 3 periodic solvers via "manufactured solutions":
//     Randomly place Stokeslets inside each particle, with random forces, assumed to be periodic like the particles
//     Evaluate flow field from Stokeslet on surface of each particle as imposed BC
//     Solve BIE for periodic spheres given BC
//     Evaluate flow field at targets exterior to spheres through BIE solution and compare with exact flow field from Stokeslets.

// Loop over copies and add consecutively, to reduce memory requirements. Perhaps do 2D planes at a time.
template <class Real> sctl::Vector<Real> exact_field(const sctl::Vector<Real>& Xtrg, const sctl::Vector<Real>& Xsrc, const sctl::Vector<Real>& sigma, const sctl::Long Ncopy, const sctl::Integer peri_mode) {
    sctl::Stokes3D_FxU ker;
  
    const sctl::Long N = Xtrg.Dim()/3;
    sctl::Vector<Real> U(N*3);
    U = 0.;
    PeriodicGeom<Real> obj;

    if (peri_mode == 1) {
        sctl::Vector<Real> Xsrc_ = obj.X_nbr_copy(Xsrc,Ncopy,peri_mode);
        sctl::Vector<Real> sigma_ = obj.vec_nbr_copy(sigma,Ncopy,peri_mode);
        ker.Eval(U,Xtrg,Xsrc_,Xsrc_,sigma_);
    } else {
        const sctl::Vector<Real> Xsrc_2D = obj.X_nbr_copy(Xsrc,Ncopy,2);
        const sctl::Vector<Real> sigma_2D = obj.vec_nbr_copy(sigma,Ncopy,2);
        sctl::Vector<Real> Xsrc_plane = Xsrc_2D;
        sctl::Vector<Real> Uloc = U;
        for (int k3=-Ncopy; k3<=Ncopy; k3++) {
            Uloc = 0.; // reset U
            Xsrc_plane = Xsrc_2D; // reset to center plane each time.
            for (int src_ind = 0; src_ind < Xsrc_2D.Dim()/3; src_ind ++) {
                Xsrc_plane[src_ind*3 + 2] += k3; // move z direction up and down.
            }
            ker.Eval(Uloc,Xtrg,Xsrc_plane,Xsrc_plane,sigma_2D);
            U += Uloc; // Add contribution to Xtrg from this layer of sources.
        }
    }

    return U;
}

template <class Real> void test(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Integer peri_mode, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, sctl::Long Ncopy, sctl::Vector<sctl::Long> level_lst, sctl::Vector<sctl::Long> m0_lst) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real tol = 1e-14;
    Real gmres_tol = 1e-12;
    const sctl::Long ElemOrder = 10;
    
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
    } else { 
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 0, 1, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 1, peri_mode, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst_nbr = std::get<0>(build_nbr);
        NormalOrient = std::get<1>(build_nbr);
    }
    const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); 
    Nptcl = ptcls_rs.Dim(); 
    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);

    // Create point charges at random locations close to particle center, by a distance of at most 0.2r.
    sctl::Long Ncharge;
    if (Nptcl < 150) {
        // two equal and opposite charges per particle
        Ncharge = 2*Nptcl;
    } else {
        // only first 150 particles get charges inside. (arbitrary, to limit true solution timing)
        Ncharge = 2*150;
    }
    // Currently one charge per particle.
    sctl::Vector<Real> Xsrc(Ncharge*3);
    sctl::Vector<Real> Stokeslet_sigma(Ncharge*3);
    srand48(2);
    for (sctl::Long i=0; i<Ncharge/2; i++) {

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
        
    }
    // sctl::Vector<Real> field_on_surf = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy, peri_mode);
    sctl::Vector<Real> field_on_surf;
    if (peri_mode == 1) {
        field_on_surf = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy, peri_mode);
    } else {
        // Setting peri_mode=1 and Ncopy = 0 calculates ker.Eval from current copy to current copy.
        field_on_surf = exact_field(X0, Xsrc, Stokeslet_sigma, 0, 1);
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst_nbr);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal, comm); // potential from elem_lst0 to proxy points
    LayerPotenOp_proxy.AddElemList(elem_lst0);
    LayerPotenOp_proxy.SetAccuracy(tol);

    // PRECONDITIONING ////////////////////////////
    sctl::Vector<sctl::Long> ptcls_pre;
    sctl::Vector<Real> ptcls_Xcs_pre, ptcls_rs_pre;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_precond = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm.Self(), ptcls_pre, ptcls_rs_pre, ptcls_Xcs_pre, geom_mode);
    sctl::SlenderElemList<Real> elem_lst_precond = std::get<0>(build_precond);
    sctl::Vector<Real> X0_precond; // target coordinates
    elem_lst_precond.GetNodeCoord(&X0_precond, nullptr, nullptr);

    StokesBIO Precond_bio(SL_scal, DL_scal, comm.Self());
    Precond_bio.SetAccuracy(tol); // set quadrature accuracy
    Precond_bio.AddElemList(elem_lst_precond);
    Precond_bio.SetTargetCoord(X0_precond);

    const auto BIO_1ptcl = [&DL_scal,&Precond_bio](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        U->SetZero();
        Precond_bio.ComputePotential(*U, sigma);
        (*U) += sigma * 0.5 * DL_scal;
    };

    sctl::Long A11size = 3*ElemOrder*FourierOrder*Nelem;
    sctl::Vector<sctl::Vector<Real>> PrecondMat(A11size);
    sctl::Vector<Real> SigmaCol_precond(A11size);
    for (sctl::Long col=0; col < A11size; col ++) {
      SigmaCol_precond = 0.;
      SigmaCol_precond[col] = 1.;
      BIO_1ptcl(PrecondMat.begin() + col,SigmaCol_precond);
    }
    sctl::Matrix<Real> A11(A11size,A11size);
    for (long col=0; col < A11size; col++) {
      for (long row = 0; row < A11size; row++) {
        A11(row,col) = PrecondMat[col][row];
      }
    }      
    sctl::Matrix<Real> Usvd, VT, S, SforInv;
    sctl::Matrix<Real> A11forSVD = sctl::Matrix<Real>(A11);
    A11forSVD.SVD(Usvd, S, VT);
    SforInv = sctl::Matrix<Real>(S);
    sctl::Matrix<Real> Sinv = SforInv.pinv(1e-16);

    // Apply A11inv to each panel of vec.
    const auto AinvApply = [&Usvd,&Sinv,&VT,&A11size, &comm](const sctl::Vector<Real>& vec) {
        sctl::Long N = vec.Dim();
        sctl::Long Nptcl = N / A11size; 
        sctl::Vector<Real> AinvVec(N);
        for (sctl::Long i=0; i<Nptcl; i++) {
            // for each particle, apply A11inv.
            sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
            sctl::Matrix<Real> AinvVecMat = VT.Transpose() * (Sinv * (Usvd.Transpose() * vecMat));
            for (sctl::Long j=0; j<A11size; j++) {
                AinvVec[i*A11size + j] = AinvVecMat(j,0);
            }
        }
        return AinvVec;
    };

    //////////////////////////////////////////////
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> A11invF = AinvApply(field_on_surf);

    PeriodicGeom<Real> trg;    
    CubeVolumeVisShifted<Real> vol_vis(20, 1.0, comm);
    // VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
    // X0 = vol_vis.GetCoord();
    sctl::Vector<Real> X0_all = vol_vis.GetCoord();
    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    const sctl::Vector<Real> X0_trg = std::get<0>(trg_tuple);
    const sctl::Vector<sctl::Long> filtered_inds = std::get<1>(trg_tuple);

    /// Proxy: LOOP through level and m0 values. 
    for (int l_ind=0; l_ind<level_lst.Dim(); l_ind++) {
        for (int m_ind=0; m_ind<m0_lst.Dim(); m_ind++) {
            sctl::Long level = level_lst[l_ind];
            sctl::Long m0 = m0_lst[m_ind];
            // if (!comm.Rank()) {
            std::cout << "level = " << level << ", m0 = " << m0 << std::endl;
            // }
            sctl::Vector<Real> X_proxy;
            if (peri_mode == 1) {
                X_proxy = Periodize1D<Real>::GetProxySurf(level,m0); // proxy points coordinates
            } else if (peri_mode == 3) {
                X_proxy = Periodize3D<Real>::GetProxySurf(); // proxy points coordinates
            } else {
                SCTL_ASSERT(false);
            }
            LayerPotenOp_proxy.ClearSetup();
            LayerPotenOp_proxy.SetTargetCoord(X_proxy);
            LayerPotenOp_proxy.Setup();
            elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
            LayerPotenOp0.ClearSetup();
            LayerPotenOp0.SetTargetCoord(X0);
            LayerPotenOp0.Setup();

            // periodized layer potential operator
            auto BIO = [&DL_scal,&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat,NormalOrient,&peri_mode,&level,&m0](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
                const sctl::Long N = sigma.Dim();

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
                        Periodize1D<Real>::EvalFarField(U_far, X0, U_proxy, level, m0);
                    } else if (peri_mode==3) {
                        // 3-periodic
                        Periodize3D<Real>::EvalFarField(U_far, X0, U_proxy);
                    } else {
                        std::cout << "2-periodic not yet implemented." << std::endl;
                        SCTL_ASSERT(false);
                    }
                    // std::cout << "size of U is " << U->Dim() << ", size of Ufar is " << U_far.Dim() << std::endl;
                    (*U) += U_far;
                } 
                // comm.Barrier();
            };

            auto BIO_precond = [&BIO,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
                sctl::Vector<Real> Uloc;
                BIO(&Uloc,sigma);
                // LEFT PRECONDITIONER: u -> A11inv*u
                (*U) = AinvApply(Uloc);
            };

            if (m0<=4 || level<=7) {
                gmres_tol = 1e-6;
            } else if (level<=9)
                gmres_tol = 1e-8;
            else {
                gmres_tol = 1e-12;
            }

            sctl::Vector<Real> sigma;
            if (l_ind==0 && m_ind==0) {
                // std::cout << "Set up Krylov preconditioner" << std::endl;
                sctl::Vector<Real> sigma_temp;
                solver(&sigma_temp, BIO_precond, A11invF, gmres_tol, -1, false, nullptr, &krylov_precond); 
            }
            solver(&sigma, BIO_precond, A11invF, gmres_tol, -1, false, nullptr, &krylov_precond);

            { // Evaluate in interior, and write visualization
                X0 = X0_trg;
                LayerPotenOp0.ClearSetup();
                LayerPotenOp0.SetTargetCoord(X0);
                LayerPotenOp0.Setup();
                // Exact solution from Xsrc and Stokeslet_sigma
                sctl::Vector<Real> field_on_trg = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy, peri_mode);
                sctl::Vector<Real> U;
                BIO(&U, sigma);

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
                    std::cout<<"Max error = "<< std::setprecision(15) << err_all[0] << ", Max relative error = " << err_all[0] / u_all[0] << std::endl;
                    // std::cout<<"Max u = "<< std::setprecision(15) << u_all[0] << std::endl;
                    // std::cout<<"Max relative error = "<< std::setprecision(15) << err_all[0] / u_all[0] << std::endl;
                }
                comm.Barrier();
            }
        }
    }
}


int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    sctl::Comm comm = sctl::Comm::World();
    long Nelem_ptcl = std::stol(argv[1]); // number of elements
    long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
    int peri_mode = std::stoi(argv[3]); // what kind of periodicity does the system have; peri_mode = j for j-periodic.
    long Nptcl = std::stol(argv[4]); // number of particles inside
    long geom_mode = std::stol(argv[5]); // =0: spheres; =1: spheroids; =3: bacteria; =4: loop.
    long Ncopy = std::stol(argv[6]); // Number of copies on each side to add to sources 
    sctl::Vector<sctl::Long> level_lst;
    // for (int i=1; i<=10; i++) { // all params
    //     level_lst.PushBack(i);
    // }
    level_lst.PushBack(10); // just for params that timed out
    level_lst.PushBack(15);
    // level_lst.PushBack(20);
    // level_lst.PushBack(30);
    sctl::Vector<sctl::Long> m0_lst;
    for (int i=16; i<20; i*=2) {
        m0_lst.PushBack(i);
    }
    // m0_lst.PushBack(20); // looks like same error as m0=16

    test<Real>(Nelem_ptcl, FourierOrder, peri_mode, comm, Nptcl, geom_mode, Ncopy, level_lst, m0_lst);
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}
