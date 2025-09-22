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

// Uniform background flow in x direction.
template <class Real> sctl::Vector<Real> bg_unif_flow(const sctl::Vector<Real>& X) {
    sctl::Vector<Real> U = X;
    const sctl::Long N = X.Dim() /3;
    U = 1.; // background flow diagonal to avoid planes of unaffected flows between periods.
    // for (sctl::Long i = 0; i < N; i++) {
    //     U[i*3+0] = 1.;
    // }
    return U;
}

/**
 * Reference solution for checking error.
 */
template <class Real> sctl::Vector<Real> u_ref(const sctl::Vector<Real>& X) {
  const sctl::Long N = X.Dim()/3;
  sctl::Vector<Real> U(N*3);
  for (sctl::Long i = 0; i < N; i++) {
    const auto x = X.begin() + i*3;
    U[i*3+0] = 1e-2 - ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
    U[i*3+1] = 0;
    U[i*3+2] = 0;
  }
  return U;
}

template <class Real> void test(sctl::Long Nelem_channel, sctl::Long FourierOrder, bool write_ref, sctl::Comm comm, sctl::Long Nptcl, sctl::Long channel_mode, sctl::Long geom_mode, sctl::Long precond_mode) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real tol = 1e-14;
    const Real gmres_tol = 1e-13;
    const sctl::Long ElemOrder = 10;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Long ptcl_ord = 4;
    if (Nptcl>0) {
        ptcls.ReInit(Nptcl);
        ptcls = ptcl_ord;
    }
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    sctl::Long peri_mode = 1;
    if (channel_mode == 0) {// straight channel
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_straight(Nelem_channel, ElemOrder, FourierOrder, 0, peri_mode, 0.2, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_straight(Nelem_channel, ElemOrder, FourierOrder, 1, peri_mode, 0.2, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst_nbr = std::get<0>(build_nbr);
        NormalOrient = std::get<1>(build_nbr); // only needed for self eval which uses nbr object.
    } else if (channel_mode == 1) {// sinusoidal with mag = 0.1
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_sinusoidal(Nelem_channel, ElemOrder, FourierOrder, 0, peri_mode, 0.2, 0.1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_sinusoidal(Nelem_channel, ElemOrder, FourierOrder, 1, peri_mode, 0.2, 0.1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst_nbr = std::get<0>(build_nbr);
        NormalOrient = std::get<1>(build_nbr);
    } else if (channel_mode == 2) {// sinusoidal with mag = 0.3, changed recently, unchecked.
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_sinusoidal(Nelem_channel, ElemOrder, FourierOrder, 0, peri_mode, 0.1, 0.3, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_sinusoidal(Nelem_channel, ElemOrder, FourierOrder, 1, peri_mode, 0.1, 0.3, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);  
        elem_lst_nbr = std::get<0>(build_nbr);
        NormalOrient = std::get<1>(build_nbr);
    } else if (channel_mode == 3) { // rmin = 0.1, rmax = 0.2
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_conv_div(Nelem_channel, ElemOrder, FourierOrder, 0, peri_mode, 0.1, 0.2, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, geom_mode);
        elem_lst0 = std::get<0>(build0);
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_conv_div(Nelem_channel, ElemOrder, FourierOrder, 1, peri_mode, 0.1, 0.2, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, geom_mode);  
        elem_lst_nbr = std::get<0>(build_nbr);
        NormalOrient = std::get<1>(build_nbr);
    } else if (channel_mode == 4) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_spiral(Nelem_channel, ElemOrder, FourierOrder, 0, peri_mode, 0.5, 0.05, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_spiral(Nelem_channel, ElemOrder, FourierOrder, 1, peri_mode, 0.5, 0.05, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);  
        elem_lst_nbr = std::get<0>(build_nbr);
        NormalOrient = std::get<1>(build_nbr);
    } else if (channel_mode == 5) { // rmin = 0.01, rmax = 0.04, set inside function.
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_trefoil(Nelem_channel, ElemOrder, FourierOrder, 0, peri_mode, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, geom_mode);
        elem_lst0 = std::get<0>(build0);
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_trefoil(Nelem_channel, ElemOrder, FourierOrder, 1, peri_mode, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, geom_mode);  
        elem_lst_nbr = std::get<0>(build_nbr);
        NormalOrient = std::get<1>(build_nbr);
    } else {
        SCTL_ASSERT(false); // not implemented
    }
    // std::cout << "Size of elem_lst_nbr is " << elem_lst_nbr.Size() << ", Size of elem_lst0 is " << elem_lst0.Size() <<std::endl;
    const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); // should be 3
    Nptcl = ptcls_rs.Dim(); // Number of particles could have changed after initializing.
    std::cout << "periodic mode is " << peri_mode << ", Nrepeat is " << Nrepeat << ", Nptcl is " << Nptcl << std::endl;

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    // const auto X_proxy = Periodize1D<Real>::GetProxySurf(); // proxy points coordinates
    sctl::Vector<Real> X_proxy;
    if (peri_mode == 1) {
        X_proxy = Periodize1D<Real>::GetProxySurf(); // proxy points coordinates
    } else if (peri_mode == 3) {
        X_proxy = Periodize3D<Real>::GetProxySurf(); // proxy points coordinates
    } else {
        SCTL_ASSERT(false);
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst_nbr);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);

    StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal, comm); // potential from elem_lst0 to proxy points
    LayerPotenOp_proxy.AddElemList(elem_lst0);
    LayerPotenOp_proxy.SetTargetCoord(X_proxy);
    LayerPotenOp_proxy.SetAccuracy(tol);

    // elem_lst0.WriteVTK("vis/Channel_ptcl_precond",X0,comm);

    // Get global index of the starting panel on this process
    sctl::Vector<sctl::Long> ElemOrderVec_temp(Nelem_channel + ptcl_ord * Nptcl);
    ElemOrderVec_temp = ElemOrder;
    sctl::Vector<sctl::Long> FourierOrderVec_temp(ElemOrderVec_temp);
    FourierOrderVec_temp = FourierOrder;
    std::tuple<sctl::Long,sctl::Long> indtpl = obj.GetGlobalIdx(ElemOrderVec_temp, FourierOrderVec_temp, comm);
    sctl::Long loc_elem_cnt = std::get<0>(indtpl);
    sctl::Long loc_elem_dsp = std::get<1>(indtpl);
    std::cout << "Rank " << comm.Rank() << " loc elem cnt = " << loc_elem_cnt << ", loc elem dsp = " << loc_elem_dsp << std::endl;

    // ======================= PRECONDITIONING : CYLINDER ====================================================
    sctl::Vector<Real> Xc_precond, eps_precond; 
    sctl::Vector<sctl::Long> ElemOrderVec_precond(1), FourierOrderVec_precond(1);
    ElemOrderVec_precond[0] = ElemOrder;
    FourierOrderVec_precond[0] = FourierOrder;
    // Determine approximate radius of channel based on channel_mode
    Real channel_radius = 0.5;
    if (channel_mode == 0 || channel_mode == 1) {
        channel_radius = 0.2;
    } else if (channel_mode == 2) {
        channel_radius = 0.1;
    } else if (channel_mode == 3) {
        channel_radius = 0.15;
    } else if (channel_mode == 4) {
        channel_radius = 0.05;
    } else if (channel_mode == 5) {
        channel_radius = 0.035;
    } else {
        SCTL_ASSERT(false);
    }
    // std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_precond = obj.build_straight(1, ElemOrder, FourierOrder, 0, peri_mode, 0.25, comm.Self(), ptcls_pre, ptcls_rs_pre, ptcls_Xcs_pre, geom_mode);
    // elem_lst_precond = std::get<0>(build_precond);

    // ALTERNATIVE: smaller panel matching channel panel length and radius.
    const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrder);
    for (sctl::Long j = 0; j < ElemOrder; j++) { // loop over panel nodes
      const Real x = (nodes[j]) / Nelem_channel; // size of precond panel should be same as one panel on pipe
      Xc_precond.PushBack(x+0.5); //  shift panel to center of unit box, arbitrary.
      Xc_precond.PushBack(0.5); 
      Xc_precond.PushBack(0.5); 
      eps_precond.PushBack(channel_radius); 
    }
    sctl::SlenderElemList<Real> elem_lst_precond(ElemOrderVec_precond, FourierOrderVec_precond, Xc_precond, eps_precond);

    sctl::Vector<Real> X0_precond; // target coordinates
    elem_lst_precond.GetNodeCoord(&X0_precond, nullptr, nullptr);

    // Alternatively, use strictly self-to-self -- no neighbors either.
    StokesBIO Precond_bio(SL_scal, DL_scal, comm.Self());
    Precond_bio.SetAccuracy(tol); // set quadrature accuracy
    Precond_bio.AddElemList(elem_lst_precond);
    Precond_bio.SetTargetCoord(X0_precond);

    const auto BIO_1panel = [&DL_scal,&Precond_bio](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        U->SetZero();
        Precond_bio.ComputePotential(*U, sigma);
        (*U) -= sigma * 0.5 * DL_scal; // for preconditioner, will always be self-to-self so always add. For panels (on channel), normal orient = 1.
    };

    sctl::Long A11size = 3*ElemOrder*FourierOrder;
    sctl::Vector<sctl::Vector<Real>> PrecondMat(A11size);
    sctl::Vector<Real> SigmaCol_precond(A11size);
    for (sctl::Long col=0; col < A11size; col ++) {
        SigmaCol_precond = 0.;
        SigmaCol_precond[col] = 1.;
        BIO_1panel(PrecondMat.begin()+col,SigmaCol_precond);
    }
    sctl::Matrix<Real> A11(A11size,A11size);
    for (sctl::Long col=0; col < A11size; col++) {
        for (sctl::Long row = 0; row < A11size; row++) {
            A11(row,col) = PrecondMat[col][row];
        }
    }
    
    sctl::Matrix<Real> Usvd, VT, S, SforInv;
    sctl::Matrix<Real> A11forSVD = sctl::Matrix<Real>(A11);
    A11forSVD.SVD(Usvd, S, VT);
    SforInv = sctl::Matrix<Real>(S);
    sctl::Matrix<Real> Sinv = SforInv.pinv(1e-16);

    // =============== PRECONDITIONING : Particle =======================================
    sctl::Vector<sctl::Long> ptcls_pre;
    sctl::Vector<Real> ptcls_rs_pre, ptcls_Xcs_pre;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_precond_ptcl = obj.many_ptcls1(ptcl_ord, ElemOrder, FourierOrder, 0, 1, comm.Self(), ptcls_pre, ptcls_rs_pre, ptcls_Xcs_pre, geom_mode); 
    sctl::SlenderElemList<Real> elem_lst_precond_ptcl = std::get<0>(build_precond_ptcl);
    sctl::Vector<Real> X0_precond_ptcl; // target coordinates
    elem_lst_precond_ptcl.GetNodeCoord(&X0_precond_ptcl, nullptr, nullptr);

    StokesBIO Precond_bio_ptcl(SL_scal, DL_scal, comm.Self());
    Precond_bio_ptcl.SetAccuracy(tol); // set quadrature accuracy
    Precond_bio_ptcl.AddElemList(elem_lst_precond_ptcl);
    Precond_bio_ptcl.SetTargetCoord(X0_precond_ptcl);
    const auto BIO_1ptcl = [&DL_scal,&Precond_bio_ptcl](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        U->SetZero();
        Precond_bio_ptcl.ComputePotential(*U, sigma);
        (*U) += sigma * 0.5 * DL_scal; // particles always have normal orient = -1; (when defined U -= normal orient * sigma)
    };

    sctl::Long A11size_ptcl = 3*ElemOrder*FourierOrder*ptcl_ord;
    sctl::Vector<sctl::Vector<Real>> PrecondMat_ptcl(A11size_ptcl);
    sctl::Vector<Real> SigmaCol_precond_ptcl(A11size_ptcl);
    for (sctl::Long col=0; col < A11size_ptcl; col ++) {
      SigmaCol_precond_ptcl = 0.;
      SigmaCol_precond_ptcl[col] = 1.;
      BIO_1ptcl(PrecondMat_ptcl.begin() + col,SigmaCol_precond_ptcl);
    }
    sctl::Matrix<Real> A11_ptcl(A11size_ptcl,A11size_ptcl);
    for (long col=0; col < A11size_ptcl; col++) {
      for (long row = 0; row < A11size_ptcl; row++) {
        A11_ptcl(row,col) = PrecondMat_ptcl[col][row];
      }
    }      
    sctl::Matrix<Real> Usvd_p, VT_p, S_p, SforInv_p;
    sctl::Matrix<Real> A11forSVD_ptcl = sctl::Matrix<Real>(A11_ptcl);
    A11forSVD_ptcl.SVD(Usvd_p, S_p, VT_p);
    SforInv_p = sctl::Matrix<Real>(S_p);
    sctl::Matrix<Real> Sinv_p = SforInv_p.pinv(1e-16);

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
    };

    // Apply A11inv to each panel of a vector.
    // Look at global panel index and determine whether belongs to a particle or the channel. ASSUMES no particles are split up among processors.
    const auto AinvApply = [&Usvd,&Sinv,&VT,&A11size,&Usvd_p,&Sinv_p,&VT_p,&A11size_ptcl,&Nelem_channel,&loc_elem_cnt,&loc_elem_dsp,&Nptcl](const sctl::Vector<Real>& vec) {
        sctl::Vector<Real> AinvVec(vec.Dim());
        if (Nptcl == 0 || (loc_elem_dsp+loc_elem_cnt) <= Nelem_channel) {
            // std::cout << "All panels" << std::endl;
            // if no particles in channel or if all panels here are on channel
            sctl::Long N = vec.Dim();
            sctl::Long Npanels = N / A11size; 
            for (sctl::Long i=0; i<Npanels; i++) {
                sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
                sctl::Matrix<Real> AinvVecMat = VT.Transpose() * (Sinv * (Usvd.Transpose() * vecMat));
                for (sctl::Long j=0; j<A11size; j++) {
                    AinvVec[i*A11size + j] = AinvVecMat(j,0);
                }
            }
        } else {
            if (loc_elem_dsp >= Nelem_channel) {
                // std::cout << "All particles" << std::endl;
                // all panels here are ptcl
                sctl::Long N = vec.Dim();
                sctl::Long Nptcls = N / A11size_ptcl; 
                for (sctl::Long i=0; i<Nptcls; i++) {
                    sctl::Matrix<Real> vecMat(A11size_ptcl,1,(sctl::Iterator<Real>) vec.begin() + i*A11size_ptcl,true);
                    sctl::Matrix<Real> AinvVecMat = VT_p.Transpose() * (Sinv_p * (Usvd_p.Transpose() * vecMat));
                    for (sctl::Long j=0; j<A11size_ptcl; j++) {
                        AinvVec[i*A11size_ptcl + j] = AinvVecMat(j,0);
                    }
                }
            } else {
                sctl::Long Npanels_here = Nelem_channel - loc_elem_dsp;
                sctl::Long Nptcls_here = loc_elem_cnt - Npanels_here;
                // DEBUG
                // std::cout << "CHECK panel-ptcl split: Npanel = " << Npanels_here << ", Nptcl = " << Nptcls_here << std::endl;
                ///////////////
                for (sctl::Long i=0; i<Npanels_here; i++) {
                    sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
                    sctl::Matrix<Real> AinvVecMat = VT.Transpose() * (Sinv * (Usvd.Transpose() * vecMat));
                    for (sctl::Long j=0; j<A11size; j++) {
                        AinvVec[i*A11size + j] = AinvVecMat(j,0);
                    }
                }
                for (sctl::Long i=0; i<Nptcls_here; i++) {
                    sctl::Matrix<Real> vecMat(A11size_ptcl,1,(sctl::Iterator<Real>) vec.begin() + i*A11size_ptcl + Npanels_here*A11size,true);
                    sctl::Matrix<Real> AinvVecMat = VT_p.Transpose() * (Sinv_p * (Usvd_p.Transpose() * vecMat));
                    for (sctl::Long j=0; j<A11size_ptcl; j++) {
                        AinvVec[Npanels_here*A11size + i*A11size_ptcl + j] = AinvVecMat(j,0);
                    }
                }
            }
        }
        return AinvVec;
    };

    const auto BIO_precond = [&BIO,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> Uloc;
        BIO(&Uloc,sigma);
        // LEFT PRECONDITIONER: u -> A11inv*u
        (*U) = AinvApply(Uloc);
    };

    // first gmres to remove timing for matrix loading, and set Krylov preconditioner.
    sctl::Vector<Real> sigma_temp;
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> A11invF = AinvApply(-bg_unif_flow(X0));

    if (precond_mode == 0) {
        // no precond
        solver(&sigma_temp, BIO, -bg_unif_flow(X0), 1e0); 
    } else if (precond_mode == 1) {
        // K no A11inv
        solver(&sigma_temp, BIO, -bg_unif_flow(X0), gmres_tol, -1, false, nullptr, &krylov_precond);
    } else if (precond_mode == 2) {
        // A11inv no K
        solver(&sigma_temp, BIO_precond, A11invF, 1e0); 
    } else {
        // A11inv and K
        solver(&sigma_temp, BIO_precond, A11invF, gmres_tol, -1, false, nullptr, &krylov_precond); 
    }
    sctl::Profile::reset();

    LayerPotenOp0.ClearSetup();
    sctl::Profile::Tic("Setup");
    LayerPotenOp0.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);

    sctl::Vector<Real> sigma;
    sctl::Profile::Tic("Solver");
    if (precond_mode == 0) {
        // no precond
        solver(&sigma, BIO, -bg_unif_flow(X0), gmres_tol); 
    } else if (precond_mode == 1) {
        // K no A11inv
        solver(&sigma, BIO, -bg_unif_flow(X0), gmres_tol, -1, false, nullptr, &krylov_precond);
    } else if (precond_mode == 2) {
        // A11inv no K
        solver(&sigma, BIO_precond, A11invF, gmres_tol);
    } else {
        // A11inv and K
        solver(&sigma, BIO_precond, A11invF, gmres_tol, -1, false, nullptr, &krylov_precond);
    }
    // solver(&sigma, BIO, -bg_flow(X0), gmres_tol);
    // solver(&sigma, BIO_precond, A11invF, gmres_tol);
    
    sctl::Profile::Toc();
    sctl::Profile::print(&comm, {"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    sctl::Profile::reset();
    comm.Barrier();
    if (!comm.Rank()) {
        std::cout << "------------------- DONE WITH SOLVE ======================" << std::endl;
    }

    X0.ReInit(6);
    X0[0] = 0.2;
    X0[1] = 0.5;
    X0[2] = 0.52;
    X0[3] = 0.75;
    X0[4] = 0.5;
    X0[5] = 0.602; 
    LayerPotenOp0.SetTargetCoord(X0);
    sctl::Vector<Real> U_Xtrg;
    BIO(&U_Xtrg, sigma);
    U_Xtrg += bg_unif_flow(X0);
    std::cout << std::setprecision(12) << "U at (0.2,0,0.02) is (" << U_Xtrg[0] << ", "<< U_Xtrg[1] << ", "<< U_Xtrg[2] << "), U at (0.75,0,0.102) is (" << U_Xtrg[3] << ", "<< U_Xtrg[4] << ", "<< U_Xtrg[5] << std::endl;

    // /*
    
    { // Evaluate in interior, and write visualization
        // std::cout << "Rank " << comm.Rank()<< " calculating target points." << std::endl;
        PeriodicGeom<Real> trg;
        const sctl::Long Nelem_trg = Nelem_channel;
        const sctl::Long FourierOrder_trg = 16;
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg;
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        // peri_mode = 1 for 1-periodic
        if (channel_mode == 0) {
            std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_straight(Nelem_trg, ElemOrder, FourierOrder_trg, 0, peri_mode, 0.15, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, geom_mode);
            elem_lst_trg = std::get<0>(build_trg);
        } else if (channel_mode == 1) {
            std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_sinusoidal(Nelem_trg, ElemOrder, FourierOrder_trg, 0, peri_mode, 0.15, 0.1, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, geom_mode);
            elem_lst_trg = std::get<0>(build_trg);
        } else if (channel_mode == 2) {
            std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_sinusoidal(Nelem_trg, ElemOrder, FourierOrder_trg, 0, peri_mode, 0.1, 0.3, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, geom_mode);
            elem_lst_trg = std::get<0>(build_trg);
        } else if (channel_mode == 3) {
            std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_conv_div(Nelem_channel, ElemOrder, FourierOrder, 0, peri_mode, 0.1, 0.1, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, 1, geom_mode);
            elem_lst_trg = std::get<0>(build_trg);
        } else if (channel_mode == 4) {
            std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_spiral(Nelem_channel, ElemOrder, FourierOrder, 0, peri_mode, 0.5, 0.05, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, geom_mode);
            elem_lst_trg = std::get<0>(build_trg);
        } else if (channel_mode == 5) {
            std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_trefoil(Nelem_channel, ElemOrder, FourierOrder, 0, peri_mode, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, 1, geom_mode);
            elem_lst_trg = std::get<0>(build_trg);
        } else {
            SCTL_ASSERT(false); // not implemented
        }

        VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        sctl::Vector<Real> X0_all = vol_vis.GetCoord(); // set new target coordinates
        // std::cout << "dim of X0_all is " << X0_all.Dim() << std::endl;
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        if (!ptcls.Dim()) {
            X0 = X0_all;
            filtered_inds = 0; // all targets valid fi no particles inside.
        } else {
            std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
            X0 = std::get<0>(trg_tuple);
            filtered_inds = std::get<1>(trg_tuple);
        }
        // std::cout << "dim of X0 is " << X0.Dim() << std::endl;
        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U += bg_unif_flow(X0);
        sctl::Vector<Real> U_vis(X0_all.Dim());
        if (!ptcls.Dim()) {
            U_vis = U;
            if (channel_mode == 0) {
                // if straight channel and no particles, use reference solution to check accuracy.
                sctl::Vector<Real> U_ref = u_ref(X0_all);
                sctl::Vector<Real> err_ref = U_ref - U_vis;
                Real max_err_ref = 0.;
                for (Real e : err_ref) max_err_ref = std::max<Real>(max_err_ref, sctl::fabs(e));
                std::cout << "Rank " << comm.Rank() << " has max err " << std::setprecision(12) << max_err_ref << std::endl;
                sctl::Vector<Real> eref_loc(1);
                eref_loc[0] = max_err_ref;
                sctl::Vector<Real> eref_all(1);
                eref_all[0] = 0.;
                comm.Allreduce((sctl::Iterator<Real>) eref_loc.begin(), (sctl::Iterator<Real>) eref_all.begin(), 1, sctl::CommOp::MAX);
                if (!comm.Rank()) {
                    std::cout << "max err from reference solution is " << std::setprecision(12) << eref_all[0] << std::endl;
                }
            }
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
        // std::cout << "after Uvis" << std::endl;
        sctl::Vector<sctl::Long> size_loc(1);
        size_loc[0] = X0_all.Dim();
        sctl::Vector<sctl::Long> size_all(1);
        comm.Allreduce((sctl::Iterator<sctl::Long>) size_loc.begin(), (sctl::Iterator<sctl::Long>) size_all.begin(), 1, sctl::CommOp::SUM);
        std::cout << "rank " << comm.Rank() << " size loc = " << size_loc[0] << ", size all is " << size_all[0] << std::endl;
        // std::string filename = "ConvDiv_U_exact_"+std::to_string(comm.Rank());
        // std::string filename = "Trefoil";
        std::string filename = "Sphere_in_Straight";
        std::string filename_out = "out/"+filename+".txt";
        std::string filename_vis = "vis/"+filename;
        
        Real max_u = 0.;
        for (const auto e : U) max_u = std::max<Real>(max_u, sctl::fabs(e));
        sctl::Vector<Real> u_loc(1);
        u_loc[0] = max_u;
        sctl::Vector<Real> u_all(1);
        u_all[0] = 0.;
        comm.Allreduce((sctl::Iterator<Real>) u_loc.begin(), (sctl::Iterator<Real>) u_all.begin(), 1, sctl::CommOp::MAX);
        if (!comm.Rank()) {
            std::cout << "max u is " << std::setprecision(12) << u_all[0] << std::endl;
        }
        
        // if (write_ref) {
        //     vol_vis.WriteVTK(filename_vis, U_vis);
        //     // sctl::Vector<Real> U_vis_all(size_all[0]);
        //     // comm.Allgather((sctl::Iterator<Real>) U_vis.begin(), size_loc[0], (sctl::Iterator<Real>) U_vis_all.begin(), size_all[0]);
        //     // if (!comm.Rank()) {
        //         // U_vis_all.Write(filename_out.c_str());
        //         // std::cout << "Rank 0 finished writing." << std::endl;
        //     // }
        //     // U.Write(filename_out.c_str());
        // } else {
        //     sctl::Vector<Real> U_ref;
        //     // if (!comm.Rank()) {
        //     //   U_ref.Read(filename_out.c_str());
        //     // }
        //     // comm.PartitionN(U_ref,size_loc[0]);
        //     U_ref.Read(filename_out.c_str());
        //     const auto err = U_vis - U_ref;
        //     double max_err = 0;
        //     for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
        //     sctl::Vector<Real> err_loc(1);
        //     err_loc[0] = max_err;
        //     sctl::Vector<Real> err_all(1);
        //     err_all[0] = 0;
        //     comm.Allreduce((sctl::Iterator<sctl::Long>) err_loc.begin(), (sctl::Iterator<sctl::Long>) err_all.begin(), 1, sctl::CommOp::MAX);
        //     if (!comm.Rank()) {
        //         std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << std::endl;
        //     }
        // }
        
    }

    // */

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
    long precond_mode = std::stol(argv[7]);
    test<Real>(Nelem_channel, FourierOrder, (write_ref==1), comm, Nptcl, channel_mode, geom_mode, precond_mode);
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}

