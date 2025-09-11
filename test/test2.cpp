#include "periodize.hpp"
#include "utils.hpp"
// Test script for calculation and timing of particle-only 1 and 3 periodic problems with background pressure flow. 

/**
 * Background flow with unit pressure gradient along X-axis.
 */
template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
    const Real pdrive = 1;
    const sctl::Long N = X.Dim()/3;
    sctl::Vector<Real> U(N*3);
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        U[i*3+0] = -pdrive * ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
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

template <class Real> void test(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Integer peri_mode, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, const Real gmres_tol, const Real tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

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
    // std::cout << "periodic mode is " << peri_mode << ", Nrepeat is " << Nrepeat << std::endl;

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    if (write_ref) {
        elem_lst0.WriteVTK("vis/"+std::to_string(Nptcl)+"spheres",X0,comm);
    }  

    sctl::Vector<Real> X_proxy;
    if (peri_mode == 1) {
        X_proxy = Periodize1D<Real>::GetProxySurf();
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

    // // =============== PRECONDITIONING =======================================
    // // Store preconditioner matrix, or make new if not present.
    // std::string precond0_file = "data/precond0_ptcl_Np"+std::to_string(Nelem)+"_Nf"+std::to_string(FourierOrder)+".mat";
    // std::string precond1_file = "data/precond1_ptcl_Np"+std::to_string(Nelem)+"_Nf"+std::to_string(FourierOrder)+".mat";
    // sctl::Matrix<Real> PrecondMat0, PrecondMat1;
    // PrecondMat0.template Read<Real>(precond0_file.c_str());

    // sctl::Long A11size;

    // comm.Barrier();
    // if (PrecondMat0.Dim(0) || PrecondMat0.Dim(1)) {
    //     // std::cout << " successfully read file." << std::endl;
    //     PrecondMat1.template Read<Real>(precond1_file.c_str());
    //     A11size = PrecondMat0.Dim(1);
    // } else {
    //     sctl::Vector<sctl::Long> ptcls_pre;
    //     sctl::Vector<Real> ptcls_Xcs_pre, ptcls_rs_pre;
    //     std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_precond = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm.Self(), ptcls_pre, ptcls_rs_pre, ptcls_Xcs_pre, geom_mode);
    //     sctl::SlenderElemList<Real> elem_lst_precond = std::get<0>(build_precond);
    //     sctl::Vector<Real> X0_precond; // target coordinates
    //     elem_lst_precond.GetNodeCoord(&X0_precond, nullptr, nullptr);
    //     StokesBIO Precond_bio(SL_scal, DL_scal, comm.Self());
    //     Precond_bio.SetAccuracy(tol); // set quadrature accuracy
    //     Precond_bio.AddElemList(elem_lst_precond);
    //     Precond_bio.SetTargetCoord(X0_precond);
    //     const auto BIO_1ptcl = [&DL_scal,&Precond_bio](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    //         U->SetZero();
    //         Precond_bio.ComputePotential(*U, sigma);
    //         (*U) += sigma * 0.5 * DL_scal;
    //     };
    //     A11size = 3*ElemOrder*FourierOrder*Nelem;
    //     sctl::Vector<sctl::Vector<Real>> PrecondMat(A11size);
    //     sctl::Vector<Real> SigmaCol_precond(A11size);
    //     for (sctl::Long col=0; col < A11size; col ++) {
    //         SigmaCol_precond = 0.;
    //         SigmaCol_precond[col] = 1.;
    //         BIO_1ptcl(PrecondMat.begin() + col,SigmaCol_precond);
    //     }
    //     sctl::Matrix<Real> A11(A11size,A11size);
    //     for (long col=0; col < A11size; col++) {
    //         for (long row = 0; row < A11size; row++) {
    //             A11(row,col) = PrecondMat[col][row];
    //         }
    //     }      
    //     sctl::Matrix<Real> Usvd, VT, S, SforInv;
    //     sctl::Matrix<Real> A11forSVD = sctl::Matrix<Real>(A11);
    //     A11forSVD.SVD(Usvd, S, VT);
    //     SforInv = sctl::Matrix<Real>(S);
    //     sctl::Matrix<Real> Sinv = SforInv.pinv(1e-16);

    //     PrecondMat0 = VT.Transpose();
    //     PrecondMat1 = Sinv * Usvd.Transpose();
    //     if (!comm.Rank()) {
    //         PrecondMat0.template Write<Real>(precond0_file.c_str());
    //         PrecondMat1.template Write<Real>(precond1_file.c_str());
    //     }
    // }

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

    // // Apply A11inv to each panel of vec.
    // const auto AinvApply = [&PrecondMat0,&PrecondMat1,&A11size, &comm](const sctl::Vector<Real>& vec) {
    //     sctl::Long N = vec.Dim();
    //     sctl::Long Nptcl = N / A11size; 
    //     sctl::Vector<Real> AinvVec(N);
    //     for (sctl::Long i=0; i<Nptcl; i++) {
    //         // for each particle, apply A11inv.
    //         sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
    //         sctl::Matrix<Real> AinvVecMat = PrecondMat1 * (PrecondMat0 * vecMat);
    //         for (sctl::Long j=0; j<A11size; j++) {
    //             AinvVec[i*A11size + j] = AinvVecMat(j,0);
    //         }
    //     }
    //     return AinvVec;
    // };

    // const auto BIO_precond = [&BIO,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    //     sctl::Vector<Real> Uloc;
    //     BIO(&Uloc,sigma);
    //     // LEFT PRECONDITIONER: u -> A11inv*u
    //     (*U) = AinvApply(Uloc);
    // };

    // sctl::Vector<Real> A11invF = AinvApply( -bg_unif_flow(X0));
    // sctl::GMRES<Real> solver(comm);
    // sctl::Vector<Real> sigma;
    // solver(&sigma, BIO_precond, A11invF, gmres_tol);

    // first gmres to remove timing for matrix loading, and set Krylov preconditioner.
    sctl::Vector<Real> sigma_temp;
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    // sctl::Vector<Real> A11invF = AinvApply(-bg_unif_flow(X0));
    // solver(&sigma_temp, BIO_precond, A11invF, 1e0);
    solver(&sigma_temp, BIO, -bg_unif_flow(X0), 1e0);
    sctl::Profile::reset();

    LayerPotenOp0.ClearSetup();
    sctl::Profile::Tic("Setup");
    LayerPotenOp0.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);
    sctl::Profile::reset();

    // sctl::Profile::Tic("Solve without KrylovPrecond");
    sctl::Profile::Tic("Solve debug high residual");
    // solver(&sigma_temp, BIO_precond, A11invF, gmres_tol, -1, false);
    solver(&sigma_temp,BIO,-bg_unif_flow(X0), gmres_tol, -1, false);
    sctl::Profile::Toc();
    sctl::Profile::print(&comm, {"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    sctl::Profile::reset();
    comm.Barrier();

    sctl::Vector<Real> sigma;
    sctl::Profile::Tic("Solver: KrylovPrecond_setup");
    // solver(&sigma, BIO_precond, A11invF, gmres_tol, -1, false, nullptr, &krylov_precond);
    solver(&sigma,BIO,-bg_unif_flow(X0), gmres_tol, -1, false, nullptr, &krylov_precond);
    // solver(&sigma, BIO_precond, A11invF, gmres_tol, -1, false);
    sctl::Profile::Toc();
    sctl::Profile::print(&comm, {"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    sctl::Profile::reset();
    comm.Barrier();

    sctl::Vector<Real> sigma1;
    sctl::Profile::Tic("Solver1");
    // PRECOND with Krylov
    // solver(&sigma1, BIO_precond, A11invF, gmres_tol, -1, false, nullptr, &krylov_precond);
    solver(&sigma1,BIO,-bg_unif_flow(X0), gmres_tol, -1, false, nullptr, &krylov_precond);
    sctl::Profile::Toc();
    sctl::Profile::print(&comm, {"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    sctl::Profile::reset();
    comm.Barrier();
    if (!comm.Rank()) {
        std::cout << "------------------- DONE WITH SOLVE ======================" << std::endl;
    }

    if (write_ref) { 
        PeriodicGeom<Real> trg;    
        CubeVolumeVisShifted<Real> vol_vis(100, 1.0, comm);
        // VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        // X0 = vol_vis.GetCoord();
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);
        // std::cout << "number of target points: " << X0.Dim() << std::endl;

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U += bg_unif_flow(X0);

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
        vol_vis.WriteVTK("vis/"+std::to_string(Nptcl)+"streamlines", U_vis); 
    } 
}


int main(int argc, char** argv) {

    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Comm comm = sctl::Comm::World();
        sctl::Profile::Enable(true);
        long Nelem_ptcl = std::stol(argv[1]); // number of elements
        long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
        int write_ref = std::stoi(argv[3]);
        int peri_mode = std::stoi(argv[4]); // what kind of periodicity does the system have; peri_mode = j for j-periodic.
        long Nptcl = std::stol(argv[5]); // number of particles inside
        long geom_mode = std::stol(argv[6]); // =0: spheres; =1: spheroids; =3: bacteria; =4: loop.
        double gmres_tol = std::stod(argv[7]);
        double tol = std::stod(argv[8]);

        test<Real>(Nelem_ptcl, FourierOrder, (write_ref==1), peri_mode, comm, Nptcl, geom_mode, gmres_tol, tol);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
