#include "periodize.hpp"
#include "utils.hpp"

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

template <class Real> void test(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Integer peri_mode, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, sctl::Long Ncopy) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    Real tol = 1e-14;
    const Real gmres_tol = 1e-7; // tolerances set up to give 6 digts of accuracy.
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
    LayerPotenOp_proxy.SetTargetCoord(X_proxy);
    LayerPotenOp_proxy.SetAccuracy(tol);

            // periodized layer potential operator
    const auto BIO = [&DL_scal,&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat,NormalOrient,&peri_mode](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
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
                Periodize1D<Real>::EvalFarField(U_far, X0, U_proxy,30,20);
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

    // =============== PRECONDITIONING =======================================
    // Store preconditioner matrix, or make new if not present.
    std::string precond0_file = "data/precond0_ptcl_Np"+std::to_string(Nelem)+"_Nf"+std::to_string(FourierOrder)+".mat";
    std::string precond1_file = "data/precond1_ptcl_Np"+std::to_string(Nelem)+"_Nf"+std::to_string(FourierOrder)+".mat";
    sctl::Matrix<Real> PrecondMat0, PrecondMat1;
    PrecondMat0.template Read<Real>(precond0_file.c_str());

    sctl::Long A11size;

    comm.Barrier();
    if (PrecondMat0.Dim(0) || PrecondMat0.Dim(1)) {
        // std::cout << " successfully read file." << std::endl;
        PrecondMat1.template Read<Real>(precond1_file.c_str());
        A11size = PrecondMat0.Dim(1);
    } else {
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
        A11size = 3*ElemOrder*FourierOrder*Nelem;
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

        PrecondMat0 = VT.Transpose();
        PrecondMat1 = Sinv * Usvd.Transpose();
        if (!comm.Rank()) {
            PrecondMat0.template Write<Real>(precond0_file.c_str());
            PrecondMat1.template Write<Real>(precond1_file.c_str());
        }
    }
    

    // Apply A11inv to each panel of vec.
    const auto AinvApply = [&PrecondMat0,&PrecondMat1,&A11size, &comm](const sctl::Vector<Real>& vec) {
        sctl::Long N = vec.Dim();
        sctl::Long Nptcl = N / A11size; 
        sctl::Vector<Real> AinvVec(N);
        for (sctl::Long i=0; i<Nptcl; i++) {
            // for each particle, apply A11inv.
            sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
            sctl::Matrix<Real> AinvVecMat = PrecondMat1 * (PrecondMat0 * vecMat);
            for (sctl::Long j=0; j<A11size; j++) {
                AinvVec[i*A11size + j] = AinvVecMat(j,0);
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
    sctl::Vector<Real> A11invF = AinvApply(field_on_surf);
    // PRECOND with Krylov
    solver(&sigma_temp, BIO_precond, A11invF, gmres_tol, -1, false, nullptr, &krylov_precond);
    sctl::Profile::reset();

    LayerPotenOp0.ClearSetup();
    sctl::Profile::Tic("Setup");
    LayerPotenOp0.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);

    sctl::Vector<Real> sigma;
    sctl::Profile::Tic("Solver");
    // PRECOND with Krylov
    solver(&sigma, BIO_precond, A11invF, gmres_tol, -1, false, nullptr, &krylov_precond);
    sctl::Profile::Toc();
    sctl::Profile::print(&comm, {"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    sctl::Profile::reset();
    comm.Barrier();
    if (!comm.Rank()) {
        std::cout << "------------------- DONE WITH SOLVE ======================" << std::endl;
    }
    if (peri_mode ==1) { 
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

        sctl::Vector<Real> field_on_trg = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy, peri_mode);
        // get max abs error
        const auto err = U - field_on_trg;
        double max_err = 0;
        Real max_u = 0.;
        for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
        for (const auto e : field_on_trg) max_u = std::max<Real>(max_u, sctl::fabs(e));
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
            std::cout<<"Max error = "<< std::setprecision(15) << err_all[0] << std::endl;
            std::cout<<"Max u = "<< std::setprecision(15) << u_all[0] << std::endl;
            std::cout<<"Max relative error = "<< std::setprecision(15) << err_all[0] / u_all[0] << std::endl;
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
        int peri_mode = std::stoi(argv[3]); // what kind of periodicity does the system have; peri_mode = j for j-periodic.
        long Nptcl = std::stol(argv[4]); // number of particles inside
        long geom_mode = std::stol(argv[5]); // =0: spheres; =1: spheroids; =3: bacteria; =4: loop.
        long Ncopy = std::stol(argv[6]);

        test<Real>(Nelem_ptcl, FourierOrder, peri_mode, comm, Nptcl, geom_mode, Ncopy);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
