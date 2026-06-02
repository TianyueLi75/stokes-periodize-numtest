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
        U[i*3+0] = - ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
        U[i*3+1] = 0;
        U[i*3+2] = 0;
    }
    return U;
}

template <class Real> void SurfaceIntegral(sctl::Vector<Real>& I, const sctl::Vector<Real>& vals, const sctl::Vector<Real>& wts) {
  const sctl::Long dof = vals.Dim() / wts.Dim();
  SCTL_ASSERT(vals.Dim() == wts.Dim() * dof);
  if (I.Dim() != dof) I.ReInit(dof);
  I = 0;
  for (sctl::Long i = 0; i < wts.Dim(); i++) {
    for (sctl::Long j = 0; j < dof; j++) {
      I[j] += vals[i*dof + j] * wts[i];
    }
  }
}

template <class Real> void AddConstVec(sctl::Vector<Real>& vals, const sctl::Vector<Real>& c0) {
  const sctl::Long dof = c0.Dim();
  const sctl::Long N = vals.Dim() / dof;
  SCTL_ASSERT(vals.Dim() == N * dof);
  for (sctl::Long i = 0; i < N; i++) {
    for (sctl::Long j = 0; j < dof; j++) {
      vals[i*dof + j] += c0[j];
    }
  }
}

// Loop over copies and add consecutively, to reduce memory requirements. Perhaps do 2D planes at a time.
template <class Real> sctl::Vector<Real> exact_field(const sctl::Vector<Real>& Xtrg, const sctl::Vector<Real>& Xsrc, const sctl::Vector<Real>& sigma, const sctl::Long Ncopy) {
    sctl::Stokes3D_FxU ker;
  
    const sctl::Long N = Xtrg.Dim()/3;
    sctl::Vector<Real> U(N*3);
    U = 0.;
    PeriodicGeom<Real> obj;

    sctl::Vector<Real> Xsrc_ = obj.X_nbr_copy(Xsrc,Ncopy,1);
    sctl::Vector<Real> sigma_ = obj.vec_nbr_copy(sigma,Ncopy,1);
    ker.Eval(U,Xtrg,Xsrc_,Xsrc_,sigma_);

    return U;
}

// First check that two 
template <class Real> void exact_field_check(sctl::Comm comm, sctl::Long Nptcl, const sctl::Long geom_mode, const sctl::Long Ncopy1, const sctl::Long Ncopy2) {
    const sctl::Long Nelem = 4;
    const sctl::Long FourierOrder = 64;
    const sctl::Long ElemOrder = 10;
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else if (Nptcl == 3) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls3(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else { 
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    }
    Nptcl = ptcls_rs.Dim(); 
    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);

    sctl::Long Ncharge;
    if (Nptcl < 150) {
        // two equal and opposite charges per particle
        Ncharge = 2*Nptcl;
    } else {
        // only first 150 particles get charges inside. (arbitrary, to limit true solution timing)
        Ncharge = 2*150;
    }
    // Currently one Stokeslet doublet per particle (for a simple net-force-zero scenario)
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

    PeriodicGeom<Real> trg;    
    CubeVolumeVisShifted<Real> vol_vis(10, 0.95, comm);
    X0 = vol_vis.GetCoord();
    
    sctl::Vector<Real> field_on_surf_1 = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy1);
    sctl::Vector<Real> field_on_surf_2 = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy2);
    sctl::Vector<Real> diff = field_on_surf_1 - field_on_surf_2;
    Real max_err = 0.;
    Real max_field1 = 0.;
    for (const auto d : diff) max_err = std::max(std::abs(d), max_err);
    for (const auto d : field_on_surf_1) max_field1 = std::max(std::abs(d), max_field1);
    Real max_rel_err = max_err / max_field1;

    std::cout << "max relative error between Ncopy1 = " << Ncopy1 << " and " << Ncopy2 << " is " << std::setprecision(15) << max_rel_err << std::endl;
}

template <class Real> void manufactured_soln_1peri(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, sctl::Long Ncopy) {

    std::cout << "Nelem = " << Nelem << ", Fourier = " << FourierOrder << std::endl;

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    Real tol = 1e-14;
    Real gmres_tol = 1e-10; 
    const sctl::Long ElemOrder = 10;
    const Real period_length = 1.;

    if (FourierOrder < 20) {
        gmres_tol = 1e-6;
    } else if (Nelem < 4) {
        gmres_tol = 1e-8;
    } else if (FourierOrder < 36) {
        gmres_tol = 1e-10;
    } else {
        gmres_tol = 1e-12;
    }

    if (FourierOrder < 20) {
        gmres_tol = 1e-6;
    } else if (Nelem < 4) {
        gmres_tol = 1e-8;
    } else if (FourierOrder < 36) {
        gmres_tol = 1e-10;
    } else {
        gmres_tol = 1e-12;
    }
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else if (Nptcl == 3) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls3(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else { 
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    }
    Nptcl = ptcls_rs.Dim(); 
    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    Real surface_area;
    sctl::Vector<Real> wts;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        // MPI
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area = sa_all[0];
        // surface_area = surface_area_[0];
    }
    std::cout << "DEBUG surface area = " << surface_area << "." << std::endl;

    // Create point charges at random locations close to particle center, by a distance of at most 0.2r.
    sctl::Long Ncharge;
    if (Nptcl < 150) {
        // two equal and opposite charges per particle
        Ncharge = 2*Nptcl;
    } else {
        // only first 150 particles get charges inside. (arbitrary, to limit true solution timing)
        Ncharge = 2*150;
    }
    // Ncharge = 2;
    // Currently one Stokeslet doublet per particle (for a simple net-force-zero scenario)
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
    // std::cout << "Computing exact field: " << std::endl;
    sctl::Vector<Real> field_on_surf = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy);
    // std::cout << "Done with exact field: " << std::endl;

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);

    // periodized layer potential operator
    const auto BIO = [&wts,&surface_area,&elem_lst0,&DL_scal,&LayerPotenOp0,&X0,NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> sigma_mean, sigma0;
        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_;
            elem_lst0.GetFarFieldDensity(sigma_, sigma);
            SurfaceIntegral(sigma_mean, sigma_, wts);
            //MPI
            sctl::Vector<Real> sa_loc = sigma_mean;
            sctl::Vector<Real> sa_all(3);
            sa_all = 0;
            comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
            comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+1, (sctl::Iterator<Real>) sa_all.begin()+1, 1, sctl::CommOp::SUM);
            comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+2, (sctl::Iterator<Real>) sa_all.begin()+2, 1, sctl::CommOp::SUM);
            sigma_mean = sa_all;
            sigma_mean *= (1/surface_area);

            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);

            // // DEBUG: check that sigma-sigma_mean has surface integral = 0:
            // sctl::Vector<Real> sigma1 = sigma_;
            // AddConstVec(sigma1, -sigma_mean);
            // sctl::Vector<Real> sigma_test_;
            // SurfaceIntegral(sigma_test_, sigma1, wts);
            // std::cout << "Surface integral of sigma - sigma bar = " << sigma_test_[0] << ", "<< sigma_test_[1] << ", " << sigma_test_[2] << ". "<< std::endl;
        
        }
        
        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma0.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer
    
        AddConstVec(*U, sigma_mean);
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
        std::cout << " successfully read file." << std::endl;
        PrecondMat1.template Read<Real>(precond1_file.c_str());
        A11size = PrecondMat0.Dim(1);
    } else {
        std::cout << " Making precond files " << std::endl;
        sctl::Vector<sctl::Long> ptcls_pre;
        sctl::Vector<Real> ptcls_Xcs_pre, ptcls_rs_pre;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_precond = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, comm.Self(), ptcls_pre, ptcls_rs_pre, ptcls_Xcs_pre, geom_mode);
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
            sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
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

    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> A11invF = AinvApply(field_on_surf);

    sctl::Vector<Real> sigma;
    solver(&sigma, BIO_precond, A11invF, gmres_tol, -1, false, nullptr, &krylov_precond);
    // solver(&sigma,BIO,field_on_surf,gmres_tol, -1, false, nullptr, &krylov_precond); 

    PeriodicGeom<Real> trg;    
    CubeVolumeVisShifted<Real> vol_vis(5, 0.95, comm);
    sctl::Vector<Real> X0_all = vol_vis.GetCoord();
    sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    X0 = std::get<0>(trg_tuple);
    filtered_inds = std::get<1>(trg_tuple);
    // std::cout<< "Number of targets before filter: " << X0_all.Dim() << ", after filter: " << X0.Dim() << std::endl;
    LayerPotenOp0.SetTargetCoord(X0);
    sctl::Vector<Real> U;
    BIO(&U, sigma);

    sctl::Vector<Real> field_on_trg = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy);
    // get max abs error
    sctl::Vector<Real> err = U - field_on_trg;

    // // DEBUGGING: just print the errors to check whether x,y,z dependence, constant, etc.
    // for (int i=0; i<err.Dim()/3; i++) {
    //     std::cout << "err: " << std::setprecision(10) << err[i*3+0] << ", " << err[i*3+1] << ", " << err[i*3+2] << ". " << std::endl;
    // }

    // Subtract mean to remove constant difference
    sctl::Vector<Real> sum_err(3);
    sum_err = 0.;
    for (sctl::Long i=0; i<err.Dim()/3; i++) {
        for (sctl::Long k=0; k<3; k++) {
            sum_err[k] += err[i*3+k];
        }
    }
    sctl::Long err_size = err.Dim()/3;
    //MPI
    sctl::Vector<Real> sum_err_loc = sum_err;
    sctl::Vector<Real> sum_err_all(3);
    sum_err_all = 0;
    comm.Allreduce((sctl::Iterator<Real>) sum_err_loc.begin(), (sctl::Iterator<Real>) sum_err_all.begin(), 1, sctl::CommOp::SUM);
    comm.Allreduce((sctl::Iterator<Real>) sum_err_loc.begin()+1, (sctl::Iterator<Real>) sum_err_all.begin()+1, 1, sctl::CommOp::SUM);
    comm.Allreduce((sctl::Iterator<Real>) sum_err_loc.begin()+2, (sctl::Iterator<Real>) sum_err_all.begin()+2, 1, sctl::CommOp::SUM);
    sum_err = sum_err_all;

    sctl::Vector<sctl::Long> err_size_loc(1);
    err_size_loc[0] = err_size;
    sctl::Vector<sctl::Long> err_size_all(1); 
    err_size_all[0] = 0;
    comm.Allreduce((sctl::Iterator<Real>) err_size_loc.begin(), (sctl::Iterator<Real>) err_size_all.begin(), 1, sctl::CommOp::SUM);
    // avg err
    sctl::Vector<Real> avg_err = sum_err / err_size_all[0];
    AddConstVec(err,-avg_err); // relative error with offset: max ((Ucalc - C) - Uexact) / Uexact, since C = Ucalc_exact - Uexact ~ E[Ucalc - Uexact]
    // std::cout << "avg err: " << avg_err[0] << ", " << avg_err[1] << ", " << avg_err[2] << std::endl;
    // for (int i=0; i<err.Dim()/3; i++) {
    //     std::cout << "err after subtracting avg err: " << std::setprecision(10) << err[i*3+0] << ", " << err[i*3+1] << ", " << err[i*3+2] << ". " << std::endl;
    // }

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
        std::cout<<"Max error = "<< std::setprecision(15) << err_all[0] << ", Max u = " << u_all[0] << ", Max relative error = " << err_all[0] / u_all[0] << std::endl;
    }   
}


int main(int argc, char** argv) {
    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        // sctl::Profile::Enable(true);
        sctl::Comm comm = sctl::Comm::World();
        long Nelem_ptcl = std::stol(argv[1]); // number of elements
        long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
        long Nptcl = std::stol(argv[3]); // number of particles inside
        long geom_mode = std::stol(argv[4]); // =0: spheres; =1: spheroids; =3: bacteria; =4: loop.
        long Ncopy = std::stol(argv[5]);

        manufactured_soln_1peri<Real>(Nelem_ptcl, FourierOrder, comm, Nptcl, geom_mode, Ncopy);
        // exact_field_check<Real>(comm, Nptcl, geom_mode, 60000, 70000);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
