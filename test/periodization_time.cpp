#include "utils.hpp"

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

template <class Real> void timing_3peri(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Comm comm, sctl::Long Nptcl, const Real tol, const Real gmres_tol) {

    // std::cout << "Nelem = " << Nelem << ", Fourier = " << FourierOrder << std::endl;

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    // const Real tol = 1e-14;
    // Real gmres_tol = 1e-10; 
    const sctl::Long ElemOrder = 10;
    const Real period_length = 1.;
    const Real pressure_drop = -1.;
    // const sctl::Long Nptcl = 25;
    const sctl::Long geom_mode = 0;
    const sctl::Long max_gmres_iter = 200;

    // if (FourierOrder < 20) {
    //     gmres_tol = 1e-6;
    // } else if (Nelem < 4) {
    //     gmres_tol = 1e-8;
    // } else if (FourierOrder < 36) {
    //     gmres_tol = 1e-10;
    // } else {
    //     gmres_tol = 1e-12;
    // }

    // if (FourierOrder < 20) {
    //     gmres_tol = 1e-6;
    // } else if (Nelem < 4) {
    //     gmres_tol = 1e-8;
    // } else if (FourierOrder < 36) {
    //     gmres_tol = 1e-10;
    // } else {
    //     gmres_tol = 1e-12;
    // }
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Vector<Real> X0surf = X0;
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

    // sctl::Profile::Tic("LPO 3 periodic");
    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length);

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
        }
        
        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma0.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer
    
        AddConstVec(*U, sigma_mean);
    };

    const auto BIO_precond = [&BIO,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> Uloc;
        BIO(&Uloc,sigma);
        // LEFT PRECONDITIONER: u -> A11inv*u
        (*U) = AinvApply(Uloc);
    };

    const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
        sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
        AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

        sctl::Vector<Real> U0;
        LayerPotenOp0.ComputeSL(U0, force_density);
        return U0;
    };

    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> A11invF = AinvApply(eval_rhs(pressure_drop)); 

    LayerPotenOp0.ClearSetup();
    sctl::Profile::Tic("3-periodic LPO Setup");
    LayerPotenOp0.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);

    sctl::Vector<Real> sigma_temp;
    solver(&sigma_temp, BIO_precond, A11invF, 1e-4, max_gmres_iter); // first run to remove high timing counts
    sctl::Profile::reset();

    sctl::Vector<Real> sigma;
    sctl::Profile::Tic("3-periodic LPO solve with preconds");
    solver(&sigma, BIO_precond, A11invF, gmres_tol, max_gmres_iter, false, nullptr, &krylov_precond);
    // solver(&sigma,BIO,field_on_surf,gmres_tol, -1, false, nullptr, &krylov_precond); 
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);
    sctl::Profile::reset();

    PeriodicGeom<Real> trg;    
    Real fact = Nptcl / 25;
    int Ngrid = (int) fact * 20;
    std::cout << "New vis grid is " << Ngrid << std::endl;
    CubeVolumeVisShifted<Real> vol_vis(Ngrid, 0.95, comm);
    sctl::Vector<Real> X0_all = vol_vis.GetCoord();
    sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    X0 = std::get<0>(trg_tuple);
    filtered_inds = std::get<1>(trg_tuple);
    std::cout << "3-periodic Number of targets after filtering: " << X0.Dim() << std::endl;

    LayerPotenOp0.SetTargetCoord(X0);

    LayerPotenOp0.ClearSetup();
    sctl::Profile::Tic("3-periodic LPO Eval Setup");
    LayerPotenOp0.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);

    sctl::Vector<Real> Utemp;
    BIO(&Utemp, sigma);
    Utemp -= eval_rhs(pressure_drop);

    const sctl::Long Ncp = 10;

    sctl::Profile::reset();
    sctl::Vector<Real> U;
    sctl::Profile::Tic("3-periodic LPO eval");
    for (int cp=0; cp < Ncp; cp++) {
        BIO(&U, sigma);
        U -= eval_rhs(pressure_drop);  
    }
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);
    sctl::Profile::reset();

    // sctl::Profile::Tic("3-periodic LPO eval_rhs");
    // sctl::Profile::Toc();
    // sctl::Profile::print(&comm);
    // sctl::Profile::reset();


    std::cout << "=========== FREE SPACE LPO =============" << std::endl;

    StokesBIO LayerPotenOp1(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp1.AddElemList(elem_lst0);
    LayerPotenOp1.SetTargetCoord(X0surf);
    LayerPotenOp1.SetAccuracy(tol);
    const auto BIO_1 = [&DL_scal,&LayerPotenOp1,NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        U->SetZero();
        LayerPotenOp1.ComputePotential(*U, sigma);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma*0.5*NormalOrient * DL_scal; // for double-layer
    };

    const auto BIO_precond_1 = [&BIO_1,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> Uloc;
        BIO_1(&Uloc,sigma);
        // LEFT PRECONDITIONER: u -> A11inv*u
        (*U) = AinvApply(Uloc);
    };

    // const auto bg_flow = [](const sctl::Vector<Real>& X) {
    //     const sctl::Long N = X.Dim()/3;
    //     sctl::Vector<Real> U(N*3);
    //     for (sctl::Long i = 0; i < N; i++) {
    //         const auto x = X.begin() + i*3;
    //         U[i*3+0] = - ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
    //         U[i*3+1] = 0;
    //         U[i*3+2] = 0;
    //     }
    //     return U;
    // };

    const auto eval_rhs_1 = [&LayerPotenOp1,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
        sctl::Vector<Real> force_density(LayerPotenOp1.Dim(0)); force_density = 0;
        AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

        sctl::Vector<Real> U0;
        LayerPotenOp1.ComputeSL(U0, force_density);
        return U0;
    };

    LayerPotenOp1.ClearSetup();
    sctl::Profile::Tic("free space LPO Setup");
    LayerPotenOp1.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);
    sctl::Profile::reset();

    sctl::GMRES<Real> solver_1(comm);
    sctl::KrylovPrecond<Real> krylov_precond_1;
    // sctl::Vector<Real> A11invF_1 = AinvApply(bg_flow(X0surf) * (pressure_drop) / (period_length));
    sctl::Vector<Real> A11invF_1 = AinvApply(eval_rhs_1(pressure_drop));

    sctl::Vector<Real> sigma_1_temp;
    solver_1(&sigma_1_temp, BIO_precond_1, A11invF_1, 1e-4, max_gmres_iter);
    sctl::Profile::reset();

    sctl::Vector<Real> sigma_1;
    sctl::Profile::Tic("free space LPO solve with preconds");
    solver_1(&sigma_1, BIO_precond_1, A11invF_1, gmres_tol, max_gmres_iter, false, nullptr, &krylov_precond_1);
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);
    sctl::Profile::reset();

    std::cout << "free space Number of targets (double checking size doesn't change): " << X0.Dim() << std::endl;

    LayerPotenOp1.SetTargetCoord(X0);

    LayerPotenOp1.ClearSetup();
    sctl::Profile::Tic("free space LPO Eval Setup");
    LayerPotenOp1.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);

    sctl::Vector<Real> U_1temp;
    BIO_1(&U_1temp, sigma_1);
    U_1temp -= eval_rhs_1(pressure_drop);

    sctl::Profile::reset();
    sctl::Profile::Tic("free space LPO eval");
    sctl::Vector<Real> U_1;
    for (int cp=0; cp<Ncp; cp++) {
        BIO_1(&U_1, sigma_1);
        U_1 -= eval_rhs_1(pressure_drop);
    }
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);
    sctl::Profile::reset();
    // sctl::Profile::Tic("free space LPO eval_rhs()");
    // // U_1 -= bg_flow(X0) * (pressure_drop) / (period_length); 
    // U_1 -= eval_rhs_1(pressure_drop);
    // sctl::Profile::Toc();
    // sctl::Profile::print(&comm);
    // sctl::Profile::reset();

}


int main(int argc, char** argv) {
    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Profile::Enable(true);
        sctl::Comm comm = sctl::Comm::World();
        long Nelem_ptcl = std::stol(argv[1]); // number of elements
        long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
        long Nptcl = std::stol(argv[3]);
        double gmres_tol = std::stod(argv[4]);
        double tol = std::stod(argv[5]);

        timing_3peri<Real>(Nelem_ptcl, FourierOrder, comm, Nptcl, tol, gmres_tol);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
