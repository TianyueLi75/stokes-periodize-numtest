/**
 * Background flow with unit pressure gradient along X-axis.
 */
template <class Real> sctl::Vector<Real> bg_flow_1peri(const sctl::Vector<Real>& X) {
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

template <class Real> sctl::Vector<Real> bg_flow_2peri(const sctl::Vector<Real>& X) {
    const sctl::Long N = X.Dim()/3;
    sctl::Vector<Real> U(N*3);
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        U[i*3+0] = - 0.5 * ((x[2]-0.5)*(x[2]-0.5));
        U[i*3+1] = 0;
        U[i*3+2] = 0;
    }
    return U;
}



template <class Real> sctl::Vector<Real> GetVslip(const sctl::Vector<Real>& ptcls_Xnsurf, const sctl::Vector<Real>& ptcls_Xcs, const sctl::Vector<Real>& ptcls_sizes, const sctl::Vector<Real>& ptcls_u0s, const sctl::Vector<Real>& ptcls_thetas, const sctl::Vector<Real>& ptcls_phis, const sctl::Vector<sctl::Long>& ptcls_ifprolate, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder) {
    const sctl::Long Nnodes_per_ptcl = Nelem * ElemOrder * FourierOrder;
    const sctl::Long Nptcls = ptcls_sizes.Dim();

    sctl::Vector<Real> Vslip(ptcls_Xnsurf.Dim());
    srand48(1);
    Vslip.SetZero();
    for (sctl::Long i=0; i<Nptcls; i++) {
        const Real a_here = ptcls_sizes[i];
        const Real u0_here = ptcls_u0s[i];
        const int if_prolate_here = ptcls_ifprolate[i];
        sctl::Vector<Real> Vslip_here(3*Nnodes_per_ptcl, (sctl::Iterator<Real>) Vslip.begin() + 3*Nnodes_per_ptcl*i, false);
        sctl::Vector<Real> center_here(3, (sctl::Iterator<Real>) ptcls_Xcs.begin() + 3*i, false);
        const Real scalar = drand48()*0.8 + 0.1; // randomly scaled slip velocity by (0.1,0.9).
        const Real theta_rotate = ptcls_thetas[i];
        const Real phi_rotate = ptcls_phis[i];
        const Real cos_theta_rotate = sctl::cos<Real>(theta_rotate);
        const Real sin_theta_rotate = sctl::sin<Real>(theta_rotate);
        const Real cos_phi_rotate = sctl::cos<Real>(phi_rotate);
        const Real sin_phi_rotate = sctl::sin<Real>(phi_rotate);

        for (sctl::Long panel=0; panel < Nelem; panel++) {
            const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrder);

            for (sctl::Long el=0; el<ElemOrder; el++) {
                const Real theta = sctl::const_pi<Real>() * (panel+nodes[el])/Nelem;
                const Real coeff = scalar * sctl::sin<Real>(theta); // sinusoidal slip magnitude, 0 at north and south poles

                for (sctl::Long fl=0; fl<FourierOrder; fl++) {
                    Real phi = 2. * sctl::const_pi<Real>() * fl / FourierOrder;
                    const sctl::Long idx = panel*ElemOrder*FourierOrder + el*FourierOrder + fl;
                    const sctl::Vector<Real> Xn_here(3, (sctl::Iterator<Real>) ptcls_Xnsurf.begin() + 3*Nnodes_per_ptcl*i + idx*3, false);

                    Real t1_unrotated = a_here * u0_here * (-sctl::sin<Real>(theta));
                    Real t2_unrotated, t3_unrotated;
                    if (if_prolate_here) {
                        t2_unrotated = a_here * sctl::sqrt<Real>(u0_here*u0_here-1) * sctl::cos<Real>(theta) * sctl::cos<Real>(phi);
                        t3_unrotated = a_here * sctl::sqrt<Real>(u0_here*u0_here-1) * sctl::cos<Real>(theta) * sctl::sin<Real>(phi);
                    } else {
                        t2_unrotated = a_here * sctl::sqrt<Real>(u0_here*u0_here+1) * sctl::cos<Real>(theta) * sctl::cos<Real>(phi);
                        t3_unrotated = a_here * sctl::sqrt<Real>(u0_here*u0_here+1) * sctl::cos<Real>(theta) * sctl::sin<Real>(phi);
                    }
                    
                    Real mag2_tang = t1_unrotated*t1_unrotated + t2_unrotated*t2_unrotated + t3_unrotated*t3_unrotated;
                    Real mag_tang = sctl::sqrt<Real>(mag2_tang);
                    // tangential sinusoidal slip, spheroid-cetnered space
                    sctl::Vector<Real> Vslip_here_unrotated(3);
                    Vslip_here_unrotated[0] = coeff * t1_unrotated / mag_tang;
                    Vslip_here_unrotated[1] = coeff * t2_unrotated / mag_tang;
                    Vslip_here_unrotated[2] = coeff * t3_unrotated / mag_tang;
                    // tagential slip, rotated
                    sctl::Vector<Real> Vslip_here_rotated(3);
                    Vslip_here_rotated[0] = cos_phi_rotate * cos_theta_rotate * Vslip_here_unrotated[0] - sin_phi_rotate * Vslip_here_unrotated[1] + cos_phi_rotate * sin_theta_rotate * Vslip_here_unrotated[2];
                    Vslip_here_rotated[1] = sin_phi_rotate * cos_theta_rotate * Vslip_here_unrotated[0] + cos_phi_rotate * Vslip_here_unrotated[1] + sin_phi_rotate * sin_theta_rotate * Vslip_here_unrotated[2];
                    Vslip_here_rotated[2] = - sin_theta_rotate * Vslip_here_unrotated[0] + cos_theta_rotate * Vslip_here_unrotated[2];

                    Real vdotn = Vslip_here_rotated[0] * Xn_here[0] + Vslip_here_rotated[1] * Xn_here[1] + Vslip_here_rotated[2] * Xn_here[2];
                    if (sctl::fabs(vdotn)>1e-8) {
                        std::cout << "ERROR: tangent dot n is nonzero: " << vdotn << std::endl;
                    }
                    Vslip_here[idx*3+0] = Vslip_here_rotated[0];
                    Vslip_here[idx*3+1] = Vslip_here_rotated[1];
                    Vslip_here[idx*3+2] = Vslip_here_rotated[2];
                }
            }
        }
    }   
    return Vslip;
}


/*
    Channel_radius for different geometries implemented in utils.cpp:
        straight: (param) 0.2
        sinusoidal: (param) 0.1
        conv div: (param, avg) 0.15
        spiral: (param) 0.05
        trefoil: (fixed) 0.035 
*/
// Note: physical length of cylinder precond panel is based on the number of elements on the whole channel, but the discretization only has one panel on the precond cylinder.
template <class Real> sctl::Long precond_channel(sctl::Matrix<Real>& PrecondMat0, sctl::Matrix<Real>& PrecondMat1, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const Real channel_radius, const Real SL_scal, const Real DL_scal, sctl::Comm comm) {
    // Store preconditioner matrix, or make new if not present.
    std::string precond0_file = "data/precond0_cyln_Np"+std::to_string(Nelem)+"_Nf"+std::to_string(FourierOrder)+".mat";
    std::string precond1_file = "data/precond1_cyln_Np"+std::to_string(Nelem)+"_Nf"+std::to_string(FourierOrder)+".mat";
    PrecondMat0.template Read<Real>(precond0_file.c_str());

    const Real tol = 1e-15; // this may be different from the tol with which some of the earlier files were made.

    sctl::Long A11size;

    comm.Barrier();
    if (PrecondMat0.Dim(0) || PrecondMat0.Dim(1)) {
        std::cout << " successfully read preconditioner file for channel (N_p = " << Nelem << ", N_f = " << FourierOrder << "): " << precond0_file << std::endl;
        PrecondMat1.template Read<Real>(precond1_file.c_str());
        A11size = PrecondMat0.Dim(1);
    } else {
        std::cout << " Making preconditioner file for channel (N_p = " << Nelem << ", N_f = " << FourierOrder << ")." << std::endl;
        sctl::Vector<Real> Xc_precond, eps_precond; 
        sctl::Vector<sctl::Long> ElemOrderVec_precond(1), FourierOrderVec_precond(1);
        ElemOrderVec_precond[0] = ElemOrder;
        FourierOrderVec_precond[0] = FourierOrder;
        const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrder);
        for (sctl::Long j = 0; j < ElemOrder; j++) { // loop over panel nodes
            const Real x = (nodes[j]) / Nelem; // size of precond panel should be same as one panel on pipe
            Xc_precond.PushBack(x+0.5); //  shift panel to center of unit box, arbitrary.
            Xc_precond.PushBack(0.5); 
            Xc_precond.PushBack(0.5); 
            eps_precond.PushBack(channel_radius); 
        }
        sctl::SlenderElemList<Real> elem_lst_precond(ElemOrderVec_precond, FourierOrderVec_precond, Xc_precond, eps_precond);
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
        A11size = 3*ElemOrder*FourierOrder;
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
        sctl::Matrix<Real> Sinv = SforInv.pinv(tol);

        PrecondMat0 = VT.Transpose();
        PrecondMat1 = Sinv * Usvd.Transpose();
        if (!comm.Rank()) {
            PrecondMat0.template Write<Real>(precond0_file.c_str());
            PrecondMat1.template Write<Real>(precond1_file.c_str());
        }
    }

    return A11size;
}

template <class Real> sctl::Long precond_ptcl(sctl::Matrix<Real>& PrecondMat0, sctl::Matrix<Real>& PrecondMat1, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const Real SL_scal, const Real DL_scal, sctl::Comm comm) {
    // Store preconditioner matrix, or make new if not present.
    std::string precond0_file = "data/precond0_ptcl_Np"+std::to_string(Nelem)+"_Nf"+std::to_string(FourierOrder)+".mat";
    std::string precond1_file = "data/precond1_ptcl_Np"+std::to_string(Nelem)+"_Nf"+std::to_string(FourierOrder)+".mat";
    PrecondMat0.template Read<Real>(precond0_file.c_str());

    const Real tol = 1e-15; // this may be different from the tol with which some of the earlier files were made.

    sctl::Long A11size;

    comm.Barrier();
    if (PrecondMat0.Dim(0) || PrecondMat0.Dim(1)) {
        std::cout << " successfully read file " << precond0_file << std::endl;
        PrecondMat1.template Read<Real>(precond1_file.c_str());
        A11size = PrecondMat0.Dim(1);
    } else {
        std::cout << " Making precond files " << std::endl;
        PeriodicGeom<Real> obj;
        sctl::Vector<sctl::Long> ptcls_pre;
        sctl::Vector<Real> ptcls_Xcs_pre, ptcls_rs_pre;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_precond = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, comm.Self(), ptcls_pre, ptcls_rs_pre, ptcls_Xcs_pre, 0);
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
        sctl::Matrix<Real> Sinv = SforInv.pinv(tol);

        PrecondMat0 = VT.Transpose();
        PrecondMat1 = Sinv * Usvd.Transpose();
        if (!comm.Rank()) {
            PrecondMat0.template Write<Real>(precond0_file.c_str());
            PrecondMat1.template Write<Real>(precond1_file.c_str());
        }
    }

    return A11size;
}

/**
 Differences between solutions possibly a constant, so remove mean from error before taking max relative errors.
*/
template <class Real> void renormalize_error(sctl::Vector<Real>& err, sctl::Comm comm) {
    sctl::Long Nnodes = err.Dim()/3;
    // Subtract mean to remove constant difference
    sctl::Vector<Real> sum_err(3);
    sum_err = 0.;
    for (sctl::Long i=0; i<Nnodes; i++) {
        for (sctl::Long k=0; k<3; k++) {
            sum_err[k] += err[i*3+k];
        }
    }
    //MPI
    sctl::Vector<Real> sum_err_loc = sum_err;
    sctl::Vector<Real> sum_err_all(3);
    sum_err_all = 0;
    comm.Allreduce((sctl::Iterator<Real>) sum_err_loc.begin(), (sctl::Iterator<Real>) sum_err_all.begin(), 1, sctl::CommOp::SUM);
    comm.Allreduce((sctl::Iterator<Real>) sum_err_loc.begin()+1, (sctl::Iterator<Real>) sum_err_all.begin()+1, 1, sctl::CommOp::SUM);
    comm.Allreduce((sctl::Iterator<Real>) sum_err_loc.begin()+2, (sctl::Iterator<Real>) sum_err_all.begin()+2, 1, sctl::CommOp::SUM);
    sum_err = sum_err_all;

    sctl::Vector<sctl::Long> Nnodes_loc(1);
    Nnodes_loc[0] = Nnodes;
    sctl::Vector<sctl::Long> Nnodes_all(1); 
    Nnodes_all[0] = 0;
    comm.Allreduce((sctl::Iterator<Real>) Nnodes_loc.begin(), (sctl::Iterator<Real>) Nnodes_all.begin(), 1, sctl::CommOp::SUM);
    // avg err
    sctl::Vector<Real> avg_err = sum_err / Nnodes_all[0];
    AddConstVec(err,-avg_err); // relative error with offset: max ((Ucalc - C) - Uexact) / Uexact, since C = Ucalc_exact - Uexact ~ E[Ucalc - Uexact]
}