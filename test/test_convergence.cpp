// export OMP_NUM_THREADS=16; time make DEBUG=0 -B bin/test1 && time mpirun -n 1 --map-by slot:pe=$OMP_NUM_THREADS ./bin/test1

#include "periodize.hpp"
#include "utils.hpp"

/**
 * Background flow with unit pressure gradient along X-axis.
 */
template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
    const Real pdrive = 10;
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

// Trefoil knot channel without particles
template <class Real> void trefoil_self_conv(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Comm comm, sctl::Long precond_mode) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real tol = 1e-10;
    const Real gmres_tol = 1e-13;
    const sctl::Long ElemOrder = 10;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::Long ptcl_ord = 1;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_trefoil(Nelem, ElemOrder, FourierOrder, 0, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, 0);
    elem_lst0 = std::get<0>(build0);
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_trefoil(Nelem, ElemOrder, FourierOrder, 1, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, 0);  
    elem_lst_nbr = std::get<0>(build_nbr);
    NormalOrient = std::get<1>(build_nbr);

    const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); // should be 3
    sctl::Long Nptcl = ptcls_rs.Dim();

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Vector<Real> X_proxy;
    X_proxy = Periodize1D<Real>::GetProxySurf(30,20); // proxy points coordinates

    if (write_ref) {
        elem_lst0.WriteVTK("vis/trefoil",X0,comm);
    }   

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst_nbr);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);

    StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal, comm); // potential from elem_lst0 to proxy points
    LayerPotenOp_proxy.AddElemList(elem_lst0);
    LayerPotenOp_proxy.SetTargetCoord(X_proxy);
    LayerPotenOp_proxy.SetAccuracy(tol);

    sctl::Matrix<Real> Usvd_p, VT_p, S_p, SforInv_p, Usvd, VT, S, SforInv;
    sctl::Matrix<Real> Sinv_p, Sinv;
    sctl::Long A11size_ptcl, A11size;

    if (precond_mode) {
        sctl::Profile::Tic("Preconditioners");
        // PRECONDITIONING PANEL: CYLINDER
        sctl::Vector<Real> Xc_precond, eps_precond; 
        sctl::Vector<sctl::Long> ElemOrderVec_precond(1), FourierOrderVec_precond(1);
        ElemOrderVec_precond[0] = ElemOrder;
        FourierOrderVec_precond[0] = FourierOrder;
        // Determine approximate radius of channel based on channel_mode
        Real channel_radius = 0.035;
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
        const auto BIO_1panel = [&DL_scal,&Precond_bio,&Nrepeat](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
            U->SetZero();
            Precond_bio.ComputePotential(*U, sigma);
            (*U) -= sigma*0.5 * DL_scal;
        };

        A11size = 3*ElemOrder*FourierOrder;
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
        sctl::Matrix<Real> A11forSVD = sctl::Matrix<Real>(A11);
        A11forSVD.SVD(Usvd, S, VT);
        SforInv = sctl::Matrix<Real>(S);
        Sinv = SforInv.pinv(1e-16);

        // =============== PRECONDITIONING : Particle =======================================
        sctl::Vector<sctl::Long> ptcls_pre;
        sctl::Vector<Real> ptcls_rs_pre, ptcls_Xcs_pre;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_precond_ptcl = obj.many_ptcls1(ptcl_ord,ElemOrder, FourierOrder, 0, 1, comm.Self(), ptcls_pre, ptcls_rs_pre, ptcls_Xcs_pre, 0); 
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

        A11size_ptcl = 3*ElemOrder*FourierOrder*ptcl_ord;
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
        
        sctl::Matrix<Real> A11forSVD_ptcl = sctl::Matrix<Real>(A11_ptcl);
        A11forSVD_ptcl.SVD(Usvd_p, S_p, VT_p);
        SforInv_p = sctl::Matrix<Real>(S_p);
        Sinv_p = SforInv_p.pinv(1e-16);

        sctl::Profile::Toc();
        sctl::Profile::print(&comm);
        sctl::Profile::reset();

    }

    // Get global index of the starting panel on this process
    sctl::Vector<sctl::Long> ElemOrderVec_temp(Nelem + ptcl_ord * Nptcl);
    ElemOrderVec_temp = ElemOrder;
    sctl::Vector<sctl::Long> FourierOrderVec_temp(ElemOrderVec_temp);
    FourierOrderVec_temp = FourierOrder;
    std::tuple<sctl::Long,sctl::Long> indtpl = obj.GetGlobalIdx(ElemOrderVec_temp, FourierOrderVec_temp, comm);
    sctl::Long loc_elem_cnt = std::get<0>(indtpl);
    sctl::Long loc_elem_dsp = std::get<1>(indtpl);

        // periodized layer potential operator
    const auto BIO = [&DL_scal,&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat,NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        const sctl::Long N = sigma.Dim();
        // std::cout << "in BIO, dim of sigma is " << N << ", output dim of LayerOp is " << LayerPotenOp0.Dim(1) << std::endl;

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
            Periodize1D<Real>::EvalFarField(U_far, X0, U_proxy,30,20);
            (*U) += U_far;
        } 
    };

    // Apply A11inv to each panel of a vector.
    // Look at global panel index and determine whether belongs to a particle or the channel. ASSUMES no particles are split up among processors.
    const auto AinvApply = [&Usvd,&Sinv,&VT,&A11size,&Usvd_p,&Sinv_p,&VT_p,&A11size_ptcl,&Nelem,&loc_elem_cnt,&loc_elem_dsp,&Nptcl](const sctl::Vector<Real>& vec) {
        sctl::Vector<Real> AinvVec(vec.Dim());
        if (Nptcl == 0 || (loc_elem_dsp+loc_elem_cnt) <= Nelem) {
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
            if (loc_elem_dsp >= Nelem) {
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
                sctl::Long Npanels_here = Nelem - loc_elem_dsp;
                sctl::Long Nptcls_here = loc_elem_cnt - Npanels_here;
                // // DEBUG
                // std::cout << "CHECK panel-ptcl split: Npanel = " << Npanels_here << ", Nptcl = " << Nptcls_here << std::endl;
                // ///////////////
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



    sctl::GMRES<Real> solver(comm);

    sctl::Vector<Real> sigma_temp;
    sctl::Profile::Tic("Setup");
    if (precond_mode) {
        solver(&sigma_temp, BIO_precond, X0, 1e0);
    } else {
        solver(&sigma_temp, BIO, X0, 1e0);
    }
    sctl::Profile::Toc();
    sctl::Profile::print(&comm,{"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    sctl::Profile::reset();

    sctl::Vector<Real> sigma;
    sctl::Profile::Tic("Solve");
    if (precond_mode) {
        sctl::Vector<Real> A11invF = AinvApply(-bg_flow(X0));
        solver(&sigma, BIO_precond, A11invF, gmres_tol);
    } else {
        solver(&sigma,BIO,-bg_flow(X0),gmres_tol);
    }
    sctl::Profile::Toc();
    sctl::Profile::print(&comm,{"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    

    { // Evaluate in interior, and write visualization
        // std::cout << "Rank " << comm.Rank()<< " calculating target points." << std::endl;
        PeriodicGeom<Real> trg;
        const sctl::Long Nelem_trg = 200;
        const sctl::Long FourierOrder_trg = 16;
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg;
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_trefoil(Nelem_trg, ElemOrder, FourierOrder_trg, 0, 1, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, 1, 0);
        elem_lst_trg = std::get<0>(build_trg);

        VolumeVis<Real> vol_vis(elem_lst_trg, comm, true); 
        X0 = vol_vis.GetCoord(); // set new target coordinates
        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U += bg_flow(X0);

        sctl::Vector<sctl::Long> size_loc(1);
        size_loc[0] = X0.Dim();
        sctl::Vector<sctl::Long> size_all(1);
        comm.Allreduce((sctl::Iterator<sctl::Long>) size_loc.begin(), (sctl::Iterator<sctl::Long>) size_all.begin(), 1, sctl::CommOp::SUM);
        // std::cout << "rank " << comm.Rank() << " size loc = " << size_loc[0] << ", size all is " << size_all[0] << std::endl;
        // std::string filename = "Trefoil_U_exact";
        std::string filename = "Trefoil_U_exact_"+std::to_string(comm.Rank());
        std::string filename_out = "out/"+filename+".txt";
        std::string filename_vis = "vis/"+filename;
        
        if (write_ref) {
            // // std::cout << "write ref" << std::endl;
            // sctl::Vector<Real> U_all(size_all[0]);
            // comm.Allgather((sctl::Iterator<Real>) U.begin(), size_loc[0], (sctl::Iterator<Real>) U_all.begin(), size_all[0]);
            // if (!comm.Rank()) {
            //     U_all.Write(filename_out.c_str());
            // }
            
            // VolumeVis<Real> vol_vis_write(elem_lst_trg, comm); 
            // X0 = vol_vis_write.GetCoord();
            // LayerPotenOp0.SetTargetCoord(X0);
            // sctl::Vector<Real> Uvis;
            // BIO(&Uvis, sigma);
            // Uvis += bg_flow(X0);
            // vol_vis_write.WriteVTK(filename_vis, Uvis); 

            U.Write(filename_out.c_str());
            
        } else {
            sctl::Vector<Real> U_ref;
            // if (!comm.Rank()) {
            //     U_ref.Read(filename_out.c_str());
            // }
            // comm.PartitionN(U_ref,size_loc[0]);
            // std::cout << "dim of U ref is " << U_ref.Dim() << ", dim of U vis is " << U_vis.Dim() << std::endl;
            U_ref.Read(filename_out.c_str());

            const auto err = U - U_ref;
            double max_err = 0;
            for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
            Real max_u = 0.;
            for (const auto e : U_ref) max_u = std::max<Real>(max_u, sctl::fabs(e));
            sctl::Vector<Real> err_loc(1);
            err_loc[0] = max_err;
            sctl::Vector<Real> err_all(1);
            err_all[0] = 0;
            comm.Allreduce((sctl::Iterator<sctl::Long>) err_loc.begin(), (sctl::Iterator<sctl::Long>) err_all.begin(), 1, sctl::CommOp::MAX);

            sctl::Vector<Real> u_loc(1);
            u_loc[0] = max_u;
            sctl::Vector<Real> u_all(1);
            u_all[0] = 0.;
            comm.Allreduce((sctl::Iterator<Real>) u_loc.begin(), (sctl::Iterator<Real>) u_all.begin(), 1, sctl::CommOp::MAX);

            if (!comm.Rank()) {
                std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << std::endl;
                std::cout<<"Max relative error = "<< std::setprecision(10) << err_all[0] / u_all[0] << std::endl;
            }
        }

    }
}

template <class Real> void particle_self_conv(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Integer peri_mode, sctl::Comm comm, const sctl::Long Nptcl) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real tol = 1e-10;
    const Real gmres_tol = 1e-10;
    const sctl::Long ElemOrder = 10;

    if (!comm.Rank()) {
        std::cout << "in particle self conv" << std::endl;
    }

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    // std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 0, 1, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, 0);
    // elem_lst0 = std::get<0>(build0);
    // std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 1, peri_mode, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, 0);
    // elem_lst_nbr = std::get<0>(build_nbr);
    // NormalOrient = std::get<1>(build_nbr);
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, 0);
    elem_lst0 = std::get<0>(build0);
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 1, peri_mode, comm, ptcls, ptcls_rs, ptcls_Xcs, 0);
    elem_lst_nbr = std::get<0>(build_nbr);
    NormalOrient = std::get<1>(build_nbr);

    const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); // should be 3
    if (!comm.Rank()) {
        std::cout << "periodic mode is " << peri_mode << ", Nrepeat is " << Nrepeat << std::endl;
        std::cout << "total number of particles is " << ptcls.Dim() << ", number of elements per processor is " << elem_lst0.Size() << std::endl;
    }
    
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

    // std::string nbr_vis = "vis/Particle";
    // sctl::Vector<Real> Xnbr;
    // elem_lst_nbr.GetNodeCoord(&Xnbr, nullptr, nullptr);
    // elem_lst_nbr.WriteVTK(nbr_vis,Xnbr,comm);

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst_nbr);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);

    StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal, comm); // potential from elem_lst0 to proxy points
    LayerPotenOp_proxy.AddElemList(elem_lst0);
    LayerPotenOp_proxy.SetTargetCoord(X_proxy);
    LayerPotenOp_proxy.SetAccuracy(tol);

    // periodized layer potential operator
    const auto BIO = [&DL_scal,&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat,&peri_mode,NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        const sctl::Long N = sigma.Dim();
        std::cout << "in BIO, dim of sigma is " << N << ", output dim of LayerOp is " << LayerPotenOp0.Dim(1) << std::endl;

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
            // Periodize1D<Real>::EvalFarField(U_far, X0, U_proxy);
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
            // std::cout << "U far size is " << U_far.Dim() << std::endl;
            (*U) += U_far;
        } 
    };

    // Solve for sigma to satisfy no-slip boundary conditions: BIO(sigma) + bg_flow = 0
    sctl::Vector<Real> sigma;
    sctl::GMRES<Real> solver(comm);
    solver(&sigma, BIO, -bg_flow(X0), gmres_tol);

    { // Evaluate in interior, and write visualization
        // std::cout << "Rank " << comm.Rank()<< " calculating target points." << std::endl;
        PeriodicGeom<Real> trg;
        CubeVolumeVisShifted<Real> vol_vis(20, 1.0, comm); // -- does not play nice with all gather, not fault of filter_target.
        // X0 = vol_vis.GetCoord();
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        std::cout << "size of X0 all is " << X0_all.Dim();
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);
        
        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U += bg_flow(X0);

        sctl::Vector<sctl::Long> size_loc(1);
        // size_loc[0] = X0.Dim();
        size_loc[0] = X0.Dim();
        sctl::Vector<sctl::Long> size_all(1);
        comm.Allreduce((sctl::Iterator<sctl::Long>) size_loc.begin(), (sctl::Iterator<sctl::Long>) size_all.begin(), 1, sctl::CommOp::SUM);
        std::cout << "rank " << comm.Rank() << " size loc = " << size_loc[0] << ", size all is " << size_all[0] << std::endl;
        std::string filename = "Particle"+std::to_string(Nptcl)+"_"+std::to_string(peri_mode)+"_peri_U_exact_"+std::to_string(comm.Rank());
        std::string filename_out = "out/"+filename+".txt";
        std::string filename_vis = "vis/"+filename;
        if (write_ref) {
            
            // // TODO: debug all gather
            // sctl::Vector<Real> U_all(size_all[0]);
            // comm.Allgather((sctl::Iterator<Real>) U.begin(), size_loc[0], (sctl::Iterator<Real>) U_all.begin(), size_all[0]);
            // // comm.Allgather((sctl::Iterator<Real>) U.begin(), size_loc[0], (sctl::Iterator<Real>) U_all.begin(), size_all[0]);
            // std::cout << U_all.Dim() << std::endl;
            // if (!comm.Rank()) {
            //     U_all.Write(filename_out.c_str());
            // }
            U.Write(filename_out.c_str());

        } else {
            sctl::Vector<Real> U_ref;
            // if (!comm.Rank()) {
            //     U_ref.Read(filename_out.c_str());
            // }
            U_ref.Read(filename_out.c_str());
            // comm.PartitionN(U_ref,size_loc[0]);
            // std::cout << "dim of U ref is " << U_ref.Dim() << ", dim of U vis is " << U_vis.Dim() << std::endl;
            const auto err = U - U_ref;
            double max_err = 0;
            for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
            Real max_u = 0.;
            for (const auto e : U_ref) max_u = std::max<Real>(max_u, sctl::fabs(e));

            sctl::Vector<Real> err_loc(1);
            err_loc[0] = max_err;
            sctl::Vector<Real> err_all(1);
            err_all[0] = 0;
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
        }
    }

}

// Conv div channel
template <class Real> void channel_self_conv(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Comm comm, sctl::Long precond_mode) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real tol = 1e-8;
    const Real gmres_tol = 1e-10;
    const sctl::Long ElemOrder = 10;

    PeriodicGeom<Real> obj;
    sctl::Long Nptcl = 50; // placeholder; will be replaced inside conv-div channel build.
    sctl::Long ptcl_ord = 1;
    sctl::Vector<sctl::Long> ptcls(Nptcl);
    ptcls = ptcl_ord;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    sctl::Long peri_mode = 1;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_conv_div(Nelem, ElemOrder, FourierOrder, 0, peri_mode, 0.1, 0.1, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, 0);
    elem_lst0 = std::get<0>(build0);
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_conv_div(Nelem, ElemOrder, FourierOrder, 1, peri_mode, 0.1, 0.1, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, 0);  
    elem_lst_nbr = std::get<0>(build_nbr);
    NormalOrient = std::get<1>(build_nbr);
  
    // std::cout << "Size of elem_lst_nbr is " << elem_lst_nbr.Size() << ", Size of elem_lst0 is " << elem_lst0.Size() <<std::endl;
    const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); // should be 3
    Nptcl = ptcls_rs.Dim(); // Number of particles could have changed after initializing.
    if (!comm.Rank()) {
        std::cout << "periodic mode is " << peri_mode << ", Nrepeat is " << Nrepeat << std::endl;
        std::cout << "total number of particles is " << Nptcl << ", number of elements per processor is " << elem_lst0.Size() << std::endl;
    }

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Vector<Real> X_proxy = Periodize1D<Real>::GetProxySurf(30,20); // proxy points coordinates

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst_nbr);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);

    StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal, comm); // potential from elem_lst0 to proxy points
    LayerPotenOp_proxy.AddElemList(elem_lst0);
    LayerPotenOp_proxy.SetTargetCoord(X_proxy);
    LayerPotenOp_proxy.SetAccuracy(tol);

    sctl::Matrix<Real> Usvd_p, VT_p, S_p, SforInv_p, Usvd, VT, S, SforInv;
    sctl::Matrix<Real> Sinv_p, Sinv;
    sctl::Long A11size_ptcl, A11size;

    if (precond_mode) {
        sctl::Profile::Tic("Preconditioners");
        // PRECONDITIONING PANEL: CYLINDER
        sctl::Vector<Real> Xc_precond, eps_precond; 
        sctl::Vector<sctl::Long> ElemOrderVec_precond(1), FourierOrderVec_precond(1);
        ElemOrderVec_precond[0] = ElemOrder;
        FourierOrderVec_precond[0] = FourierOrder;
        // Determine approximate radius of channel based on channel_mode
        Real channel_radius = 0.035;
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
        const auto BIO_1panel = [&DL_scal,&Precond_bio,&Nrepeat](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
            U->SetZero();
            Precond_bio.ComputePotential(*U, sigma);
            (*U) -= sigma*0.5 * DL_scal;
        };

        A11size = 3*ElemOrder*FourierOrder;
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
        sctl::Matrix<Real> A11forSVD = sctl::Matrix<Real>(A11);
        A11forSVD.SVD(Usvd, S, VT);
        SforInv = sctl::Matrix<Real>(S);
        Sinv = SforInv.pinv(1e-16);

        // =============== PRECONDITIONING : Particle =======================================
        sctl::Vector<sctl::Long> ptcls_pre;
        sctl::Vector<Real> ptcls_rs_pre, ptcls_Xcs_pre;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_precond_ptcl = obj.many_ptcls1(ptcl_ord,ElemOrder, FourierOrder, 0, 1, comm.Self(), ptcls_pre, ptcls_rs_pre, ptcls_Xcs_pre, 0); 
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

        A11size_ptcl = 3*ElemOrder*FourierOrder*ptcl_ord;
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
        
        sctl::Matrix<Real> A11forSVD_ptcl = sctl::Matrix<Real>(A11_ptcl);
        A11forSVD_ptcl.SVD(Usvd_p, S_p, VT_p);
        SforInv_p = sctl::Matrix<Real>(S_p);
        Sinv_p = SforInv_p.pinv(1e-16);

        sctl::Profile::Toc();
        sctl::Profile::print(&comm);
        sctl::Profile::reset();

    }

    // Get global index of the starting panel on this process
    sctl::Vector<sctl::Long> ElemOrderVec_temp(Nelem + ptcl_ord * Nptcl);
    ElemOrderVec_temp = ElemOrder;
    sctl::Vector<sctl::Long> FourierOrderVec_temp(ElemOrderVec_temp);
    FourierOrderVec_temp = FourierOrder;
    std::tuple<sctl::Long,sctl::Long> indtpl = obj.GetGlobalIdx(ElemOrderVec_temp, FourierOrderVec_temp, comm);
    sctl::Long loc_elem_cnt = std::get<0>(indtpl);
    sctl::Long loc_elem_dsp = std::get<1>(indtpl);

    // periodized layer potential operator
    const auto BIO = [&DL_scal,&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat,NormalOrient,&peri_mode, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
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
            Periodize1D<Real>::EvalFarField(U_far, X0, U_proxy, 30, 20);
            (*U) += U_far;
        } 
    };

    // Apply A11inv to each panel of a vector.
    // Look at global panel index and determine whether belongs to a particle or the channel. ASSUMES no particles are split up among processors.
    const auto AinvApply = [&Usvd,&Sinv,&VT,&A11size,&Usvd_p,&Sinv_p,&VT_p,&A11size_ptcl,&Nelem,&loc_elem_cnt,&loc_elem_dsp,&Nptcl](const sctl::Vector<Real>& vec) {
        sctl::Vector<Real> AinvVec(vec.Dim());
        if (Nptcl == 0 || (loc_elem_dsp+loc_elem_cnt) <= Nelem) {
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
            if (loc_elem_dsp >= Nelem) {
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
                sctl::Long Npanels_here = Nelem - loc_elem_dsp;
                sctl::Long Nptcls_here = loc_elem_cnt - Npanels_here;
                // // DEBUG
                // std::cout << "CHECK panel-ptcl split: Npanel = " << Npanels_here << ", Nptcl = " << Nptcls_here << std::endl;
                // ///////////////
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

    sctl::GMRES<Real> solver(comm);

    sctl::Vector<Real> sigma_temp;
    sctl::Profile::Tic("Setup");
    if (precond_mode) {
        solver(&sigma_temp, BIO_precond, X0, 1e0);
    } else {
        solver(&sigma_temp, BIO, X0, 1e0);
    }
    sctl::Profile::Toc();
    sctl::Profile::print(&comm,{"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    sctl::Profile::reset();

    sctl::Vector<Real> sigma;
    sctl::Profile::Tic("Solve");
    if (precond_mode) {
        sctl::Vector<Real> A11invF = AinvApply(-bg_flow(X0));
        solver(&sigma, BIO_precond, A11invF, gmres_tol);
    } else {
        solver(&sigma,BIO,-bg_flow(X0),gmres_tol);
    }
    sctl::Profile::Toc();
    sctl::Profile::print(&comm,{"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});

    { // Evaluate in interior, and write visualization
        // std::cout << "Rank " << comm.Rank()<< " calculating target points." << std::endl;
        PeriodicGeom<Real> trg;
        const sctl::Long Nelem_trg = 16;
        const sctl::Long FourierOrder_trg = 16;
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg; // dim = 0 so no particles are first generated
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_conv_div(Nelem_trg, ElemOrder, FourierOrder_trg, 0, peri_mode, 0.1, 0.1, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, ptcl_ord, 0);
        elem_lst_trg = std::get<0>(build_trg);

        VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        // X0 = vol_vis.GetCoord();
        sctl::Vector<Real> X0_all = vol_vis.GetCoord(); // set new target coordinates
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U += bg_flow(X0);

        sctl::Vector<sctl::Long> size_loc(1);
        size_loc[0] = X0.Dim();
        sctl::Vector<sctl::Long> size_all(1);
        comm.Allreduce((sctl::Iterator<sctl::Long>) size_loc.begin(), (sctl::Iterator<sctl::Long>) size_all.begin(), 1, sctl::CommOp::SUM);
        // std::cout << "rank " << comm.Rank() << " size loc = " << size_loc[0] << ", size all is " << size_all[0] << std::endl;

        // std::string filename = "Conv_div.txt";
        std::string filename = "ConvDiv_U_exact_"+std::to_string(comm.Rank());
        std::string filename_out = "out/"+filename+".txt";
        std::string filename_vis = "vis/"+filename;
        if (write_ref) {
            // // try with all gather, when no filter target is employed
            // sctl::Vector<Real> U_all(size_all[0]);
            // comm.Allgather((sctl::Iterator<Real>) U.begin(), size_loc[0], (sctl::Iterator<Real>) U_all.begin(), size_all[0]);
            // std::cout << U_all.Dim() << std::endl;
            // if (!comm.Rank()) {
            //     U_all.Write(filename_out.c_str());
            // }

            U.Write(filename_out.c_str());

        } else {
            sctl::Vector<Real> U_ref;
            // if (!comm.Rank()) {
            //     U_ref.Read(filename_out.c_str());
            // }
            U_ref.Read(filename_out.c_str());
            // comm.PartitionN(U_ref,size_loc[0]);
            const auto err = U - U_ref;
            double max_err = 0;
            for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
            Real max_u = 0.;
            for (const auto e : U_ref) max_u = std::max<Real>(max_u, sctl::fabs(e));
            sctl::Vector<Real> err_loc(1);
            err_loc[0] = max_err;
            sctl::Vector<Real> err_all(1);
            err_all[0] = 0;
            comm.Allreduce((sctl::Iterator<sctl::Long>) err_loc.begin(), (sctl::Iterator<sctl::Long>) err_all.begin(), 1, sctl::CommOp::MAX);

            sctl::Vector<Real> u_loc(1);
            u_loc[0] = max_u;
            sctl::Vector<Real> u_all(1);
            u_all[0] = 0.;
            comm.Allreduce((sctl::Iterator<Real>) u_loc.begin(), (sctl::Iterator<Real>) u_all.begin(), 1, sctl::CommOp::MAX);

            if (!comm.Rank()) {
                std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << std::endl;
                std::cout<<"Max relative error = "<< std::setprecision(10) << err_all[0] / u_all[0] << std::endl;
            }
        }

  }
}

// Conv div channel
template <class Real> void trefoil_ptcl_self_conv(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Comm comm, sctl::Long precond_mode) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real tol = 1e-8;
    const Real gmres_tol = 1e-13;
    const sctl::Long ElemOrder = 10;

    PeriodicGeom<Real> obj;
    sctl::Long Nptcl = 50; // placeholder; will be replaced inside conv-div channel build.
    sctl::Long ptcl_ord = 1;
    sctl::Vector<sctl::Long> ptcls(Nptcl);
    ptcls = ptcl_ord;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    sctl::Long peri_mode = 1;
    // std::cout << "right before build trefoils" << std::endl;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_trefoil(Nelem, ElemOrder, FourierOrder, 0, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, 0);
    elem_lst0 = std::get<0>(build0);
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.build_trefoil(Nelem, ElemOrder, FourierOrder, 1, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, 0);  
    elem_lst_nbr = std::get<0>(build_nbr);
    NormalOrient = std::get<1>(build_nbr);
  
    // std::cout << "Size of elem_lst_nbr is " << elem_lst_nbr.Size() << ", Size of elem_lst0 is " << elem_lst0.Size() <<std::endl;
    const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); // should be 3
    Nptcl = ptcls_rs.Dim(); // Number of particles could have changed after initializing.
    if (!comm.Rank()) {
        std::cout << "periodic mode is " << peri_mode << ", Nrepeat is " << Nrepeat << std::endl;
        std::cout << "total number of particles is " << Nptcl << ", number of elements per processor is " << elem_lst0.Size() << std::endl;
    }

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Vector<Real> X_proxy = Periodize1D<Real>::GetProxySurf(30,20); // proxy points coordinates

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst_nbr);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);

    StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal, comm); // potential from elem_lst0 to proxy points
    LayerPotenOp_proxy.AddElemList(elem_lst0);
    LayerPotenOp_proxy.SetTargetCoord(X_proxy);
    LayerPotenOp_proxy.SetAccuracy(tol);

    sctl::Matrix<Real> Usvd_p, VT_p, S_p, SforInv_p, Usvd, VT, S, SforInv;
    sctl::Matrix<Real> Sinv_p, Sinv;
    sctl::Long A11size_ptcl, A11size;

    if (precond_mode) {
        sctl::Profile::Tic("Preconditioners");
        // PRECONDITIONING PANEL: CYLINDER
        sctl::Vector<Real> Xc_precond, eps_precond; 
        sctl::Vector<sctl::Long> ElemOrderVec_precond(1), FourierOrderVec_precond(1);
        ElemOrderVec_precond[0] = ElemOrder;
        FourierOrderVec_precond[0] = FourierOrder;
        // Determine approximate radius of channel based on channel_mode
        Real channel_radius = 0.035;
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
        const auto BIO_1panel = [&DL_scal,&Precond_bio,&Nrepeat](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
            U->SetZero();
            Precond_bio.ComputePotential(*U, sigma);
            (*U) -= sigma*0.5 * DL_scal;
        };

        A11size = 3*ElemOrder*FourierOrder;
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
        sctl::Matrix<Real> A11forSVD = sctl::Matrix<Real>(A11);
        A11forSVD.SVD(Usvd, S, VT);
        SforInv = sctl::Matrix<Real>(S);
        Sinv = SforInv.pinv(1e-16);

        // =============== PRECONDITIONING : Particle =======================================
        sctl::Vector<sctl::Long> ptcls_pre;
        sctl::Vector<Real> ptcls_rs_pre, ptcls_Xcs_pre;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_precond_ptcl = obj.many_ptcls1(ptcl_ord,ElemOrder, FourierOrder, 0, 1, comm.Self(), ptcls_pre, ptcls_rs_pre, ptcls_Xcs_pre, 0); 
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

        A11size_ptcl = 3*ElemOrder*FourierOrder*ptcl_ord;
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
        
        sctl::Matrix<Real> A11forSVD_ptcl = sctl::Matrix<Real>(A11_ptcl);
        A11forSVD_ptcl.SVD(Usvd_p, S_p, VT_p);
        SforInv_p = sctl::Matrix<Real>(S_p);
        Sinv_p = SforInv_p.pinv(1e-16);

        sctl::Profile::Toc();
        sctl::Profile::print(&comm);
        sctl::Profile::reset();

    }

    // Get global index of the starting panel on this process
    sctl::Vector<sctl::Long> ElemOrderVec_temp(Nelem + ptcl_ord * Nptcl);
    ElemOrderVec_temp = ElemOrder;
    sctl::Vector<sctl::Long> FourierOrderVec_temp(ElemOrderVec_temp);
    FourierOrderVec_temp = FourierOrder;
    std::tuple<sctl::Long,sctl::Long> indtpl = obj.GetGlobalIdx(ElemOrderVec_temp, FourierOrderVec_temp, comm);
    sctl::Long loc_elem_cnt = std::get<0>(indtpl);
    sctl::Long loc_elem_dsp = std::get<1>(indtpl);

    // periodized layer potential operator
    const auto BIO = [&DL_scal,&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat,NormalOrient,&peri_mode, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
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
            Periodize1D<Real>::EvalFarField(U_far, X0, U_proxy, 30, 20);
            (*U) += U_far;
        } 
    };

    // Apply A11inv to each panel of a vector.
    // Look at global panel index and determine whether belongs to a particle or the channel. ASSUMES no particles are split up among processors.
    const auto AinvApply = [&Usvd,&Sinv,&VT,&A11size,&Usvd_p,&Sinv_p,&VT_p,&A11size_ptcl,&Nelem,&loc_elem_cnt,&loc_elem_dsp,&Nptcl](const sctl::Vector<Real>& vec) {
        sctl::Vector<Real> AinvVec(vec.Dim());
        if (Nptcl == 0 || (loc_elem_dsp+loc_elem_cnt) <= Nelem) {
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
            if (loc_elem_dsp >= Nelem) {
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
                sctl::Long Npanels_here = Nelem - loc_elem_dsp;
                sctl::Long Nptcls_here = loc_elem_cnt - Npanels_here;
                // // DEBUG
                // std::cout << "CHECK panel-ptcl split: Npanel = " << Npanels_here << ", Nptcl = " << Nptcls_here << std::endl;
                // ///////////////
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

    sctl::GMRES<Real> solver(comm);

    sctl::Vector<Real> sigma_temp;
    sctl::Profile::Tic("Setup");
    if (precond_mode) {
        solver(&sigma_temp, BIO_precond, X0, 1e0);
    } else {
        solver(&sigma_temp, BIO, X0, 1e0);
    }
    sctl::Profile::Toc();
    sctl::Profile::print(&comm,{"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    sctl::Profile::reset();

    sctl::Vector<Real> sigma;
    sctl::Profile::Tic("Solve");
    if (precond_mode) {
        sctl::Vector<Real> A11invF = AinvApply(-bg_flow(X0));
        solver(&sigma, BIO_precond, A11invF, gmres_tol);
    } else {
        solver(&sigma,BIO,-bg_flow(X0),gmres_tol);
    }
    sctl::Profile::Toc();
    sctl::Profile::print(&comm,{"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});

    { // Evaluate in interior, and write visualization
        // std::cout << "Rank " << comm.Rank()<< " calculating target points." << std::endl;
        PeriodicGeom<Real> trg;
        const sctl::Long Nelem_trg = 200;
        const sctl::Long FourierOrder_trg = 16;
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg;
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_trefoil(Nelem_trg, ElemOrder, FourierOrder_trg, 0, 1, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, 1, 0);
        elem_lst_trg = std::get<0>(build_trg);

        VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        // X0 = vol_vis.GetCoord();
        sctl::Vector<Real> X0_all = vol_vis.GetCoord(); // set new target coordinates
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U += bg_flow(X0);

        sctl::Vector<sctl::Long> size_loc(1);
        size_loc[0] = X0.Dim();
        sctl::Vector<sctl::Long> size_all(1);
        comm.Allreduce((sctl::Iterator<sctl::Long>) size_loc.begin(), (sctl::Iterator<sctl::Long>) size_all.begin(), 1, sctl::CommOp::SUM);
        // std::cout << "rank " << comm.Rank() << " size loc = " << size_loc[0] << ", size all is " << size_all[0] << std::endl;

        std::string filename = "Trefoil_ptcl_U_exact_"+std::to_string(comm.Rank());
        std::string filename_out = "out/"+filename+".txt";
        std::string filename_vis = "vis/"+filename;
        if (write_ref) {
            U.Write(filename_out.c_str());

        } else {
            sctl::Vector<Real> U_ref;
            // if (!comm.Rank()) {
            //     U_ref.Read(filename_out.c_str());
            // }
            U_ref.Read(filename_out.c_str());
            // comm.PartitionN(U_ref,size_loc[0]);
            const auto err = U - U_ref;
            double max_err = 0;
            for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
            Real max_u = 0.;
            for (const auto e : U_ref) max_u = std::max<Real>(max_u, sctl::fabs(e));
            sctl::Vector<Real> err_loc(1);
            err_loc[0] = max_err;
            sctl::Vector<Real> err_all(1);
            err_all[0] = 0;
            comm.Allreduce((sctl::Iterator<sctl::Long>) err_loc.begin(), (sctl::Iterator<sctl::Long>) err_all.begin(), 1, sctl::CommOp::MAX);

            sctl::Vector<Real> u_loc(1);
            u_loc[0] = max_u;
            sctl::Vector<Real> u_all(1);
            u_all[0] = 0.;
            comm.Allreduce((sctl::Iterator<Real>) u_loc.begin(), (sctl::Iterator<Real>) u_all.begin(), 1, sctl::CommOp::MAX);

            if (!comm.Rank()) {
                std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << std::endl;
                std::cout<<"Max relative error = "<< std::setprecision(10) << err_all[0] / u_all[0] << std::endl;
            }
        }
    }
}

int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;
//   std::cout << "In main" << std::endl;

  {
    // sctl::Profile::Enable(true);
    sctl::Comm comm = sctl::Comm::World();
    long test_mode = std::stol(argv[1]); // =0 for trefoil, =1 for 1-particle; =2 for conv div
    long peri_mode = std::stol(argv[2]); // 1- or 3- periodic
    long precond_mode = std::stol(argv[3]); // =0 for no precond, =1 for panel+particle precond.
    // std::cout << "Test mode: " << test_mode << ", peri mode: " << peri_mode << ", precond_mode: " << precond_mode << std::endl;

    sctl::Vector<sctl::Long> Nelem_lst;
    if (test_mode==1) { // particle, doesn't need many panels -- but if using Nptcl = 25 will need more panels
        for (int i=1; i<10; i += 2) {
            Nelem_lst.PushBack(2*i);
        }
    } else if (test_mode == 0 || test_mode == 3) { // trefoil needs more panels
        Nelem_lst.PushBack(200);
        Nelem_lst.PushBack(400);
        Nelem_lst.PushBack(800);
    } else {
        // for conv div channel
        Nelem_lst.PushBack(6);
        Nelem_lst.PushBack(12);
        Nelem_lst.PushBack(24);
    } 
    
    sctl::Vector<sctl::Long> FourierOrder_lst;
    if (test_mode == 2) {
        FourierOrder_lst.PushBack(16);
        FourierOrder_lst.PushBack(32);
        FourierOrder_lst.PushBack(64);
    } else {
        // FourierOrder_lst.PushBack(16);
        FourierOrder_lst.PushBack(32);
        FourierOrder_lst.PushBack(64);
        FourierOrder_lst.PushBack(96); 
    }

    sctl::Profile::Enable(true);
    // trefoil_self_conv<Real>(200, 32, false, comm, 0);
    // particle_self_conv<Real>(1, 16, 1, true, comm, 1);
    // channel_self_conv<Real>(4, 16, false, comm);
    // trefoil_ptcl_self_conv<Real>(200, 64, true, comm, 0);
    // trefoil_ptcl_self_conv<Real>(200, 64, true, comm, 1);
    
    sctl::Long Nelem, FourierOrder;
    for (int i=Nelem_lst.Dim()-1; i>=0; i--) {
        for (int j=FourierOrder_lst.Dim()-1; j>=0; j--) {
            Nelem = Nelem_lst[i];
            FourierOrder = FourierOrder_lst[j];
            if (!comm.Rank()) {
                std::cout << "Nelem = " << Nelem << ", FourierOrder = " << FourierOrder << "; " << std::endl;
            }
            if (i==Nelem_lst.Dim()-1 && j == FourierOrder_lst.Dim()-1) {
                if (test_mode==0) {
                    trefoil_self_conv<Real>(Nelem, FourierOrder, true, comm, precond_mode);
                } else if (test_mode == 1) {
                    particle_self_conv<Real>(Nelem, FourierOrder, true, peri_mode, comm, 25);
                } else if (test_mode==2) {
                    channel_self_conv<Real>(Nelem, FourierOrder, true, comm, precond_mode); 
                } else {
                    trefoil_ptcl_self_conv<Real>(Nelem, FourierOrder, true, comm, precond_mode);
                }
                // continue;
            } else {
                if (test_mode==0) {
                    trefoil_self_conv<Real>(Nelem, FourierOrder, false, comm, precond_mode);
                } else if (test_mode == 1){
                    particle_self_conv<Real>(Nelem, FourierOrder, false, peri_mode, comm, 25); 
                } else if (test_mode==2) {
                    channel_self_conv<Real>(Nelem, FourierOrder, false, comm, precond_mode); 
                } else {
                    trefoil_ptcl_self_conv<Real>(Nelem, FourierOrder, false, comm, precond_mode);
                }
                
            }
        }
    }
    
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}

