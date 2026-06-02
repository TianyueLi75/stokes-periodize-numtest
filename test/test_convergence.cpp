// export OMP_NUM_THREADS=16; time make DEBUG=0 -B bin/test1 && time mpirun -n 1 --map-by slot:pe=$OMP_NUM_THREADS ./bin/test1

#include "periodize.hpp"
#include "utils.hpp"
#include "planeNaive.hpp"

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
        std::cout << " successfully read file " << precond0_file << std::endl;
        PrecondMat1.template Read<Real>(precond1_file.c_str());
        A11size = PrecondMat0.Dim(1);
    } else {
        std::cout << " Making precond files " << std::endl;
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
    // std::cout << "avg err: " << avg_err[0] << ", " << avg_err[1] << ", " << avg_err[2] << std::endl;
    // for (int i=0; i<err.Dim()/3; i++) {
    //     std::cout << "err after subtracting avg err: " << std::setprecision(10) << err[i*3+0] << ", " << err[i*3+1] << ", " << err[i*3+2] << ". " << std::endl;
    // }
}

/**
 * Build element list for a straight channel with a sphere inside.
 */
template <class Real> sctl::SlenderElemList<Real> build_elem_lst(const sctl::Long Nelem_sphere, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Vector<Real>* NormalOrient_ptr = nullptr) {
    sctl::Vector<Real> Xc, eps, orient;
    sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
    sctl::Long Nelem_channel = 20;
    for (sctl::Long i = 0; i < Nelem_channel; i++) {
        ElemOrderVec.PushBack(ElemOrder);
        FourierOrderVec.PushBack(FourierOrder);
        const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
        for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
            const Real x = (i+nodes[j])/Nelem_channel;
            Xc.PushBack(x);
            Xc.PushBack(0.5);
            Xc.PushBack(0.5);
            eps.PushBack(0.2);

            orient.PushBack(0);
            orient.PushBack(0);
            orient.PushBack(1);
        }
    }

    for (sctl::Long i = 0; i < Nelem_sphere; i++) { // add a sphere
        ElemOrderVec.PushBack(ElemOrder);
        FourierOrderVec.PushBack(FourierOrder);
        const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
        for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
            const Real r = 0.1;
            const Real theta = sctl::const_pi<Real>() * (i+nodes[j])/Nelem_sphere;
            Xc.PushBack(0.5+r*sctl::cos<Real>(theta));
            Xc.PushBack(0.5);
            Xc.PushBack(0.5);
            eps.PushBack(r*sctl::sin<Real>(theta));

            orient.PushBack(0);
            orient.PushBack(0);
            orient.PushBack(1);
        }
    }

    sctl::SlenderElemList<Real> elem_lst0(ElemOrderVec, FourierOrderVec, Xc, eps, orient);

    if (NormalOrient_ptr != nullptr) {
        NormalOrient_ptr->ReInit(0);
        constexpr sctl::Integer COORD_DIM = 3;
        sctl::Vector<sctl::Long> elem_wise_node_cnt;
        elem_lst0.GetNodeCoord(nullptr, nullptr, &elem_wise_node_cnt);
        for (sctl::Long i = 0; i < elem_wise_node_cnt.Dim(); i++) {
            for (sctl::Long j = 0; j < elem_wise_node_cnt[i]*COORD_DIM; j++) {
                NormalOrient_ptr->PushBack(i < Nelem_channel ? -1 : 1);
            }
        }
    }

    return elem_lst0;
}

// self convergence test for a straight channel with a sphere
template <class Real> void channel_sphere_self_conv(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Comm comm, const Real tol, const Real gmres_tol) {
    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real pressure_drop = -1.0;
    const Real period_length = 1;

    const sctl::Long gmres_max_iter = 200;
    const sctl::Long ElemOrder = 10;

    sctl::Vector<Real> NormalOrient; // normal orientation (1 if normal into fluid, else -1)
    const auto elem_lst0 = build_elem_lst(Nelem, ElemOrder, FourierOrder, &NormalOrient); // geometry in the unit box [0,1]^3
    // Hardcode the properties of the centered sphere.
    sctl::Vector<sctl::Long> ptcls(1);
    ptcls = Nelem;
    sctl::Vector<Real> ptcls_rs(1);
    ptcls_rs = 0.1;
    sctl::Vector<Real> ptcls_Xcs(3);
    ptcls_Xcs = 0.5;

    Real surface_area;
    sctl::Vector<Real> X0, wts;
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Vector<Real> X0surf = X0;
    if (write_ref) {
        elem_lst0.WriteVTK("vis/channel_sphere",X0surf,comm);
    }  
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

    StokesBIO<Real> LayerPotenOp0(SL_scal, DL_scal, comm);
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetTargetCoord(X0surf);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,&NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
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
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) += sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    sctl::Vector<Real> U, sigma;
    sctl::GMRES<Real> solver(comm);
    solver(&sigma, BIO, bg_flow(X0surf) * (pressure_drop/period_length), gmres_tol, gmres_max_iter);

    { // Evaluate in interior, and write visualization
        PeriodicGeom<Real> trg;
        VolumeVis<Real> vol_vis(elem_lst0, comm);
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);
        LayerPotenOp0.SetTargetCoord(X0);
        BIO(&U, sigma);
        U -= bg_flow(X0) * (pressure_drop/period_length);

        std::string filename = "Channel_sphere_1peri_U_exact_"+std::to_string(comm.Rank());
        std::string filename_out = "out/"+filename+".txt";
        std::string filename_vis = "vis/"+filename;
        if (write_ref) {
            
            U.Write(filename_out.c_str());
            // Create array of velocity for all target points, including filtered out ones.
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
            // Write visualization to VTK
            vol_vis.WriteVTK(filename_vis, U_vis);

        } else {

            sctl::Vector<Real> U_ref;
            U_ref.Read(filename_out.c_str());
            sctl::Vector<Real> err = U - U_ref;
            renormalize_error(err,comm);
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
                std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << ", Max u = " << u_all[0] << ", Max relative error = " << err_all[0] / u_all[0] << std::endl;
            }
        }
    }
}

// Trefoil knot channel without particles
template <class Real> void trefoil_self_conv(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Comm comm, const Real tol, const Real gmres_tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;
    const sctl::Long ElemOrder = 10;
    const sctl::Long gmres_max_iter = 200;

    const Real pressure_drop = -1.0;
    const Real period_length = 1;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::Long ptcl_ord = 1;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_trefoil(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, 0);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);    
    // sctl::Long Nptcl = ptcls_rs.Dim();

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);

    if (write_ref) {
        elem_lst0.WriteVTK("vis/trefoil",X0,comm);
    }   

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

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);

    // // Preconditioning
    // sctl::Matrix<Real> PrecondMat0, PrecondMat1, PrecondMat0_ptcl, PrecondMat1_ptcl;
    // sctl::Long A11size_ptcl, A11size;
    // Real channel_radius = 0.035; // for trefoil geometry
    // A11size = precond_channel(PrecondMat0, PrecondMat1, Nelem, ElemOrder, FourierOrder, channel_radius, SL_scal, DL_scal, comm);
    // A11size_ptcl = precond_ptcl(PrecondMat0_ptcl, PrecondMat1_ptcl, ptcl_ord, ElemOrder, FourierOrder, SL_scal, DL_scal, comm);

    // // Get global index of the starting panel on this process
    // sctl::Vector<sctl::Long> ElemOrderVec_temp(Nelem + ptcl_ord * Nptcl);
    // ElemOrderVec_temp = ElemOrder;
    // sctl::Vector<sctl::Long> FourierOrderVec_temp(ElemOrderVec_temp);
    // FourierOrderVec_temp = FourierOrder;
    // std::tuple<sctl::Long,sctl::Long> indtpl = obj.GetGlobalIdx(ElemOrderVec_temp, FourierOrderVec_temp, comm);
    // sctl::Long loc_elem_cnt = std::get<0>(indtpl);
    // sctl::Long loc_elem_dsp = std::get<1>(indtpl);

    const auto BIO = [&wts,&surface_area,&elem_lst0,&DL_scal,&LayerPotenOp0,NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
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

    // // Apply A11inv to each panel of a vector.
    // // Look at global panel index and determine whether belongs to a particle or the channel. ASSUMES no particles are split up among processors.
    // const auto AinvApply = [&PrecondMat0,&PrecondMat1,&A11size,&PrecondMat0_ptcl,&PrecondMat1_ptcl,&A11size_ptcl,&Nelem,&loc_elem_cnt,&loc_elem_dsp,&Nptcl](const sctl::Vector<Real>& vec) {
    //     sctl::Vector<Real> AinvVec(vec.Dim());
    //     if (Nptcl == 0 || (loc_elem_dsp+loc_elem_cnt) <= Nelem) {
    //         // std::cout << "All panels" << std::endl;
    //         // if no particles in channel or if all panels here are on channel
    //         sctl::Long N = vec.Dim();
    //         sctl::Long Npanels = N / A11size; 
    //         for (sctl::Long i=0; i<Npanels; i++) {
    //             sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
    //             sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
    //             for (sctl::Long j=0; j<A11size; j++) {
    //                 AinvVec[i*A11size + j] = AinvVecMat(j,0);
    //             }
    //         }
    //     } else {
    //         if (loc_elem_dsp >= Nelem) {
    //             std::cout << "All particles" << std::endl;
    //             // all panels here are ptcl
    //             sctl::Long N = vec.Dim();
    //             sctl::Long Nptcls = N / A11size_ptcl; 
    //             for (sctl::Long i=0; i<Nptcls; i++) {
    //                 sctl::Matrix<Real> vecMat(A11size_ptcl,1,(sctl::Iterator<Real>) vec.begin() + i*A11size_ptcl,true);
    //                 sctl::Matrix<Real> AinvVecMat = PrecondMat0_ptcl * (PrecondMat1_ptcl * vecMat);
    //                 for (sctl::Long j=0; j<A11size_ptcl; j++) {
    //                     AinvVec[i*A11size_ptcl + j] = AinvVecMat(j,0);
    //                 }
    //             }
    //         } else {
    //             sctl::Long Npanels_here = Nelem - loc_elem_dsp;
    //             sctl::Long Nptcls_here = loc_elem_cnt - Npanels_here;
    //             // // DEBUG
    //             std::cout << "CHECK panel-ptcl split: Npanel = " << Npanels_here << ", Nptcl = " << Nptcls_here << std::endl;
    //             // ///////////////
    //             for (sctl::Long i=0; i<Npanels_here; i++) {
    //                 sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
    //                 sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
    //                 for (sctl::Long j=0; j<A11size; j++) {
    //                     AinvVec[i*A11size + j] = AinvVecMat(j,0);
    //                 }
    //             }
    //             for (sctl::Long i=0; i<Nptcls_here; i++) {
    //                 sctl::Matrix<Real> vecMat(A11size_ptcl,1,(sctl::Iterator<Real>) vec.begin() + i*A11size_ptcl + Npanels_here*A11size,true);
    //                 sctl::Matrix<Real> AinvVecMat = PrecondMat0_ptcl * (PrecondMat1_ptcl * vecMat);
    //                 for (sctl::Long j=0; j<A11size_ptcl; j++) {
    //                     AinvVec[Npanels_here*A11size + i*A11size_ptcl + j] = AinvVecMat(j,0);
    //                 }
    //             }
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

    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov;
    sctl::Vector<Real> sigma;
    sctl::Vector<Real> rhs = bg_flow(X0) * (pressure_drop/period_length);
    // sctl::Vector<Real> A11invF = AinvApply(rhs);
    // solver(&sigma, BIO_precond, A11invF, gmres_tol, gmres_max_iter);
    solver(&sigma, BIO, rhs, gmres_tol, gmres_max_iter, false, nullptr, &krylov);

    { // Evaluate in interior, and write visualization
        PeriodicGeom<Real> trg;
        const sctl::Long Nelem_trg = 200;
        const sctl::Long FourierOrder_trg = 4;
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg;
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_trefoil(Nelem_trg, ElemOrder, FourierOrder_trg, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, 1, 0);
        elem_lst_trg = std::get<0>(build_trg);

        VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        X0 = vol_vis.GetCoord(); // set new target coordinates
        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U -= bg_flow(X0)* (pressure_drop/period_length);

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

            U.Write(filename_out.c_str());
            vol_vis.WriteVTK(filename_vis, U);
            
        } else {

            sctl::Vector<Real> U_ref;
            U_ref.Read(filename_out.c_str());

            sctl::Vector<Real> err = U - U_ref;
            // std::cout << "===== error before renorm: ======" << std::endl;
            // for (int i=0; i<err.Dim()/3; i++) {
            //     std::cout << std::setprecision(8) << err[i*3+0] << ", " << err[i*3+1] << ", " << err[i*3+2] << "; " << std::endl;
            // }
            renormalize_error(err,comm);
            // std::cout << "===== error after renorm: ======" << std::endl;
            // for (int i=0; i<err.Dim()/3; i++) {
            //     std::cout << std::setprecision(8) << err[i*3+0] << ", " << err[i*3+1] << ", " << err[i*3+2] << "; " << std::endl;
            // }
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
                std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << ", Max u = " << u_all[0] << ", Max relative error = " << err_all[0] / u_all[0] << std::endl;
            }

            // Try average error instead of max? (can probably use omp_scan.. )
            double avg_err = 0.;
            for (const auto e : err) avg_err += sctl::fabs(e);
            sctl::Vector<Real> avg_err_loc(1);
            avg_err_loc[0] = avg_err;
            sctl::Vector<Real> avg_err_all(1);
            avg_err_all[0] = 0;
            comm.Allreduce((sctl::Iterator<Real>) avg_err_loc.begin(), (sctl::Iterator<Real>) avg_err_all.begin(), 1, sctl::CommOp::SUM);
            std::cout << "on MPI process " << comm.Rank() << ", sum of err locally is " << avg_err << ", total error is " << avg_err_all[0] << std::endl;

            long size_err = err.Dim();
            sctl::Vector<sctl::Long> size_err_loc(1);
            size_err_loc[0] = size_err;
            sctl::Vector<sctl::Long> size_err_all(1);
            size_err_all[0] = 0;
            comm.Allreduce((sctl::Iterator<sctl::Long>) size_err_loc.begin(), (sctl::Iterator<sctl::Long>) size_err_all.begin(), 1, sctl::CommOp::SUM);

            std::cout << "size of error on process" << comm.Rank() << " is " << size_err << ", after collection total size is " << size_err_all[0] << std::endl;

            avg_err = avg_err_all[0] / (1.0*size_err_all[0]);

            if (!comm.Rank()) {
                std::cout<<"Average error = "<< std::setprecision(10) << avg_err << ", average relative error (divide by max u above) = " << avg_err / u_all[0] << std::endl;
            }
        }

    }
}

template <class Real> void particle_self_conv(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Integer peri_mode, sctl::Comm comm, const sctl::Long Nptcl, const Real tol, const Real gmres_tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const sctl::Long ElemOrder = 10;
    const sctl::Long geom_mode = 0;
    const Real period_length = 1.;
    const sctl::Long gmres_max_iter = 100;
    const Real pressure_drop = -1.;

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
    // Nptcl = ptcls_rs.Dim(); // many_ptcls will not change number of particles.

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
    // std::cout << "DEBUG surface area = " << surface_area << "." << std::endl;

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    
    if (peri_mode == 1) {
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);
    } else if (peri_mode == 2) {
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XY, period_length);
    } else if (peri_mode == 3) {
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length);
    } else {
        SCTL_ASSERT(false);
    }

    const auto BIO = [&wts,&surface_area,&elem_lst0,&DL_scal,&LayerPotenOp0,NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
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

    // sctl::Matrix<Real> PrecondMat0, PrecondMat1;
    // sctl::Long A11size;
    // A11size = precond_ptcl(PrecondMat0, PrecondMat1, Nelem, ElemOrder, FourierOrder, SL_scal, DL_scal, comm);

    // // Apply A11inv to each panel of vec.
    // const auto AinvApply = [&PrecondMat0,&PrecondMat1,&A11size, &comm](const sctl::Vector<Real>& vec) {
    //     sctl::Long N = vec.Dim();
    //     sctl::Long Nptcl = N / A11size; 
    //     sctl::Vector<Real> AinvVec(N);
    //     for (sctl::Long i=0; i<Nptcl; i++) {
    //         // for each particle, apply A11inv.
    //         sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
    //         sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
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

    const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
        sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
        AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

        sctl::Vector<Real> U0;
        LayerPotenOp0.ComputeSL(U0, force_density);
        return U0;
    };

    // Solve for sigma to satisfy no-slip boundary conditions: BIO(sigma) + bg_flow = 0
    sctl::Vector<Real> sigma;
    sctl::GMRES<Real> solver(comm, false);
    sctl::KrylovPrecond<Real> krylov;
    // sctl::Vector<Real> rhs = bg_flow(X0)*(pressure_drop/period_length);
    sctl::Vector<Real> rhs;
    if (peri_mode == 1 || peri_mode == 2) {
        rhs = bg_flow(X0)*(pressure_drop/period_length);
    } else {
        rhs = eval_rhs(pressure_drop);
    }
    
    solver(&sigma, BIO, rhs, gmres_tol, gmres_max_iter, false, nullptr, &krylov);
    // solver(&sigma, BIO_precond, AinvApply(rhs), gmres_tol, gmres_max_iter);

    { // Evaluate in interior, and write visualization
        PeriodicGeom<Real> trg;
        CubeVolumeVisShifted<Real> vol_vis(20, 0.9, comm); 
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);
        
        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);

        if (peri_mode == 1 || peri_mode == 2) {
            rhs = bg_flow(X0) * (pressure_drop/period_length);
            U -= rhs;
        } else {
            rhs = eval_rhs(pressure_drop);
            U -= rhs;
        }

        sctl::Vector<sctl::Long> size_loc(1);
        size_loc[0] = X0.Dim();
        sctl::Vector<sctl::Long> size_all(1);
        comm.Allreduce((sctl::Iterator<sctl::Long>) size_loc.begin(), (sctl::Iterator<sctl::Long>) size_all.begin(), 1, sctl::CommOp::SUM);
        // std::cout << "rank " << comm.Rank() << " size loc = " << size_loc[0] << ", size all is " << size_all[0] << std::endl;
        std::string filename = "Particle"+std::to_string(Nptcl)+"_"+std::to_string(peri_mode)+"_peri_U_exact_"+std::to_string(comm.Rank());
        std::string filename_out = "out/"+filename+".txt";
        std::string filename_vis = "vis/"+filename;
        if (write_ref) {
            
            U.Write(filename_out.c_str());
            // Create array of velocity for all target points, including filtered out ones.
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
            // Write visualization to VTK
            vol_vis.WriteVTK(filename_vis, U_vis);

        } else {

            sctl::Vector<Real> U_ref;
            U_ref.Read(filename_out.c_str());
            sctl::Vector<Real> err = U - U_ref;
            renormalize_error(err,comm);
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
                std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << ", Max u = " << u_all[0] << ", Max relative error = " << err_all[0] / u_all[0] << std::endl;
            }
        }
    }

}

// Conv div channel
template <class Real> void convdiv_self_conv(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Comm comm, const Real tol, const Real gmres_tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const sctl::Long ElemOrder = 10;
    const Real pressure_drop = -1.0;
    const Real period_length = 1;
    const sctl::Long geom_mode = 1;
    const sctl::Long gmres_max_iter = 400;

    PeriodicGeom<Real> obj;
    sctl::Long Nptcl = 3; // correct, but will be replaced inside conv-div channel build.
    sctl::Long ptcl_ord = 48; // Be careful with this if using particle preconditioner, since it's possible that not all panels on each particle are on the same MPI process.
    sctl::Vector<sctl::Long> ptcls(Nptcl);
    ptcls = ptcl_ord;
    // /////// DEBUG no particle in conv div, shoudl converge
    // sctl::Long Nptcl = 0;
    // sctl::Long ptcl_ord = 1;
    // sctl::Vector<sctl::Long> ptcls;
    // ////////////
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient, ptcls_thetas, ptcls_phis;
    // sctl::Long peri_mode = 1;
    Real channel_r1 = 0.1;
    Real channel_r2 = 0.2;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.build_conv_div_sph(Nelem, ElemOrder, FourierOrder, channel_r1, channel_r2, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    ptcls_thetas = std::get<2>(build0);
    ptcls_phis = std::get<3>(build0);
    Nptcl = ptcls_rs.Dim();

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);

    if (write_ref) {
        elem_lst0.WriteVTK("vis/channel", X0);
    }

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
    if (!comm.Rank()) {
        std::cout << "DEBUG surface area = " << surface_area << "." << std::endl;
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);

    // sctl::Matrix<Real> PrecondMat0, PrecondMat1, PrecondMat0_ptcl, PrecondMat1_ptcl;
    // sctl::Long A11size_ptcl, A11size;

    // Real channel_radius = 0.15; // for convdiv geometry
    // A11size = precond_channel(PrecondMat0, PrecondMat1, Nelem, ElemOrder, FourierOrder, channel_radius, SL_scal, DL_scal, comm);
    // A11size_ptcl = precond_ptcl(PrecondMat0_ptcl, PrecondMat1_ptcl, ptcl_ord, ElemOrder, FourierOrder, SL_scal, DL_scal, comm);

    // // Get global index of the starting panel on this process
    // sctl::Vector<sctl::Long> ElemOrderVec_temp(Nelem + ptcl_ord * Nptcl);
    // ElemOrderVec_temp = ElemOrder;
    // sctl::Vector<sctl::Long> FourierOrderVec_temp(ElemOrderVec_temp);
    // FourierOrderVec_temp = FourierOrder;
    // std::tuple<sctl::Long,sctl::Long> indtpl = obj.GetGlobalIdx(ElemOrderVec_temp, FourierOrderVec_temp, comm);
    // sctl::Long loc_elem_cnt = std::get<0>(indtpl);
    // sctl::Long loc_elem_dsp = std::get<1>(indtpl);

    const auto BIO = [&wts,&surface_area,&elem_lst0,&DL_scal,&LayerPotenOp0,NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
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

    // // Apply A11inv to each panel of a vector.
    // // Look at global panel index and determine whether belongs to a particle or the channel. ASSUMES no particles are split up among processors.
    // const auto AinvApply = [&PrecondMat0,&PrecondMat1,&A11size,&PrecondMat0_ptcl,&PrecondMat1_ptcl,&A11size_ptcl,&Nelem,&loc_elem_cnt,&loc_elem_dsp,&Nptcl](const sctl::Vector<Real>& vec) {
    //     sctl::Vector<Real> AinvVec(vec.Dim());
    //     if (Nptcl == 0 || (loc_elem_dsp+loc_elem_cnt) <= Nelem) {
    //         // std::cout << "All panels" << std::endl;
    //         // if no particles in channel or if all panels here are on channel
    //         sctl::Long N = vec.Dim();
    //         sctl::Long Npanels = N / A11size; 
    //         for (sctl::Long i=0; i<Npanels; i++) {
    //             sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
    //             sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
    //             for (sctl::Long j=0; j<A11size; j++) {
    //                 AinvVec[i*A11size + j] = AinvVecMat(j,0);
    //             }
    //         }
    //     } else {
    //         if (loc_elem_dsp >= Nelem) {
    //             std::cout << "All particles" << std::endl;
    //             // all panels here are ptcl
    //             sctl::Long N = vec.Dim();
    //             sctl::Long Nptcls = N / A11size_ptcl; 
    //             for (sctl::Long i=0; i<Nptcls; i++) {
    //                 sctl::Matrix<Real> vecMat(A11size_ptcl,1,(sctl::Iterator<Real>) vec.begin() + i*A11size_ptcl,true);
    //                 sctl::Matrix<Real> AinvVecMat = PrecondMat0_ptcl * (PrecondMat1_ptcl * vecMat);
    //                 for (sctl::Long j=0; j<A11size_ptcl; j++) {
    //                     AinvVec[i*A11size_ptcl + j] = AinvVecMat(j,0);
    //                 }
    //             }
    //         } else {
    //             sctl::Long Npanels_here = Nelem - loc_elem_dsp;
    //             sctl::Long Nptcls_here = loc_elem_cnt - Npanels_here;
    //             // // DEBUG
    //             std::cout << "CHECK panel-ptcl split: Npanel = " << Npanels_here << ", Nptcl = " << Nptcls_here << std::endl;
    //             // ///////////////
    //             for (sctl::Long i=0; i<Npanels_here; i++) {
    //                 sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
    //                 sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
    //                 for (sctl::Long j=0; j<A11size; j++) {
    //                     AinvVec[i*A11size + j] = AinvVecMat(j,0);
    //                 }
    //             }
    //             for (sctl::Long i=0; i<Nptcls_here; i++) {
    //                 sctl::Matrix<Real> vecMat(A11size_ptcl,1,(sctl::Iterator<Real>) vec.begin() + i*A11size_ptcl + Npanels_here*A11size,true);
    //                 sctl::Matrix<Real> AinvVecMat = PrecondMat0_ptcl * (PrecondMat1_ptcl * vecMat);
    //                 for (sctl::Long j=0; j<A11size_ptcl; j++) {
    //                     AinvVec[Npanels_here*A11size + i*A11size_ptcl + j] = AinvVecMat(j,0);
    //                 }
    //             }
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

    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    sctl::KrylovPrecond<Real> krylov;
    sctl::Vector<Real> rhs = bg_flow(X0) * (pressure_drop/period_length);
    // sctl::Vector<Real> A11invF = AinvApply(rhs);
    // solver(&sigma, BIO_precond, A11invF, gmres_tol, gmres_max_iter);
    solver(&sigma, BIO, rhs, gmres_tol, gmres_max_iter, false, nullptr, &krylov);

    { // Evaluate in interior, and write visualization
        PeriodicGeom<Real> trg;
        const sctl::Long Nelem_trg = 16;
        const sctl::Long FourierOrder_trg = 16;
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg; // dim = 0 so no particles are first generated
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build_trg = trg.build_conv_div_sph(Nelem_trg, ElemOrder, FourierOrder_trg, channel_r1, channel_r2, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, ptcl_ord);
        elem_lst_trg = std::get<0>(build_trg);

        VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        sctl::Vector<Real> X0_all = vol_vis.GetCoord(); // set new target coordinates
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target_rotated(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode, ptcls_thetas, ptcls_phis);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U -= bg_flow(X0) * (pressure_drop / period_length);

        // sctl::Vector<sctl::Long> size_loc(1);
        // size_loc[0] = X0.Dim();
        // sctl::Vector<sctl::Long> size_all(1);
        // comm.Allreduce((sctl::Iterator<sctl::Long>) size_loc.begin(), (sctl::Iterator<sctl::Long>) size_all.begin(), 1, sctl::CommOp::SUM);
        // std::cout << "rank " << comm.Rank() << " size loc = " << size_loc[0] << ", size all is " << size_all[0] << std::endl;

        // std::string filename = "Conv_div.txt";
        std::string filename = "ConvDiv_U_exact_"+std::to_string(comm.Rank());
        std::string filename_out = "out/"+filename+".txt";
        std::string filename_vis = "vis/"+filename;
        if (write_ref) {

            U.Write(filename_out.c_str());
            // Create array of velocity for all target points, including filtered out ones.
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
            // Write visualization to VTK
            vol_vis.WriteVTK(filename_vis, U_vis);

        } else {
            sctl::Vector<Real> U_ref;
            U_ref.Read(filename_out.c_str());
            sctl::Vector<Real> err = U - U_ref;
            renormalize_error(err,comm);
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
                std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << ", Max u = " << u_all[0] << ", Max relative error = " << err_all[0] / u_all[0] << std::endl;
            }

            // Try average error instead of max? (can probably use omp_scan.. )
            double avg_err = 0.;
            for (const auto e : err) avg_err += sctl::fabs(e);
            sctl::Vector<Real> avg_err_loc(1);
            avg_err_loc[0] = avg_err;
            sctl::Vector<Real> avg_err_all(1);
            avg_err_all[0] = 0;
            comm.Allreduce((sctl::Iterator<Real>) avg_err_loc.begin(), (sctl::Iterator<Real>) avg_err_all.begin(), 1, sctl::CommOp::SUM);
            std::cout << "on MPI process " << comm.Rank() << ", sum of err locally is " << avg_err << ", total error is " << avg_err_all[0] << std::endl;

            long size_err = err.Dim();
            sctl::Vector<sctl::Long> size_err_loc(1);
            size_err_loc[0] = size_err;
            sctl::Vector<sctl::Long> size_err_all(1);
            size_err_all[0] = 0;
            comm.Allreduce((sctl::Iterator<sctl::Long>) size_err_loc.begin(), (sctl::Iterator<sctl::Long>) size_err_all.begin(), 1, sctl::CommOp::SUM);

            std::cout << "size of error on process" << comm.Rank() << " is " << size_err << ", after collection total size is " << size_err_all[0] << std::endl;

            avg_err = avg_err_all[0] / (1.0*size_err_all[0]);

            if (!comm.Rank()) {
                std::cout<<"Average error = "<< std::setprecision(10) << avg_err << ", average relative error (divide by max u above) = " << avg_err / u_all[0] << std::endl;
            }

            // Create array of velocity for all target points, including filtered out ones.
            sctl::Vector<Real> err_vis(X0_all.Dim());
            err_vis = 0.;
            sctl::Long X1_ptr = 0;
            for (sctl::Long i=0; i<X0_all.Dim()/3; i++) {
                if (filtered_inds[i] == 0) {
                    err_vis[i*3] = err[X1_ptr*3];
                    err_vis[i*3+1] = err[X1_ptr*3+1];
                    err_vis[i*3+2] = err[X1_ptr*3+2];
                    X1_ptr += 1;
                }
            }
            // Write visualization to VTK
            vol_vis.WriteVTK("vis/ConvDiv_U_err", err_vis);
        }

    } 
}

template <class Real> void trefoil_ptcl_self_conv(sctl::Long Nelem, sctl::Long FourierOrder, const bool write_ref, sctl::Comm comm, const Real tol, const Real gmres_tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    // const Real tol = 1e-8;
    // const Real gmres_tol = 1e-13;
    const sctl::Long ElemOrder = 10;
    const Real pressure_drop = -1.0;
    const Real period_length = 1;
    const sctl::Long geom_mode = 0;
    const sctl::Long gmres_max_iter = 400;

    PeriodicGeom<Real> obj;
    sctl::Long Nptcl = 50; // placeholder; will be replaced inside conv-div channel build.
    sctl::Long ptcl_ord = 1;
    sctl::Vector<sctl::Long> ptcls(Nptcl);
    ptcls = ptcl_ord;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    // sctl::Long peri_mode = 1;
    // std::cout << "right before build trefoils" << std::endl;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_trefoil(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, geom_mode);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    Nptcl = ptcls_rs.Dim(); 

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    if (write_ref) {
        elem_lst0.WriteVTK("vis/trefoil_ptcl", X0);
    }

    Real surface_area;
    sctl::Vector<Real> wts;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        // surface_area = surface_area_[0];
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area = sa_all[0];
        if (!comm.Rank()) {
            std::cout << "Total surface area of channel + particles is " << surface_area  << std::endl;
        }
    }
    std::cout << "DEBUG surface area = " << surface_area << "." << std::endl;

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);

    // sctl::Matrix<Real> PrecondMat0, PrecondMat1, PrecondMat0_ptcl, PrecondMat1_ptcl;
    // sctl::Long A11size_ptcl, A11size;

    // Real channel_radius = 0.035;
    // A11size = precond_channel(PrecondMat0, PrecondMat1, Nelem, ElemOrder, FourierOrder, channel_radius, SL_scal, DL_scal, comm);
    // A11size_ptcl = precond_ptcl(PrecondMat0_ptcl, PrecondMat1_ptcl, ptcl_ord, ElemOrder, FourierOrder, SL_scal, DL_scal, comm);

    // // Get global index of the starting panel on this process
    // sctl::Vector<sctl::Long> ElemOrderVec_temp(Nelem + ptcl_ord * Nptcl);
    // ElemOrderVec_temp = ElemOrder;
    // sctl::Vector<sctl::Long> FourierOrderVec_temp(ElemOrderVec_temp);
    // FourierOrderVec_temp = FourierOrder;
    // std::tuple<sctl::Long,sctl::Long> indtpl = obj.GetGlobalIdx(ElemOrderVec_temp, FourierOrderVec_temp, comm);
    // sctl::Long loc_elem_cnt = std::get<0>(indtpl);
    // sctl::Long loc_elem_dsp = std::get<1>(indtpl);

    const auto BIO = [&wts,&surface_area,&elem_lst0,&DL_scal,&LayerPotenOp0,NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
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

            // DEBUG: check that sigma-sigma_mean has surface integral = 0:
            sctl::Vector<Real> sigma1 = sigma_;
            AddConstVec(sigma1, -sigma_mean);
            sctl::Vector<Real> sigma_test_;
            SurfaceIntegral(sigma_test_, sigma1, wts);
            std::cout << "Surface integral of sigma - sigma bar = " << sigma_test_[0] << ", "<< sigma_test_[1] << ", " << sigma_test_[2] << ". "<< std::endl;
        
        }
        
        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma0.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    // // Apply A11inv to each panel of a vector.
    // // Look at global panel index and determine whether belongs to a particle or the channel. ASSUMES no particles are split up among processors.
    // const auto AinvApply = [&PrecondMat0,&PrecondMat1,&A11size,&PrecondMat0_ptcl,&PrecondMat1_ptcl,&A11size_ptcl,&Nelem,&loc_elem_cnt,&loc_elem_dsp,&Nptcl](const sctl::Vector<Real>& vec) {
    //     sctl::Vector<Real> AinvVec(vec.Dim());
    //     if (Nptcl == 0 || (loc_elem_dsp+loc_elem_cnt) <= Nelem) {
    //         // std::cout << "All panels" << std::endl;
    //         // if no particles in channel or if all panels here are on channel
    //         sctl::Long N = vec.Dim();
    //         sctl::Long Npanels = N / A11size; 
    //         for (sctl::Long i=0; i<Npanels; i++) {
    //             sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
    //             sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
    //             for (sctl::Long j=0; j<A11size; j++) {
    //                 AinvVec[i*A11size + j] = AinvVecMat(j,0);
    //             }
    //         }
    //     } else {
    //         if (loc_elem_dsp >= Nelem) {
    //             std::cout << "All particles" << std::endl;
    //             // all panels here are ptcl
    //             sctl::Long N = vec.Dim();
    //             sctl::Long Nptcls = N / A11size_ptcl; 
    //             for (sctl::Long i=0; i<Nptcls; i++) {
    //                 sctl::Matrix<Real> vecMat(A11size_ptcl,1,(sctl::Iterator<Real>) vec.begin() + i*A11size_ptcl,true);
    //                 sctl::Matrix<Real> AinvVecMat = PrecondMat0_ptcl * (PrecondMat1_ptcl * vecMat);
    //                 for (sctl::Long j=0; j<A11size_ptcl; j++) {
    //                     AinvVec[i*A11size_ptcl + j] = AinvVecMat(j,0);
    //                 }
    //             }
    //         } else {
    //             sctl::Long Npanels_here = Nelem - loc_elem_dsp;
    //             sctl::Long Nptcls_here = loc_elem_cnt - Npanels_here;
    //             // // DEBUG
    //             std::cout << "CHECK panel-ptcl split: Npanel = " << Npanels_here << ", Nptcl = " << Nptcls_here << std::endl;
    //             // ///////////////
    //             for (sctl::Long i=0; i<Npanels_here; i++) {
    //                 sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
    //                 sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
    //                 for (sctl::Long j=0; j<A11size; j++) {
    //                     AinvVec[i*A11size + j] = AinvVecMat(j,0);
    //                 }
    //             }
    //             for (sctl::Long i=0; i<Nptcls_here; i++) {
    //                 sctl::Matrix<Real> vecMat(A11size_ptcl,1,(sctl::Iterator<Real>) vec.begin() + i*A11size_ptcl + Npanels_here*A11size,true);
    //                 sctl::Matrix<Real> AinvVecMat = PrecondMat0_ptcl * (PrecondMat1_ptcl * vecMat);
    //                 for (sctl::Long j=0; j<A11size_ptcl; j++) {
    //                     AinvVec[Npanels_here*A11size + i*A11size_ptcl + j] = AinvVecMat(j,0);
    //                 }
    //             }
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

    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    sctl::KrylovPrecond<Real> krylov;
    sctl::Vector<Real> rhs = bg_flow(X0)*pressure_drop / period_length;
    // sctl::Vector<Real> A11invF = AinvApply(bg_flow(X0) * pressure_drop / period_length);
    // solver(&sigma, BIO_precond, A11invF, gmres_tol, gmres_max_iter);
    solver(&sigma, BIO, rhs, gmres_tol, gmres_max_iter, false, nullptr, &krylov);
    
    { // Evaluate in interior, and write visualization
        PeriodicGeom<Real> trg;
        const sctl::Long Nelem_trg = 200;
        const sctl::Long FourierOrder_trg = 16;
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg;
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_trefoil(Nelem_trg, ElemOrder, FourierOrder_trg, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, 1, 0);
        elem_lst_trg = std::get<0>(build_trg);

        VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        sctl::Vector<Real> X0_all = vol_vis.GetCoord(); // set new target coordinates
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U -= bg_flow(X0) * pressure_drop / period_length;

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
            // Create array of velocity for all target points, including filtered out ones.
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
            // Write visualization to VTK
            vol_vis.WriteVTK(filename_vis, U_vis);

        } else {
            sctl::Vector<Real> U_ref;
            U_ref.Read(filename_out.c_str());
            sctl::Vector<Real> err = U - U_ref;
            renormalize_error(err,comm);
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
                std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << ", Max u = " << u_all[0] << ", Max relative error = " << err_all[0] / u_all[0] << std::endl;
            }
        }
    }
}

template <class Real> void plane_ptcl_self_conv(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, const sctl::Long Nptcl, const sctl::Long geom_mode, sctl::Comm comm, const Real tol, const Real gmres_tol) {
    SCTL_ASSERT(Nptcl == 1 || Nptcl == 3);

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const sctl::Long ElemOrder = 10;

    const Real period_length = 1.;
    const sctl::Long gmres_max_iter = 100;
    const Real pressure_drop = -1.;
    // const sctl::Long peri_mode = 2;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient, ptcls_thetas, ptcls_phis;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else if (geom_mode == 0) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls3(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode); // checking with three spheres
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.many_loops3(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
        ptcls_thetas = std::get<2>(build0);
        ptcls_phis = std::get<3>(build0);
    }

    sctl::Vector<Real> X0_ptcl; // target coordinates
    elem_lst0.GetNodeCoord(&X0_ptcl, nullptr, nullptr);
    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); 
    LayerPotenOp0.AddElemList(elem_lst0,"1");
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XY, period_length);
    
    // Plane
    // TODO: put plane on only Rank 0 MPI?
    sctl::Vector<Real> X0_wall;
    const sctl::Long gl_order = 49;
    const sctl::Long Nelem_x = 2;
    const sctl::Long Nelem_y = 2;
    const Real z_offset = 0.01;
    sctl::PlaneIntegral<Real> plane(gl_order, Nelem_x, Nelem_y, z_offset);
    
    plane.GetNodeCoord(&X0_wall, nullptr, nullptr);
    LayerPotenOp0.AddElemList(plane,"2"); // replacing all LPO2

    sctl::Vector<Real> X0;
    X0.ReInit(X0_ptcl.Dim() + X0_wall.Dim());
    for (int j=0; j<X0_ptcl.Dim(); j++) {
        X0[j] = X0_ptcl[j];
    }
    for (int j=0; j<X0_wall.Dim(); j++) {
        X0[j+X0_ptcl.Dim()] = X0_wall[j];
    }
    // Add plane normal orient as well 
    sctl::Vector<Real> NormalOrient_(NormalOrient.Dim() + X0_wall.Dim());
    NormalOrient_ = -1.; // Normal orient = -1 (-sign below) means all normals point into fluid (exterior problem)
    NormalOrient_.Swap(NormalOrient);
    LayerPotenOp0.SetTargetCoord(X0);

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
    Real surface_area_wall;
    sctl::Vector<Real> wts_wall;
    {
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        plane.GetFarFieldNodes(X, Xn, wts_wall, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts_wall*0+1, wts_wall);
        // surface_area_wall = surface_area_[0];
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area_wall = sa_all[0];
        std::cout << "DEBUG wall surface area, computed to be " << surface_area_wall << std::endl;
    }
    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&wts_wall,&surface_area_wall,&plane,&LayerPotenOp0,&DL_scal,&X0_ptcl,&X0_wall,&NormalOrient, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        
        sctl::Vector<Real> ptcl_dens(X0_ptcl.Dim(), (sctl::Iterator<Real>) sigma.begin(), true);
        sctl::Vector<Real> wall_dens(X0_wall.Dim(), (sctl::Iterator<Real>) sigma.begin()+ptcl_dens.Dim(), true);
            
        sctl::Vector<Real> sigma_mean, sigma0;

        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_mean_ptcl, sigma_mean_wall;
            sctl::Vector<Real> sigma_ptcl_;
            elem_lst0.GetFarFieldDensity(sigma_ptcl_, ptcl_dens);
            SurfaceIntegral(sigma_mean_ptcl, sigma_ptcl_, wts);
            sctl::Vector<Real> sigma_wall_ = wall_dens;
            SurfaceIntegral(sigma_mean_wall, sigma_wall_, wts_wall);
            // MPI
            sctl::Vector<Real> sa_loc(6);
            for (int i=0; i<3; i++) {
                sa_loc[i] = sigma_mean_ptcl[i];
            }
            for (int i=0; i<3; i++) {
                sa_loc[i+3] = sigma_mean_wall[i];
            }
            sctl::Vector<Real> sa_all(6);
            sa_all = 0.;
            for (int i=0; i<3; i++) {
                comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+i, (sctl::Iterator<Real>) sa_all.begin()+i, 1, sctl::CommOp::SUM);
            }
            for (int i=0; i<3; i++) {
                comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+i+3, (sctl::Iterator<Real>) sa_all.begin()+i+3, 1, sctl::CommOp::SUM);
            }
            // sigma_mean = sigma_mean_ptcl + sigma_mean_wall; 
            sigma_mean.ReInit(sigma_mean_ptcl.Dim());
            for (int i=0; i<3; i++) {
                sigma_mean[i] = sa_all[i] + sa_all[i+3]; // adding total sigma in {x,y,z} coord separately
            }
            sigma_mean *= (1/(surface_area + surface_area_wall));
            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);

        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma0.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    const auto bg_flow_2peri = [](const sctl::Vector<Real>& X) {
        const sctl::Long N = X.Dim()/3;
        sctl::Vector<Real> U(N*3);
        for (sctl::Long i = 0; i < N; i++) {
            const auto x = X.begin() + i*3;
            U[i*3+0] = - 0.5 * ((x[2]-0.5)*(x[2]-0.5)); // 2-periodic flow between plates.
            U[i*3+1] = 0.;
            U[i*3+2] = 0.;
        }
        return U;
    };

    // Solve for sigma to satisfy no-slip boundary conditions: BIO(sigma) + bg_flow = 0
    sctl::Vector<Real> sigma;
    sctl::GMRES<Real> solver(comm, false);
    sctl::KrylovPrecond<Real> krylov;
    sctl::Vector<Real> rhs = bg_flow_2peri(X0)*(pressure_drop/period_length);
    solver(&sigma, BIO, rhs, gmres_tol, gmres_max_iter, false, nullptr, &krylov);

    { // Evaluate in interior, and write visualization
        PeriodicGeom<Real> trg;
        CubeVolumeVisShifted<Real> vol_vis(10, 0.8, comm); // leave space away from planes to avoid close eval errors
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        if (geom_mode == 0) {
            std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
            X0 = std::get<0>(trg_tuple);
            filtered_inds = std::get<1>(trg_tuple);
        } else {
            std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target_rotated(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode, ptcls_thetas, ptcls_phis);
            X0 = std::get<0>(trg_tuple);
            filtered_inds = std::get<1>(trg_tuple);
        }
        
        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U -= bg_flow_2peri(X0) * (pressure_drop/period_length);

        sctl::Vector<sctl::Long> size_loc(1);
        size_loc[0] = X0.Dim();
        sctl::Vector<sctl::Long> size_all(1);
        comm.Allreduce((sctl::Iterator<sctl::Long>) size_loc.begin(), (sctl::Iterator<sctl::Long>) size_all.begin(), 1, sctl::CommOp::SUM);
        std::string filename = "Plane_ptcl_2_peri_U_exact_"+std::to_string(comm.Rank());
        std::string filename_out = "out/"+filename+".txt";
        std::string filename_vis = "vis/"+filename;
        if (write_ref) {
            
            U.Write(filename_out.c_str());
            // Create array of velocity for all target points, including filtered out ones.
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
            // Write visualization to VTK
            vol_vis.WriteVTK(filename_vis, U_vis);

        } else {

            sctl::Vector<Real> U_ref;
            U_ref.Read(filename_out.c_str());
            sctl::Vector<Real> err = U - U_ref;
            renormalize_error(err,comm);
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
                std::cout<<"Max error = "<< std::setprecision(10) << err_all[0] << ", Max u = " << u_all[0] << ", Max relative error = " << err_all[0] / u_all[0] << std::endl;
            }

            // Try average error instead of max? (can probably use omp_scan.. )
            double avg_err = 0.;
            for (const auto e : err) avg_err += sctl::fabs(e);
            sctl::Vector<Real> avg_err_loc(1);
            avg_err_loc[0] = avg_err;
            sctl::Vector<Real> avg_err_all(1);
            avg_err_all[0] = 0;
            comm.Allreduce((sctl::Iterator<Real>) avg_err_loc.begin(), (sctl::Iterator<Real>) avg_err_all.begin(), 1, sctl::CommOp::SUM);
            // std::cout << "on MPI process " << comm.Rank() << ", sum of err locally is " << avg_err << ", total error is " << avg_err_all[0] << std::endl;

            long size_err = err.Dim();
            sctl::Vector<sctl::Long> size_err_loc(1);
            size_err_loc[0] = size_err;
            sctl::Vector<sctl::Long> size_err_all(1);
            size_err_all[0] = 0;
            comm.Allreduce((sctl::Iterator<sctl::Long>) size_err_loc.begin(), (sctl::Iterator<sctl::Long>) size_err_all.begin(), 1, sctl::CommOp::SUM);

            // std::cout << "size of error on process" << comm.Rank() << " is " << size_err << ", after collection total size is " << size_err_all[0] << std::endl;

            avg_err = avg_err_all[0] / (1.0*size_err_all[0]);

            if (!comm.Rank()) {
                std::cout<<"Average error = "<< std::setprecision(10) << avg_err << ", average relative error (divide by max u above) = " << avg_err / u_all[0] << std::endl;
            }

            // Create array of velocity for all target points, including filtered out ones.
            sctl::Vector<Real> err_vis(X0_all.Dim());
            err_vis = 0.;
            sctl::Long X1_ptr = 0;
            for (sctl::Long i=0; i<X0_all.Dim()/3; i++) {
                if (filtered_inds[i] == 0) {
                    err_vis[i*3] = err[X1_ptr*3];
                    err_vis[i*3+1] = err[X1_ptr*3+1];
                    err_vis[i*3+2] = err[X1_ptr*3+2];
                    X1_ptr += 1;
                }
            }
            // Write visualization to VTK
            // vol_vis.WriteVTK("vis/Plane_ptcl_2_peri_U_err_", err_vis);
        }
    }

}

int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    // sctl::Profile::Enable(true);
    sctl::Comm comm = sctl::Comm::World();
    long test_mode = std::stol(argv[1]); // =0 for trefoil, =1 for particle; =2 for conv div, =3 for trefoil with particle. =4 for 2peri plane with particle(s); =5 for 1peri with sphere
    long peri_mode = std::stol(argv[2]); // 1-, 2-, or 3- periodic
    long Nptcl = std::stol(argv[3]); // Number of particles, for tests 1 and 4
    long geom_mode = std::stol(argv[4]); // type of particles, for test 4.
    std::cout << "Test mode: " << test_mode << "; peri mode (if test 1 for particles): " << peri_mode << "; Nptcl (if test 1 or test 4): " << Nptcl << "; geom mode (if test 4): " << geom_mode << std::endl;

    sctl::Vector<sctl::Long> Nelem_lst, FourierOrder_lst;

    if (test_mode == 0) { // trefoil empty
        Nelem_lst.PushBack(400);
        Nelem_lst.PushBack(600);
        Nelem_lst.PushBack(800);

        // FourierOrder_lst.PushBack(16);
        // FourierOrder_lst.PushBack(32);
        // FourierOrder_lst.PushBack(64);
        FourierOrder_lst.PushBack(80);
        FourierOrder_lst.PushBack(96); 
    } else if (test_mode == 1 || test_mode == 5) { // particle self conv or channle with one particle
        if (Nptcl < 10) {
            for (int i=1; i<6; i += 2) {
                Nelem_lst.PushBack(2*i); 
            }
            Nelem_lst.PushBack(18); // Np = 2, 6, 10, 18

            FourierOrder_lst.PushBack(4);
            FourierOrder_lst.PushBack(16);
            FourierOrder_lst.PushBack(32);
            FourierOrder_lst.PushBack(64);
            FourierOrder_lst.PushBack(80);
            FourierOrder_lst.PushBack(96);
        } else { // larger parameters for the more dense system of 25 particles
            for (int i=1; i<11; i += 2) {
                Nelem_lst.PushBack(2*i);
            }

            // FourierOrder_lst.PushBack(16); // cannot achieve desired tolerance
            FourierOrder_lst.PushBack(32);
            FourierOrder_lst.PushBack(64);
            FourierOrder_lst.PushBack(80);
            FourierOrder_lst.PushBack(96);
        }

        
    } else if (test_mode == 2) { // convdiv with particles
        // Nelem_lst.PushBack(25);
        // Nelem_lst.PushBack(50);
        Nelem_lst.PushBack(100);

        // FourierOrder_lst.PushBack(16);
        // FourierOrder_lst.PushBack(32);
        // FourierOrder_lst.PushBack(64);
        FourierOrder_lst.PushBack(80);
        FourierOrder_lst.PushBack(96); 
    } else if (test_mode == 3) { // trefoil with particles
        Nelem_lst.PushBack(400);
        Nelem_lst.PushBack(600);
        Nelem_lst.PushBack(800);

        // FourierOrder_lst.PushBack(16);
        // FourierOrder_lst.PushBack(32);
        // FourierOrder_lst.PushBack(64);
        FourierOrder_lst.PushBack(80);
        FourierOrder_lst.PushBack(96); 
    } else if (test_mode == 4) { // plane with particles
        if (geom_mode == 0) {
            for (int i=1; i<10; i += 2) {
            // for (int i=1; i<4; i += 2) {
                Nelem_lst.PushBack(2*i);
            }

            FourierOrder_lst.PushBack(16);
            FourierOrder_lst.PushBack(32);
            FourierOrder_lst.PushBack(64);
            FourierOrder_lst.PushBack(80);
            FourierOrder_lst.PushBack(96);
        } else {
            Nelem_lst.PushBack(8);
            Nelem_lst.PushBack(16);
            Nelem_lst.PushBack(32);
            Nelem_lst.PushBack(80);

            FourierOrder_lst.PushBack(16);
            FourierOrder_lst.PushBack(32);
            FourierOrder_lst.PushBack(64);
            FourierOrder_lst.PushBack(80);
            FourierOrder_lst.PushBack(96); 
        }        
    }
    
    sctl::Long Nelem, FourierOrder;
    for (int i=Nelem_lst.Dim()-1; i>=0; i--) {
        for (int j=FourierOrder_lst.Dim()-1; j>=0; j--) {
            Nelem = Nelem_lst[i];
            FourierOrder = FourierOrder_lst[j];
            if (!comm.Rank()) {
                std::cout << "Nelem = " << Nelem << ", FourierOrder = " << FourierOrder << "; " << std::endl;
            }
            if (i==Nelem_lst.Dim()-1 && j == FourierOrder_lst.Dim()-1) {
                Real tol = 1e-14;
                Real gmres_tol = 1e-14;
                if (test_mode==0) {
                    trefoil_self_conv<Real>(Nelem, FourierOrder, true, comm, tol, gmres_tol);
                } else if (test_mode == 1) {
                    particle_self_conv<Real>(Nelem, FourierOrder, true, peri_mode, comm, Nptcl, tol, gmres_tol);
                } else if (test_mode==2) {
                    convdiv_self_conv<Real>(Nelem, FourierOrder, true, comm, tol, gmres_tol); 
                } else if (test_mode==3) {
                    trefoil_ptcl_self_conv<Real>(Nelem, FourierOrder, true, comm, tol, gmres_tol);
                } else if (test_mode==4) {
                    plane_ptcl_self_conv<Real>(Nelem, FourierOrder, true, Nptcl, geom_mode, comm, tol, gmres_tol);
                } else if (test_mode==5) {
                    channel_sphere_self_conv<Real>(Nelem, FourierOrder, true, comm, tol, gmres_tol);
                } else {
                    SCTL_ASSERT(false);
                }
                // continue;
            } else {
                Real tol = 1e-14;
                Real gmres_tol;
                if (i < 3) {
                    gmres_tol = 1e-10;
                } else {
                    gmres_tol = 1e-12;
                }
                if (test_mode==0) {
                    trefoil_self_conv<Real>(Nelem, FourierOrder, false, comm, tol, gmres_tol);
                } else if (test_mode == 1){
                    particle_self_conv<Real>(Nelem, FourierOrder, false, peri_mode, comm, Nptcl, tol, gmres_tol); 
                } else if (test_mode==2) {
                    convdiv_self_conv<Real>(Nelem, FourierOrder, false, comm, tol, gmres_tol); 
                } else if (test_mode==3) {
                    trefoil_ptcl_self_conv<Real>(Nelem, FourierOrder, false, comm, tol, gmres_tol);
                } else if (test_mode==4) {
                    plane_ptcl_self_conv<Real>(Nelem, FourierOrder, false, Nptcl, geom_mode, comm, tol, gmres_tol);
                } else if (test_mode==5) {
                    channel_sphere_self_conv<Real>(Nelem, FourierOrder, false, comm, tol, gmres_tol);
                } else {
                    SCTL_ASSERT(false);
                }
                
            }
        }
    }
    
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}

