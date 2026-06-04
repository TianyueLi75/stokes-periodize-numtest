/*
    Self-convergence tests on singly-, doubly-, and triply- periodic examples.
*/

// Boundary integral operators
#include "stokes_bio.hpp" 

// Geometry for tests
#include "planeNaive.hpp"
#include "utils_geom.hpp"

// Other util functions
#include "utils_tests.cpp" 

// Visualization
#include "utils_vis.hpp" 

// self convergence test for a straight channel with a sphere
// If using MPI, will need to make sure particle panels will not be split among processes.
template <class Real> void channel_sphere_self_conv(
    const sctl::Long Nelem, 
    const sctl::Long FourierOrder, 
    const bool write_ref, 
    sctl::Comm comm, 
    const Real gmres_tol, 
    const Real tol) 
    {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real pressure_drop = -1.0;
    const Real period_length = 1;

    const sctl::Long ElemOrder = 10;
    const sctl::Long geom_mode = 0;
    const sctl::Long Nptcl = 1;
    const sctl::Long Nelem_channel = 20;
    const Real channel_radius = 0.2;

    const sctl::Long gmres_max_iter = 100;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Long ptcl_ord = Nelem;
    if (Nptcl>0) {
        ptcls.ReInit(Nptcl);
        ptcls = ptcl_ord;
    }
    sctl::Vector<Real> ptcls_Xcs, ptcls_rs, NormalOrient;
    sctl::SlenderElemList<Real> elem_lst0;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_straight(Nelem_channel, ElemOrder, FourierOrder, channel_radius, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);

    sctl::Vector<Real> X0, wts;
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Vector<Real> X0surf = X0;
    if (write_ref) {
        elem_lst0.WriteVTK("vis/channel_sphere",X0surf,comm);
    }  
    Real surface_area;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area = sa_all[0];
    }

    // Get block-preconditioner on channel and particle geometry
    sctl::Matrix<Real> PrecondMat0, PrecondMat1, PrecondMat0_ptcl, PrecondMat1_ptcl;
    sctl::Long A11size_ptcl, A11size;
    A11size = precond_channel(PrecondMat0, PrecondMat1, Nelem_channel, ElemOrder, FourierOrder, channel_radius, SL_scal, DL_scal, comm);
    A11size_ptcl = precond_ptcl(PrecondMat0_ptcl, PrecondMat1_ptcl, ptcl_ord, ElemOrder, FourierOrder, SL_scal, DL_scal, comm);

    // Get global index of the starting panel on this process
    sctl::Vector<sctl::Long> ElemOrderVec_temp(Nelem_channel + ptcl_ord * Nptcl);
    ElemOrderVec_temp = ElemOrder;
    sctl::Vector<sctl::Long> FourierOrderVec_temp(ElemOrderVec_temp);
    FourierOrderVec_temp = FourierOrder;
    std::tuple<sctl::Long,sctl::Long> indtpl = obj.GetGlobalIdx(ElemOrderVec_temp, FourierOrderVec_temp, comm);
    sctl::Long loc_elem_cnt = std::get<0>(indtpl);
    sctl::Long loc_elem_dsp = std::get<1>(indtpl);

    // Apply A11inv to each panel of a vector.
    // Look at global panel index and determine whether belongs to a particle or the channel. 
    // ASSUMES no particles are split up among processors; also assumes channel panels are filled in before particles.
    const auto AinvApply = [&PrecondMat0,&PrecondMat1,&A11size,&PrecondMat0_ptcl,&PrecondMat1_ptcl,&A11size_ptcl,&Nelem_channel,&ptcl_ord,&loc_elem_cnt,&loc_elem_dsp,&Nptcl](const sctl::Vector<Real>& vec) {
        sctl::Vector<Real> AinvVec(vec.Dim());
        if (Nptcl == 0 || (loc_elem_dsp+loc_elem_cnt) <= Nelem_channel) { // if no particles in channel or if all panels here are on channel
            // std::cout << "All panels on this MPI process" << std::endl;
            sctl::Long N = vec.Dim();
            sctl::Long Npanels = N / A11size; 
            for (sctl::Long i=0; i<Npanels; i++) {
                sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
                sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
                for (sctl::Long j=0; j<A11size; j++) {
                    AinvVec[i*A11size + j] = AinvVecMat(j,0);
                }
            }
        } else {
            if (loc_elem_dsp >= Nelem_channel) { // all panels here belong to particles
                // std::cout << "All particles on this MPI process" << std::endl;
                sctl::Long N = vec.Dim();
                sctl::Long Nptcls = N / A11size_ptcl; 
                for (sctl::Long i=0; i<Nptcls; i++) {
                    sctl::Matrix<Real> vecMat(A11size_ptcl,1,(sctl::Iterator<Real>) vec.begin() + i*A11size_ptcl,true);
                    sctl::Matrix<Real> AinvVecMat = PrecondMat0_ptcl * (PrecondMat1_ptcl * vecMat);
                    for (sctl::Long j=0; j<A11size_ptcl; j++) {
                        AinvVec[i*A11size_ptcl + j] = AinvVecMat(j,0);
                    }
                }
            } else {
                sctl::Long Npanels_here = Nelem_channel - loc_elem_dsp;
                sctl::Long Nptcls_here = (int) (loc_elem_cnt - Npanels_here) / ptcl_ord;
                // std::cout << "There are " << Npanels_here << " panels and " << Nptcls_here << " particles on this MPI process." << std::endl;
                for (sctl::Long i=0; i<Npanels_here; i++) {
                    sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
                    sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
                    for (sctl::Long j=0; j<A11size; j++) {
                        AinvVec[i*A11size + j] = AinvVecMat(j,0);
                    }
                }
                for (sctl::Long i=0; i<Nptcls_here; i++) {
                    sctl::Matrix<Real> vecMat(A11size_ptcl,1,(sctl::Iterator<Real>) vec.begin() + i*A11size_ptcl + Npanels_here*A11size,true);
                    sctl::Matrix<Real> AinvVecMat = PrecondMat0_ptcl * (PrecondMat1_ptcl * vecMat);
                    for (sctl::Long j=0; j<A11size_ptcl; j++) {
                        AinvVec[Npanels_here*A11size + i*A11size_ptcl + j] = AinvVecMat(j,0);
                    }
                }
            }
        }
        return AinvVec;
    };

    StokesBIO<Real> LayerPotenOp0(SL_scal, DL_scal, comm);
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetTargetCoord(X0surf);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);

    // peridozizedlayer potential operator
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
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    // Left preconditioning using block-preconditioner u -> A11inv*u
    const auto BIO_precond = [&BIO,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> Uloc;
        BIO(&Uloc,sigma);
        (*U) = AinvApply(Uloc);
    };

    // Set up B.C. on right hand side
    const sctl::Vector<Real> RHS = bg_flow_1peri(X0surf) * (pressure_drop/period_length);
    // Left precondition 
    sctl::Vector<Real> A11invF = AinvApply(RHS);

    sctl::Vector<Real> U, sigma;
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    solver(&sigma, BIO_precond, A11invF, gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);

    { // Evaluate in interior, and write visualization
        PeriodicGeom<Real> trg;
        VolumeVis<Real> vol_vis(elem_lst0, comm);
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        // Filter out targets not in fluid domain
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        LayerPotenOp0.SetTargetCoord(X0);
        BIO(&U, sigma);
        U -= bg_flow_1peri(X0) * (pressure_drop/period_length);

        std::string filename = "Channel_sphere_1peri_U_exact_"+std::to_string(comm.Rank());
        std::string filename_out = "out/"+filename+".txt";
        std::string filename_vis = "vis/"+filename;
        if (write_ref) {
            // Save solution to file as reference solution.
            U.Write(filename_out.c_str());

            // Create array of velocity for all target points, including filtered out ones, for visualization
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
            vol_vis.WriteVTK(filename_vis, U_vis);

        } else {
            // Grab reference solution
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

template <class Real> void particle_self_conv(
    const sctl::Long Nelem, 
    const sctl::Long FourierOrder, 
    const bool write_ref, 
    sctl::Comm comm,
    sctl::Long Nptcl, 
    const sctl::Integer peri_mode, 
    const Real gmres_tol, 
    const Real tol) 
    {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const sctl::Long ElemOrder = 10;
    const sctl::Long geom_mode = 0;
    
    const Real period_length = 1.;
    const Real pressure_drop = -1.;

    const sctl::Long gmres_max_iter = 200;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs, ptcls_rs, NormalOrient;
    sctl::SlenderElemList<Real> elem_lst0;
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

    sctl::Vector<Real> X0; 
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    Real surface_area;
    sctl::Vector<Real> wts;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area = sa_all[0];
    }

    // Get block-preconditioner on particle geometry
    sctl::Matrix<Real> PrecondMat0, PrecondMat1;
    sctl::Long A11size = precond_ptcl(PrecondMat0, PrecondMat1, Nelem, ElemOrder, FourierOrder, SL_scal, DL_scal, comm);
    
    // Apply A11inv to each panel of <vec>, assuming <vec> contains whole particles (i.e. vec.Dim() = A11size * Nptcl)
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
        }
        
        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma0.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    const auto BIO_precond = [&BIO,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> Uloc;
        BIO(&Uloc,sigma);
        (*U) = AinvApply(Uloc);
    };

    const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
        sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
        AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

        sctl::Vector<Real> U0;
        LayerPotenOp0.ComputeSL(U0, force_density);
        return U0;
    };

    sctl::Vector<Real> sigma;
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> rhs;
    if (peri_mode == 1) {
        rhs = bg_flow_1peri(X0)*(pressure_drop/period_length);
    } else if (peri_mode == 2) {
        rhs = bg_flow_2peri(X0)*(pressure_drop/period_length);
    } else if (peri_mode == 3) {
        rhs = eval_rhs(pressure_drop);
    } else {
        std::cout << "Only singly-, doubly-, or triply-periodic solvers supported." << std::endl;
        SCTL_ASSERT(false);
    }
    
    solver(&sigma, BIO_precond, AinvApply(rhs), gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);

    { // Evaluate in interior, and write visualization
        PeriodicGeom<Real> trg;
        CubeVolumeVisShifted<Real> vol_vis(20, 0.9, comm); 
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        // Filter out targets inside particles
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);
        
        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);

        if (peri_mode == 1) {
            rhs = bg_flow_1peri(X0) * (pressure_drop/period_length);
        } else if (peri_mode == 2) {
            rhs = bg_flow_2peri(X0) * (pressure_drop/period_length);
        } else {
            rhs = eval_rhs(pressure_drop);
        }
        U -= rhs;

        sctl::Vector<sctl::Long> size_loc(1);
        size_loc[0] = X0.Dim();
        sctl::Vector<sctl::Long> size_all(1);
        comm.Allreduce((sctl::Iterator<sctl::Long>) size_loc.begin(), (sctl::Iterator<sctl::Long>) size_all.begin(), 1, sctl::CommOp::SUM);
        std::string filename = "Particle"+std::to_string(Nptcl)+"_"+std::to_string(peri_mode)+"_peri_U_exact_"+std::to_string(comm.Rank());
        std::string filename_out = "out/"+filename+".txt";
        std::string filename_vis = "vis/"+filename;
        if (write_ref) {
            // Save solution as reference.
            U.Write(filename_out.c_str());

            // Create array of velocity for all target points, including filtered out ones, for visualization
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
            vol_vis.WriteVTK(filename_vis, U_vis);

        } else {
            // Grab reference solution
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

template <class Real> void plane_ptcl_self_conv(
    const sctl::Long Nelem, 
    const sctl::Long FourierOrder, 
    const bool write_ref,
    sctl::Comm comm,
    const sctl::Long Nptcl, 
    const sctl::Long geom_mode, 
    const Real gmres_tol, 
    const Real tol) 
    {

    SCTL_ASSERT(comm.Size() == 1); // This code does not support multi-process. See examples.cpp for parallelized plane-bound flow code.

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const sctl::Long ElemOrder = 10;

    const Real period_length = 1.;
    const Real pressure_drop = -1.;

    const sctl::Long gmres_max_iter = 200;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs, ptcls_rs, NormalOrient, ptcls_thetas, ptcls_phis;
    sctl::SlenderElemList<Real> elem_lst0;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else if (geom_mode == 0) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls3(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
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
    sctl::Vector<Real> X0_wall;
    const sctl::Long gl_order = 49;
    const sctl::Long Nelem_x = 2;
    const sctl::Long Nelem_y = 2;
    const Real z_offset = 0.01; // shift planes slightly away from z=0,1 to avoid overlapping the unit cube
    sctl::PlaneIntegral<Real> plane(gl_order, Nelem_x, Nelem_y, z_offset);
    
    plane.GetNodeCoord(&X0_wall, nullptr, nullptr);
    LayerPotenOp0.AddElemList(plane,"2"); 

    sctl::Vector<Real> X0;
    X0.ReInit(X0_ptcl.Dim() + X0_wall.Dim());
    for (int j=0; j<X0_ptcl.Dim(); j++) {
        X0[j] = X0_ptcl[j];
    }
    for (int j=0; j<X0_wall.Dim(); j++) {
        X0[j+X0_ptcl.Dim()] = X0_wall[j];
    }
    sctl::Vector<Real> NormalOrient_(NormalOrient.Dim() + X0_wall.Dim());
    NormalOrient_ = -1.; 
    NormalOrient_.Swap(NormalOrient);

    LayerPotenOp0.SetTargetCoord(X0);

    Real surface_area;
    sctl::Vector<Real> wts;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area = sa_all[0];
    }
    Real surface_area_wall;
    sctl::Vector<Real> wts_wall;
    {
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        plane.GetFarFieldNodes(X, Xn, wts_wall, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts_wall*0+1, wts_wall);
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area_wall = sa_all[0];
    }

    // periodized boundary integral operator
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

    // Solve for sigma to satisfy no-slip boundary conditions: BIO(sigma) + bg_flow = 0
    sctl::Vector<Real> sigma;
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> rhs = bg_flow_2peri(X0)*(pressure_drop/period_length);
    solver(&sigma, BIO, rhs, gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);

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

int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    // sctl::Profile::Enable(true);
    sctl::Comm comm = sctl::Comm::World();
    long test_mode = std::stol(argv[1]); // =0 for 1-peri channel-sphere, =1 for 2-peri plane-particles; =2 for triply-periodic spheres
    long peri_mode = std::stol(argv[2]); // peri_mode = j for j-periodic
    long Nptcl = std::stol(argv[3]); // Number of particles
    long geom_mode = std::stol(argv[4]); // type of particles, for test 4. sphere: 0, sheroids: 1, loops: 2
    std::cout << "Test mode: " << test_mode << "; peri mode: " << peri_mode << "; Nptcl: " << Nptcl << "; geom mode: " << geom_mode << std::endl;

    sctl::Vector<sctl::Long> Nelem_lst, FourierOrder_lst;

    if (test_mode == 0) {
        for (int i=1; i<6; i += 2) {
            Nelem_lst.PushBack(2*i);
        }
        // Nelem_lst.PushBack(18);

        FourierOrder_lst.PushBack(4);
        FourierOrder_lst.PushBack(16);
        FourierOrder_lst.PushBack(32);
        FourierOrder_lst.PushBack(64);
        // FourierOrder_lst.PushBack(96);
        
    } else if (test_mode == 1) {
        if (geom_mode == 0) {
            for (int i=1; i<10; i += 2) {
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
    } else if (test_mode == 2) {
        if (Nptcl < 10) {
            for (int i=1; i<4; i += 2) {
                Nelem_lst.PushBack(2*i);
            }
            FourierOrder_lst.PushBack(4);
            FourierOrder_lst.PushBack(16);
            FourierOrder_lst.PushBack(32);
            // FourierOrder_lst.PushBack(64);
        } else { // larger parameters for the more dense system of 25 particles
            for (int i=1; i<11; i += 2) {
                Nelem_lst.PushBack(2*i);
            }
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
                    channel_sphere_self_conv<Real>(Nelem, FourierOrder, true, comm, gmres_tol, tol);
                } else if (test_mode == 1) {
                    plane_ptcl_self_conv<Real>(Nelem, FourierOrder, true, comm, Nptcl, geom_mode, gmres_tol, tol);
                } else if (test_mode==2) {
                    particle_self_conv<Real>(Nelem, FourierOrder, true, comm, Nptcl, peri_mode, gmres_tol, tol); 
                } else {
                    SCTL_ASSERT(false);
                }
            } else {
                Real tol = 1e-14;
                Real gmres_tol;
                if (i < 3) {
                    gmres_tol = 1e-10;
                } else {
                    gmres_tol = 1e-12;
                }
                if (test_mode==0) {
                    channel_sphere_self_conv<Real>(Nelem, FourierOrder, false, comm, gmres_tol, tol);
                } else if (test_mode == 1) {
                    plane_ptcl_self_conv<Real>(Nelem, FourierOrder, false, comm, Nptcl, geom_mode, gmres_tol, tol);
                } else if (test_mode==2) {
                    particle_self_conv<Real>(Nelem, FourierOrder, false, comm, Nptcl, peri_mode, gmres_tol, tol); 
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

