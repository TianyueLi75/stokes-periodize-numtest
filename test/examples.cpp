// =============================================================================
// examples.cpp
//
// Worked examples demonstrating the solver on complex bounded periodic geometries, 
// with volumetric flow fields written to vis/ for rendering (e.g. in ParaView).
//
// Usage:
//   make examples
//   mpirun -n <Nproc> --map-by numa:pe=$OMP_NUM_THREADS ./bin/examples <mode>
//
// Arguments:
//   mode  0 = singly-periodic converging-diverging channel seeded with
//             rotated spheroids (slip boundary condition on the particles)
//        1 = doubly-periodic flow between two flat walls past toroidal loops,
//             driven by a background pressure drop
//   Example:   ./bin/examples 0
//
// Method:
//   The combined-field BIE is solved with the surface-mean projection over the
//   union of the bounding surface (channel wall or flat planes) and the interior
//   particles; flat walls are represented by the PlaneIntegral element list.
//   After solving, the velocity is evaluated on an interior target grid with
//   points inside the particles filtered out, and written as VTK.
// =============================================================================

// Boundary integral operators
#include "stokes_bio.hpp" 

// Geometry for tests
#include "planeNaive.hpp"
#include "utils_geom.hpp"

// Other util functions
#include "utils_tests.cpp" 

// Visualization
#include "utils_vis.hpp" 

// Singly-periodic converging-diverging channel with interior spheroids under a slip BC;
// writes the geometry, slip BC, and solution flow to vis/.
template <class Real> void channel_with_particle(sctl::Comm comm) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    // Set quadrature parameters
    const sctl::Long Nelem_channel = 62;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 48;
    const Real tol = 1e-8; // quadrature tolerance

    // Set geometry parameters
    const Real period_length = 1; // length of periodic box.
    const Real channel_r1 = 0.025; 
    const Real channel_r2 = 0.1;
    const sctl::Long Nspheroids_start = 2560; // Use the unit cube geom of <Nspheroids_start> suspension then keep only spheroids inside channel.
    // Set GMRES parameters
    const Real gmres_tol = 1e-7;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls(Nspheroids_start); 
    sctl::Long ptcl_ord = 2; // number of panels on each particle (same fourier order as channel)
    ptcls = ptcl_ord;
    // Arrays to hold particle properties.
    sctl::Vector<Real> ptcls_Xcs, ptcls_rs, ptcls_u0s, ptcls_thetas, ptcls_phis, NormalOrient; 
    sctl::Vector<sctl::Long> ptcls_ifprolate;
    sctl::SlenderElemList<Real> elem_lst0; 
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.build_conv_div(Nelem_channel, ElemOrder, FourierOrder, channel_r1, channel_r2, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcls_u0s, ptcls_ifprolate, ptcl_ord, Nspheroids_start);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    ptcls_thetas = std::get<2>(build0); 
    ptcls_phis = std::get<3>(build0);
    sctl::Long Nptcl = ptcls.Dim();
    std::cout << "done with forming conv div channel, number of particles is " << ptcls.Dim() << std::endl;

    sctl::Vector<Real> X0, Xnsurf; // target coordinates
    elem_lst0.GetNodeCoord(&X0, &Xnsurf, nullptr);
    elem_lst0.WriteVTK("vis/ConvDiv_geometry", Xnsurf, comm);

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

    // Get block-preconditioner on channel and particle geometry
    sctl::Matrix<Real> PrecondMat0, PrecondMat1, PrecondMat0_ptcl, PrecondMat1_ptcl;
    sctl::Long A11size_ptcl, A11size;
    A11size = precond_channel(PrecondMat0, PrecondMat1, Nelem_channel, ElemOrder, FourierOrder, channel_r1, SL_scal, DL_scal, comm);
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
            std::cout << "All panels on this MPI process" << std::endl;
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
                std::cout << "All particles on this MPI process" << std::endl;
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
                std::cout << "There are " << Npanels_here << " panels and " << Nptcls_here << " particles on this MPI process." << std::endl;
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

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); 
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
        }
        
        // Compute the periodic solution using sigma0
        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma0.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer
    
        // Add back sigma_mean
        AddConstVec(*U, sigma_mean);
    };
    
    // Left preconditioning using block-preconditioner u -> A11inv*u
    const auto BIO_precond = [&BIO,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> Uloc;
        BIO(&Uloc,sigma);
        (*U) = AinvApply(Uloc);
    };

    sctl::GMRES<Real> solver(comm);
    // sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> sigma;

    // Slip velocity on particles.
    sctl::Long Nnodes_channel = Nelem_channel * ElemOrder * FourierOrder;
    sctl::Long Nnodes_per_ptcl = ptcl_ord * ElemOrder * FourierOrder;
    sctl::Long Nelem_tot = Nelem_channel + ptcl_ord * ptcls.Dim();
    sctl::Long Nelem_this_mpi = static_cast<int>(Nelem_tot / comm.Size());
    sctl::Long rmd = Nelem_tot - Nelem_this_mpi * comm.Size();
    sctl::Long elem_dsp_this_mpi = Nelem_this_mpi * comm.Rank();
    if (comm.Rank() < rmd) {
        Nelem_this_mpi += 1;
        elem_dsp_this_mpi += comm.Rank();
    } else {
        elem_dsp_this_mpi += rmd;
    }
    // Assumes: all channel panels on MPI rank 0; 
    //          no particle is split between two processes
    sctl::Long Nptcl_this_mpi, Nptcl_dsp_this_mpi;
    if (!comm.Rank()) {
        Nptcl_this_mpi = X0.Dim() - 3*Nnodes_channel;
        Nptcl_this_mpi = static_cast<int>(Nptcl_this_mpi / 3 / Nnodes_per_ptcl);
        Nptcl_dsp_this_mpi = 0;
    } else {
        Nptcl_this_mpi = static_cast<int>(X0.Dim() / Nnodes_per_ptcl / 3);
        Nptcl_dsp_this_mpi = static_cast<int>((elem_dsp_this_mpi - Nelem_channel) / ptcl_ord);
    }
    sctl::Vector<Real> rhs(X0.Dim());
    rhs.SetZero();
    sctl::Vector<Real> vslip_ptcl;
    if (!comm.Rank()) {
        const sctl::Vector<Real> Xnptcl(3 * Nnodes_per_ptcl * Nptcl_this_mpi, (sctl::Iterator<Real>) Xnsurf.begin() + 3*Nnodes_channel, false);
        sctl::Vector<Real> ptcls_Xcs_here(3*Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_Xcs.begin() + 3*Nptcl_dsp_this_mpi, true);
        sctl::Vector<Real> ptcls_sizes_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_rs.begin() + Nptcl_dsp_this_mpi, true);
        sctl::Vector<Real> ptcls_u0s_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_u0s.begin() + Nptcl_dsp_this_mpi, true);
        sctl::Vector<Real> ptcls_thetas_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_thetas.begin() + Nptcl_dsp_this_mpi, true);
        sctl::Vector<Real> ptcls_phis_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_phis.begin() + Nptcl_dsp_this_mpi, true);
        sctl::Vector<sctl::Long> ptcls_ifprolate_here(Nptcl_this_mpi, (sctl::Iterator<sctl::Long>) ptcls_ifprolate.begin() + Nptcl_dsp_this_mpi, true);
        vslip_ptcl = GetVslip(Xnptcl, ptcls_Xcs_here, ptcls_sizes_here, ptcls_u0s_here, ptcls_thetas_here, ptcls_phis_here, ptcls_ifprolate_here, ptcl_ord, ElemOrder, FourierOrder);
        for (int j=3*Nnodes_channel; j<rhs.Dim(); j++) {
            rhs[j] = vslip_ptcl[j-3*Nnodes_channel];
        }
    } else {
        const sctl::Vector<Real> Xptcl = X0;
        sctl::Vector<Real> ptcls_Xcs_here(3*Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_Xcs.begin() + 3*Nptcl_dsp_this_mpi, true);
        sctl::Vector<Real> ptcls_sizes_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_rs.begin() + Nptcl_dsp_this_mpi, true);
        sctl::Vector<Real> ptcls_u0s_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_u0s.begin() + Nptcl_dsp_this_mpi, true);
        sctl::Vector<Real> ptcls_thetas_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_thetas.begin() + Nptcl_dsp_this_mpi, true);
        sctl::Vector<Real> ptcls_phis_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_phis.begin() + Nptcl_dsp_this_mpi, true);
        sctl::Vector<sctl::Long> ptcls_ifprolate_here(Nptcl_this_mpi, (sctl::Iterator<sctl::Long>) ptcls_ifprolate.begin() + Nptcl_dsp_this_mpi, true);
        vslip_ptcl = GetVslip(Xnsurf, ptcls_Xcs_here, ptcls_sizes_here, ptcls_u0s_here, ptcls_thetas_here, ptcls_phis_here, ptcls_ifprolate_here, ptcl_ord, ElemOrder, FourierOrder);
        rhs = vslip_ptcl;
    }
    elem_lst0.WriteVTK("vis/ConvDiv_slip_BC", rhs, comm);
    // Check that u dot n is always zero
    for (int ii=0; ii<X0.Dim()/3; ii++) {
        const sctl::Vector<Real> vslip_here(3, (sctl::Iterator<Real>) rhs.begin() + ii*3, false);
        const sctl::Vector<Real> xn_here(3, (sctl::Iterator<Real>) Xnsurf.begin() + ii*3, false);
        const Real udotn = vslip_here[0]*xn_here[0] + vslip_here[1]*xn_here[1] + vslip_here[2]*xn_here[2];
        if (sctl::fabs(udotn) > 1e-8) {
            std::cout << "ERROR: u dot n at node " << ii << " is nonzero: " << udotn << ". Normal is " << xn_here[0] << ", " << xn_here[1] << ", " << xn_here[2] << ", vslip is " << vslip_here[0] << ", " << vslip_here[1] << ", " << vslip_here[2] << std::endl;
            return;
        }
    }

    sctl::Vector<Real> A11invF = AinvApply(rhs);
    solver(&sigma, BIO_precond, A11invF, gmres_tol);

    { 
        // Create uniform grid in cylindrical coordinates for targets inside channel
        PeriodicGeom<Real> trg;
        const sctl::Long Nelem_trg = 80;
        const sctl::Long FourierOrder_trg = 64;
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg; 
        sctl::Vector<Real> ptcls_Xcs_trg, ptcls_rs_trg, ptcls_u0s_trg;
        sctl::Vector<sctl::Long> ptcls_ifprolate_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build_trg = trg.build_conv_div(Nelem_trg, ElemOrder, FourierOrder_trg, channel_r1, channel_r2, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, ptcls_u0s_trg, ptcls_ifprolate_trg, ptcl_ord);
        elem_lst_trg = std::get<0>(build_trg);

        VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        sctl::Vector<Real> X0_all = vol_vis.GetCoord(); 
        // Filter out targets inside particles
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_spheroids_rotated(X0_all, ptcls_rs, ptcls_u0s, ptcls_Xcs, ptcls_ifprolate, ptcls_thetas, ptcls_phis);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        // Evaluate solution flow at targets
        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);

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
        vol_vis.WriteVTK("vis/ConvDiv_slip_U", U_vis); 
    }
}

// Doubly-periodic array of toroidal loops between two infinite flat walls, driven by a
// background pressure drop; writes the solution flow to vis/.
template <class Real> void planes_with_loops(sctl::Comm comm) {

    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    // Set quadrature parameters
    const sctl::Long Nelem = 40; // N_p on loops
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 96; // N_f on loops
    const Real tol = 1e-14; // quadrature tolerance
    const sctl::Long gl_order = 49; // *Cheb order on each panel of plane
    const sctl::Long Nelem_x = 2; // Number of panels in x .. 
    const sctl::Long Nelem_y = 2; // .. and y directions
    const Real z_offset = 0.005; // planes located at z = <z_offset> and z = 1 - <z_offset>
    
    // Set geometry parameters
    const Real pressure_drop = -1.0; // background pressure driven flow.
    const Real period_length = 1; // length of periodic box.
    
    // Set GMRES parameters
    const Real gmres_tol = 1e-8;
    const sctl::Long gmres_max_iter = 250;

    // Generate random sized loops from data file of system of spheres in unit cube
    const sctl::Long Nptcls = 25;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_rs, ptcls_Xcs, ptcls_major_rs, ptcls_minor_rs, ptcls_thetas, ptcls_phis;
    PeriodicGeom<Real> obj;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_loops2(Nelem, ElemOrder, FourierOrder, comm, Nptcls, ptcls, ptcls_rs, ptcls_Xcs, ptcls_major_rs, ptcls_minor_rs, ptcls_thetas, ptcls_phis);
    sctl::SlenderElemList<Real> elem_lst0 = std::get<0>(build0);
    sctl::Vector<Real> NormalOrient = std::get<1>(build0);

    sctl::Vector<Real> X0_ptcl;
    elem_lst0.GetNodeCoord(&X0_ptcl, nullptr, nullptr);
    elem_lst0.WriteVTK("vis/PlanesPtcl_ptcls", X0_ptcl, comm);

    // Create plane object on all process, but only use it on Rnak 0
    sctl::PlaneIntegral<Real> plane(gl_order, Nelem_x, Nelem_y, z_offset);

    sctl::Vector<Real> X0_wall;
    if (!comm.Rank()) {
        plane.GetNodeCoord(&X0_wall, nullptr, nullptr);
    } else {
        X0_wall.ReInit(0);
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); 
    if (!comm.Rank()) {
        LayerPotenOp0.AddElemList(plane,"2"); 
    }
    LayerPotenOp0.AddElemList(elem_lst0,"1");
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XY, period_length);
    
    sctl::Vector<Real> X0;
    if (!comm.Rank()) {
        X0.ReInit(X0_ptcl.Dim() + X0_wall.Dim());
        for (int j=0; j<X0_ptcl.Dim(); j++) {
            X0[j] = X0_ptcl[j];
        }
        for (int j=0; j<X0_wall.Dim(); j++) {
            X0[j+X0_ptcl.Dim()] = X0_wall[j];
        }
    } else {
        X0 = X0_ptcl;
    }
    LayerPotenOp0.SetTargetCoord(X0);

    sctl::Vector<Real> NormalOrient_(NormalOrient.Dim() + X0_wall.Dim());
    NormalOrient_ = -1.; 
    NormalOrient_.Swap(NormalOrient);

    Real surface_area_wall, total_surface_area;
    sctl::Vector<Real> wts, wts_wall;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        if (!comm.Rank()) {
            // Add plane area
            sctl::Vector<Real> X_wall, Xn_wall, dist_far_wall, surface_area_wall_;
            sctl::Vector<sctl::Long> element_wise_node_cnt_wall;
            plane.GetFarFieldNodes(X_wall, Xn_wall, wts_wall, dist_far_wall, element_wise_node_cnt_wall, 1);
            SurfaceIntegral(surface_area_wall_, wts_wall*0+1, wts_wall);
            surface_area_wall = surface_area_wall_[0];
            surface_area_[0] += surface_area_wall;
        }
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0]; // rank 0 process will contribute the planes area.
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        total_surface_area = sa_all[0];
    }

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&total_surface_area,&elem_lst0,&wts_wall,&LayerPotenOp0,&DL_scal,&X0_ptcl,&X0_wall,&NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {

        sctl::Vector<Real> sigma_mean, sigma0;
        sctl::Vector<Real> sa_loc;
        if (!comm.Rank()) {
            sctl::Vector<Real> ptcl_dens(X0_ptcl.Dim(), (sctl::Iterator<Real>) sigma.begin(), true);
            sctl::Vector<Real> wall_dens(X0_wall.Dim(), (sctl::Iterator<Real>) sigma.begin()+X0_ptcl.Dim(), true);
            sctl::Vector<Real> sigma_mean_ptcl, sigma_mean_wall;
            sctl::Vector<Real> sigma_ptcl_;
            elem_lst0.GetFarFieldDensity(sigma_ptcl_, ptcl_dens);
            SurfaceIntegral(sigma_mean_ptcl, sigma_ptcl_, wts);
            sctl::Vector<Real> sigma_wall_ = wall_dens;
            SurfaceIntegral(sigma_mean_wall, sigma_wall_, wts_wall);
            sigma_mean = sigma_mean_ptcl + sigma_mean_wall; 
            sa_loc = sigma_mean;
        } else {
            sctl::Vector<Real> sigma_ptcl_;
            elem_lst0.GetFarFieldDensity(sigma_ptcl_, sigma);
            SurfaceIntegral(sigma_mean, sigma_ptcl_, wts);
            sa_loc = sigma_mean;
        }
        sctl::Vector<Real> sa_all(3);
        sa_all = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+1, (sctl::Iterator<Real>) sa_all.begin()+1, 1, sctl::CommOp::SUM);
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+2, (sctl::Iterator<Real>) sa_all.begin()+2, 1, sctl::CommOp::SUM);
        sigma_mean = sa_all;
        sigma_mean *= 1./total_surface_area;
        sigma0 = sigma;
        AddConstVec(sigma0, -sigma_mean);

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma0.Dim()){
            (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer
        } 

        AddConstVec(*U, sigma_mean);
    };

    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    sctl::KrylovPrecond<Real> krylov_precond;
    solver(&sigma,BIO, bg_flow_2peri(X0) * (pressure_drop/period_length), gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);

    {
        PeriodicGeom<Real> trg;    
        CubeVolumeVisShifted<Real> vol_vis(60, 0.9, comm);
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_loops_rotated(X0_all, ptcls_major_rs, ptcls_minor_rs, ptcls_Xcs, ptcls_thetas, ptcls_phis);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U(X0.Dim());
        BIO(&U, sigma);
        U -= bg_flow_2peri(X0) * (pressure_drop/period_length);
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
        vol_vis.WriteVTK("vis/PlanesPtcl_U", U_vis);
    }
}


int main(int argc, char** argv) {

    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Comm comm = sctl::Comm::World();
        long mode = std::stol(argv[1]);

        if (mode == 0) {
            channel_with_particle<Real>(comm);
        } else if (mode == 1) {
            planes_with_loops<Real>(comm);
        } else {
            std::cout << "Only supports test mode 0: channel-bound singly-periodic flow past obstacles, or mode 1: plane-bound doubly periodic flow through loops." << std::endl;
            SCTL_ASSERT(false);
        }
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
