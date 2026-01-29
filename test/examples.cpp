#include "periodize.hpp"
#include "utils.hpp"
#include "planeNaive.hpp"

// Test script for calculation and timing of 1, 2, and 3 periodic problems with background pressure flow. 

/**
    Takes the surface integral, populate into <I>, using weights <wts>   
    Supporting functions for imposing net-force-zero densities during gmres solve.
*/
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

/**
    Assuming <vals> is a dimension-fast-node-slow list of values, add <c0> to each node (if dimension matches).
    Supporting functions for imposing net-force-zero densities during gmres solve.
*/
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

/**
    Example of a converging-divering channel with some spheroids on the interior. Visualizations are stored in vis/ folder.
*/
template <class Real> void channel_with_particle(sctl::Comm comm) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;
    // Set quadrature parameters
    const sctl::Long Nelem_channel = 100;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 64;
    const sctl::Long geom_mode = 1; 
    const Real tol = 1e-8; // quadrature tolerance
    // Set geometry parameters
    const Real pressure_drop = -1.0; // background pressure driven flow.
    const Real period_length = 1; // length of periodic box.
    const Real channel_r1 = 0.1; 
    const Real channel_r2 = 0.2;
    // Set GMRES parameters
    const Real gmres_tol = 1e-8;
    const sctl::Long gmres_max_iter = 400;

    // Create channel + particles object: "elem_lst0".
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls(1); // length of ptcls will get updated through the conv-div channel constructor. (this is not the case for all channel constructors.)
    sctl::Long ptcl_ord = 4; // number of panels on each particle (same fourier order as channel)
    ptcls = ptcl_ord;
    // Arrays to hold particle...
    sctl::Vector<Real> ptcls_Xcs; // ...center positions
    sctl::Vector<Real> ptcls_rs; // ...radii
    sctl::SlenderElemList<Real> elem_lst0; 
    sctl::Vector<Real> NormalOrient; // -1 if outward normal on surface points into fluid, +1 otherwise. (since BIO defined with a -1/2*sigma already)
    sctl::Vector<Real> ptcls_thetas, ptcls_phis;
    sctl::Long peri_mode = 1;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.build_conv_div(Nelem_channel, ElemOrder, FourierOrder, channel_r1, channel_r2, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    ptcls_thetas = std::get<2>(build0);
    ptcls_phis = std::get<3>(build0);

    // TODO: change this to conv.div channel but with a couple of differently shaped particles

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    elem_lst0.WriteVTK("vis/ConvDiv_geometry", X0);
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

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); 
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,&NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> sigma_mean, sigma0;
        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_;
            elem_lst0.GetFarFieldDensity(sigma_, sigma);
            SurfaceIntegral(sigma_mean, sigma_, wts);
            // MPI for total surface area
            sctl::Vector<Real> sa_loc = sigma_mean;
            sctl::Vector<Real> sa_all(3);
            sa_all = 0;
            comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
            comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+1, (sctl::Iterator<Real>) sa_all.begin()+1, 1, sctl::CommOp::SUM);
            comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+2, (sctl::Iterator<Real>) sa_all.begin()+2, 1, sctl::CommOp::SUM);
            sigma_mean = sa_all;
            sigma_mean *= (1./surface_area);

            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;

    // background flow with unit pressure drop 
    const auto bg_flow = [](const sctl::Vector<Real>& X) {
        const sctl::Long N = X.Dim()/3;
        sctl::Vector<Real> U(N*3);
        for (sctl::Long i = 0; i < N; i++) {
            const auto x = X.begin() + i*3;
            U[i*3+0] = - ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5)) / 4;
            U[i*3+1] = 0;
            U[i*3+2] = 0;
        }
        return U;
    };

    sctl::Vector<Real> sigma;
    solver(&sigma, BIO, bg_flow(X0) * (pressure_drop/period_length), gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);

    // Visualization
    { 
        // Create a new conv-div channel with no particle inside to sample target points from.
        PeriodicGeom<Real> trg;
        const sctl::Long Nelem_trg = 16;
        const sctl::Long FourierOrder_trg = 16;
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg; // place holder arrays for particle locations, but no initial length tells constructor not to include any particles.
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build_trg = trg.build_conv_div(Nelem_trg, ElemOrder, FourierOrder_trg, channel_r1, channel_r2, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, ptcl_ord);
        elem_lst_trg = std::get<0>(build_trg);

        // Generate interior sample points inside the channel, then filter out targets that landed within a particle (of actual geometry).
        VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        sctl::Vector<Real> X0_all = vol_vis.GetCoord(); // set new target coordinates
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target_rotated(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode, ptcls_thetas, ptcls_phis);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        // Evaluate solution flow at targets
        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U -= bg_flow(X0) * (pressure_drop / period_length);

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
        vol_vis.WriteVTK("vis/ConvDiv_U", U_vis); 
    }
}

/**
    Example of a 2 periodic array of several loops sandwiched between two infinite planes.
    // TODO: this function is not parallelized 
*/

template <class Real> void planes_with_particle(sctl::Comm comm) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;
    // Set quadrature parameters
    const sctl::Long Nelem = 4;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 64;

    const sctl::Long geom_mode = 3; 

    const Real tol = 1e-10; // quadrature tolerance
    const sctl::Long gl_order = 49; // *Cheb order on each panel of plane
    const sctl::Long Nelem_x = 2; // Number of panels in x .. 
    const sctl::Long Nelem_y = 2; // .. and y directions
    const Real z_offset = 0.01; // planes located at z = <z_offset> and z = 1 - <z_offset>
    // Set geometry parameters
    const Real pressure_drop = -1.0; // background pressure driven flow.
    const Real period_length = 1; // length of periodic box.
    // Set GMRES parameters
    const Real gmres_tol = 1e-8;
    const sctl::Long gmres_max_iter = 150;

    // Make SlenderElemList object for one particle in the center of the unit box.
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient, ptcls_thetas, ptcls_phis;
    // std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    // std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.many_spheroids3(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs);
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.many_loops3(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    ptcls_thetas = std::get<2>(build0);
    ptcls_phis = std::get<3>(build0);
    sctl::Vector<Real> X0_ptcl;
    elem_lst0.GetNodeCoord(&X0_ptcl, nullptr, nullptr);
    elem_lst0.WriteVTK("vis/planes_ptcl_geometry", X0_ptcl);

    // Plane (note: not parallelized)
    sctl::PlaneIntegral<Real> plane(gl_order, Nelem_x, Nelem_y, z_offset);
    sctl::Vector<Real> X0_wall;
    plane.GetNodeCoord(&X0_wall, nullptr, nullptr);
    // plane.WriteVTK("vis/planes_geometry", X0_wall, comm);

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); 
    LayerPotenOp0.AddElemList(elem_lst0,"1");
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XY, period_length);
    LayerPotenOp0.AddElemList(plane,"2"); 

    // Collect all surface points
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
        surface_area = surface_area_[0];
    }
    Real surface_area_wall;
    sctl::Vector<Real> wts_wall;
    {
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        plane.GetFarFieldNodes(X, Xn, wts_wall, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts_wall*0+1, wts_wall);
        surface_area_wall = surface_area_[0];
    }
    

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&wts_wall,&surface_area_wall,&plane,&LayerPotenOp0,&DL_scal,&X0_ptcl,&X0_wall,&NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        
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
            sigma_mean = sigma_mean_ptcl + sigma_mean_wall; 
            sigma_mean *= (1/(surface_area + surface_area_wall));
            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma0.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    const auto bg_flow = [](const sctl::Vector<Real>& X) {
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

    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    solver(&sigma,BIO, bg_flow(X0) * (pressure_drop/period_length), gmres_tol, gmres_max_iter);

    // Evaluation
    {
        // Create a uniform grid in the unit cube below the planes, then filter out points inside the particle
        PeriodicGeom<Real> trg;    
        CubeVolumeVisShifted<Real> vol_vis(10, 0.9, comm);
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        // std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target_rotated(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode, ptcls_thetas, ptcls_phis);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U(X0.Dim());
        BIO(&U, sigma);
        U -= bg_flow(X0) * (pressure_drop/period_length);
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
        vol_vis.WriteVTK("vis/planes_U", U_vis);
    }

    // {
    //     sctl::Vector<Real> X0_all(3);
    //     X0_all[0] = 1;
    //     X0_all[0] = 1;
    //     X0_all[0] = 1;
    //     PeriodicGeom<Real> trg;    
    //     std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target_rotated(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode, ptcls_thetas, ptcls_phis);
    // }
    
}

/**
    Example of a 3 periodic array of 25 spheres under a background pressure drop.
*/
template <class Real> void particle_3peri(sctl::Comm comm) {

   // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;
    // Set quadrature parameters
    const sctl::Long Nelem = 4;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 64;
    const sctl::Long geom_mode = 0; // designate particles to be spheres.
    const sctl::Long Nptcl = 25;
    const Real tol = 1e-10; // quadrature tolerance
    // Set geometry parameters
    const Real pressure_drop = -1.0; // background pressure driven flow.
    const Real period_length = 1; // length of periodic box.
    // Set GMRES parameters
    const Real gmres_tol = 1e-8;
    const sctl::Long gmres_max_iter = 150;
    
    // Set up SlenderElemList object for spheres
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    elem_lst0.WriteVTK("vis/25spheres_geometry", X0); 

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
    }

    // // DEBUG MPI Surface area:
    // if (!comm.Rank()) {
    //     Real surfA_manual = 0.;
    //     for (sctl::Long i=0; i<Nptcl; i++) {
    //         Real r_i = ptcls_rs[i];
    //         surfA_manual += 4.*sctl::const_pi<Real>() * r_i * r_i;
    //     }
    //     std::cout << "Surface area computed for a total of " << Nptcl << " spheres is " << surface_area << "; manual calculation gives " << surfA_manual << std::endl;
    // }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); 
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length);

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

            // DEBUG: check that sigma-sigma_mean has surface integral = 0:
            sctl::Vector<Real> sigma1 = sigma_;
            AddConstVec(sigma1, -sigma_mean);
            sctl::Vector<Real> sigma_test_;
            SurfaceIntegral(sigma_test_, sigma1, wts);
            std::cout << "Surface integral of sigma - sigma bar = " << sigma_test_[0] << ", "<< sigma_test_[1] << ", " << sigma_test_[2] << ". "<< std::endl;
        
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
        sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
        AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

        sctl::Vector<Real> U0;
        LayerPotenOp0.ComputeSL(U0, force_density);
        return U0;
    };

    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    LayerPotenOp0.SetTargetCoord(X0);
    solver(&sigma, BIO, eval_rhs(pressure_drop), gmres_tol, gmres_max_iter);

    {
        PeriodicGeom<Real> trg;    
        CubeVolumeVisShifted<Real> vol_vis(50, 0.95, comm);
        // Filter out target points inside spheres 
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);
        U -= eval_rhs(pressure_drop);

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
        vol_vis.WriteVTK("vis/25spheres_U", U_vis); 
    }

}


int main(int argc, char** argv) {

    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Comm comm = sctl::Comm::World();
        channel_with_particle<Real>(comm);
        // planes_with_particle<Real>(comm);
        // particle_3peri<Real>(comm);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
