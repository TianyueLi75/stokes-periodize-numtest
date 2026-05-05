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

/**
    Example of a converging-divering channel with some spheroids on the interior. Visualizations are stored in vis/examples/ folder.
*/
template <class Real> void channel_with_particle(sctl::Comm comm) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;
    // Set quadrature parameters
    const sctl::Long Nelem_channel = 176;
    // const sctl::Long Nelem_channel = 264; // only for regularized spheroids -- 77 contained in geom, not 62 anymore. 
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 48;
    const Real tol = 1e-10; // quadrature tolerance
    // Set geometry parameters
    const Real period_length = 1; // length of periodic box.
    const Real channel_r1 = 0.025; 
    const Real channel_r2 = 0.1;
    // Set GMRES parameters
    const Real gmres_tol = 1e-8;
    const sctl::Long gmres_max_iter = 200;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls(1); 
    sctl::Long ptcl_ord = 8; // number of panels on each particle (same fourier order as channel)
    ptcls = ptcl_ord;
    // Arrays to hold particle...
    sctl::Vector<Real> ptcls_Xcs, ptcls_rs, ptcls_u0s, ptcls_thetas, ptcls_phis, NormalOrient; // -1 if outward normal on surface points into fluid, +1 otherwise. (since BIO defined with a -1/2*sigma already)
    sctl::Vector<sctl::Long> ptcls_ifprolate;
    sctl::SlenderElemList<Real> elem_lst0; 
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.build_conv_div(Nelem_channel, ElemOrder, FourierOrder, channel_r1, channel_r2, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcls_u0s, ptcls_ifprolate, ptcl_ord, 3);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    ptcls_thetas = std::get<2>(build0); 
    ptcls_phis = std::get<3>(build0);
    std::cout << "done with forming conv div channel, number of particles is " << ptcls.Dim() << std::endl;

    sctl::Vector<Real> X0, Xnsurf; // target coordinates
    elem_lst0.GetNodeCoord(&X0, &Xnsurf, nullptr);
    elem_lst0.WriteVTK("vis/examples/ConvDiv_geometry", Xnsurf, comm);

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
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); 
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length); 

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,&NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> sigma_mean, sigma0;
        { 
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
    sctl::Vector<Real> sigma;

    // Slip velocity on particles.
    sctl::Long Nnodes_channel = Nelem_channel * ElemOrder * FourierOrder;
    sctl::Long Nnodes_per_ptcl = ptcl_ord * ElemOrder * FourierOrder;
    // MPI
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
    // HARD ASSUMPTIONS FOR NOW: all channel panels on MPI rank 0; no particles split between two processes
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
    elem_lst0.WriteVTK("vis/examples/ConvDiv_slipBC", rhs, comm);
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

    solver(&sigma, BIO, rhs, gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);

    { 
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
        vol_vis.WriteVTK("vis/examples/ConvDiv_slipU", U_vis); 
    }
}

/**
    Example of a 2 periodic array of several loops sandwiched between two infinite planes.
*/
template <class Real> void planes_with_loops(sctl::Comm comm) {

    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;
    // Set quadrature parameters
    const sctl::Long Nelem = 40;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 96;

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

    // System of loops
    const sctl::Long Nptcl = 25;
    std::string data_filename = "data/sphere_data_"+std::to_string(Nptcl)+"_larger.txt";
    sctl::Matrix<Real> Xc_from_file(Nptcl,4);
    std::ifstream infile(data_filename);
    if (!infile) {
        std::cerr << "Error opening file " << data_filename << std::endl;
        SCTL_ASSERT(false);
    }
    for (sctl::Long row=0; row < Nptcl; row++) {
        for (sctl::Long col=0; col < 4; col++) {
        if (!(infile >> Xc_from_file(row,col))) {
            std::cerr << "not enough entries in data file" << std::endl;
        }
        }
    }
    sctl::Vector<Real> ptcls_Xcs, ptcls_rs, ptcls_major_rs, ptcls_minor_rs, ptcls_thetas, ptcls_phis;
    sctl::Vector<sctl::Long> ptcls;
    ptcls_Xcs.ReInit(0);
    ptcls_rs.ReInit(0);
    ptcls.ReInit(0);
    ptcls_major_rs.ReInit(0);
    ptcls_minor_rs.ReInit(0);
    ptcls_thetas.ReInit(0);
    ptcls_phis.ReInit(0);
    srand48(2);
    for (sctl::Long i=0; i<Nptcl; i++) {
        ptcls_Xcs.PushBack(Xc_from_file(i,0));
        ptcls_Xcs.PushBack(Xc_from_file(i,1));
        ptcls_Xcs.PushBack(Xc_from_file(i,2));
        ptcls_rs.PushBack(Xc_from_file(i,3));
        ptcls.PushBack(Nelem);

        const Real minor_r = 0.03 + (drand48()-0.5) * 0.005; // 0.05 +- 0.0025
        ptcls_minor_rs.PushBack(minor_r);
        ptcls_major_rs.PushBack(ptcls_rs[i] - minor_r);

        const Real theta_rotate = drand48() * sctl::const_pi<Real>() * 2.;
        const Real phi_rotate = drand48() * sctl::const_pi<Real>();
        ptcls_thetas.PushBack(theta_rotate);
        ptcls_phis.PushBack(phi_rotate);
    }
    PeriodicGeom<Real> obj;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.loops_system(ptcls, ElemOrder, FourierOrder, ptcls_Xcs, ptcls_major_rs, ptcls_minor_rs, ptcls_thetas, ptcls_phis, comm);
    sctl::SlenderElemList<Real> elem_lst0 = std::get<0>(build0);
    sctl::Vector<Real> NormalOrient = std::get<1>(build0);

    sctl::Vector<Real> X0_ptcl;
    elem_lst0.GetNodeCoord(&X0_ptcl, nullptr, nullptr);
    elem_lst0.WriteVTK("vis/examples/PlanesPtcl_ptcls", X0_ptcl, comm);

    // Create plane object on all process, but only use it on Rnak 0
    sctl::PlaneIntegral<Real> plane(gl_order, Nelem_x, Nelem_y, z_offset);

    sctl::Vector<Real> X0_wall;
    if (!comm.Rank()) {
        plane.GetNodeCoord(&X0_wall, nullptr, nullptr);
        plane.WriteVTK("vis/examples/PlanesPtcl_wall", X0_wall, comm);
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
    NormalOrient_ = -1.; // Normal orient = -1 (-sign below) means all normals point into fluid (exterior problem)
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
        vol_vis.WriteVTK("vis/examples/PlanesPtcl_U", U_vis);
    }
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
    const sctl::Long geom_mode = 0; 
    const sctl::Long Nptcl = 25;
    const Real tol = 1e-14; // quadrature tolerance
    // Set geometry parameters
    const Real pressure_drop = -1.0; // background pressure driven flow.
    const Real period_length = 1; // length of periodic box.
    // Set GMRES parameters
    const Real gmres_tol = 1e-12;
    const sctl::Long gmres_max_iter = 150;
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs, ptcls_rs, NormalOrient;
    sctl::SlenderElemList<Real> elem_lst0;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    elem_lst0.WriteVTK("vis/examples/Ptcl3peri_geometry", X0, comm); 

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

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); 
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length);

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,&NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> sigma_mean, sigma0;
        { 
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
        CubeVolumeVisShifted<Real> vol_vis(60, 0.95, comm);
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
        vol_vis.WriteVTK("vis/examples/Ptcl3peri_U", U_vis); 
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
        } else if (mode == 2) {
            particle_3peri<Real>(comm);
        } 
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
