#include "periodize.hpp"
#include "utils_geom.hpp"
#include "bio_operator.hpp"

template <class Real> sctl::Vector<Real> GetVslip_spheres(const sctl::Vector<Real>& ptcls_Xnsurf, const sctl::Vector<Real>& ptcls_Xcs, const sctl::Vector<Real>& ptcls_thetas, const sctl::Vector<Real>& ptcls_phis, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder) {
    const sctl::Long Nnodes_per_ptcl = Nelem * ElemOrder * FourierOrder;
    const sctl::Long Nptcls = ptcls_Xcs.Dim()/3;
    // std::cout << "DEBUG in vslip: Nnodes per ptcl is " << Nnodes_per_ptcl << ", Nptcls is " << Nptcls << std::endl;
    
    sctl::Vector<Real> Vslip(ptcls_Xnsurf.Dim());
    srand48(1);
    Vslip.SetZero();
    for (sctl::Long i=0; i<Nptcls; i++) {
        sctl::Vector<Real> Vslip_here(3*Nnodes_per_ptcl, (sctl::Iterator<Real>) Vslip.begin() + 3*Nnodes_per_ptcl*i, false);
        sctl::Vector<Real> center_here(3, (sctl::Iterator<Real>) ptcls_Xcs.begin() + 3*i, false);
        // const Real scalar = drand48()*0.8 + 0.1; // randomly scaled slip velocity by (0.1,0.9).
        const Real scalar = 1.; // no change in slip magnitude on different spheres/spheroids
        Real theta_rotate, phi_rotate;
        if (ptcls_thetas.Dim()>0) {
            theta_rotate = ptcls_thetas[i];
            phi_rotate = ptcls_phis[i];
        } else {
            theta_rotate = 0.;
            phi_rotate = 0.;
        }
        const Real cos_theta_rotate = sctl::cos<Real>(theta_rotate);
        const Real sin_theta_rotate = sctl::sin<Real>(theta_rotate);
        const Real cos_phi_rotate = sctl::cos<Real>(phi_rotate);
        const Real sin_phi_rotate = sctl::sin<Real>(phi_rotate);

        for (sctl::Long panel=0; panel < Nelem; panel++) {
            const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrder);

            for (sctl::Long el=0; el<ElemOrder; el++) {
                const Real theta = sctl::const_pi<Real>() * (panel+nodes[el])/Nelem;
                const Real coeff = scalar * sctl::sin<Real>(theta); // sinusoidal slip magnitude, 0 at north and south poles
                // const Real coeff = 1.;

                for (sctl::Long fl=0; fl<FourierOrder; fl++) {
                    Real phi = 2. * sctl::const_pi<Real>() * fl / FourierOrder;
                    const sctl::Long idx = panel*ElemOrder*FourierOrder + el*FourierOrder + fl;
                    const sctl::Vector<Real> Xn_here(3, (sctl::Iterator<Real>) ptcls_Xnsurf.begin() + 3*Nnodes_per_ptcl*i + idx*3, false);

                    Real t1_unrotated = -sctl::sin<Real>(theta);
                    Real t2_unrotated = sctl::cos<Real>(theta) * sctl::cos<Real>(phi);
                    Real t3_unrotated = sctl::cos<Real>(theta) * sctl::sin<Real>(phi);
                    
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
                    
                    // std::cout << "Xn here is (" << Xn_here[0] << "," << Xn_here[1] << ", " << Xn_here[2] << "), vslip unrotated = " << Vslip_here_unrotated[0] << ", " << Vslip_here_unrotated[1] << "," << Vslip_here_unrotated[2] << std::endl;
                    if (sctl::fabs(vdotn)>1e-6) {
                        std::cout << "tangent dot n is nonzero: " << vdotn << std::endl;
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

template <class Real> sctl::Vector<Real> GetVslip_spheroids(const sctl::Vector<Real>& ptcls_Xnsurf, const sctl::Vector<Real>& ptcls_Xcs, const sctl::Vector<Real>& ptcls_sizes, const sctl::Vector<Real>& ptcls_u0s, const sctl::Vector<Real>& ptcls_thetas, const sctl::Vector<Real>& ptcls_phis, const sctl::Vector<sctl::Long>& ptcls_ifprolate, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder) {
    const sctl::Long Nnodes_per_ptcl = Nelem * ElemOrder * FourierOrder;
    const sctl::Long Nptcls = ptcls_sizes.Dim();
    std::cout << "DEBUG in vslip: Nnodes per ptcl is " << Nnodes_per_ptcl << ", Nptcls is " << Nptcls << std::endl;
    
    sctl::Vector<Real> Vslip(ptcls_Xnsurf.Dim());
    srand48(1);
    Vslip.SetZero();
    for (sctl::Long i=0; i<Nptcls; i++) {
        const Real a_here = ptcls_sizes[i];
        const Real u0_here = ptcls_u0s[i];
        const int if_prolate_here = ptcls_ifprolate[i];
        sctl::Vector<Real> Vslip_here(3*Nnodes_per_ptcl, (sctl::Iterator<Real>) Vslip.begin() + 3*Nnodes_per_ptcl*i, false);
        sctl::Vector<Real> center_here(3, (sctl::Iterator<Real>) ptcls_Xcs.begin() + 3*i, false);
        // const Real scalar = drand48()*0.8 + 0.1; // randomly scaled slip velocity by (0.1,0.9).
        const Real scalar = 1.; // no change in slip magnitude on different spheres/spheroids
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
                // const Real coeff = 1.;

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
                    
                    // std::cout << "Xn here is (" << Xn_here[0] << "," << Xn_here[1] << ", " << Xn_here[2] << "), vslip unrotated = " << Vslip_here_unrotated[0] << ", " << Vslip_here_unrotated[1] << "," << Vslip_here_unrotated[2] << std::endl;
                    if (sctl::fabs(vdotn)>1e-6) {
                        std::cout << "tangent dot n is nonzero: " << vdotn << std::endl;
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
    Given arrays of properties for spheroids, return the SlenderELem List representing this setup.
*/
template <class Real> sctl::SlenderElemList<Real> spheroid_system(                                                           
                                                                const sctl::Long Nelem_ptcl, 
                                                                const sctl::Long ElemOrder, 
                                                                const sctl::Long FourierOrder, 
                                                                const sctl::Vector<Real> Xcenter_lst, 
                                                                const sctl::Vector<sctl::Long> if_prolate_lst, // true if prolate
                                                                const sctl::Vector<Real> u0_lst, 
                                                                const sctl::Vector<Real> r_lst, 
                                                                const sctl::Vector<Real> theta_lst, 
                                                                const sctl::Vector<Real> phi_lst, 
                                                                const sctl::Comm comm) 
{
    const sctl::Long Nptcls = u0_lst.Dim();
    sctl::Vector<sctl::Long> ElemOrderVec(Nptcls * Nelem_ptcl);
    sctl::Vector<sctl::Long> FourierOrderVec(Nptcls * Nelem_ptcl);
    ElemOrderVec = ElemOrder;
    FourierOrderVec = FourierOrder;
    
    const auto prolate_geom = [](Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& cylindrical_radius, const Real polar_angle, const Real size, const Real u0) {
        SCTL_ASSERT(u0 > 1.); 
        x = size * u0 * sctl::cos<Real>(polar_angle);
        y = 0.;
        z = 0.;
        cylindrical_radius = size * sctl::sqrt<Real>(u0*u0 - 1) * sctl::sin<Real>(polar_angle);
        ex = 0.;
        ey = 1.;
        ez = 0.;
    };

    const auto oblate_geom = [](Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& cylindrical_radius, const Real polar_angle, const Real size, const Real u0) {
        x = size * u0 * sctl::cos<Real>(polar_angle); 
        y = 0.;
        z = 0.;
        cylindrical_radius = size * sctl::sqrt<Real>(u0*u0 + 1) * sctl::sin<Real>(polar_angle);
        ex = 0.;
        ey = 1.;
        ez = 0.;
    };

    sctl::Vector<Real> Xc, eps, orient;
    for (sctl::Long pid = 0; pid < Nptcls; pid ++) {
        const Real ptcl_size = r_lst[pid];
        const Real ptcl_u0 = u0_lst[pid];
        const Real theta_rotate = theta_lst[pid];
        const Real phi_rotate = phi_lst[pid];
        const sctl::Long if_prolate_here = if_prolate_lst[pid];
        const Real cos_theta_rotate = sctl::cos<Real>(theta_rotate);
        const Real sin_theta_rotate = sctl::sin<Real>(theta_rotate);
        const Real cos_phi_rotate = sctl::cos<Real>(phi_rotate);
        const Real sin_phi_rotate = sctl::sin<Real>(phi_rotate);
        for (sctl::Long i=0; i < Nelem_ptcl; i++) {
            const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
            for (sctl::Long j=0; j<ElemOrderVec[i]; j++) {
                const Real theta = sctl::const_pi<Real>() * (i+nodes[j])/Nelem_ptcl;
                Real x, y, z, ex, ey, ez, eps_j;
                if (if_prolate_here) {
                    prolate_geom(x,y,z,ex,ey,ez,eps_j,theta,ptcl_size,ptcl_u0);
                } else {
                    oblate_geom(x,y,z,ex,ey,ez,eps_j,theta,ptcl_size,ptcl_u0);
                }
                Real x_rotated = cos_theta_rotate * cos_phi_rotate * x - sin_phi_rotate * y + sin_theta_rotate * cos_phi_rotate * z;
                Real y_rotated = cos_theta_rotate * sin_phi_rotate * x + cos_phi_rotate * y + sin_theta_rotate * sin_phi_rotate * z;
                Real z_rotated = -sin_theta_rotate * x + cos_theta_rotate * z;
                Real ex_rotated = cos_theta_rotate * cos_phi_rotate * ex - sin_phi_rotate * ey + sin_theta_rotate * cos_phi_rotate * ez;
                Real ey_rotated = cos_theta_rotate * sin_phi_rotate * ex + cos_phi_rotate * ey + sin_theta_rotate * sin_phi_rotate * ez;
                Real ez_rotated = -sin_theta_rotate * ex + cos_theta_rotate * ez;

                Xc.PushBack(Xcenter_lst[pid*3+0]+x_rotated);
                Xc.PushBack(Xcenter_lst[pid*3+1]+y_rotated);
                Xc.PushBack(Xcenter_lst[pid*3+2]+z_rotated);
                eps.PushBack(eps_j);
                orient.PushBack(ex_rotated);
                orient.PushBack(ey_rotated);
                orient.PushBack(ez_rotated);
            }
        }
    }

    const auto init_elem_lst = [](sctl::SlenderElemList<Real>& elem_lst, const sctl::Vector<sctl::Long>& ElemOrderVec, const sctl::Vector<sctl::Long>& FourierOrderVec, const sctl::Vector<Real>& X, const sctl::Vector<Real>& R, const sctl::Vector<Real>& OrientVec, const sctl::Comm comm) {
        const sctl::Long Nelem_tot = ElemOrderVec.Dim();
        sctl::Long loc_elem_cnt, loc_elem_dsp;
        if (Nelem_tot) { // Set loc_elem_cnt, loc_elem_dsp 
            // node_dsp.ReInit(Nelem_tot);
            // sctl::Vector<sctl::Long> node_cnt(Nelem_tot);
            sctl::Vector<sctl::Long> node_cnt(Nelem_tot), node_dsp(Nelem_tot); node_dsp = 0;
            for (sctl::Long i = 0; i < Nelem_tot; i++) {
            node_cnt[i] = ElemOrderVec[i] * FourierOrderVec[i] * FourierOrderVec[i];
            }
            sctl::omp_par::scan(node_cnt.begin(), node_dsp.begin(), Nelem_tot);
            const sctl::Long Nnodes = node_cnt[Nelem_tot-1] + node_dsp[Nelem_tot-1];

            const sctl::Long Np = comm.Size();
            const sctl::Long rank = comm.Rank();
            sctl::Long a = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+0)/Np) - node_dsp.begin();
            sctl::Long b = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+1)/Np) - node_dsp.begin();
            if (rank == Np - 1) b = Nelem_tot;
            if (rank == 0) a = 0;
            loc_elem_cnt = b-a;
            loc_elem_dsp = a;
        } else {
            loc_elem_cnt = 0;
            loc_elem_dsp = 0;
        }

        const sctl::Vector<sctl::Long> LocElemOrder(loc_elem_cnt, (sctl::Iterator<sctl::Long>)ElemOrderVec.begin() + loc_elem_dsp, false);
        const sctl::Vector<sctl::Long> LocFourierOrder(loc_elem_cnt, (sctl::Iterator<sctl::Long>)FourierOrderVec.begin() + loc_elem_dsp, false);

        sctl::Long dsp = 0, cnt = 0;
        for (sctl::Long i = 0; i < loc_elem_dsp; i++) dsp += ElemOrderVec[i];
        for (sctl::Long i = 0; i < loc_elem_cnt; i++) cnt += ElemOrderVec[loc_elem_dsp+i];
        const sctl::Vector<Real> X_(cnt*3, (sctl::Iterator<Real>)X.begin() + dsp*3, false);
        const sctl::Vector<Real> R_(cnt, (sctl::Iterator<Real>)R.begin() + dsp, false);
        const sctl::Vector<Real> OrientVec_(cnt*3, (sctl::Iterator<Real>)OrientVec.begin() + dsp*3, false);

        elem_lst.template Init<Real>(LocElemOrder, LocFourierOrder, X_, R_, OrientVec_);  
    };

    sctl::SlenderElemList<Real> elem_lst;
    init_elem_lst(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, orient, comm);
  
    return elem_lst;

}



template <class Real> void test_function(bool if_DL, bool if_noslip, sctl::Comm comm) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    Real DL_scal = 0;
    if (if_DL) {
        DL_scal = 1.;
    }
    // Set quadrature parameters
    const sctl::Long Nelem_channel = 100;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 72;
    const Real tol = 1e-12; // quadrature tolerance
    // Set geometry parameters
    const Real pressure_drop = -1.0; // background pressure driven flow.
    const Real period_length = 1; // length of periodic box.
    const Real channel_r1 = 0.025; 
    const Real channel_r2 = 0.1;
    // Set GMRES parameters
    const Real gmres_tol = 1e-8;
    const sctl::Long gmres_max_iter = 1000;

    // Create channel + particles object: "elem_lst0".
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls(1); // length of ptcls will get updated through the conv-div channel constructor. (this is not the case for all channel constructors.)
    sctl::Long ptcl_ord = 4; 
    ptcls = ptcl_ord;
    // Arrays to hold particle...
    sctl::Vector<Real> ptcls_Xcs; // ...center positions
    sctl::Vector<Real> ptcls_rs; // ...radii
    sctl::Vector<Real> ptcls_u0s; // .. ~ aspect ratio
    sctl::Vector<Real> ptcls_thetas, ptcls_phis;
    sctl::Vector<sctl::Long> ptcls_ifprolate;
    sctl::SlenderElemList<Real> elem_lst0; 
    sctl::Vector<Real> NormalOrient; // -1 if outward normal on surface points into fluid, +1 otherwise. (since BIO defined with a -1/2*sigma already)
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.build_conv_div(Nelem_channel, ElemOrder, FourierOrder, channel_r1, channel_r2, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcls_u0s, ptcls_ifprolate, ptcl_ord);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    ptcls_thetas = std::get<2>(build0); 
    ptcls_phis = std::get<3>(build0);
    std::cout << "done with forming conv div channel, number of particles is " << ptcls.Dim() << std::endl;

    sctl::Vector<Real> X0, Xnsurf; // target coordinates
    elem_lst0.GetNodeCoord(&X0, &Xnsurf, nullptr);
    // elem_lst0.WriteVTK("vis/ConvDiv_geometry_4", Xnsurf, comm);
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
        // std::cout << "DEBUG, in BIO" << std::endl;
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
            // std::cout << "DEBUG: Sigma mean= " << sigma_mean[0] << ", " << sigma_mean[1] << ", " << sigma_mean[2] << std::endl;

            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        // std::cout << "DEBUG: After compute potentials, before DL add sigma" << std::endl;
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> sigma;

    sctl::Vector<Real> rhs(X0.Dim());
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
    if (if_noslip) {
        rhs = bg_flow(X0)* (pressure_drop/period_length);
    } else {
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
        rhs.SetZero();
        sctl::Vector<Real> vslip_ptcl;
        if (!comm.Rank()) {
            const sctl::Vector<Real> Xnptcl(3 * Nnodes_per_ptcl * Nptcl_this_mpi, (sctl::Iterator<Real>) Xnsurf.begin() + 3*Nnodes_channel, false);
            // std::cout << "starting point of Xnptcl: " << 3*Nnodes_channel << std::endl;
            // Input only the relevant Xc, theta, phi
            sctl::Vector<Real> ptcls_Xcs_here(3*Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_Xcs.begin() + 3*Nptcl_dsp_this_mpi, true);
            sctl::Vector<Real> ptcls_sizes_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_rs.begin() + Nptcl_dsp_this_mpi, true);
            sctl::Vector<Real> ptcls_u0s_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_u0s.begin() + Nptcl_dsp_this_mpi, true);
            sctl::Vector<Real> ptcls_thetas_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_thetas.begin() + Nptcl_dsp_this_mpi, true);
            sctl::Vector<Real> ptcls_phis_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_phis.begin() + Nptcl_dsp_this_mpi, true);
            sctl::Vector<sctl::Long> ptcls_ifprolate_here(Nptcl_this_mpi, (sctl::Iterator<sctl::Long>) ptcls_ifprolate.begin() + Nptcl_dsp_this_mpi, true);
            // vslip_ptcl = GetVslip_spheres(Xnptcl, ptcls_Xcs_here, ptcls_thetas_here, ptcls_phis_here, ptcl_ord, ElemOrder, FourierOrder);
            vslip_ptcl = GetVslip_spheroids(Xnptcl, ptcls_Xcs_here, ptcls_sizes_here, ptcls_u0s_here, ptcls_thetas_here, ptcls_phis_here, ptcls_ifprolate_here, ptcl_ord, ElemOrder, FourierOrder);
            // Populate only the non-channel part of rhs
            for (int j=3*Nnodes_channel; j<rhs.Dim(); j++) {
                rhs[j] = vslip_ptcl[j-3*Nnodes_channel];
            }
        } else {
            const sctl::Vector<Real> Xptcl = X0;
            sctl::Long Nptcl_this_rank = static_cast<int>(X0.Dim() / Nnodes_per_ptcl / 3);
            // Input only the relevant Xc, theta, phi
            sctl::Vector<Real> ptcls_Xcs_here(3*Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_Xcs.begin() + 3*Nptcl_dsp_this_mpi, true);
            sctl::Vector<Real> ptcls_sizes_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_rs.begin() + Nptcl_dsp_this_mpi, true);
            sctl::Vector<Real> ptcls_u0s_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_u0s.begin() + Nptcl_dsp_this_mpi, true);
            sctl::Vector<Real> ptcls_thetas_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_thetas.begin() + Nptcl_dsp_this_mpi, true);
            sctl::Vector<Real> ptcls_phis_here(Nptcl_this_mpi, (sctl::Iterator<Real>) ptcls_phis.begin() + Nptcl_dsp_this_mpi, true);
            sctl::Vector<sctl::Long> ptcls_ifprolate_here(Nptcl_this_mpi, (sctl::Iterator<sctl::Long>) ptcls_ifprolate.begin() + Nptcl_dsp_this_mpi, true);
            // vslip_ptcl = GetVslip_spheres(Xnsurf, ptcls_Xcs_here, ptcls_thetas_here, ptcls_phis_here, ptcl_ord, ElemOrder, FourierOrder);
            vslip_ptcl = GetVslip_spheroids(Xnsurf, ptcls_Xcs_here, ptcls_sizes_here, ptcls_u0s_here, ptcls_thetas_here, ptcls_phis_here, ptcls_ifprolate_here, ptcl_ord, ElemOrder, FourierOrder);
            rhs = vslip_ptcl;
        }
        // std::cout << "X0 size on rank " << comm.Rank() << " is " << X0.Dim() << ", rhs size is " << rhs.Dim() << std::endl;
        // DEBUG
        // std::cout << "done with vslip, writing to VTK" << std::endl;
        // elem_lst0.WriteVTK("vis/channel_vslip_4", rhs, comm);
        // Check that u dot n is always zero
        for (int ii=0; ii<X0.Dim()/3; ii++) {
            const sctl::Vector<Real> vslip_here(3, (sctl::Iterator<Real>) rhs.begin() + ii*3, false);
            const sctl::Vector<Real> xn_here(3, (sctl::Iterator<Real>) Xnsurf.begin() + ii*3, false);
            const Real udotn = vslip_here[0]*xn_here[0] + vslip_here[1]*xn_here[1] + vslip_here[2]*xn_here[2];
            if (sctl::fabs(udotn) > 1e-6) {
                // std::cout << "u dot n at node " << ii << " is nonzero: " << udotn << ". For context, Nelem_channel = " << Nelem_channel << ", Nelem_ptcl = " << ptcl_ord << ", FourierORder = " << FourierOrder << std::endl;
                std::cout << "u dot n at node " << ii << " is nonzero: " << udotn << ". Normal is " << xn_here[0] << ", " << xn_here[1] << ", " << xn_here[2] << ", vslip is " << vslip_here[0] << ", " << vslip_here[1] << ", " << vslip_here[2] << std::endl;
                return;
            }
        }
    }

    solver(&sigma, BIO, rhs, gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);
    // elem_lst0.WriteVTK("vis/ConvDiv_density_4", sigma, comm);

    // Evaluation
    {
        // Create a new conv-div channel with no particle inside to sample target points from.
        PeriodicGeom<Real> trg;
        const sctl::Long Nelem_trg = 80;
        const sctl::Long FourierOrder_trg = 64;
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg; // place holder arrays for particle locations, but no initial length tells constructor not to include any particles.
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        sctl::Vector<Real> ptcls_u0s_trg;
        sctl::Vector<sctl::Long> ptcls_ifprolate_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build_trg = trg.build_conv_div(Nelem_trg, ElemOrder, FourierOrder_trg, channel_r1, channel_r2, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, ptcls_u0s_trg, ptcls_ifprolate_trg, ptcl_ord);
        elem_lst_trg = std::get<0>(build_trg);

        // Generate interior sample points inside the channel, then filter out targets that landed within a particle (of actual geometry).
        VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        sctl::Vector<Real> X0_all = vol_vis.GetCoord(); // set new target coordinates
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        // std::cout << "size of target points before filter: " << X0_all.Dim()<< std::endl;
        // std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0); // ChocoFeb2026: filter spheres, vslip debug
        // std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_spheroids(X0_all, ptcls_rs, ptcls_u0s, ptcls_Xcs, ptcls_ifprolate);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_spheroids_rotated(X0_all, ptcls_rs, ptcls_u0s, ptcls_Xcs, ptcls_ifprolate, ptcls_thetas, ptcls_phis);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);
        // std::cout << "size of target points after filter: " << X0.Dim()<< std::endl;

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U(X0.Dim());
        BIO(&U, sigma);
        if (if_noslip) {
            U -= bg_flow(X0) * (pressure_drop/period_length);
        }
        sctl::Vector<Real> U_vis(X0_all.Dim());
        // U_vis = 0.;
        U_vis = std::nan(" ");
        sctl::Long X1_ptr = 0;
        for (sctl::Long i=0; i<X0_all.Dim()/3; i++) {
            if (filtered_inds[i] == 0) {
                U_vis[i*3] = U[X1_ptr*3];
                U_vis[i*3+1] = U[X1_ptr*3+1];
                U_vis[i*3+2] = U[X1_ptr*3+2];
                X1_ptr += 1;
            }
        }
        vol_vis.WriteVTK("vis/ConvDiv_U_4_nan", U_vis);
    }
}

template <class Real> void test_nochannel_spheres(bool if_DL, bool if_noslip, sctl::Comm comm) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    Real DL_scal = 0;
    if (if_DL) {
        DL_scal = 1.;
    }
    // Set quadrature parameters
    const sctl::Long Nelem = 4;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 64;
    const Real tol = 1e-10; // quadrature tolerance
    // Set geometry parameters
    const Real pressure_drop = -1.0; // background pressure driven flow.
    const Real period_length = 1; // length of periodic box.
    // Spheres
    const sctl::Long geom_mode = 0;
    const sctl::Long Nptcl = 25; 
    // Set GMRES parameters
    const Real gmres_tol = 1e-6;
    const sctl::Long gmres_max_iter = 200;

    // Set up SlenderElemList object for spheres
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient, ptcls_thetas, ptcls_phis;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);

    sctl::Vector<Real> X0, Xnsurf; // target coordinates
    elem_lst0.GetNodeCoord(&X0, &Xnsurf, nullptr);
    elem_lst0.WriteVTK("vis/ConvDiv_geometry_3", Xnsurf, comm);
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
        // std::cout << "DEBUG, in BIO" << std::endl;
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
            // std::cout << "DEBUG: Sigma mean= " << sigma_mean[0] << ", " << sigma_mean[1] << ", " << sigma_mean[2] << std::endl;

            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        // std::cout << "DEBUG: After compute potentials, before DL add sigma" << std::endl;
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> sigma;

    sctl::Vector<Real> rhs(X0.Dim());
    if (if_noslip) {
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
        rhs = bg_flow(X0);
    } else {
        rhs = GetVslip_spheres(Xnsurf, ptcls_Xcs, ptcls_thetas, ptcls_phis, Nelem, ElemOrder, FourierOrder);
        elem_lst0.WriteVTK("vis/channel_vslip_3", rhs, comm);
        // Check that u dot n is always zero
        for (int ii=0; ii<X0.Dim()/3; ii++) {
            const sctl::Vector<Real> vslip_here(3, (sctl::Iterator<Real>) rhs.begin() + ii*3, false);
            const sctl::Vector<Real> xn_here(3, (sctl::Iterator<Real>) Xnsurf.begin() + ii*3, false);
            const Real udotn = vslip_here[0]*xn_here[0] + vslip_here[1]*xn_here[1] + vslip_here[2]*xn_here[2];
            if (sctl::fabs(udotn) > 1e-6) {
                // std::cout << "u dot n at node " << ii << " is nonzero: " << udotn << ". For context, Nelem_channel = " << Nelem_channel << ", Nelem_ptcl = " << ptcl_ord << ", FourierORder = " << FourierOrder << std::endl;
                std::cout << "u dot n at node " << ii << " is nonzero: " << udotn << ". Normal is " << xn_here[0] << ", " << xn_here[1] << ", " << xn_here[2] << ", vslip is " << vslip_here[0] << ", " << vslip_here[1] << ", " << vslip_here[2] << std::endl;
                return;
            }
        }
    }

    solver(&sigma, BIO, rhs, gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);

}

template <class Real> void test_debug_channel(bool if_DL, bool if_noslip, sctl::Comm comm) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    Real DL_scal = 0;
    if (if_DL) {
        DL_scal = 1.;
    }
    // Set quadrature parameters
    const sctl::Long Nelem_channel = 100;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 96;
    const Real tol = 1e-10; // quadrature tolerance
    // Set geometry parameters
    const Real pressure_drop = -1.0; // background pressure driven flow.
    const Real period_length = 1; // length of periodic box.
    const Real channel_r1 = 0.025; 
    const Real channel_r2 = 0.1;
    // Set GMRES parameters
    const Real gmres_tol = 1e-6;
    const sctl::Long gmres_max_iter = 400;

    // Create channel + particles object: "elem_lst0".
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls(1); // length of ptcls will get updated through the conv-div channel constructor. (this is not the case for all channel constructors.)
    sctl::Long ptcl_ord = 1; 
    ptcls = ptcl_ord;
    // Arrays to hold particle...
    sctl::Vector<Real> ptcls_Xcs; // ...center positions
    sctl::Vector<Real> ptcls_rs; // ...radii
    sctl::Vector<Real> ptcls_u0s; // .. ~ aspect ratio
    sctl::Vector<Real> ptcls_thetas, ptcls_phis;
    sctl::Vector<sctl::Long> ptcls_ifprolate;
    sctl::SlenderElemList<Real> elem_lst0; 
    sctl::Vector<Real> NormalOrient; // -1 if outward normal on surface points into fluid, +1 otherwise. (since BIO defined with a -1/2*sigma already)
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.build_conv_div(Nelem_channel, ElemOrder, FourierOrder, channel_r1, channel_r2, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcls_u0s, ptcls_ifprolate, ptcl_ord);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    ptcls_thetas.ReInit(ptcls_rs.Dim());
    ptcls_phis.ReInit(ptcls_rs.Dim());
    ptcls_thetas.SetZero();
    ptcls_phis.SetZero();
    std::cout << "done with forming conv div channel, number of particles is " << ptcls.Dim() << std::endl;

    sctl::Vector<Real> X0, Xnsurf; // target coordinates
    elem_lst0.GetNodeCoord(&X0, &Xnsurf, nullptr);
    // elem_lst0.WriteVTK("vis/ConvDiv_geometry_3", Xnsurf, comm);
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
        // std::cout << "DEBUG, in BIO" << std::endl;
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
            // std::cout << "DEBUG: Sigma mean= " << sigma_mean[0] << ", " << sigma_mean[1] << ", " << sigma_mean[2] << std::endl;

            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        // std::cout << "DEBUG: After compute potentials, before DL add sigma" << std::endl;
        if (DL_scal && U->Dim() == sigma.Dim()){
            // std::cout << "in DL jump condition" << std::endl;
            (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer
        } 

        AddConstVec(*U, sigma_mean);
    };

    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> sigma;

    sctl::Vector<Real> rhs(X0.Dim());
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
    if (if_noslip) {
        rhs = bg_flow(X0) * (pressure_drop/period_length);
    } else {
        rhs = GetVslip_spheres(Xnsurf, ptcls_Xcs, ptcls_thetas, ptcls_phis, ptcl_ord, ElemOrder, FourierOrder);
        elem_lst0.WriteVTK("vis/channel_vslip_3", rhs, comm);
        // Check that u dot n is always zero
        for (int ii=0; ii<X0.Dim()/3; ii++) {
            const sctl::Vector<Real> vslip_here(3, (sctl::Iterator<Real>) rhs.begin() + ii*3, false);
            const sctl::Vector<Real> xn_here(3, (sctl::Iterator<Real>) Xnsurf.begin() + ii*3, false);
            const Real udotn = vslip_here[0]*xn_here[0] + vslip_here[1]*xn_here[1] + vslip_here[2]*xn_here[2];
            if (sctl::fabs(udotn) > 1e-6) {
                // std::cout << "u dot n at node " << ii << " is nonzero: " << udotn << ". For context, Nelem_channel = " << Nelem_channel << ", Nelem_ptcl = " << ptcl_ord << ", FourierORder = " << FourierOrder << std::endl;
                std::cout << "u dot n at node " << ii << " is nonzero: " << udotn << ". Normal is " << xn_here[0] << ", " << xn_here[1] << ", " << xn_here[2] << ", vslip is " << vslip_here[0] << ", " << vslip_here[1] << ", " << vslip_here[2] << std::endl;
                return;
            }
        }
    }

    solver(&sigma, BIO, rhs, gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);
    elem_lst0.WriteVTK("vis/ConvDiv_density_3", sigma, comm);

    // Evaluation
    {
        // // Create a uniform grid in the unit cube below the planes, then filter out points inside the particle
        // PeriodicGeom<Real> trg;    
        // CubeVolumeVisShifted<Real> vol_vis(60, 0.9, comm);
        // sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        // sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        // std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
        // X0 = std::get<0>(trg_tuple);
        // filtered_inds = std::get<1>(trg_tuple);

        // Create a new conv-div channel with no particle inside to sample target points from.
        PeriodicGeom<Real> trg;
        const sctl::Long Nelem_trg = 80;
        const sctl::Long FourierOrder_trg = 64;
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg; // place holder arrays for particle locations, but no initial length tells constructor not to include any particles.
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        sctl::Vector<Real> ptcls_u0s_trg;
        sctl::Vector<sctl::Long> ptcls_ifprolate_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build_trg = trg.build_conv_div(Nelem_trg, ElemOrder, FourierOrder_trg, channel_r1, channel_r2, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, ptcls_u0s_trg, ptcls_ifprolate_trg, ptcl_ord);
        elem_lst_trg = std::get<0>(build_trg);

        // Generate interior sample points inside the channel, then filter out targets that landed within a particle (of actual geometry).
        VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        sctl::Vector<Real> X0_all = vol_vis.GetCoord(); // set new target coordinates
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        // std::cout << "size of target points before filter: " << X0_all.Dim()<< std::endl;
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0); // ChocoFeb2026: filter spheres, vslip debug
        // std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_spheroids(X0_all, ptcls_rs, ptcls_u0s, ptcls_Xcs, ptcls_ifprolate);
        // std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_spheroids_rotated(X0_all, ptcls_rs, ptcls_u0s, ptcls_Xcs, ptcls_ifprolate, ptcls_thetas, ptcls_phis);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);
        // std::cout << "size of target points after filter: " << X0.Dim()<< std::endl;

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U(X0.Dim());
        BIO(&U, sigma);
        if (if_noslip) {
            U -= bg_flow(X0) * (pressure_drop/period_length);
        }
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
        vol_vis.WriteVTK("vis/ConvDiv_U_3", U_vis);
    }

}

template <class Real> void test_nochannel_spheroids(bool if_DL, bool if_noslip, sctl::Comm comm) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    Real DL_scal = 0;
    if (if_DL) {
        DL_scal = 1.;
    }
    // Set quadrature parameters
    const sctl::Long Nelem = 4;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 64;
    const Real tol = 1e-10; // quadrature tolerance
    // Set geometry parameters
    const Real pressure_drop = -1.0; // background pressure driven flow.
    const Real period_length = 1; // length of periodic box.
    // Set GMRES parameters
    const Real gmres_tol = 1e-6;
    const sctl::Long gmres_max_iter = 200;

    // System of spheroids
    const sctl::Long N = 20;
    sctl::Vector<Real> ptcls_Xcs, ptcls_u0s, ptcls_rs, ptcls_thetas, ptcls_phis;
    sctl::Vector<sctl::Long> ptcls_ifprolate;
    std::string filename_Xcenter, filename_u0, filename_size, filename_ifprolate;
    filename_Xcenter = "data/bie_spheroids/Xcenter_"+std::to_string(N)+".txt";
    filename_u0 = "data/bie_spheroids/u0_"+std::to_string(N)+".txt";
    filename_size = "data/bie_spheroids/size_"+std::to_string(N)+".txt";
    filename_ifprolate = "data/bie_spheroids/ifprolate_"+std::to_string(N)+".txt";

    const auto read_file_double = [](const std::string filename, sctl::Vector<Real>& Out, const sctl::Long expected_length) {
        Out.ReInit(expected_length);
        Out.SetZero();
        
        std::ifstream infile(filename);
        if (!infile) {
            std::cerr << "Error opening file " << filename << std::endl;
            SCTL_ASSERT(false);
        }

        for (sctl::Long i=0; i<expected_length; i++ ) {
            if (!(infile >> Out[i])) {
                std::cerr << "File cut short at " << i << ", before expected length " << expected_length << std::endl;
            }
        }
    };
    const auto read_file_int = [](const std::string filename, sctl::Vector<sctl::Long>& Out, const sctl::Long expected_length) {
        Out.ReInit(expected_length);
        Out.SetZero();
        
        std::ifstream infile(filename);
        if (!infile) {
            std::cerr << "Error opening file " << filename << std::endl;
            SCTL_ASSERT(false);
        }

        for (sctl::Long i=0; i<expected_length; i++ ) {
            if (!(infile >> Out[i])) {
                std::cerr << "File cut short at " << i << ", before expected length " << expected_length << std::endl;
            }
        }
    };

    read_file_double(filename_Xcenter, ptcls_Xcs, N * 3);
    read_file_double(filename_u0, ptcls_u0s, N);
    read_file_double(filename_size, ptcls_rs, N);
    read_file_int(filename_ifprolate, ptcls_ifprolate, N);

    sctl::Vector<sctl::Long> ptcls(N);
    ptcls = Nelem;

    srand48(2);
    for (sctl::Long i=0; i<N; i++) {
        const Real theta_rotate = drand48() * sctl::const_pi<Real>() * 2.;
        const Real phi_rotate = drand48() * sctl::const_pi<Real>();
        ptcls_thetas.PushBack(theta_rotate);
        ptcls_phis.PushBack(phi_rotate);
    }
    sctl::SlenderElemList<Real> elem_lst0 = spheroid_system(Nelem, ElemOrder, FourierOrder, ptcls_Xcs, ptcls_ifprolate, ptcls_u0s, ptcls_rs, ptcls_thetas, ptcls_phis, comm);

    sctl::Vector<Real> X0, Xnsurf; // target coordinates
    elem_lst0.GetNodeCoord(&X0, &Xnsurf, nullptr);
    elem_lst0.WriteVTK("vis/ConvDiv_geometry_3", Xnsurf, comm);
    sctl::Vector<Real> NormalOrient(X0.Dim());
    NormalOrient = -1.;

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
        // std::cout << "DEBUG, in BIO" << std::endl;
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
            // std::cout << "DEBUG: Sigma mean= " << sigma_mean[0] << ", " << sigma_mean[1] << ", " << sigma_mean[2] << std::endl;

            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        // std::cout << "DEBUG: After compute potentials, before DL add sigma" << std::endl;
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> sigma;

    sctl::Vector<Real> rhs(X0.Dim());
    if (if_noslip) {
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
        rhs = bg_flow(X0);
    } else {
        rhs = GetVslip_spheroids(Xnsurf, ptcls_Xcs, ptcls_rs, ptcls_u0s, ptcls_thetas, ptcls_phis, ptcls_ifprolate, Nelem, ElemOrder, FourierOrder);
        
        elem_lst0.WriteVTK("vis/channel_vslip_3", rhs, comm);
        // Check that u dot n is always zero
        for (int ii=0; ii<X0.Dim()/3; ii++) {
            const sctl::Vector<Real> vslip_here(3, (sctl::Iterator<Real>) rhs.begin() + ii*3, false);
            const sctl::Vector<Real> xn_here(3, (sctl::Iterator<Real>) Xnsurf.begin() + ii*3, false);
            const Real udotn = vslip_here[0]*xn_here[0] + vslip_here[1]*xn_here[1] + vslip_here[2]*xn_here[2];
            if (sctl::fabs(udotn) > 1e-6) {
                // std::cout << "u dot n at node " << ii << " is nonzero: " << udotn << ". For context, Nelem_channel = " << Nelem_channel << ", Nelem_ptcl = " << ptcl_ord << ", FourierORder = " << FourierOrder << std::endl;
                std::cout << "u dot n at node " << ii << " is nonzero: " << udotn << ". Normal is " << xn_here[0] << ", " << xn_here[1] << ", " << xn_here[2] << ", vslip is " << vslip_here[0] << ", " << vslip_here[1] << ", " << vslip_here[2] << std::endl;
                return;
            }
        }
    }

    solver(&sigma, BIO, rhs, gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);

}

int main(int argc, char** argv) {

    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Comm comm = sctl::Comm::World();
        long if_dl = std::stol(argv[1]);
        long if_noslip = std::stol(argv[2]);
        test_function<Real>((if_dl==1), (if_noslip==1), comm);
        // test_nochannel_spheres<Real>((if_dl==1), (if_noslip==1), comm);
        // test_debug_channel<Real>((if_dl==1), (if_noslip==1), comm);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}