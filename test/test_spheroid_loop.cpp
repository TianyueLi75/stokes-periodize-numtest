#include "utils_geom.hpp"
#include "utils_tests.cpp"

// template <class Real> sctl::SlenderElemList<Real> spheroid_system(                                                           
//                                                                 const sctl::Long Nelem_ptcl, 
//                                                                 const sctl::Long ElemOrder, 
//                                                                 const sctl::Long FourierOrder, 
//                                                                 const sctl::Vector<Real> Xcenter_lst, 
//                                                                 const sctl::Vector<sctl::Long> if_prolate_lst, // true if prolate
//                                                                 const sctl::Vector<Real> u0_lst, 
//                                                                 const sctl::Vector<Real> r_lst, 
//                                                                 const sctl::Vector<Real> theta_lst, 
//                                                                 const sctl::Vector<Real> phi_lst, 
//                                                                 const sctl::Comm comm) 
// {
//     const sctl::Long Nptcls = u0_lst.Dim();
//     sctl::Vector<sctl::Long> ElemOrderVec(Nptcls * Nelem_ptcl);
//     sctl::Vector<sctl::Long> FourierOrderVec(Nptcls * Nelem_ptcl);
//     ElemOrderVec = ElemOrder;
//     FourierOrderVec = FourierOrder;
    
//     const auto prolate_geom = [](Real& x, Real& y, Real& z, Real& cylindrical_radius, const Real polar_angle, const Real size, const Real u0) {
//         SCTL_ASSERT(u0 > 1.); 
//         x = size * u0 * sctl::cos<Real>(polar_angle);
//         y = 0.;
//         z = 0.;
//         cylindrical_radius = size * sctl::sqrt<Real>(u0*u0 - 1) * sctl::sin<Real>(polar_angle);
//     };

//     const auto oblate_geom = [](Real& x, Real& y, Real& z, Real& cylindrical_radius, const Real polar_angle, const Real size, const Real u0) {
//         x = size * u0 * sctl::cos<Real>(polar_angle); 
//         y = 0.;
//         z = 0.;
//         cylindrical_radius = size * sctl::sqrt<Real>(u0*u0 + 1) * sctl::sin<Real>(polar_angle);
//     };

//     sctl::Vector<Real> Xc, eps;
//     for (sctl::Long pid = 0; pid < Nptcls; pid ++) {
//         const Real ptcl_size = r_lst[pid];
//         const Real ptcl_u0 = u0_lst[pid];
//         const Real theta_rotate = theta_lst[pid];
//         const Real phi_rotate = phi_lst[pid];
//         const sctl::Long if_prolate_here = if_prolate_lst[pid];
//         const Real cos_theta_rotate = sctl::cos<Real>(theta_rotate);
//         const Real sin_theta_rotate = sctl::sin<Real>(theta_rotate);
//         const Real cos_phi_rotate = sctl::cos<Real>(phi_rotate);
//         const Real sin_phi_rotate = sctl::sin<Real>(phi_rotate);
//         for (sctl::Long i=0; i < Nelem_ptcl; i++) {
//             const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
//             for (sctl::Long j=0; j<ElemOrderVec[i]; j++) {
//                 const Real theta = sctl::const_pi<Real>() * (i+nodes[j])/Nelem_ptcl;
//                 Real x, y, z, eps_j;
//                 if (if_prolate_here) {
//                     prolate_geom(x,y,z,eps_j,theta,ptcl_size,ptcl_u0);
//                 } else {
//                     oblate_geom(x,y,z,eps_j,theta,ptcl_size,ptcl_u0);
//                 }
//                 // Rz(theta) * Ry(phi) rotation
//                 Real x_rotated = cos_phi_rotate * cos_theta_rotate * x - sin_phi_rotate * y + cos_phi_rotate * sin_theta_rotate * z;
//                 Real y_rotated = sin_phi_rotate * cos_theta_rotate * x + cos_phi_rotate * y + sin_phi_rotate * sin_theta_rotate * z;
//                 Real z_rotated = - sin_theta_rotate * x + cos_theta_rotate * z;

//                 // std::cout << "rotated centerline node at " << x_rotated<<", " << y_rotated<<", " << z_rotated << std::endl;

//                 Xc.PushBack(Xcenter_lst[pid*3+0]+x_rotated);
//                 Xc.PushBack(Xcenter_lst[pid*3+1]+y_rotated);
//                 Xc.PushBack(Xcenter_lst[pid*3+2]+z_rotated);
//                 eps.PushBack(eps_j);
//             }
//         }
//     }

//     const auto init_elem_lst = [](sctl::SlenderElemList<Real>& elem_lst, const sctl::Vector<sctl::Long>& ElemOrderVec, const sctl::Vector<sctl::Long>& FourierOrderVec, const sctl::Vector<Real>& X, const sctl::Vector<Real>& R, const sctl::Comm comm) {
//         const sctl::Long Nelem_tot = ElemOrderVec.Dim();
//         sctl::Long loc_elem_cnt, loc_elem_dsp;
//         if (Nelem_tot) { // Set loc_elem_cnt, loc_elem_dsp 
//             // node_dsp.ReInit(Nelem_tot);
//             // sctl::Vector<sctl::Long> node_cnt(Nelem_tot);
//             sctl::Vector<sctl::Long> node_cnt(Nelem_tot), node_dsp(Nelem_tot); node_dsp = 0;
//             for (sctl::Long i = 0; i < Nelem_tot; i++) {
//             node_cnt[i] = ElemOrderVec[i] * FourierOrderVec[i] * FourierOrderVec[i];
//             }
//             sctl::omp_par::scan(node_cnt.begin(), node_dsp.begin(), Nelem_tot);
//             const sctl::Long Nnodes = node_cnt[Nelem_tot-1] + node_dsp[Nelem_tot-1];

//             const sctl::Long Np = comm.Size();
//             const sctl::Long rank = comm.Rank();
//             sctl::Long a = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+0)/Np) - node_dsp.begin();
//             sctl::Long b = std::lower_bound(node_dsp.begin(),  node_dsp.end(), Nnodes*(rank+1)/Np) - node_dsp.begin();
//             if (rank == Np - 1) b = Nelem_tot;
//             if (rank == 0) a = 0;
//             loc_elem_cnt = b-a;
//             loc_elem_dsp = a;
//         } else {
//             loc_elem_cnt = 0;
//             loc_elem_dsp = 0;
//         }

//         const sctl::Vector<sctl::Long> LocElemOrder(loc_elem_cnt, (sctl::Iterator<sctl::Long>)ElemOrderVec.begin() + loc_elem_dsp, false);
//         const sctl::Vector<sctl::Long> LocFourierOrder(loc_elem_cnt, (sctl::Iterator<sctl::Long>)FourierOrderVec.begin() + loc_elem_dsp, false);

//         sctl::Long dsp = 0, cnt = 0;
//         for (sctl::Long i = 0; i < loc_elem_dsp; i++) dsp += ElemOrderVec[i];
//         for (sctl::Long i = 0; i < loc_elem_cnt; i++) cnt += ElemOrderVec[loc_elem_dsp+i];
//         const sctl::Vector<Real> X_(cnt*3, (sctl::Iterator<Real>)X.begin() + dsp*3, false);
//         const sctl::Vector<Real> R_(cnt, (sctl::Iterator<Real>)R.begin() + dsp, false);

//         elem_lst.template Init<Real>(LocElemOrder, LocFourierOrder, X_, R_);  
//     };

//     sctl::SlenderElemList<Real> elem_lst;
//     init_elem_lst(elem_lst, ElemOrderVec, FourierOrderVec, Xc, eps, comm);
  
//     return elem_lst;

// }

template <class Real> void read_file(const std::string filename, sctl::Vector<Real>& Out, const sctl::Long expected_length) {
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
}

// #if 0
// // Legacy helper retained for comparison. Prefer PeriodicGeom<Real>::outside_spheroid_rotated.
// template <class Real> bool outside_spheroid(const Real x1, const Real x2, const Real x3, const Real pXc1, const Real pXc2, const Real pXc3, const Real a, const Real u0, const Real theta, const Real phi, int if_prolate) {
//     const Real v1 = x1-pXc1;
//     const Real v2 = x2-pXc2;
//     const Real v3 = x3-pXc3;
//     const Real buffer = 1.25;
    
//     // Counter rotate target -- R^{-1} = R^T
//     const Real cos_theta_rotate = sctl::cos<Real>(theta);
//     const Real sin_theta_rotate = sctl::sin<Real>(theta);
//     const Real cos_phi_rotate = sctl::cos<Real>(phi);
//     const Real sin_phi_rotate = sctl::sin<Real>(phi);
//     Real v1_rotated = cos_theta_rotate * cos_phi_rotate * v1 + sin_phi_rotate * cos_theta_rotate * v2 - sin_theta_rotate * v3;
//     Real v2_rotated = - sin_phi_rotate * v1 + cos_phi_rotate * v2;
//     Real v3_rotated = cos_phi_rotate * sin_theta_rotate * v1 + sin_phi_rotate * sin_theta_rotate * v2 + cos_theta_rotate * v3; 

//     Real A;
//     if (if_prolate) {
//         A = a * sctl::sqrt(u0*u0-1);
//     } else {
//         A = a * sctl::sqrt(u0*u0+1);
//     }
//     const Real C = a * u0;
//     const Real A2inv = 1. / (A*A); // TODO: numerical instability?
//     const Real C2inv = 1. / (C*C); 
//     const Real x1sq = v1_rotated*v1_rotated; // should be same as v1^2, etc
//     const Real x2sq = v2_rotated*v2_rotated;
//     const Real x3sq = v3_rotated*v3_rotated;
//     if (x1sq * C2inv + (x2sq+x3sq) * A2inv < buffer) { 
//         return false;
//     } else {
//         return true;
//     }
// }
// #endif

// // TEST Vslip function on spheroid
// template <class Real> sctl::Vector<Real> GetVslip(const sctl::Vector<Real>& ptcls_Xnsurf, const sctl::Vector<Real>& ptcls_Xcs, const sctl::Vector<Real>& ptcls_sizes, const sctl::Vector<Real>& ptcls_u0s, const sctl::Vector<Real>& ptcls_thetas, const sctl::Vector<Real>& ptcls_phis, const sctl::Vector<sctl::Long>& ptcls_ifprolate, const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder) {
//     const sctl::Long Nnodes_per_ptcl = Nelem * ElemOrder * FourierOrder;
//     const sctl::Long Nptcls = ptcls_sizes.Dim();
    
//     sctl::Vector<Real> Vslip(ptcls_Xnsurf.Dim());
//     srand48(1);
//     Vslip.SetZero();
//     for (sctl::Long i=0; i<Nptcls; i++) {
//         const Real a_here = ptcls_sizes[i];
//         const Real u0_here = ptcls_u0s[i];
//         const int if_prolate_here = ptcls_ifprolate[i];
//         sctl::Vector<Real> Vslip_here(3*Nnodes_per_ptcl, (sctl::Iterator<Real>) Vslip.begin() + 3*Nnodes_per_ptcl*i, false);
//         sctl::Vector<Real> center_here(3, (sctl::Iterator<Real>) ptcls_Xcs.begin() + 3*i, false);
//         // const Real scalar = drand48()*0.8 + 0.1; // randomly scaled slip velocity by (0.1,0.9).
//         const Real scalar = 1.; // no change in slip magnitude on different spheres/spheroids
//         const Real theta_rotate = ptcls_thetas[i];
//         const Real phi_rotate = ptcls_phis[i];
//         const Real cos_theta_rotate = sctl::cos<Real>(theta_rotate);
//         const Real sin_theta_rotate = sctl::sin<Real>(theta_rotate);
//         const Real cos_phi_rotate = sctl::cos<Real>(phi_rotate);
//         const Real sin_phi_rotate = sctl::sin<Real>(phi_rotate);

//         for (sctl::Long panel=0; panel < Nelem; panel++) {
//             const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrder);

//             for (sctl::Long el=0; el<ElemOrder; el++) {
//                 const Real theta = sctl::const_pi<Real>() * (panel+nodes[el])/Nelem;
//                 // const Real coeff = scalar * sctl::sin<Real>(theta); // sinusoidal slip magnitude, 0 at north and south poles
//                 const Real coeff = 1.;

//                 for (sctl::Long fl=0; fl<FourierOrder; fl++) {
//                     Real phi = 2. * sctl::const_pi<Real>() * fl / FourierOrder;
//                     if (theta_rotate == 0) {
//                         if (phi_rotate != 0) {
//                             phi = phi + sctl::const_pi<Real>()/2.;
//                         }
//                     }
//                     const sctl::Long idx = panel*ElemOrder*FourierOrder + el*FourierOrder + fl;
//                     const sctl::Vector<Real> Xn_here(3, (sctl::Iterator<Real>) ptcls_Xnsurf.begin() + idx*3, false);

//                     Real t1_unrotated = a_here * u0_here * (-sctl::sin<Real>(theta));
//                     Real t2_unrotated, t3_unrotated;
//                     if (if_prolate_here) {
//                         t2_unrotated = a_here * sctl::sqrt<Real>(u0_here*u0_here-1) * sctl::cos<Real>(theta) * sctl::cos<Real>(phi);
//                         t3_unrotated = a_here * sctl::sqrt<Real>(u0_here*u0_here-1) * sctl::cos<Real>(theta) * sctl::sin<Real>(phi);
//                     } else {
//                         t2_unrotated = a_here * sctl::sqrt<Real>(u0_here*u0_here+1) * sctl::cos<Real>(theta) * sctl::cos<Real>(phi);
//                         t3_unrotated = a_here * sctl::sqrt<Real>(u0_here*u0_here+1) * sctl::cos<Real>(theta) * sctl::sin<Real>(phi);
//                     }
                    
//                     Real mag2_tang = t1_unrotated*t1_unrotated + t2_unrotated*t2_unrotated + t3_unrotated*t3_unrotated;
//                     Real mag_tang = sctl::sqrt<Real>(mag2_tang);
//                     // tangential sinusoidal slip, spheroid-cetnered space
//                     sctl::Vector<Real> Vslip_here_unrotated(3);
//                     Vslip_here_unrotated[0] = coeff * t1_unrotated / mag_tang;
//                     Vslip_here_unrotated[1] = coeff * t2_unrotated / mag_tang;
//                     Vslip_here_unrotated[2] = coeff * t3_unrotated / mag_tang;
//                     // tagential slip, rotated
//                     sctl::Vector<Real> Vslip_here_rotated(3);
//                     Vslip_here_rotated[0] = cos_phi_rotate * cos_theta_rotate * Vslip_here_unrotated[0] - sin_phi_rotate * Vslip_here_unrotated[1] + cos_phi_rotate * sin_theta_rotate * Vslip_here_unrotated[2];
//                     Vslip_here_rotated[1] = sin_phi_rotate * cos_theta_rotate * Vslip_here_unrotated[0] + cos_phi_rotate * Vslip_here_unrotated[1] + sin_phi_rotate * sin_theta_rotate * Vslip_here_unrotated[2];
//                     Vslip_here_rotated[2] = - sin_theta_rotate * Vslip_here_unrotated[0] + cos_theta_rotate * Vslip_here_unrotated[2];

//                     Real vdotn = Vslip_here_rotated[0] * Xn_here[0] + Vslip_here_rotated[1] * Xn_here[1] + Vslip_here_rotated[2] * Xn_here[2];
//                     Vslip_here[idx*3+0] = Vslip_here_rotated[0] - vdotn * Xn_here[0];
//                     Vslip_here[idx*3+1] = Vslip_here_rotated[1] - vdotn * Xn_here[1];
//                     Vslip_here[idx*3+2] = Vslip_here_rotated[2] - vdotn * Xn_here[2];

                    
//                 }
//             }

//         }
//     }   
//     return Vslip;
// }


// Manufactured solution on spheroid or loop in free space to check their Layer potential evaluations at far and near.
// Note this does not have parallelization or periodization.

template <class Real> void function(const sctl::Long Nelem, const sctl::Long FourierOrder, const sctl::Long geom_mode, const Real gmres_tol, const Real tol, sctl::Comm comm) {
    // Set parameters
    const sctl::Long ElemOrder = 10;
    const sctl::Long gmres_max_iter = 100;

    // SCTL_ASSERT((geom_mode==1) || (geom_mode==3)); // only look at spheroids or loops

    // Create spheroid in the center of unit box (but not periodized)
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    sctl::Vector<Real> X0, Xn0;
    elem_lst0.GetNodeCoord(&X0, &Xn0, nullptr);
    // plot surface and normal
    if (geom_mode == 0) {
        elem_lst0.WriteVTK("vis/debug/sphere", Xn0); 
    } else if (geom_mode==1) {
        elem_lst0.WriteVTK("vis/debug/spheroid", Xn0); 
    } else if (geom_mode==3) {
        elem_lst0.WriteVTK("vis/debug/loop", Xn0); 
    }
    // Set layer potential on ptcl
    StokesBIO LPO(1.,1.,comm);
    LPO.AddElemList(elem_lst0);
    LPO.SetAccuracy(tol);
    // Define the boundary-integral operator: (I/2 + D + S)[sigma]
    MeanCorrectedStokesBIOOperator<Real> BIO(LPO, NormalOrient, 1.0, comm);

    // Old lambda-based BIO block retained for reference:
    // const auto BIO = [&LPO,&NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    //     U->SetZero();
    //     LPO.ComputePotential(*U, sigma);
    //     if (U->Dim() == sigma.Dim()) {
    //         (*U) -= sigma*0.5*NormalOrient; // for double-layer
    //     }
    // };

    // Since taking targets as a larger version of surface, will erroneously enter into the DL self eval loop if using BIO.
    const auto BIO_eval = [&LPO,&NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        U->SetZero();
        LPO.ComputePotential(*U, sigma);
    };

    // Set Stokeslet inside surface
    sctl::Vector<Real> Xsrc(3);
    sctl::Vector<Real> Fsrc(3);
    if (geom_mode == 1 || geom_mode == 0) {
        // Somewhat random placement of stokeslet to break symmetry, but close enough to center of spheroid.
        Xsrc[0] = 0.486;
        Xsrc[1] = 0.508;
        Xsrc[2] = 0.5026;
        Fsrc[0] = 0.4;
        Fsrc[1] = -0.5;
        Fsrc[2] = 0.8;
        // Xsrc[4] = 0.486;
        // Xsrc[5] = 0.508;
        // Xsrc[6] = 0.5026;
        // Fsrc[4] = -0.4;
        // Fsrc[5] = 0.5;
        // Fsrc[6] = -0.8;
    } else if (geom_mode == 3) {
        // May need to move after looking at visualization of the particle.
        Xsrc[0] = 0.386;
        Xsrc[1] = 0.508;
        Xsrc[2] = 0.5;
        Fsrc[0] = 0.4;
        Fsrc[1] = -0.5;
        Fsrc[2] = 0.8;
        // Xsrc[4] = 0.386;
        // Xsrc[5] = 0.508;
        // Xsrc[6] = 0.5;
        // Fsrc[4] = -0.4;
        // Fsrc[5] = 0.5;
        // Fsrc[6] = -0.8;
    }

    // exact solution
    const auto exact_U = [&Xsrc, &Fsrc](const sctl::Vector<Real> Xt) {
        sctl::Vector<Real> Ut(Xt.Dim());
        Ut.SetZero();
        
        sctl::Stokes3D_FxU ker;
        ker.Eval(Ut, Xt, Xsrc, Xsrc, Fsrc);
        return Ut;
    };

    // BC on particle
    sctl::Vector<Real> Usurf = exact_U(X0);
    // std::cout << "size of Usurf is " << Usurf.Dim() << "; first node U; " << Usurf[0] << ", " << Usurf[1] << ", " << Usurf[2] << std::endl;
    // if (geom_mode == 0) {
    //     elem_lst0.WriteVTK("vis/debug/sphereUexact", Usurf);
    // } else if (geom_mode == 1) {
    //     elem_lst0.WriteVTK("vis/debug/spheroid_Uexact", Usurf);
    // } else if (geom_mode == 3) {
    //     elem_lst0.WriteVTK("vis/debug/loop_Uexact", Usurf);
    // }
    
    // Solve Stokes BVP
    sctl::GMRES<Real> solver(comm, false);
    sctl::Vector<Real> sigma;
    LPO.SetTargetCoord(X0);
    solver(&sigma, BIO, Usurf, gmres_tol, gmres_max_iter);

    // Set target points, one layer far, one layer close.
    Real dist_far = 0.35;
    sctl::Vector<Real> Xtrg = dist_far * Xn0 + X0;
    LPO.SetTargetCoord(Xtrg);
    sctl::Vector<Real> eval_U(Xtrg.Dim());
    BIO_eval(&eval_U, sigma);
    // if (geom_mode == 0) {
    //     elem_lst0.WriteVTK("vis/debug/sphereUeval", eval_U);
    // } else if (geom_mode == 1) {
    //     elem_lst0.WriteVTK("vis/debug/spheroid_Ueval", eval_U);
    // } else if (geom_mode == 3) {
    //     elem_lst0.WriteVTK("vis/debug/loop_Ueval", eval_U);
    // }
    sctl::Vector<Real> exp_U = exact_U(Xtrg);
    // if (geom_mode == 0) {
    //     elem_lst0.WriteVTK("vis/debug/sphereUetrg", exp_U);
    // } else if (geom_mode == 1) {
    //     elem_lst0.WriteVTK("vis/debug/spheroid_Utrg", exp_U);
    // } else if (geom_mode == 3) {
    //     elem_lst0.WriteVTK("vis/debug/loop_Utrg", exp_U);
    // }
    sctl::Vector<Real> err = eval_U - exp_U;
    // for (int i=0; i<err.Dim(); i++) {
    //     std::cout << err[i] << std::endl;
    // }
    Real max_err = 0.;
    for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
    Real max_u = 0.;
    for (const auto e : exp_U) max_u = std::max<Real>(max_u, sctl::fabs(e));
    std::cout << "at distance " << dist_far << " away from surface, max error is " << std::setprecision(10) << max_err << ", max u is " << max_u << ", rel error is " << max_err / max_u << std::endl;

    Real dist_close = 0.01;
    Xtrg = dist_close * Xn0 + X0;
    LPO.SetTargetCoord(Xtrg);
    BIO_eval(&eval_U, sigma);
    exp_U = exact_U(Xtrg);
    err = eval_U - exp_U;
    // for (int i=0; i<err.Dim(); i++) {
    //     std::cout << err[i] << std::endl;
    // }
    max_err = 0.;
    for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
    max_u = 0.;
    for (const auto e : exp_U) max_u = std::max<Real>(max_u, sctl::fabs(e));
    std::cout << "at distance " << dist_close << " away from surface, max error is " << std::setprecision(10) << max_err << ", max u is " << max_u << ", rel error is " << max_err / max_u << std::endl;


}

template <class Real> void test_spheroid_tangents(const sctl::Long Nelem, const sctl::Long FourierOrder, sctl::Comm comm) {
    // Set parameters
    const sctl::Long ElemOrder = 10;
    sctl::Long N = 1;

    sctl::Vector<Real> Xcenter_all, u0_all, r_all, theta_all, phi_all;
    sctl::Vector<sctl::Long> if_prolate_all;
    Xcenter_all.PushBack(0.3);
    Xcenter_all.PushBack(0.4);
    Xcenter_all.PushBack(0.5);
    u0_all.PushBack(1.1);
    r_all.PushBack(1./1.1);
    theta_all.PushBack(sctl::const_pi<Real>()/4.);
    phi_all.PushBack(sctl::const_pi<Real>()/8.);
    // theta_all.PushBack(0);
    // phi_all.PushBack(0);
    if_prolate_all.PushBack(1);

    // Set up surface
    sctl::SlenderElemList<Real> elem_lst0 = spheroid_system<Real>(Nelem, ElemOrder, FourierOrder, Xcenter_all, if_prolate_all, u0_all, r_all, theta_all, phi_all, comm);
    sctl::Vector<Real> X0, Xn0;
    elem_lst0.GetNodeCoord(&X0, &Xn0, nullptr);
    elem_lst0.WriteVTK("vis/debug_spheroid_normal", Xn0, comm);

    sctl::Vector<Real> vslip = GetVslip(Xn0, Xcenter_all, r_all, u0_all, theta_all, phi_all, if_prolate_all, Nelem, ElemOrder, FourierOrder);
    elem_lst0.WriteVTK("vis/debug_spheroid_slip", vslip, comm);

    // Check that u dot n is always zero
    for (int ii=0; ii<X0.Dim()/3; ii++) {
        const sctl::Vector<Real> vslip_here(3, (sctl::Iterator<Real>) vslip.begin() + ii*3, false);
        const sctl::Vector<Real> xn_here(3, (sctl::Iterator<Real>) Xn0.begin() + ii*3, false);
        const Real udotn = vslip_here[0]*xn_here[0] + vslip_here[1]*xn_here[1] + vslip_here[2]*xn_here[2];
        if (sctl::fabs(udotn) > 1e-5) {
            std::cout << "udotn = " << udotn << ", vslip = " << vslip_here[0] << ", "<< vslip_here[1] << ", " << vslip_here[2] << ", xn = " << xn_here[0] << ", " << xn_here[1] << ", " << xn_here[2] << std::endl;
            // std::cout << "u dot n at node " << ii << " is nonzero: " << udotn << ". For context, Nelem_ptcl = " << Nelem << ", FourierOrder = " << FourierOrder << std::endl;
        }
    }

}


int main(int argc, char** argv) {
    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;
    // using Real=sctl::QuadReal;

    {
        sctl::Comm comm = sctl::Comm::World();
        long Nelem = std::stol(argv[1]); // =0 for trefoil, =1 for 1-particle; =2 for conv div, =3 for trefoil with particle. TODO: =4 for 2peri plane with 1 particle
        long FourierOrder = std::stol(argv[2]);
        double gmres_tol = std::stod(argv[3]);
        double tol = std::stod(argv[4]);

        // std::cout << "Running sphere manufactured solutions for Nelem = " << Nelem << ", Fourier = " << FourierOrder << ";" << std::endl;
        // function<Real>(Nelem, FourierOrder, 0, gmres_tol, tol, comm);
        // std::cout << "Running spheroid manufactured solutions for Nelem = " << Nelem << ", Fourier = " << FourierOrder << ";" << std::endl;
        // function<Real>(Nelem, FourierOrder, 1, gmres_tol, tol, comm);
        // std::cout << "Running loop manufactured solutions for Nelem = " << Nelem << ", Fourier = " << FourierOrder << ";" << std::endl;
        // function<Real>(Nelem, FourierOrder, 3, gmres_tol, tol, comm);

        test_spheroid_tangents<Real>(Nelem, FourierOrder, comm);

    }
    sctl::Comm::MPI_Finalize();
    return 0;
}