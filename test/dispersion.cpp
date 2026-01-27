// Solves an empty periodic channel problem, then save files of time evolution of Lagrangian particles following the solution flow.
// Using preconditioner.

#include "periodize.hpp"
#include "utils.hpp"

/**
 * Background flow with unit pressure gradient along X-axis.
 */
template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
    const Real dp = 15;
    const sctl::Long N = X.Dim()/3;
    sctl::Vector<Real> U(N*3);
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        U[i*3+0] = -dp*((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
        U[i*3+1] = 0;
        U[i*3+2] = 0;
    }
    return U;
}

/**
 * Reference solution for checking error.
 */
template <class Real> sctl::Vector<Real> u_ref(const sctl::Vector<Real>& X) {
  // TODO: check change with variable dp.
  const Real dp = -50;
  const sctl::Long N = X.Dim()/3;
  sctl::Vector<Real> U(N*3);
  for (sctl::Long i = 0; i < N; i++) {
    const auto x = X.begin() + i*3;
    U[i*3+0] = dp*(1e-2 - ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4);
    U[i*3+1] = 0;
    U[i*3+2] = 0;
  }
  return U;
}

template <class Real> std::tuple<bool, Real, Real, Real> in_trefoil(Real a, Real b, Real c) {
    Real r_min = 0.01;
    Real r_max = 0.04;

    if (a>1+1e-5 || a < -1e-5) { // shift x to within [0,1].
        a = a - std::floor(a);
    }

    auto get_r = [&r_min,&r_max](const Real& x) {
        Real angle = sctl::const_pi<Real>() * (16.*x - 28./3.); // =8*(t-pi/6), t = (x-0.5)*2pi
        return r_min + (r_max - r_min) * (0.5 * sctl::sin<Real>(angle) + 0.5);
    };

    auto get_xyz = [](const Real& x) {
        const Real xminus = x-0.5;
        const Real x4pi = 4.*sctl::const_pi<Real>()*xminus;
        const Real x8pi = 2.*x4pi;
        const Real xminus2 = xminus * xminus;
        const Real xminus5 = xminus2 * xminus2 * xminus;
        Real xcoeff = xminus2 * 4. - 1.;
        xcoeff = xcoeff / 5.;
        Real x_ = 0.5 * xminus * sctl::cos<Real>(x4pi) + 8. * xminus5 + 0.5;
        Real y_ = sctl::sin<Real>(x4pi) * xcoeff + 0.5;
        Real z_ = sctl::sin<Real>(x8pi) * xcoeff + 0.5;
        // std::cout << "inside getxyz, x = " << x_ << ", y = " << y_ << ", z = " << z_ << std::endl;

        return std::make_tuple(x_,y_,z_);
    };

    Real min_dist2 = 10.;
    Real closest_x = 0.;
    const int N = 4000; // resolution of the sampling
    for (int i = 0; i <= N; i++) {
        Real x = (Real)i / N; // TODO: account for distributed memory for x \in (a,b) instead of (0,1).
        // TODO: maybe look further if close to another panel, or just look +- 5 panels...
        std::tuple<Real,Real,Real> xchere = get_xyz(x);
        Real cx = std::get<0>(xchere);
        Real cy = std::get<1>(xchere);
        Real cz = std::get<2>(xchere);

        Real dx = cx - a;
        Real dy = cy - b;
        Real dz = cz - c;

        Real dist2 = dx*dx + dy*dy + dz*dz;

        if (dist2 < min_dist2) {
            min_dist2 = dist2;
            closest_x = x;
        }
    }

    Real r = get_r(closest_x);
    bool is_in_trefoil = (min_dist2 <= r*r);
    std::tuple<Real,Real,Real> closest_xyz = get_xyz(closest_x);
    Real xc = std::get<0>(closest_xyz);
    Real yc = std::get<1>(closest_xyz);
    Real zc = std::get<2>(closest_xyz);
    return std::make_tuple(is_in_trefoil, xc, yc, zc);
}

template <class Real> void trefoil_dispersion(sctl::Long Nelem_channel, sctl::Long FourierOrder, sctl::Comm comm) {
    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real tol = 1e-14;
    const Real gmres_tol = 1e-9;
    const sctl::Long ElemOrder = 10;
    const Real period_length = 1.;
    const Real pressure_drop = -1.;
    const sctl::Long gmres_max_iter = 400;

    PeriodicGeom<Real> obj;
    sctl::Long Nptcl = 0;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient;
    // set parameters (though unused)
    sctl::Long peri_mode = 1;
    sctl::Long geom_mode = 0;
    sctl::Long ptcl_ord = 1;

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_trefoil(Nelem_channel, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, geom_mode);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    Nptcl = ptcls_rs.Dim();

    sctl::Vector<Real> X0;
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    elem_lst0.WriteVTK("vis/Trefoil_dispersion_geometry",X0,comm);
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
        // if (!comm.Rank()) {
        //     std::cout << "Total surface area of channel + particles is " << surface_area  << std::endl;
        // }
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length)

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

    std::string sigma_file = "out/trefoil_density_"+std::to_string(Nelem_channel)+"_"+std::to_string(FourierOrder)+"_"+std::to_string(comm.Rank())+".txt";
    sctl::Vector<Real> sigma;
    sigma.Read(sigma_file.c_str());

    if (!sigma.Dim()) {
        std::cout << "rank " << comm.Rank() << "couldn't read density, making new." << std::endl;
        // // ======================= PRECONDITIONING : CYLINDER ====================================================
        // sctl::Vector<Real> Xc_precond, eps_precond; 
        // sctl::Vector<sctl::Long> ElemOrderVec_precond(1), FourierOrderVec_precond(1);
        // ElemOrderVec_precond[0] = ElemOrder;
        // FourierOrderVec_precond[0] = FourierOrder;
        // // Determine approximate radius of channel based on channel_mode
        // Real channel_radius = 0.035;

        // // ALTERNATIVE: smaller panel matching channel panel length and radius.
        // const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrder);
        // for (sctl::Long j = 0; j < ElemOrder; j++) { // loop over panel nodes
        //     const Real x = (nodes[j]) / Nelem_channel; // size of precond panel should be same as one panel on pipe
        //     Xc_precond.PushBack(x+0.5); //  shift panel to center of unit box, arbitrary.
        //     Xc_precond.PushBack(0.5); 
        //     Xc_precond.PushBack(0.5); 
        //     eps_precond.PushBack(channel_radius); 
        // }
        // sctl::SlenderElemList<Real> elem_lst_precond(ElemOrderVec_precond, FourierOrderVec_precond, Xc_precond, eps_precond);

        // sctl::Vector<Real> X0_precond; // target coordinates
        // elem_lst_precond.GetNodeCoord(&X0_precond, nullptr, nullptr);

        // StokesBIO Precond_bio(SL_scal, DL_scal, comm.Self());
        // Precond_bio.SetAccuracy(tol); // set quadrature accuracy
        // Precond_bio.AddElemList(elem_lst_precond);
        // Precond_bio.SetTargetCoord(X0_precond);

        // const auto BIO_1panel = [&DL_scal,&Precond_bio](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        //     U->SetZero();
        //     Precond_bio.ComputePotential(*U, sigma);
        //     (*U) -= sigma * 0.5 * DL_scal; // for preconditioner, will always be self-to-self so always add. For panels (on channel), normal orient = 1.
        // };

        // sctl::Long A11size = 3*ElemOrder*FourierOrder;
        // sctl::Vector<sctl::Vector<Real>> PrecondMat(A11size);
        // sctl::Vector<Real> SigmaCol_precond(A11size);
        // for (sctl::Long col=0; col < A11size; col ++) {
        //     SigmaCol_precond = 0.;
        //     SigmaCol_precond[col] = 1.;
        //     BIO_1panel(PrecondMat.begin()+col,SigmaCol_precond);
        // }
        // sctl::Matrix<Real> A11(A11size,A11size);
        // for (sctl::Long col=0; col < A11size; col++) {
        //     for (sctl::Long row = 0; row < A11size; row++) {
        //         A11(row,col) = PrecondMat[col][row];
        //     }
        // }

        // sctl::Matrix<Real> Usvd, VT, S, SforInv;
        // sctl::Matrix<Real> A11forSVD = sctl::Matrix<Real>(A11);
        // A11forSVD.SVD(Usvd, S, VT);
        // SforInv = sctl::Matrix<Real>(S);
        // sctl::Matrix<Real> Sinv = SforInv.pinv(1e-16);

        // // Apply A11inv to each panel of a vector.
        // const auto AinvApply = [&Usvd,&Sinv,&VT,&A11size](const sctl::Vector<Real>& vec) {
        //     sctl::Vector<Real> AinvVec(vec.Dim());
        //     sctl::Long N = vec.Dim();
        //     sctl::Long Npanels = N / A11size; 
        //     for (sctl::Long i=0; i<Npanels; i++) {
        //         sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
        //         sctl::Matrix<Real> AinvVecMat = VT.Transpose() * (Sinv * (Usvd.Transpose() * vecMat));
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

        sctl::GMRES<Real> solver(comm);
        sctl::KrylovPrecond<Real> krylov;
        solver(&sigma, BIO, -bg_flow(X0)*pressure_drop/period_length, gmres_tol, gmres_max_iter, false, nullptr, &krylov);
        // sctl::Vector<Real> A11invF = AinvApply(-bg_flow(X0));
        // solver(&sigma, BIO_precond, A11invF, gmres_tol);

        std::cout << "rank " << comm.Rank() << "done with solve, writing sigma to file named " << sigma_file << std::endl;
        sigma.Write(sigma_file.c_str());
    } 

   
    { 
        PeriodicGeom<Real> trg;
        sctl::Long Nelem_trg=200;  
        const sctl::Long FourierOrder_trg = 8; // not used
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg;
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_trefoil(Nelem_trg, ElemOrder, FourierOrder_trg, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, 1, geom_mode);
        elem_lst_trg = std::get<0>(build_trg);
        // Form targets at Ngroups cross sections, divided evenly among processes
        XsectionVis<Real> XsectVis(elem_lst_trg, comm);
        X0 = XsectVis.GetCoord();
        sctl::Vector<Real> U0;
        LayerPotenOp0.SetTargetCoord(X0);
        BIO(&U0, sigma);
        U0 += bg_flow(X0);
        XsectVis.WriteVTK("vis/XsectionVis_t0",U0);

        Real T = 100000.;
        sctl::Long Nt = 10000;
        Real dt = T / Nt; 

        // time loop
        for (sctl::Long tind = 1; tind <= Nt; tind++) {
            if (!comm.Rank()) {
                std::cout << "time step " << tind << std::endl;
            }
            // Calculate velocity at current location
            sctl::Vector<Real> U;
            if (tind > 1) {
                LayerPotenOp0.SetTargetCoord(X0);
                BIO(&U, sigma);
                U += bg_flow(X0);
            } else {
                U = U0;
            }
            for (sctl::Long xind=0; xind<X0.Dim()/3; xind++) {
                const Real current_x = X0[xind*3+0]+dt*U[xind*3+0];
                const Real current_y = X0[xind*3+1]+dt*U[xind*3+1];
                const Real current_z = X0[xind*3+2]+dt*U[xind*3+2];
                auto [is_in_trefoil, xc, yc, zc] = in_trefoil(current_x,current_y,current_z);
                if (is_in_trefoil) {
                    if (current_x > 1+1e-5 || current_x < -1e-5) { 
                        X0[xind*3+0] = current_x - std::floor(current_x);
                        X0[xind*3+1] = current_y;
                        X0[xind*3+2] = current_z;
                    } else {
                        X0[xind*3+0] = current_x;
                        X0[xind*3+1] = current_y;
                        X0[xind*3+2] = current_z;
                    }
                } else {
                    // if new point would be out of trefoil, negate normal vector for repulsive force
                    Real nx = xc - current_x;
                    Real ny = yc - current_y;
                    Real nz = zc - current_z;
                    Real n2 = nx*nx + ny*ny + nz*nz;
                    Real udotn = U[xind*3+0] * nx + U[xind*3+1] * ny + U[xind*3+2] * nz;
                    udotn = udotn / n2; // normalized n_vec
                    // u_n = udotn * <nx, ny, nz> normal direction velocity
                    // buffer by negating u_n: u_new = u - 2*u_n -- most still end up outside trefoil.
                    // buffer 2: negating u_n but also shrink by 1/100: u_new = u-1.01*u_n
                    Real U_x = U[xind*3+0] - 1.01*udotn*nx;
                    Real U_y = U[xind*3+1] - 1.01*udotn*ny;
                    Real U_z = U[xind*3+2] - 1.01*udotn*nz;
                    Real new_x = X0[xind*3+0]+dt*U_x;
                    Real new_y = X0[xind*3+1]+dt*U_y;
                    Real new_z = X0[xind*3+2]+dt*U_z;
                    auto [is_in_trefoil2, xc2, yc2, zc2] = in_trefoil(new_x, new_y, new_z);
                    if (!is_in_trefoil2) {
                        // TODO: better handling / distinguishing in this case.
                        // std::cout << " even after buffer, still ends up outside trefoil. reverting back to origional spot." << std::endl;
                    } else {
                        // std::cout << "buffer 2 worked successfully." << std::endl;
                        if (new_x > 1+1e-5 || new_x < -1e-5) { 
                            X0[xind*3+0] = new_x - std::floor(new_x);
                            X0[xind*3+1] = new_y;
                            X0[xind*3+2] = new_z;
                        } else {
                            X0[xind*3+0] = new_x;
                            X0[xind*3+1] = new_y;
                            X0[xind*3+2] = new_z;
                        }
                    }
                }
                
            }
            XsectVis.SetCoord(X0);
            if (tind % 1000 == 0) {
                XsectVis.WriteVTK("vis/XsectionVis_t"+std::to_string(tind),U);
            }
        }

    }


}

int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    sctl::Comm comm = sctl::Comm::World();
    long Nelem_channel = std::stol(argv[1]); // number of elements
    long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes

    trefoil_dispersion<Real>(Nelem_channel, FourierOrder, comm);

    
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}

