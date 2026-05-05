// Solves an empty periodic channel problem, then save files of time evolution of Lagrangian particles following the solution flow.
// Using preconditioner.

#include "utils.hpp"

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

        return std::make_tuple(x_,y_,z_);
    };

    Real min_dist2 = 10.;
    Real closest_x = 0.;
    const int N = 4000; // resolution of the sampling
    for (int i = 0; i <= N; i++) {
        Real x = (Real)i / N; // TODO: account for distributed memory for x \in (a,b) instead of (0,1).
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

template <class Real> void trefoil_dispersion(sctl::Long Nelem_channel, sctl::Long FourierOrder, sctl::Comm comm) {
    
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real tol = 1e-14;
    const Real gmres_tol = 1e-9;
    const sctl::Long ElemOrder = 10;
    const Real period_length = 1.;
    const Real pressure_drop = -15.;
    const sctl::Long gmres_max_iter = 400;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient;
    sctl::Long geom_mode = 0;
    sctl::Long ptcl_ord = 1;

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_trefoil(Nelem_channel, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, geom_mode);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);

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

    std::string sigma_file = "out/trefoil_density_"+std::to_string(Nelem_channel)+"_"+std::to_string(FourierOrder)+"_"+std::to_string(comm.Rank())+".txt";
    sctl::Vector<Real> sigma;
    sigma.Read(sigma_file.c_str());

    if (!sigma.Dim()) {
        std::cout << "rank " << comm.Rank() << "couldn't read density, making new." << std::endl;
        // Preconditioning
        sctl::Matrix<Real> PrecondMat0, PrecondMat1;
        sctl::Long A11size;
        Real channel_radius = 0.035; // for trefoil geometry
        A11size = precond_channel(PrecondMat0, PrecondMat1, Nelem_channel, ElemOrder, FourierOrder, channel_radius, SL_scal, DL_scal, comm);
        // Apply A11inv to each panel of a vector.
        const auto AinvApply = [&PrecondMat0,&PrecondMat1,&A11size](const sctl::Vector<Real>& vec) {
            sctl::Vector<Real> AinvVec(vec.Dim());
            sctl::Long N = vec.Dim();
            sctl::Long Npanels = N / A11size; 
            for (sctl::Long i=0; i<Npanels; i++) {
                sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
                sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
                for (sctl::Long j=0; j<A11size; j++) {
                    AinvVec[i*A11size + j] = AinvVecMat(j,0);
                }
            }
            return AinvVec;
        };
        // Left diagonal preconditioner
        const auto BIO_precond = [&BIO,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
            sctl::Vector<Real> Uloc;
            BIO(&Uloc,sigma);
            (*U) = AinvApply(Uloc);
        };

        sctl::GMRES<Real> solver(comm);
        sctl::KrylovPrecond<Real> krylov;
        sctl::Vector<Real> A11invF = AinvApply(bg_flow(X0)*pressure_drop/period_length);
        solver(&sigma, BIO_precond, A11invF, gmres_tol, gmres_max_iter, false, nullptr, &krylov);

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
        U0 -= bg_flow(X0)*pressure_drop/period_length;
        XsectVis.WriteVTK("vis/disperison/XsectionVis_t0",U0);

        Real T = 100000.;
        sctl::Long Nt = 10000;
        Real dt = T / Nt; 
        for (sctl::Long tind = 1; tind <= Nt; tind++) {
            sctl::Vector<Real> U;
            if (tind > 1) {
                LayerPotenOp0.SetTargetCoord(X0);
                BIO(&U, sigma);
                U -= bg_flow(X0)*pressure_drop/period_length;
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
                    udotn = udotn / n2; 
                    // buffer: u_new = u - 1.01*u_n
                    Real U_x = U[xind*3+0] - 1.01*udotn*nx;
                    Real U_y = U[xind*3+1] - 1.01*udotn*ny;
                    Real U_z = U[xind*3+2] - 1.01*udotn*nz;
                    Real new_x = X0[xind*3+0]+dt*U_x;
                    Real new_y = X0[xind*3+1]+dt*U_y;
                    Real new_z = X0[xind*3+2]+dt*U_z;
                    auto [is_in_trefoil2, xc2, yc2, zc2] = in_trefoil(new_x, new_y, new_z);
                    if (is_in_trefoil2) { // double check new position back in trefoil, otherwise stays put.
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
                XsectVis.WriteVTK("vis/dispersion/XsectionVis_t"+std::to_string(tind),U);
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

