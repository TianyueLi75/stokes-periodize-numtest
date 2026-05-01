#include "periodize.hpp"
#include "utils.hpp"
#include "planeNaive.hpp"

// Test script for calculation and timing of 1, 2, and 3 periodic problems with background pressure flow. 

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

template <class Real> void test1peri_channel(sctl::Long Nelem_channel, sctl::Long FourierOrder, sctl::Comm comm, sctl::Long Nptcl, Real gmres_tol, Real tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;
    const sctl::Long ElemOrder = 10;
    const Real pressure_drop = -1.0;
    const Real period_length = 1;
    const sctl::Long geom_mode = 0;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Long ptcl_ord = 4;
    if (Nptcl>0) {
        ptcls.ReInit(Nptcl);
        ptcls = ptcl_ord;
    }
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient, ptcls_thetas, ptcls_phis;
    sctl::Long peri_mode = 1;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.build_conv_div_sph(Nelem_channel, ElemOrder, FourierOrder, 0.1, 0.2, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);
    ptcls_thetas = std::get<2>(build0);
    ptcls_phis = std::get<3>(build0);
    Nptcl = ptcls_rs.Dim(); // Number of particles could have changed after initializing.

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    // elem_lst0.WriteVTK("vis/convdiv_ptcl", X0);
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

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
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
            // sigma_mean *= (1/surface_area);
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

            // // DEBUG: check that sigma-sigma_mean has surface integral = 0:
            // sctl::Vector<Real> sigma1 = sigma_;
            // AddConstVec(sigma1, -sigma_mean);
            // sctl::Vector<Real> sigma_test_;
            // SurfaceIntegral(sigma_test_, sigma1, wts);
            // std::cout << "Surface integral of sigma - sigma bar = " << sigma_test_[0] << ", "  << sigma_test_[1] << ", "  << sigma_test_[2] << std::endl;
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    // first gmres to remove timing for matrix loading, and set Krylov preconditioner.
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;

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
    solver(&sigma, BIO, bg_flow(X0) * (pressure_drop/period_length), gmres_tol, -1, false, nullptr, &krylov_precond);

    Real channel_radius = 0.15;
    sctl::Long Ntrg_side = 5;
    Real side_len = channel_radius * sctl::sqrt<Real>(2);
    Real gap = side_len/(Ntrg_side+1);
    X0.ReInit(Ntrg_side * Ntrg_side * 3);
    // X symmetry
    sctl::Vector<Real> X1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long yind=0; yind<Ntrg_side; yind++) {
        for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
            sctl::Long Ntrg_nodeind = yind*Ntrg_side + zind;
            X0[Ntrg_nodeind * 3 + 0] = 0.;
            X0[Ntrg_nodeind * 3 + 1] = 0.5-side_len/2.+(yind+1)*gap;
            X0[Ntrg_nodeind * 3 + 2] = 0.5-side_len/2.+(zind+1)*gap;
            X1[Ntrg_nodeind * 3 + 0] = 1.;
            X1[Ntrg_nodeind * 3 + 1] = 0.5-side_len/2.+(yind+1)*gap;
            X1[Ntrg_nodeind * 3 + 2] = 0.5-side_len/2.+(zind+1)*gap;
        }
    }
    sctl::Vector<Real> UX0(X0.Dim());
    LayerPotenOp0.SetTargetCoord(X0);
    BIO(&UX0,sigma);
    sctl::Vector<Real> UX1(X1.Dim());
    LayerPotenOp0.SetTargetCoord(X1);
    BIO(&UX1,sigma);
    UX0 -= bg_flow(X0) * (pressure_drop/period_length);
    UX1 -= bg_flow(X0) * (pressure_drop/period_length);
    std::cout << "============ X periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffX = UX0-UX1;
    for (int i=0; i<UdiffX.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffX[i*3+0]<< ", " << UdiffX[i*3+1]<< ", " << UdiffX[i*3+2]<< ". " << std::endl;
    }
}

template <class Real> void test2peri_plane(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Comm comm, const Real gmres_tol, const Real tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real pressure_drop = -10.0;
    const Real period_length = 1;

    const sctl::Long ElemOrder = 10;
    // const sctl::Long geom_mode = 0;
    const sctl::Long geom_mode = 1; // spheroids
    // const sctl::Long Nptcl = 1; // Other spherical configurations are not checked to remain between z_offsets to stay clear of the plates, so only check periodicity for one particle.
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient, ptcls_thetas, ptcls_phis;
    if (geom_mode == 0) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls3(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else if (geom_mode == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.many_spheroids3(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
        ptcls_thetas = std::get<2>(build0);
        ptcls_phis = std::get<3>(build0);
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
    elem_lst0.WriteVTK("vis/plane-ptcl-geometry", X0_ptcl, comm);
    
    // Plane
    // TODO: put plane on only Rank 0 MPI?
    sctl::Vector<Real> X0_wall;
    const sctl::Long gl_order = 49;
    const sctl::Long Nelem_x = 2;
    const sctl::Long Nelem_y = 2;
    const Real z_offset = 0.01;
    sctl::PlaneIntegral<Real> plane(gl_order, Nelem_x, Nelem_y, z_offset);
    
    plane.GetNodeCoord(&X0_wall, nullptr, nullptr);
    LayerPotenOp0.AddElemList(plane,"2");
    plane.WriteVTK("vis/plane-geometry", X0_wall, comm);

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
        // surface_area = surface_area_[0];
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
        // surface_area_wall = surface_area_[0];
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area_wall = sa_all[0];
        // std::cout << "DEBUG wall surface area, computed to be " << surface_area_wall << std::endl;
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


    // // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    // const auto BIO_eval = [&wts,&surface_area,&elem_lst0,&wts_wall,&surface_area_wall,&plane,&LayerPotenOp0,&DL_scal,&X0_ptcl,&X0_wall,&NormalOrient, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        
    //     sctl::Vector<Real> ptcl_dens(X0_ptcl.Dim(), (sctl::Iterator<Real>) sigma.begin(), true);
    //     sctl::Vector<Real> wall_dens(X0_wall.Dim(), (sctl::Iterator<Real>) sigma.begin()+ptcl_dens.Dim(), true);
            
    //     sctl::Vector<Real> sigma_mean, sigma0;

    //     { // compute sigma_mean and sigma0 = sigma - sigma_mean
    //         sctl::Vector<Real> sigma_mean_ptcl, sigma_mean_wall;
    //         sctl::Vector<Real> sigma_ptcl_;
    //         elem_lst0.GetFarFieldDensity(sigma_ptcl_, ptcl_dens);
    //         SurfaceIntegral(sigma_mean_ptcl, sigma_ptcl_, wts);
    //         sctl::Vector<Real> sigma_wall_ = wall_dens;
    //         SurfaceIntegral(sigma_mean_wall, sigma_wall_, wts_wall);
    //         // MPI
    //         sctl::Vector<Real> sa_loc(6);
    //         for (int i=0; i<3; i++) {
    //             sa_loc[i] = sigma_mean_ptcl[i];
    //         }
    //         for (int i=0; i<3; i++) {
    //             sa_loc[i+3] = sigma_mean_wall[i];
    //         }
    //         sctl::Vector<Real> sa_all(6);
    //         sa_all = 0.;
    //         for (int i=0; i<3; i++) {
    //             comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+i, (sctl::Iterator<Real>) sa_all.begin()+i, 1, sctl::CommOp::SUM);
    //         }
    //         for (int i=0; i<3; i++) {
    //             comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+i+3, (sctl::Iterator<Real>) sa_all.begin()+i+3, 1, sctl::CommOp::SUM);
    //         }
    //         // sigma_mean = sigma_mean_ptcl + sigma_mean_wall; 
    //         sigma_mean.ReInit(sigma_mean_ptcl.Dim());
    //         for (int i=0; i<3; i++) {
    //             sigma_mean[i] = sa_all[i] + sa_all[i+3]; // adding total sigma in {x,y,z} coord separately
    //         }
    //         sigma_mean *= (1/(surface_area + surface_area_wall));
    //         sigma0 = sigma;
    //         AddConstVec(sigma0, -sigma_mean);

    //     }

    //     U->SetZero();
    //     LayerPotenOp0.ComputePotential(*U, sigma0);

    //     AddConstVec(*U, sigma_mean);
    // };

    /*
        // Getting singular values
        sctl::Long Nsrc = X0_ptcl.Dim()/3 + X0_wall.Dim()/3;
        sctl::Vector<Real> sigma_eye(Nsrc*3);
        sctl::Vector<sctl::Vector<Real>> LPOvecvec(Nsrc*3);
        for (sctl::Long i=0; i<Nsrc; i++) {
            for (sctl::Long k=0; k<3; k++) {
                sigma_eye.SetZero();
                sigma_eye[i*3+k] = 1.;
                // std::cout << "Node number is = " << i << ", dimension = " << k << std::endl;
                BIO(LPOvecvec.begin()+i*3+k, sigma_eye);
            }
        }
        // SVD
        sctl::Matrix<Real> LPOmat(Nsrc*3,Nsrc*3);
        for (long i=0; i < Nsrc*3; i++) {
            for (long j = 0; j < Nsrc*3; j++) {
                LPOmat(j,i) = LPOvecvec[i][j];
            }
        }      
        sctl::Matrix<Real> Usvd_p, VT_p, S_p;
        sctl::Matrix<Real> LPOforSVD = sctl::Matrix<Real>(LPOmat);
        LPOforSVD.SVD(Usvd_p, S_p, VT_p);
        std::cout << "debug by printing the last 10 matrix singular values:" << std::endl;
        std::cout << "shape of S_p is " << S_p.Dim(0) << ", " << S_p.Dim(1) << std::endl;
        for (long i=S_p.Dim(0)-10; i<S_p.Dim(0); i++) {
            std::cout << S_p(i,i) << std::endl;
        }
    */

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
    solver(&sigma,BIO, bg_flow(X0) * (pressure_drop/period_length), gmres_tol);

    // /*
        sctl::Long Ntrg_side = 5;
        Real gap = 1./(Ntrg_side+5); // TRY: taking targets farther from plates
        X0.ReInit(Ntrg_side * Ntrg_side * 3);
        // X symmetry
        sctl::Vector<Real> X1(Ntrg_side * Ntrg_side * 3);
        for (sctl::Long yind=0; yind<Ntrg_side; yind++) {
            for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
                sctl::Long Ntrg_nodeind = yind*Ntrg_side + zind;
                X0[Ntrg_nodeind * 3 + 0] = 0.;
                X0[Ntrg_nodeind * 3 + 1] = (yind+3)*gap; // Shift to start further from the plates.
                X0[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
                X1[Ntrg_nodeind * 3 + 0] = 1.;
                X1[Ntrg_nodeind * 3 + 1] = (yind+3)*gap;
                X1[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
                
            }
        }
        sctl::Vector<Real> UX0(X0.Dim());
        LayerPotenOp0.SetTargetCoord(X0);
        // LayerPotenOp2.SetTargetCoord(X0);
        BIO(&UX0,sigma);
        sctl::Vector<Real> UX1(X1.Dim());
        LayerPotenOp0.SetTargetCoord(X1);
        // LayerPotenOp2.SetTargetCoord(X1);
        BIO(&UX1,sigma);
        UX0 -= bg_flow(X0) * (pressure_drop/period_length);
        UX1 -= bg_flow(X1) * (pressure_drop/period_length);
        std::cout << "============ X periodicity =================" << std::endl;
        sctl::Vector<Real> UdiffX = UX0-UX1;
        for (int i=0; i<UdiffX.Dim()/3; i++) {
            std::cout << std::setprecision(8) << UdiffX[i*3+0]<< ", " << UdiffX[i*3+1]<< ", " << UdiffX[i*3+2]<< ". " << std::endl;
        }

        // Y symmetry
        sctl::Vector<Real> Y0(Ntrg_side * Ntrg_side * 3);
        sctl::Vector<Real> Y1(Ntrg_side * Ntrg_side * 3);
        for (sctl::Long xind=0; xind<Ntrg_side; xind++) {
            for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
                sctl::Long Ntrg_nodeind = xind*Ntrg_side + zind;
                Y0[Ntrg_nodeind * 3 + 0] = (xind+3)*gap; // TRY: Changed here as well.
                Y0[Ntrg_nodeind * 3 + 1] = 0.;
                Y0[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
                Y1[Ntrg_nodeind * 3 + 0] = (xind+3)*gap;
                Y1[Ntrg_nodeind * 3 + 1] = 1.;
                Y1[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
            }
        }
        sctl::Vector<Real> UY0(Y0.Dim());
        LayerPotenOp0.SetTargetCoord(Y0);
        // LayerPotenOp2.SetTargetCoord(Y0);
        BIO(&UY0, sigma); 
        sctl::Vector<Real> UY1(Y1.Dim());
        LayerPotenOp0.SetTargetCoord(Y1);
        // LayerPotenOp2.SetTargetCoord(Y1);
        BIO(&UY1, sigma); 
        std::cout << "============ Y periodicity =================" << std::endl;
        sctl::Vector<Real> UdiffY = UY0-UY1;
        for (int i=0; i<UdiffY.Dim()/3; i++) {
            std::cout << std::setprecision(8) << UdiffY[i*3+0]<< ", " << UdiffY[i*3+1]<< ", " << UdiffY[i*3+2]<< ". " << std::endl;
        }
    // */

    // /*
    {
        PeriodicGeom<Real> trg;    
        CubeVolumeVisShifted<Real> vol_vis(10, 0.9, comm);
        // X0 = vol_vis.GetCoord();
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple;
        if (geom_mode > 0) {
            trg_tuple = trg.filter_target_rotated(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode, ptcls_thetas, ptcls_phis);
        } else {
            trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        }
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        LayerPotenOp0.SetTargetCoord(X0);
        // LayerPotenOp2.SetTargetCoord(X0);
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
        vol_vis.WriteVTK("vis/plane-U-2p-Ubg", U_vis);
    }
    // */
    
}

template <class Real> void test3peri(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Comm comm, sctl::Long Nptcl, const Real gmres_tol, const Real tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real pressure_drop = -1.0;
    const Real period_length = 1;

    const sctl::Long ElemOrder = 10;
    const sctl::Long geom_mode = 0;
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
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
    Nptcl = ptcls_rs.Dim(); 

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

    // DEBUG MPI Surface area:
    if (!comm.Rank()) {
        Real surfA_manual = 0.;
        for (sctl::Long i=0; i<Nptcl; i++) {
            Real r_i = ptcls_rs[i];
            surfA_manual += 4.*sctl::const_pi<Real>() * r_i * r_i;
        }
        std::cout << "Surface area computed for a total of " << Nptcl << " spheres is " << surface_area << "; manual calculation gives " << surfA_manual << std::endl;
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0,"1");
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

    const auto eval_rhs2 = [](const sctl::Vector<Real>& X) {
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

    // first gmres to remove timing for matrix loading, and set Krylov preconditioner.
    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    LayerPotenOp0.SetTargetCoord(X0);
    solver(&sigma, BIO, eval_rhs(pressure_drop), gmres_tol);
    // solver(&sigma, BIO, eval_rhs2(X0), gmres_tol);

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
        // U -= eval_rhs2(X0);

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
        vol_vis.WriteVTK("vis/"+std::to_string(Nptcl)+"spheres_U", U_vis); 
    }

    sctl::Long Ntrg_side = 5;
    Real gap = 1./(Ntrg_side+5);
    X0.ReInit(Ntrg_side * Ntrg_side * 3);
    // X symmetry
    sctl::Vector<Real> X1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long yind=0; yind<Ntrg_side; yind++) {
        for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
            sctl::Long Ntrg_nodeind = yind*Ntrg_side + zind;
            X0[Ntrg_nodeind * 3 + 0] = 0.;
            X0[Ntrg_nodeind * 3 + 1] = (yind+3)*gap;
            X0[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
            X1[Ntrg_nodeind * 3 + 0] = 0.9999999999;
            X1[Ntrg_nodeind * 3 + 1] = (yind+3)*gap;
            X1[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
        }
    }
    sctl::Vector<Real> UX0(X0.Dim());
    LayerPotenOp0.SetTargetCoord(X0);
    BIO(&UX0,sigma);
    UX0 -= eval_rhs(pressure_drop);
    sctl::Vector<Real> UX1(X1.Dim());
    LayerPotenOp0.SetTargetCoord(X1);
    BIO(&UX1,sigma);
    UX1 -= eval_rhs(pressure_drop);
    std::cout << "============ X periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffX = UX0-UX1;
    for (int i=0; i<UdiffX.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffX[i*3+0]<< ", " << UdiffX[i*3+1]<< ", " << UdiffX[i*3+2]<< ". " << std::endl;
    }
    // Y symmetry
    sctl::Vector<Real> Y0(Ntrg_side * Ntrg_side * 3);
    sctl::Vector<Real> Y1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long xind=0; xind<Ntrg_side; xind++) {
        for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
            sctl::Long Ntrg_nodeind = xind*Ntrg_side + zind;
            Y0[Ntrg_nodeind * 3 + 0] = (xind+3)*gap;
            Y0[Ntrg_nodeind * 3 + 1] = 0.0;
            Y0[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
            Y1[Ntrg_nodeind * 3 + 0] = (xind+3)*gap;
            Y1[Ntrg_nodeind * 3 + 1] = 0.9999999999;
            Y1[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
        }
    }
    sctl::Vector<Real> UY0(Y0.Dim());
    LayerPotenOp0.SetTargetCoord(Y0);
    BIO(&UY0,sigma);
    UY0 -= eval_rhs(pressure_drop);
    sctl::Vector<Real> UY1(Y1.Dim());
    LayerPotenOp0.SetTargetCoord(Y1);
    BIO(&UY1,sigma);
    UY1 -= eval_rhs(pressure_drop);
    std::cout << "============ Y periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffY = UY0-UY1;
    for (int i=0; i<UdiffY.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffY[i*3+0]<< ", " << UdiffY[i*3+1]<< ", " << UdiffY[i*3+2]<< ". " << std::endl;
    }
    // Z symmetry
    sctl::Vector<Real> Z0(Ntrg_side * Ntrg_side * 3);
    sctl::Vector<Real> Z1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long xind=0; xind<Ntrg_side; xind++) {
        for (sctl::Long yind=0; yind<Ntrg_side; yind++) {
            sctl::Long Ntrg_nodeind = xind*Ntrg_side + yind;
            Z0[Ntrg_nodeind * 3 + 0] = (xind+3)*gap;
            Z0[Ntrg_nodeind * 3 + 1] = (yind+3)*gap;
            Z0[Ntrg_nodeind * 3 + 2] = 0.0;
            Z1[Ntrg_nodeind * 3 + 0] = (xind+3)*gap;
            Z1[Ntrg_nodeind * 3 + 1] = (yind+3)*gap;
            Z1[Ntrg_nodeind * 3 + 2] = 0.9999999999;
        }
    }
    sctl::Vector<Real> UZ0(Z0.Dim());
    LayerPotenOp0.SetTargetCoord(Z0);
    BIO(&UZ0,sigma);
    UZ0 -= eval_rhs(pressure_drop);
    sctl::Vector<Real> UZ1(Z1.Dim());
    LayerPotenOp0.SetTargetCoord(Z1);
    BIO(&UZ1,sigma);
    UZ1 -= eval_rhs(pressure_drop);
    std::cout << "============ Z periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffZ = UZ0-UZ1;
    for (int i=0; i<UdiffZ.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffZ[i*3+0]<< ", " << UdiffZ[i*3+1]<< ", " << UdiffZ[i*3+2]<< ". " << std::endl;
    }
}


int main(int argc, char** argv) {

    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Comm comm = sctl::Comm::World();
        // sctl::Profile::Enable(true);
        long Nelem_ptcl = std::stol(argv[1]); // number of elements
        long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
        int peri_mode = std::stoi(argv[3]); // what kind of periodicity does the system have; peri_mode = j for j-periodic.
        long Nptcl = std::stol(argv[4]); // number of particles inside
        double gmres_tol = std::stod(argv[5]);
        double tol = std::stod(argv[6]);

        if (peri_mode==1) {
            test1peri_channel<Real>(Nelem_ptcl, FourierOrder, comm, Nptcl, gmres_tol, tol);
        } else if (peri_mode == 2) {
            test2peri_plane<Real>(Nelem_ptcl, FourierOrder, comm, gmres_tol, tol);
        } else {
            test3peri<Real>(Nelem_ptcl, FourierOrder, comm, Nptcl, gmres_tol, tol);
        }
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
