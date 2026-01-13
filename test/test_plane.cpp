#include "periodize.hpp"
#include "utils.hpp"
#include "planeNaive.hpp"

template <class Real> sctl::Vector<Real> bg_pres_flow(const sctl::Vector<Real>& X) {
    const Real dpdx = -1.;
    const sctl::Long N = X.Dim()/3;
    sctl::Vector<Real> U(N*3);
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        U[i*3+0] = dpdx * 0.5 * ((x[2]-0.5)*(x[2]-0.5)); // 2-periodic flow between plates.
        U[i*3+1] = 0.;
        U[i*3+2] = 0.;
    }
    return U;
}

// Uniform background flow in x direction.
template <class Real> sctl::Vector<Real> bg_unif_flow(const sctl::Vector<Real>& X) {
    sctl::Vector<Real> U = X;
    const sctl::Long N = X.Dim() /3;
    // U = 1.; // background flow diagonal to avoid planes of unaffected flows between periods.
    for (sctl::Long i = 0; i < N; i++) {
        U[i*3+0] = 1.; 
        U[i*3+1] = 0.5; 
        U[i*3+2] = 0.; 
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


template <class Real> void test(sctl::Long gl_order, sctl::Long Nelem_xy, Real z_offset, Real trg_size, sctl::Comm comm, const Real gmres_tol, const Real tol) {

    // Main operator being DL, no extra handling for self-eval
    const Real SL_scal = 0.0;
    const Real DL_scal = 1.0;

    sctl::Long Nelem_x = Nelem_xy;
    sctl::Long Nelem_y = Nelem_xy;
    sctl::PlaneIntegral<Real> plane(gl_order, Nelem_x, Nelem_y, z_offset);
    sctl::Vector<Real> X0_src,Xn_src;
    plane.GetNodeCoord(&X0_src, &Xn_src, nullptr);
    sctl::Vector<Real> X0 = X0_src;

    StokesBIO LayerPotenOp1(SL_scal, DL_scal, comm); 
    LayerPotenOp1.AddElemList(plane);
    LayerPotenOp1.SetTargetCoord(X0);
    LayerPotenOp1.SetAccuracy(tol);
    LayerPotenOp1.SetPeriodicity(sctl::Periodicity::XY, 1.0);
    sctl::Vector<Real> NormalOrient(X0.Dim());
    NormalOrient = -1.; 

    // Added SL operator -- need singularity subtraction for self-eval.
    StokesBIO LayerPotenOp2(1.0, 0.0, comm);
    LayerPotenOp2.AddElemList(plane);
    LayerPotenOp2.SetAccuracy(tol);
    LayerPotenOp2.SetPeriodicity(sctl::Periodicity::XY, 1.0);

    Real surface_area_wall;
    sctl::Vector<Real> wts_wall;
    { // get wts and surface area
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

    // periodized layer potential operator
    const auto BIO = [&wts_wall, &surface_area_wall, &plane, &comm, DL_scal,&LayerPotenOp1,NormalOrient,&LayerPotenOp2,&X0_src](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> sigma_mean, sigma0;
        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_mean_wall;
            sctl::Vector<Real> sigma_wall_ = sigma;
            SurfaceIntegral(sigma_mean_wall, sigma_wall_, wts_wall);
            // MPI
            sctl::Vector<Real> sa_loc(3);
            for (int i=0; i<3; i++) {
                sa_loc[i] = sigma_mean_wall[i];
            }
            sctl::Vector<Real> sa_all(3);
            sa_all = 0.;
            for (int i=0; i<3; i++) {
                comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+i, (sctl::Iterator<Real>) sa_all.begin()+i, 1, sctl::CommOp::SUM);
            }
            sigma_mean = sa_all / surface_area_wall;
            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);
        }
        
        U->SetZero();
        LayerPotenOp1.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer
        // DEBUG SL
        // U->SetZero();

        // single layer singularity subtraction
        {
            sctl::Long Nentries_plane = X0_src.Dim()/2; 
            sctl::Vector<Real> U2(3); 
            // Top plane first.
            for (int i=0; i<X0_src.Dim()/6; i++) {
                // Get sigma, x for i-th node
                sctl::Vector<Real> sigma_i(3, (sctl::Iterator<Real>) sigma0.begin() + i*3, true);
                sctl::Vector<Real> x_i(3, (sctl::Iterator<Real>) X0_src.begin() + i*3, true);
                // Change sigma_top_0 to sigma_top_0 - sigma_i
                sctl::Vector<Real> sigma_topsub = sigma0;
                sctl::Vector<Real> sigma_top(Nentries_plane, (sctl::Iterator<Real>) sigma_topsub.begin(), false); 
                AddConstVec(sigma_top, -sigma_i);
                // Evaluate SL1[sigma_top - sigma_i](x_i) + SL2[sigma_bot](x_i), subtracted (sigma_i) * (int G dS) = 0
                LayerPotenOp2.SetTargetCoord(x_i);
                LayerPotenOp2.ComputePotential(U2, sigma_topsub); 
                // Add SL vel to DL vel
                for (int k=0; k<3; k++) {
                    (*U)[i*3+k] += U2[k];
                }
            }
            // Bottom plane
            for (int i=0; i<X0_src.Dim()/6; i++) {
                sctl::Vector<Real> sigma_i(3, (sctl::Iterator<Real>) sigma0.begin() + Nentries_plane + i*3, true);
                sctl::Vector<Real> x_i(3, (sctl::Iterator<Real>) X0_src.begin() + Nentries_plane + i*3, true);
                // Change sigma_top_0 to sigma_top_0 - sigma_i
                sctl::Vector<Real> sigma_botsub = sigma0;
                sctl::Vector<Real> sigma_bot(Nentries_plane, (sctl::Iterator<Real>) sigma_botsub.begin() + Nentries_plane, false); 
                AddConstVec(sigma_bot, -sigma_i);
                LayerPotenOp2.SetTargetCoord(x_i);
                LayerPotenOp2.ComputePotential(U2, sigma_botsub); 
                for (int k=0; k<3; k++) {
                    (*U)[i*3+k+Nentries_plane] += U2[k];
                }
            }
        }
        
        AddConstVec(*U, sigma_mean);
    };

    // Eval: all far.
    const auto BIO_eval = [&wts_wall, &surface_area_wall, &comm, &LayerPotenOp1, &LayerPotenOp2](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {        
        sctl::Vector<Real> sigma_mean, sigma0;
        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_mean_wall;
            sctl::Vector<Real> sigma_wall_ = sigma;
            SurfaceIntegral(sigma_mean_wall, sigma_wall_, wts_wall);
            // MPI
            sctl::Vector<Real> sa_loc(3);
            for (int i=0; i<3; i++) {
                sa_loc[i] = sigma_mean_wall[i];
            }
            sctl::Vector<Real> sa_all(3);
            sa_all = 0.;
            for (int i=0; i<3; i++) {
                comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+i, (sctl::Iterator<Real>) sa_all.begin()+i, 1, sctl::CommOp::SUM);
            }
            sigma_mean = sa_all / surface_area_wall;
            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);
        }

        // std::cout << "DEBUG: sigma mean is " << sigma_mean[0] << ", " << sigma_mean[1] << ", " << sigma_mean[2] << std::endl;
        
        U->SetZero();
        LayerPotenOp1.ComputePotential(*U, sigma0); 
        // U->SetZero(); // DEBUG SL
        sctl::Vector<Real> U2(U->Dim());
        LayerPotenOp2.ComputePotential(U2, sigma0);
        (*U) += U2;

        AddConstVec(*U, sigma_mean);
        
    };

    const auto bg_shear_flow = [](const sctl::Vector<Real>& Xtrg) {
        sctl::Vector<Real> Utrg(Xtrg.Dim());
        Utrg.SetZero();
        sctl::Long Ntrg = Xtrg.Dim()/3;
        for (sctl::Long i=0; i<Ntrg/2; i++) {
            Utrg[i*3+0] = 1.;
            Utrg[i*3+1] = 0.5;
            // shear flow on top plane: [1,0.5,0]; bottom plate fixed.
        }
        return Utrg;
    };

    // sctl::GMRES<Real> solver(comm);
    // sctl::Vector<Real> sigma;
    // sctl::Vector<Real> Utrg = bg_shear_flow(X0);
    // solver(&sigma,BIO, Utrg, gmres_tol);
    // plane.WriteVTK("vis/plane_density", sigma, comm);

    // Getting singular values
    sctl::Long Nsrc = X0_src.Dim()/3;
    sctl::Vector<Real> sigma_eye(Nsrc*3);
    sctl::Vector<sctl::Vector<Real>> LPOvecvec(Nsrc*3);
    for (sctl::Long i=0; i<Nsrc; i++) {
        for (sctl::Long k=0; k<3; k++) {
            sigma_eye.SetZero();
            sigma_eye[i*3+k] = 1.;
            std::cout << "Node number is = " << i << ", dimension = " << k << std::endl;
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
    std::cout << "debug by printing matrix singular values:" << std::endl;
    std::cout << "shape of S_p is " << S_p.Dim(0) << ", " << S_p.Dim(1) << std::endl;
    for (long i=0; i<S_p.Dim(0); i++) {
        std::cout << S_p(i,i) << std::endl;
    }
    
    /*
    // Eval and vis at traget box
    {
        PeriodicGeom<Real> trg;    
        CubeVolumeVisShifted<Real> vol_vis(10, trg_size, comm);
        X0 = vol_vis.GetCoord();
        LayerPotenOp1.SetTargetCoord(X0);
        LayerPotenOp2.SetTargetCoord(X0);
        sctl::Vector<Real> U(X0.Dim());
        BIO_eval(&U, sigma);
        vol_vis.WriteVTK("vis/plane_nobg", U); 
        // U += bg_pres_flow(X0);
        // vol_vis.WriteVTK("vis/planeSL_withbg", U); 
    }
    */      

    /*
    // difference at x=1 vs 0, y=1 vs 0 check periodicity.
    {
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
                X1[Ntrg_nodeind * 3 + 0] = 1.;
                X1[Ntrg_nodeind * 3 + 1] = (yind+3)*gap;
                X1[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
                
            }
        }
        sctl::Vector<Real> UX0(X0.Dim());
        LayerPotenOp1.SetTargetCoord(X0);
        LayerPotenOp2.SetTargetCoord(X0);
        BIO_eval(&UX0,sigma);
        sctl::Vector<Real> UX1(X1.Dim());
        LayerPotenOp1.SetTargetCoord(X1);
        LayerPotenOp2.SetTargetCoord(X1);
        BIO_eval(&UX1,sigma);
        // UX0 += bg_pres_flow(X0);
        // UX1 += bg_pres_flow(X1);
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
                Y1[Ntrg_nodeind * 3 + 1] = 1.;
                Y1[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
            }
        }
        sctl::Vector<Real> UY0(Y0.Dim());
        LayerPotenOp1.SetTargetCoord(Y0);
        LayerPotenOp2.SetTargetCoord(Y0);
        BIO_eval(&UY0, sigma); 
        sctl::Vector<Real> UY1(Y1.Dim());
        LayerPotenOp1.SetTargetCoord(Y1);
        LayerPotenOp2.SetTargetCoord(Y1);
        BIO_eval(&UY1, sigma); 
        std::cout << "============ Y periodicity =================" << std::endl;
        sctl::Vector<Real> UdiffY = UY0-UY1;
        for (int i=0; i<UdiffY.Dim()/3; i++) {
            std::cout << std::setprecision(8) << UdiffY[i*3+0]<< ", " << UdiffY[i*3+1]<< ", " << UdiffY[i*3+2]<< ". " << std::endl;
        }
    } 
    */
}


int main(int argc, char** argv) {

    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Comm comm = sctl::Comm::World();
        // sctl::Profile::Enable(true);
        long gl_order = std::stol(argv[1]); // number of elements
        long Nelem_xy = std::stol(argv[2]);  // number of Fourier nodes
        double z_offset = std::stod(argv[3]);
        double trg_size = std::stod(argv[4]);
        double gmres_tol = std::stod(argv[5]);
        double tol = std::stod(argv[6]);

        test<Real>(gl_order, Nelem_xy, z_offset, trg_size, comm, gmres_tol, tol);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
