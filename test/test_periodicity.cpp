/*
    Tests scripts in all three periodicity on the solution's periodicity.
*/

// Boundary integral operators and their periodization
#include "stokes_bio.hpp" 
#include "bio_operator.hpp" 

// Geometry for tests
#include "planeNaive.hpp"
#include "utils_geom.hpp"

// Other util functions
#include "utils_tests.cpp" 

// Visualization
#include "utils_vis.hpp" 


template <class Real> void test1peri_channel(
    const sctl::Long Nelem_channel, 
    const sctl::Long FourierOrder, 
    sctl::Comm comm, 
    sctl::Long Nptcl, 
    const Real gmres_tol, 
    const Real tol) 
    {

    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;
    const sctl::Long ElemOrder = 10;
    const Real pressure_drop = -1.0;
    const Real period_length = 1;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Long ptcl_ord = 4;
    if (Nptcl>0) {
        ptcls.ReInit(Nptcl);
        ptcls = ptcl_ord;
    }
    sctl::Vector<Real> ptcls_Xcs, ptcls_rs, NormalOrient, ptcls_thetas, ptcls_phis;
    sctl::SlenderElemList<Real> elem_lst0;
    // std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>,sctl::Vector<Real>,sctl::Vector<Real>> build0 = obj.build_conv_div_sph(Nelem_channel, ElemOrder, FourierOrder, 0.1, 0.2, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord);
    // TODO: try with easy set up to check periodicity.
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_straight(Nelem_channel, ElemOrder, FourierOrder, 0.3, comm, ptcls, ptcls_rs, ptcls_Xcs, 0);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);

    sctl::Vector<Real> X0; 
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    // Real surface_area;
    // sctl::Vector<Real> wts;
    // { // get wts and surface area
    //     sctl::Vector<Real> X, Xn, dist_far, surface_area_;
    //     sctl::Vector<sctl::Long> element_wise_node_cnt;
    //     elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
    //     SurfaceIntegral(surface_area_, wts*0+1, wts);
    //     // surface_area = surface_area_[0];
    //     sctl::Vector<Real> sa_loc(1);
    //     sa_loc[0] = surface_area_[0];
    //     sctl::Vector<Real> sa_all(1);
    //     sa_all[0] = 0;
    //     comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
    //     surface_area = sa_all[0];
    // }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); 
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    MeanCorrectedStokesBIOOperator<Real, sctl::SlenderElemList<Real>> BIO(LayerPotenOp0, NormalOrient, DL_scal, comm);
    BIO.AddSurface(elem_lst0);

    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;

    sctl::Vector<Real> sigma;
    solver(&sigma, BIO, bg_flow_1peri(X0) * (pressure_drop/period_length), gmres_tol, -1, false, nullptr, &krylov_precond);

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
            X1[Ntrg_nodeind * 3 + 0] = 0.99999999;
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
    UX0 -= bg_flow_1peri(X0) * (pressure_drop/period_length);
    UX1 -= bg_flow_1peri(X1) * (pressure_drop/period_length);
    std::cout << "============ X periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffX = UX0-UX1;
    for (int i=0; i<UdiffX.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffX[i*3+0]<< ", " << UdiffX[i*3+1]<< ", " << UdiffX[i*3+2]<< ". " << std::endl;
    }
}

// Runs periodicity test on sets of 3 particles (spheres, spheroids, or loops).
template <class Real> void test2peri_plane(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Comm comm, const Real gmres_tol, const Real tol) {

    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;
    const Real pressure_drop = -10.0;
    const Real period_length = 1;
    const sctl::Long ElemOrder = 10;
    const sctl::Long geom_mode = 1; 

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs, ptcls_rs, NormalOrient, ptcls_thetas, ptcls_phis;
    sctl::SlenderElemList<Real> elem_lst0;
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
    
    // Plane
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
    
    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&wts_wall,&surface_area_wall,&plane,&LayerPotenOp0,&DL_scal,&X0_ptcl,&X0_wall,&NormalOrient, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        
        sctl::Vector<Real> ptcl_dens(X0_ptcl.Dim(), (sctl::Iterator<Real>) sigma.begin(), true);
        sctl::Vector<Real> wall_dens(X0_wall.Dim(), (sctl::Iterator<Real>) sigma.begin()+ptcl_dens.Dim(), true);
            
        sctl::Vector<Real> sigma_mean, sigma0;

        { 
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
    // MeanCorrectedStokesBIOOperator<Real, sctl::SlenderElemList> BIO(LayerPotenOp0, NormalOrient, DL_scal, comm);
    // BIO.AddSurface(elem_lst0);

    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    solver(&sigma,BIO, bg_flow_2peri(X0) * (pressure_drop/period_length), gmres_tol);

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
            X1[Ntrg_nodeind * 3 + 0] = 0.99999999;
            X1[Ntrg_nodeind * 3 + 1] = (yind+3)*gap;
            X1[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
            
        }
    }
    sctl::Vector<Real> UX0(X0.Dim());
    LayerPotenOp0.SetTargetCoord(X0);
    BIO(&UX0,sigma);
    sctl::Vector<Real> UX1(X1.Dim());
    LayerPotenOp0.SetTargetCoord(X1);
    BIO(&UX1,sigma);
    UX0 -= bg_flow_2peri(X0) * (pressure_drop/period_length);
    UX1 -= bg_flow_2peri(X1) * (pressure_drop/period_length);
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
            Y1[Ntrg_nodeind * 3 + 1] = 0.99999999;
            Y1[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
        }
    }
    sctl::Vector<Real> UY0(Y0.Dim());
    LayerPotenOp0.SetTargetCoord(Y0);
    BIO(&UY0, sigma); 
    sctl::Vector<Real> UY1(Y1.Dim());
    LayerPotenOp0.SetTargetCoord(Y1);
    BIO(&UY1, sigma); 
    std::cout << "============ Y periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffY = UY0-UY1;
    for (int i=0; i<UdiffY.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffY[i*3+0]<< ", " << UdiffY[i*3+1]<< ", " << UdiffY[i*3+2]<< ". " << std::endl;
    }    
}

template <class Real> void test3peri(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Comm comm, sctl::Long Nptcl, const Real gmres_tol, const Real tol) {

    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real pressure_drop = -1.0;
    const Real period_length = 1;

    const sctl::Long ElemOrder = 10;
    const sctl::Long geom_mode = 0;
    
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
        // MPI
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area = sa_all[0];
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0,"1");
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length);

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    MeanCorrectedStokesBIOOperator<Real, sctl::SlenderElemList<Real>> BIO(LayerPotenOp0, NormalOrient, DL_scal, comm);
    BIO.AddSurface(elem_lst0);

    const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
        sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
        AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

        sctl::Vector<Real> U0;
        LayerPotenOp0.ComputeSL(U0, force_density);
        return U0;
    };

    // first gmres to remove timing for matrix loading, and set Krylov preconditioner.
    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    LayerPotenOp0.SetTargetCoord(X0);
    solver(&sigma, BIO, eval_rhs(pressure_drop), gmres_tol);

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
