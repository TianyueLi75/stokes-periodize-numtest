/*
    Compare the setup and eval time when using the periodized solver (skip precomputation cost)
    to time when using the free-space solver
    to get the periodization overhead.
    Example uses triply-periodized operator.
*/

// Boundary integral operators
#include "stokes_bio.hpp" 

// Geometry for tests
#include "utils_geom.hpp"

// Other util functions
#include "utils_tests.cpp" 

// Visualization
#include "utils_vis.hpp" 

template <class Real> void periodization_time(
    const sctl::Long Nelem, 
    const sctl::Long FourierOrder, 
    sctl::Comm comm, 
    sctl::Long Nptcl, 
    const Real gmres_tol, 
    const Real tol) 
    {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const sctl::Long ElemOrder = 10;
    const sctl::Long geom_mode = 0;

    const Real period_length = 1.;
    const Real pressure_drop = -1.;
    
    const sctl::Long gmres_max_iter = 200;
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs, ptcls_rs, NormalOrient;
    sctl::SlenderElemList<Real> elem_lst0;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else { 
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    }
    Nptcl = ptcls_rs.Dim();

    // Get surface collocation nodes
    sctl::Vector<Real> X0; 
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    const sctl::Vector<Real> X0surf = X0;
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

    // Get block-preconditioner on particle geometry
    sctl::Matrix<Real> PrecondMat0, PrecondMat1;
    sctl::Long A11size = precond_ptcl(PrecondMat0, PrecondMat1, Nelem, ElemOrder, FourierOrder, SL_scal, DL_scal, comm);

    // Apply A11inv to each panel of vec.
    const auto AinvApply = [&PrecondMat0,&PrecondMat1,&A11size, &comm](const sctl::Vector<Real>& vec) {
        sctl::Long N = vec.Dim();
        sctl::Long Nptcl = N / A11size; 
        sctl::Vector<Real> AinvVec(N);
        for (sctl::Long i=0; i<Nptcl; i++) {
            // for each particle, apply A11inv.
            sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
            sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
            for (sctl::Long j=0; j<A11size; j++) {
                AinvVec[i*A11size + j] = AinvVecMat(j,0);
            }
        }
        return AinvVec;
    };

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); 
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0surf);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length);

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
        // LEFT PRECONDITIONER: u -> A11inv*u
        (*U) = AinvApply(Uloc);
    };

    // Right-hand-side for triply-periodic system to handle background pressure drop
    const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
        sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
        AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

        sctl::Vector<Real> U0;
        LayerPotenOp0.ComputeSL(U0, force_density);
        return U0;
    };

    sctl::Vector<Real> A11invF = AinvApply(eval_rhs(pressure_drop)); 

    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    // first gmres to remove timing for matrix loading, and set Krylov preconditioner.
    sctl::Vector<Real> sigma_setup;
    solver(&sigma_setup, BIO_precond, A11invF, 1e-2, gmres_max_iter);
    sctl::Profile::reset();

    LayerPotenOp0.ClearSetup();
    sctl::Profile::Tic("3-periodic LPO Setup");
    LayerPotenOp0.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);

    sctl::Vector<Real> sigma_temp;
    solver(&sigma_temp, BIO_precond, A11invF, 1e-4, gmres_max_iter); // first run to remove high timing counts
    sctl::Profile::reset();
    comm.Barrier();

    sctl::Vector<Real> sigma;
    sctl::Profile::Tic("3-periodic LPO solve with preconds");
    solver(&sigma, BIO_precond, A11invF, gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);
    sctl::Profile::reset();

    PeriodicGeom<Real> trg;    
    sctl::Long Ngrid; 
    // Grids selected to reflect size of system: {Nptcl, Ngrid} = {25, 36}, {50, 46}, {100, 57}, {200, 72}
    if (Nptcl == 25) {
        Ngrid = 36;
    } else if (Nptcl == 50) {
        Ngrid = 46;
    } else if (Nptcl == 100) {
        Ngrid = 57;
    } else if (Nptcl == 200) {
        Ngrid = 72;
    } else {
        // For other configs, just use grid of 80^3
        Ngrid = 80;
    }
    CubeVolumeVisShifted<Real> vol_vis(Ngrid, 0.95, comm);
    sctl::Vector<Real> X0_all = vol_vis.GetCoord();
    sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    X0 = std::get<0>(trg_tuple);
    filtered_inds = std::get<1>(trg_tuple);

    LayerPotenOp0.SetTargetCoord(X0);

    LayerPotenOp0.ClearSetup();
    sctl::Profile::Tic("3-periodic LPO Eval Setup");
    LayerPotenOp0.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);

    // First run to remove extra compiling and setup time
    sctl::Vector<Real> Utemp;
    BIO(&Utemp, sigma);
    Utemp -= eval_rhs(pressure_drop);

    sctl::Profile::reset();
    sctl::Vector<Real> U;
    sctl::Profile::Tic("3-periodic LPO eval");
    for (int loop=0; loop < 10; loop++) {
        BIO(&U, sigma);
        U -= eval_rhs(pressure_drop);  
    }
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);
    sctl::Profile::reset();


    std::cout << "=========== FREE SPACE LPO =============" << std::endl;
    StokesBIO LayerPotenOp1(SL_scal, DL_scal, comm);
    LayerPotenOp1.AddElemList(elem_lst0);
    LayerPotenOp1.SetTargetCoord(X0surf);
    LayerPotenOp1.SetAccuracy(tol);

    const auto BIO_free = [&DL_scal,&LayerPotenOp1,NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        U->SetZero();
        LayerPotenOp1.ComputePotential(*U, sigma);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma*0.5*NormalOrient * DL_scal; // for double-layer
    };

    const auto BIO_precond_free = [&BIO_free,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> Uloc;
        BIO_free(&Uloc,sigma);
        (*U) = AinvApply(Uloc);
    };

    const auto eval_rhs_free = [&LayerPotenOp1,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
        sctl::Vector<Real> force_density(LayerPotenOp1.Dim(0)); force_density = 0;
        AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

        sctl::Vector<Real> U0;
        LayerPotenOp1.ComputeSL(U0, force_density);
        return U0;
    };

    sctl::GMRES<Real> solver_free(comm);
    sctl::KrylovPrecond<Real> krylov_precond_free;
    sctl::Vector<Real> A11invF_free = AinvApply(eval_rhs_free(pressure_drop));

    sctl::Vector<Real> sigma_setup_free;
    solver_free(&sigma_setup_free, BIO_precond_free, A11invF_free, 1e-2, gmres_max_iter);
    sctl::Profile::reset();
    LayerPotenOp1.ClearSetup();
    sctl::Profile::Tic("free space LPO Setup");
    LayerPotenOp1.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);
    sctl::Profile::reset();

    sctl::Vector<Real> sigma_free_temp;
    solver_free(&sigma_free_temp, BIO_precond_free, A11invF_free, 1e-4, gmres_max_iter);
    sctl::Profile::reset();

    sctl::Vector<Real> sigma_free;
    sctl::Profile::Tic("free space LPO solve with preconds");
    solver_free(&sigma_free, BIO_precond_free, A11invF_free, gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond_free);
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);
    sctl::Profile::reset();

    LayerPotenOp1.SetTargetCoord(X0); // Same target grid from periodic case

    LayerPotenOp1.ClearSetup();
    sctl::Profile::Tic("free space LPO Eval Setup");
    LayerPotenOp1.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);

    sctl::Vector<Real> U_free_temp;
    BIO_free(&U_free_temp, sigma_free);
    U_free_temp -= eval_rhs_free(pressure_drop);

    sctl::Profile::reset();
    sctl::Profile::Tic("free space LPO eval");
    sctl::Vector<Real> U_free;
    for (int loop=0; loop<10; loop++) {
        BIO_free(&U_free, sigma_free);
        U_free -= eval_rhs_free(pressure_drop);
    }
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);
    sctl::Profile::reset();

}


int main(int argc, char** argv) {
    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Profile::Enable(true);
        sctl::Comm comm = sctl::Comm::World();
        long Nelem = std::stol(argv[1]); // N_p
        long FourierOrder = std::stol(argv[2]);  // N_f
        long Nptcl = std::stol(argv[3]);
        double gmres_tol = std::stod(argv[4]);
        double tol = std::stod(argv[5]);

        periodization_time<Real>(Nelem, FourierOrder, comm, Nptcl, gmres_tol, tol);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}