// =============================================================================
// test2_ptcl_conv.cpp
//
// Manufactured-solution accuracy test for a singly-periodic spherical
// suspension. Produces the digits-of-accuracy table (tab. 3) and the (N_p, N_f)
// parameter sweep figure (fig. 7) in the accompanying paper.
//
// Usage:
//   make test_manufactured_soln
//   mpirun -n <Nproc> --map-by numa:pe=$OMP_NUM_THREADS ./bin/test_manufactured_soln \
//          <N_p> <N_f> <Nptcl> <geom_mode> <Ncopy>
//
// Arguments:
//   N_p        panels per particle
//   N_f        azimuthal Fourier modes per particle
//   Nptcl      number of particles
//   geom_mode  particle shape: 0 = sphere, 1 = spheroid, 2 = loop
//   Ncopy      half-width of the lattice truncation for the reference solution
//   Example:   ./bin/test_manufactured_soln 6 32 25 0 60000
//
// Method:
//   A Stokes doublet (a pair of equal-and-opposite Stokeslets) is placed inside
//   each particle, giving zero net force per period as required for absolute
//   convergence of the lattice sum. The reference field u_e is the periodic
//   array of doublets truncated to 2*Ncopy+1 cells, summed directly. Its trace
//   on the particle boundaries is imposed as the Dirichlet data; the combined-
//   field BIE is solved and the solver field is compared with u_e on an interior
//   grid (mean removed) to report the maximum relative error.
//
//   The accuracy of the truncated direct sum can be verified through 
//   exact_field_check().
// =============================================================================

// Boundary integral operators
#include "stokes_bio.hpp" 

// Geometry for tests
#include "utils_geom.hpp"

// Other util functions
#include "utils_tests.cpp" 

// Visualization
#include "utils_vis.hpp" 

// Direct lattice sum at Xtrg of the singly-periodic array of Stokeslets (Xsrc, sigma),
// truncated to 2*Ncopy+1 cells.
template <class Real> sctl::Vector<Real> exact_field(const sctl::Vector<Real>& Xtrg, const sctl::Vector<Real>& Xsrc, const sctl::Vector<Real>& sigma, const sctl::Long Ncopy) {
    sctl::Stokes3D_FxU ker;
  
    const sctl::Long N = Xtrg.Dim()/3;
    sctl::Vector<Real> U(N*3);
    U = 0.;
    PeriodicGeom<Real> obj;

    sctl::Vector<Real> Xsrc_ = obj.X_nbr_copy(Xsrc,Ncopy,1);
    sctl::Vector<Real> sigma_ = obj.vec_nbr_copy(sigma,Ncopy,1);
    ker.Eval(U,Xtrg,Xsrc_,Xsrc_,sigma_);

    return U;
}

// Estimate the truncation error of the reference lattice sum by comparing the Ncopy1
// and Ncopy2 truncations (Ncopy2 > Ncopy1) on a uniform interior grid.
template <class Real> void exact_field_check(
    sctl::Comm comm, 
    sctl::Long Nptcl, 
    const sctl::Long geom_mode, 
    const sctl::Long Ncopy1, 
    const sctl::Long Ncopy2) 
    {
    const sctl::Long Nelem = 4;
    const sctl::Long FourierOrder = 64;
    const sctl::Long ElemOrder = 10;
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
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

    sctl::Long Ncharge;
    if (Nptcl < 150) {
        // two equal and opposite charges per particle
        Ncharge = 2*Nptcl;
    } else {
        // only first 150 particles get charges inside. (arbitrary, to limit true solution timing)
        Ncharge = 2*150;
    }
    // Currently one Stokeslet doublet per particle (for a simple net-force-zero scenario)
    sctl::Vector<Real> Xsrc(Ncharge*3);
    sctl::Vector<Real> Stokeslet_sigma(Ncharge*3);
    srand48(2);
    for (sctl::Long i=0; i<Ncharge/2; i++) {

        const Real disp = 0.2 * ptcls_rs[i];
        const Real disp_y = disp * drand48();
        const Real disp_z = disp * drand48();
        const Real rand_mag = drand48()-0.5;
        Xsrc[i*6+0] = ptcls_Xcs[i*3+0];
        Xsrc[i*6+1] = ptcls_Xcs[i*3+1] + disp_y;
        Xsrc[i*6+2] = ptcls_Xcs[i*3+2] + disp_z;
        Stokeslet_sigma[i*6+0] = 0.; 
        Stokeslet_sigma[i*6+1] = -rand_mag * disp_y; 
        Stokeslet_sigma[i*6+2] = -rand_mag * disp_z; 
        Xsrc[i*6+3] = ptcls_Xcs[i*3+0];
        Xsrc[i*6+4] = ptcls_Xcs[i*3+1] - disp_y;
        Xsrc[i*6+5] = ptcls_Xcs[i*3+2] - disp_z;
        Stokeslet_sigma[i*6+3] = 0.; 
        Stokeslet_sigma[i*6+4] = rand_mag * disp_y; 
        Stokeslet_sigma[i*6+5] = rand_mag * disp_z; 
        
    }

    PeriodicGeom<Real> trg;    
    CubeVolumeVisShifted<Real> vol_vis(10, 0.95, comm);
    X0 = vol_vis.GetCoord();
    
    sctl::Vector<Real> field_on_surf_1 = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy1);
    sctl::Vector<Real> field_on_surf_2 = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy2);
    sctl::Vector<Real> diff = field_on_surf_1 - field_on_surf_2;
    Real max_err = 0.;
    Real max_field1 = 0.;
    for (const auto d : diff) max_err = std::max(std::abs(d), max_err);
    for (const auto d : field_on_surf_1) max_field1 = std::max(std::abs(d), max_field1);
    Real max_rel_err = max_err / max_field1;

    std::cout << "max relative error between Ncopy1 = " << Ncopy1 << " and " << Ncopy2 << " is " << std::setprecision(15) << max_rel_err << std::endl;
}

// Manufactured-solution test: impose the reference lattice-sum field on the particle
// boundaries, solve the exterior Dirichlet problem, and report the interior error.
template <class Real> void manufactured_soln_1peri(
    const sctl::Long Nelem, 
    const sctl::Long FourierOrder, 
    sctl::Comm comm,
    sctl::Long Nptcl, 
    const sctl::Long geom_mode, 
    const sctl::Long Ncopy) 
    {

    std::cout << "Running singly-periodic manufactured solutions test on N_p = " << Nelem << ", N_f = " << FourierOrder << std::endl;

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    Real tol = 1e-14;
    Real gmres_tol = 1e-10; 
    const sctl::Long gmres_max_iter = 200;
    const sctl::Long ElemOrder = 10;
    const Real period_length = 1.;

    if (FourierOrder < 20) {
        gmres_tol = 1e-6;
    } else if (Nelem < 4) {
        gmres_tol = 1e-8;
    } else if (FourierOrder < 36) {
        gmres_tol = 1e-10;
    } else {
        gmres_tol = 1e-12;
    }
    
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

    sctl::Vector<Real> X0; // First solve: target coordinates are surface nodes
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    Real surface_area;
    sctl::Vector<Real> wts;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        // Collect surface area of all geometry across MPI proc
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area = sa_all[0];
    }

    // Create point charges at random locations close to particle center, by a distance of at most 0.2r.
    sctl::Long Ncharge;
    if (Nptcl < 150) {
        // two equal and opposite charges per particle
        Ncharge = 2*Nptcl;
    } else {
        // only first 150 particles get charges inside. (arbitrary, to limit true solution timing)
        Ncharge = 2*150;
    }
    // Place one Stokeslet doublet per particle (for a simple net-force-zero scenario)
    sctl::Vector<Real> Xsrc(Ncharge*3);
    sctl::Vector<Real> Stokeslet_sigma(Ncharge*3);
    srand48(2);
    for (sctl::Long i=0; i<Ncharge/2; i++) {

        const Real disp = 0.2 * ptcls_rs[i];
        const Real disp_y = disp * drand48();
        const Real disp_z = disp * drand48();
        const Real rand_mag = drand48()-0.5;
        Xsrc[i*6+0] = ptcls_Xcs[i*3+0];
        Xsrc[i*6+1] = ptcls_Xcs[i*3+1] + disp_y;
        Xsrc[i*6+2] = ptcls_Xcs[i*3+2] + disp_z;
        Stokeslet_sigma[i*6+0] = 0.; 
        Stokeslet_sigma[i*6+1] = -rand_mag * disp_y; 
        Stokeslet_sigma[i*6+2] = -rand_mag * disp_z; 
        Xsrc[i*6+3] = ptcls_Xcs[i*3+0];
        Xsrc[i*6+4] = ptcls_Xcs[i*3+1] - disp_y;
        Xsrc[i*6+5] = ptcls_Xcs[i*3+2] - disp_z;
        Stokeslet_sigma[i*6+3] = 0.; 
        Stokeslet_sigma[i*6+4] = rand_mag * disp_y; 
        Stokeslet_sigma[i*6+5] = rand_mag * disp_z; 
        
    }
    sctl::Vector<Real> field_on_surf = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy);

    // Get block-preconditioner on particle geometry
    sctl::Matrix<Real> PrecondMat0, PrecondMat1;
    sctl::Long A11size = precond_ptcl(PrecondMat0, PrecondMat1, Nelem, ElemOrder, FourierOrder, SL_scal, DL_scal, comm);
    
    // Apply A11inv to each panel of <vec>, assuming <vec> contains whole particles (i.e. vec.Dim() = A11size * Nptcl)
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
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);

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
        (*U) = AinvApply(Uloc);
    };
    // Left precondition on RHS
    sctl::Vector<Real> A11invF = AinvApply(field_on_surf);

    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> sigma;
    solver(&sigma, BIO_precond, A11invF, gmres_tol, gmres_max_iter, false, nullptr, &krylov_precond);

    // Evaluate solution at target grid in unit cube
    PeriodicGeom<Real> trg;    
    CubeVolumeVisShifted<Real> vol_vis(5, 0.95, comm);
    sctl::Vector<Real> X0_all = vol_vis.GetCoord();
    // Filter out targets inside particles (not in fluid domain)
    sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
    std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
    X0 = std::get<0>(trg_tuple);
    filtered_inds = std::get<1>(trg_tuple);
    // Evaluate periodic solution at new targets
    LayerPotenOp0.SetTargetCoord(X0);
    sctl::Vector<Real> U;
    BIO(&U, sigma);
    // Compare with reference solution
    sctl::Vector<Real> field_on_trg = exact_field(X0, Xsrc, Stokeslet_sigma, Ncopy);
    // get max abs error
    sctl::Vector<Real> err = U - field_on_trg;

    // Subtract mean to remove constant difference
    renormalize_error(err, comm);
    
    double max_err = 0;
    Real max_u = 0.;
    for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
    for (const auto e : field_on_trg) max_u = std::max<Real>(max_u, sctl::fabs(e));

    sctl::Vector<Real> err_loc(1);
    err_loc[0] = max_err;
    sctl::Vector<Real> err_all(1);
    err_all[0] = 0;
    comm.Allreduce((sctl::Iterator<Real>) err_loc.begin(), (sctl::Iterator<Real>) err_all.begin(), 1, sctl::CommOp::MAX);

    sctl::Vector<Real> u_loc(1);
    u_loc[0] = max_u;
    sctl::Vector<Real> u_all(1);
    u_all[0] = 0.;
    comm.Allreduce((sctl::Iterator<Real>) u_loc.begin(), (sctl::Iterator<Real>) u_all.begin(), 1, sctl::CommOp::MAX);

    if (!comm.Rank()) {
        std::cout<<"Max error = "<< std::setprecision(15) << err_all[0] << ", Max u = " << u_all[0] << ", Max relative error = " << err_all[0] / u_all[0] << std::endl;
    }   
}


int main(int argc, char** argv) {
    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        // sctl::Profile::Enable(true);
        sctl::Comm comm = sctl::Comm::World();
        long Nelem_ptcl = std::stol(argv[1]); // N_p
        long FourierOrder = std::stol(argv[2]);  // N_f
        long Nptcl = std::stol(argv[3]); // number of particles in suspension
        long geom_mode = std::stol(argv[4]); // type of particles in suspension: spheres: 0; spheroids: 1; loop: 2.
        long Ncopy = std::stol(argv[5]); // Number of periodic box summed for reference solution

        manufactured_soln_1peri<Real>(Nelem_ptcl, FourierOrder, comm, Nptcl, geom_mode, Ncopy);
        // exact_field_check<Real>(comm, Nptcl, geom_mode, 60000, 70000);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
