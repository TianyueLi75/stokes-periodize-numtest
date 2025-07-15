#include "periodize.hpp"
#include "utils.hpp"
// Test script for calculation and timing of particle-only 1 and 3 periodic problems with background pressure flow. 

/**
 * Background flow with unit pressure gradient along X-axis.
 */
template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
    const Real pdrive = 1;
    const sctl::Long N = X.Dim()/3;
    sctl::Vector<Real> U(N*3);
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        U[i*3+0] = -pdrive * ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
        U[i*3+1] = 0;
        U[i*3+2] = 0;
    }
    return U;
}

template <class Real> void test(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Integer peri_mode, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    // const Real tol = 1e-15;
    const Real tol = 1e-8;
    const Real gmres_tol = 1e-12;
    const sctl::Long ElemOrder = 10;
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 1, peri_mode, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst_nbr = std::get<0>(build_nbr);
        NormalOrient = std::get<1>(build_nbr);
    } else { 
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 0, 1, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_nbr = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 1, peri_mode, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst_nbr = std::get<0>(build_nbr);
        NormalOrient = std::get<1>(build_nbr);
    }
    const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); 
    Nptcl = ptcls_rs.Dim(); 
    std::cout << "periodic mode is " << peri_mode << ", Nrepeat is " << Nrepeat << std::endl;

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Vector<Real> X_proxy;
    if (peri_mode == 1) {
        X_proxy = Periodize1D<Real>::GetProxySurf(); // proxy points coordinates
    } else if (peri_mode == 3) {
        X_proxy = Periodize3D<Real>::GetProxySurf(); // proxy points coordinates
    } else {
        SCTL_ASSERT(false);
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst_nbr);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);

    StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal, comm); // potential from elem_lst0 to proxy points
    LayerPotenOp_proxy.AddElemList(elem_lst0);
    LayerPotenOp_proxy.SetTargetCoord(X_proxy);
    LayerPotenOp_proxy.SetAccuracy(tol);

    // periodized layer potential operator
    const auto BIO = [&DL_scal,&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat,NormalOrient,&peri_mode, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        const sctl::Long N = sigma.Dim();
        // std::cout << "in BIO, dim of sigma is " << N << std::endl;

        sctl::Vector<Real> sigma_nbr(Nrepeat*N); // repeat sigma Nrepeat times
        for (sctl::Long k = 0; k < Nrepeat; k++) {
            for (sctl::Long i = 0; i < N; i++) {
                sigma_nbr[k*N+i] = sigma[i];
            }
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma_nbr);
        if (DL_scal && U->Dim() == N) {
            (*U) -= sigma*0.5*NormalOrient * DL_scal;
        }

        { // Add far-field
            sctl::Vector<Real> U_proxy, U_far;
            LayerPotenOp_proxy.ComputePotential(U_proxy, sigma);
            if (peri_mode==1) {
                // 1-periodic
                Periodize1D<Real>::EvalFarField(U_far, X0, U_proxy);
            } else if (peri_mode==3) {
                // 3-periodic
                Periodize3D<Real>::EvalFarField(U_far, X0, U_proxy);
            } else {
                std::cout << "2-periodic not yet implemented." << std::endl;
                SCTL_ASSERT(false);
            }
            
            (*U) += U_far;
        } 
        // comm.Barrier();
    };

    // first gmres to remove timing for matrix loading.
    sctl::Vector<Real> sigma_temp;
    sctl::GMRES<Real> solver(comm);
    // sctl::KrylovPrecond<Real> krylov_precond;
    // solver(&sigma_temp, BIO, -bg_flow(X0), 1e-2, -1, false, nullptr, &krylov_precond);
    solver(&sigma_temp, BIO, -bg_flow(X0), 1e0);
    sctl::Profile::reset();

    LayerPotenOp0.ClearSetup();
    sctl::Profile::Tic("Setup");
    LayerPotenOp0.Setup();
    sctl::Profile::Toc();
    sctl::Profile::print(&comm);

    sctl::Vector<Real> sigma;
    sctl::Profile::Tic("Solver");
    solver(&sigma, BIO, -bg_flow(X0), gmres_tol);
    // solver(&sigma, BIO, -bg_flow(X0), gmres_tol, -1, false, nullptr, &krylov_precond);
    sctl::Profile::Toc();
    sctl::Profile::print(&comm, {"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    sctl::Profile::reset();

}


int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    sctl::Profile::Enable(true);
    sctl::Comm comm = sctl::Comm::World();
    long Nelem_ptcl = std::stol(argv[1]); // number of elements
    long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
    int peri_mode = std::stoi(argv[3]); // what kind of periodicity does the system have; peri_mode = j for j-periodic.
    long Nptcl = std::stol(argv[4]); // number of particles inside
    long geom_mode = std::stol(argv[5]); // =0: spheres; =1: spheroids; =3: bacteria; =4: loop.

    test<Real>(Nelem_ptcl, FourierOrder, peri_mode, comm, Nptcl, geom_mode);
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}
