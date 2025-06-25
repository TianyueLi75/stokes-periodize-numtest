#include "periodize.hpp"
#include "utils.hpp"

/**
 * Background flow with unit pressure gradient along X-axis.
 */
template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
  const sctl::Long N = X.Dim()/3;
  sctl::Vector<Real> U(N*3);
  for (sctl::Long i = 0; i < N; i++) {
    const auto x = X.begin() + i*3;
    U[i*3+0] = -((x[1])*(x[1]) + (x[2])*(x[2]))/4;
    U[i*3+1] = 0;
    U[i*3+2] = 0;
  }
  return U;
}

template <class Real> void test(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Long Npairs, sctl::Integer peri_mode, sctl::Comm comm) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real tol = 1e-15;
    const Real gmres_tol = 1e-10;
    const sctl::Long ElemOrder = 10;
    const sctl::Long Nrepeat = Npairs*2+1;
    // const Real dist = 0.4; // arbitrary distance between equator of spheres, distance between centers is dx = 2*r + dist
    const Real radi = 0.3;
    const sctl::Long Nptcls = (peri_mode == 1) ? Nrepeat : Nrepeat*Nrepeat*Nrepeat;
    std::cout << "Nptcl = " << Nptcls << std::endl;

    sctl::Vector<sctl::Long> ptcls(Nptcls);
    ptcls = Nelem;
    sctl::Vector<Real> ptcls_rs(Nptcls);
    ptcls_rs = radi;
    
    PeriodicGeom<Real> obj;
    sctl::Vector<Real> center_Xc(3); // center of center particle
    center_Xc = 0.; // centered at the origin.
    sctl::Vector<Real> ptcls_Xcs = obj.X_nbr_copy(center_Xc, Npairs, peri_mode); // automatically has dist = 1.
    sctl::SlenderElemList<Real> elem_lst = obj.free_ptcls(ElemOrder,FourierOrder,comm,ptcls,ptcls_rs,ptcls_Xcs);

    std::string filename_vis = "vis/free_space_"+std::to_string(Nptcls)+"_spheres";
    sctl::Vector<Real> X0, Xnr; // target coordinates
    elem_lst.GetNodeCoord(&X0, &Xnr, nullptr);
    elem_lst.WriteVTK(filename_vis,Xnr,comm);

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);

    // std::cout << "Rank " << comm.Rank() << " dim of Op is " << LayerPotenOp0.Dim(1) << ", " << LayerPotenOp0.Dim(0) << std::endl;

    // layer potential operator -- no periodicity assumed.
    const auto BIO = [&LayerPotenOp0](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        // std::cout << "in BIO, dim sigma is" << sigma.Dim() <<", dim of U is " << U->Dim() << std::endl;
        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma);
        (*U) += sigma*0.5;
    };

    // Solve for sigma to satisfy no-slip boundary conditions: BIO(sigma) + bg_flow = 0
    sctl::Vector<Real> sigma;
    sctl::GMRES<Real> solver(comm);
    solver(&sigma, BIO, -bg_flow(X0), gmres_tol);

    { // Evaluate in interior, and write visualization
        filename_vis = filename_vis + "_U";

        // only look at cube containing 1 pair of copies.
        const Real cube_len = 0.5;
        sctl::CubeVolumeVis<Real> cube(50, cube_len, comm); 
        
        // X0 = cube.GetCoord();
        PeriodicGeom<Real> trg;
        sctl::Vector<Real> X0_all = cube.GetCoord();
        // std::cout << "X0all dim = " << X0_all.Dim() << std::endl;
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);
        // std::cout << "X0 dim = " << X0.Dim() << std::endl;

        LayerPotenOp0.SetTargetCoord(X0);// Set target coordinates for the operator
        sctl::Vector<Real> U;                            // Vector for storing the velocity field
        LayerPotenOp0.ComputePotential(U, sigma);     // Evaluate the velocity field
        U += bg_flow(X0);

        // cube.WriteVTK(filename_vis, U);    // Write to a VTK file
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
        cube.WriteVTK(filename_vis, U_vis);
    }
}


int main(int argc, char** argv) {
    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Comm comm = sctl::Comm::World();
        long Nelem = std::stol(argv[1]); // number of elements
        long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
        long Npairs = std::stol(argv[3]); // # of particle repeats in all direction.
        long peri_mode = std::stol(argv[4]); // 1- or 3- periodic geometry
        test<Real>(Nelem, FourierOrder, Npairs, peri_mode, comm);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
