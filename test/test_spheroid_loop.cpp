#include "utils.hpp"

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
    const auto BIO = [&LPO,&NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        U->SetZero();
        LPO.ComputePotential(*U, sigma);
        if (U->Dim() == sigma.Dim()) {
            (*U) -= sigma*0.5*NormalOrient; // for double-layer
        }
    };

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



int main(int argc, char** argv) {
    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Comm comm = sctl::Comm::World();
        long Nelem = std::stol(argv[1]); // =0 for trefoil, =1 for 1-particle; =2 for conv div, =3 for trefoil with particle. TODO: =4 for 2peri plane with 1 particle
        long FourierOrder = std::stol(argv[2]);
        double gmres_tol = std::stod(argv[3]);
        double tol = std::stod(argv[4]);

        std::cout << "Running sphere manufactured solutions for Nelem = " << Nelem << ", Fourier = " << FourierOrder << ";" << std::endl;
        function<Real>(Nelem, FourierOrder, 0, gmres_tol, tol, comm);
        std::cout << "Running spheroid manufactured solutions for Nelem = " << Nelem << ", Fourier = " << FourierOrder << ";" << std::endl;
        function<Real>(Nelem, FourierOrder, 1, gmres_tol, tol, comm);
        std::cout << "Running loop manufactured solutions for Nelem = " << Nelem << ", Fourier = " << FourierOrder << ";" << std::endl;
        function<Real>(Nelem, FourierOrder, 3, gmres_tol, tol, comm);

    }
    sctl::Comm::MPI_Finalize();
    return 0;
}