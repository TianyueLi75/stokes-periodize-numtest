// Straight channel with a sphere inside, 3-periodic and 1-periodic formulations.
// To test 3-periodic new errors in periodicity after most recent update

#include "periodize.hpp"
#include "utils.hpp"
#include "bio_operator.hpp"

/**
 * Background flow with unit pressure gradient along X-axis.
 */
template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
    const sctl::Long N = X.Dim()/3;
    sctl::Vector<Real> U(N*3);
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        U[i*3+0] = -((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
        U[i*3+1] = 0;
        U[i*3+2] = 0;
    }
    return U;
}

/**
 * Build element list for a straight channel with a sphere inside.
 */
template <class Real> sctl::SlenderElemList<Real> build_elem_lst(const sctl::Long Nelem_channel, const sctl::Long ElemOrder, const sctl::Long FourierOrder, sctl::Vector<Real>* NormalOrient_ptr = nullptr) {
    sctl::Vector<Real> Xc, eps, orient;
    sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
    for (sctl::Long i = 0; i < Nelem_channel; i++) {
        ElemOrderVec.PushBack(ElemOrder);
        FourierOrderVec.PushBack(FourierOrder);
        const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
        for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
            const Real x = (i+nodes[j])/Nelem_channel;
            Xc.PushBack(x);
            Xc.PushBack(0.5);
            Xc.PushBack(0.5);
            eps.PushBack(0.2);

            orient.PushBack(0);
            orient.PushBack(0);
            orient.PushBack(1);
        }
    }

    const sctl::Long Nelem_sphere = 4;
    for (sctl::Long i = 0; i < Nelem_sphere; i++) { // add a sphere
        ElemOrderVec.PushBack(ElemOrder);
        FourierOrderVec.PushBack(FourierOrder);
        const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
        for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
            const Real r = 0.1;
            const Real theta = sctl::const_pi<Real>() * (i+nodes[j])/Nelem_sphere;
            Xc.PushBack(0.5+r*sctl::cos<Real>(theta));
            Xc.PushBack(0.5);
            Xc.PushBack(0.5);
            eps.PushBack(r*sctl::sin<Real>(theta));

            orient.PushBack(0);
            orient.PushBack(0);
            orient.PushBack(1);
        }
    }

    sctl::SlenderElemList<Real> elem_lst0(ElemOrderVec, FourierOrderVec, Xc, eps, orient);

    if (NormalOrient_ptr != nullptr) {
        NormalOrient_ptr->ReInit(0);
        constexpr sctl::Integer COORD_DIM = 3;
        sctl::Vector<sctl::Long> elem_wise_node_cnt;
        elem_lst0.GetNodeCoord(nullptr, nullptr, &elem_wise_node_cnt);
        for (sctl::Long i = 0; i < elem_wise_node_cnt.Dim(); i++) {
            for (sctl::Long j = 0; j < elem_wise_node_cnt[i]*COORD_DIM; j++) {
                NormalOrient_ptr->PushBack(i < Nelem_channel ? -1 : 1);
            }
        }
    }

    return elem_lst0;
}

template <class Real> void test(const sctl::Comm& comm) {
    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real pressure_drop = -1.0;
    const Real period_length = 1;

    const Real tol = 1e-8;
    const Real gmres_tol = 1e-8;
    const sctl::Long gmres_max_iter = 200;
    const sctl::Long Nelem_channel = 4;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 28;

    sctl::Vector<Real> NormalOrient; // normal orientation (1 if normal into fluid, else -1)
    const auto elem_lst0 = build_elem_lst(Nelem_channel, ElemOrder, FourierOrder, &NormalOrient); // geometry in the unit box [0,1]^3

    Real surface_area;
    sctl::Vector<Real> X0, wts;
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Vector<Real> X0surf = X0;
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


    StokesBIO<Real> LayerPotenOp0(SL_scal, DL_scal, comm);
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetAccuracy(tol);

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    MeanCorrectedStokesBIOOperator<Real> BIO(LayerPotenOp0, NormalOrient, DL_scal, comm, +1.0);
    BIO.AddSurface(elem_lst0, wts, surface_area);

    // Old lambda-based BIO block retained for reference:
    // const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,&NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    //     sctl::Vector<Real> sigma_mean, sigma0;
    //     { // compute sigma_mean and sigma0 = sigma - sigma_mean
    //         sctl::Vector<Real> sigma_;
    //         elem_lst0.GetFarFieldDensity(sigma_, sigma);
    //         SurfaceIntegral(sigma_mean, sigma_, wts);
    //         //MPI
    //         sctl::Vector<Real> sa_loc = sigma_mean;
    //         sctl::Vector<Real> sa_all(3);
    //         sa_all = 0;
    //         comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
    //         comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+1, (sctl::Iterator<Real>) sa_all.begin()+1, 1, sctl::CommOp::SUM);
    //         comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+2, (sctl::Iterator<Real>) sa_all.begin()+2, 1, sctl::CommOp::SUM);
    //         sigma_mean = sa_all;
    //         sigma_mean *= (1/surface_area);
    //
    //         sigma0 = sigma;
    //         AddConstVec(sigma0, -sigma_mean);
    //
    //         // DEBUG: check that sigma-sigma_mean has surface integral = 0:
    //         sctl::Vector<Real> sigma1 = sigma_;
    //         AddConstVec(sigma1, -sigma_mean);
    //         sctl::Vector<Real> sigma_test_;
    //         SurfaceIntegral(sigma_test_, sigma1, wts);
    //         std::cout << "Surface integral of sigma - sigma bar = " << sigma_test_[0] << ", "  << sigma_test_[1] << ", "  << sigma_test_[2] << std::endl;
    //     }
    //
    //     U->SetZero();
    //     LayerPotenOp0.ComputePotential(*U, sigma0);
    //     if (DL_scal && U->Dim() == sigma.Dim()) (*U) += sigma0*0.5*NormalOrient * DL_scal; // for double-layer
    //
    //     AddConstVec(*U, sigma_mean);
    // };

    const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length,wts,&elem_lst0](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
        sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
        AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

        sctl::Vector<Real> U0;
        LayerPotenOp0.ComputeSL(U0, force_density);
        return U0;
    };
  
    sctl::Vector<Real> Ubg3p, Up3p;
    { // three-periodic
        LayerPotenOp0.SetTargetCoord(X0surf);
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length);

        sctl::Vector<Real> sigma;
        sctl::GMRES<Real> solver(comm);
        solver(&sigma, BIO, eval_rhs(pressure_drop), gmres_tol, gmres_max_iter);
        // elem_lst0.WriteVTK("vis/Up-3p", eval_rhs(pressure_drop), comm);
        // elem_lst0.WriteVTK("vis/sigma-3p", sigma, comm);
        // std::cout << "Eval rhs in 3p" << std::endl;
        // sctl::Vector<Real> rhs_test =  eval_rhs(pressure_drop);
        // for (int i=0; i<rhs_test.Dim()/3; i++) {
        //     std::cout << rhs_test[i*3+0] << ", " << rhs_test[i*3+1] << ", " << rhs_test[i*3+2] << std::endl;
        // }

        LayerPotenOp0.SetTargetCoord(X0surf);
        sctl::Vector<Real> sigma_ubg;
        solver(&sigma_ubg, BIO, bg_flow(X0surf)*(pressure_drop/period_length), gmres_tol, gmres_max_iter);
        // elem_lst0.WriteVTK("vis/Ubg-3p", bg_flow(X0surf) * (pressure_drop/period_length), comm);

    
        { // Evaluate in interior, and write visualization
            VolumeVis<Real> vol_vis(elem_lst0, comm);
            const auto& X0 = vol_vis.GetCoord(); // set new target coordinates
            LayerPotenOp0.SetTargetCoord(X0);
            // sctl::Vector<Real> U;
            BIO(&Ubg3p, sigma_ubg);
            // vol_vis.WriteVTK("vis/U-3p-Ubg-b4", Ubg3p);
            Ubg3p -= bg_flow(X0) * (pressure_drop/period_length);
            // vol_vis.WriteVTK("vis/U-3p-Ubg", Ubg3p);
        }

        { // Evaluate in interior, and write visualization
            VolumeVis<Real> vol_vis(elem_lst0, comm);
            const auto& X0 = vol_vis.GetCoord(); // set new target coordinates
            LayerPotenOp0.SetTargetCoord(X0);
            // sctl::Vector<Real> U;
            BIO(&Up3p, sigma);
            // vol_vis.WriteVTK("vis/U-3p-b4", Up3p);
            Up3p -= eval_rhs(pressure_drop);
            // vol_vis.WriteVTK("vis/U-3p", Up3p);
        }

        // Compare diff magnitude
        std::cout << "3 peri Magnitude of difference between U from background flow and U from force densities: " << std::endl;
        sctl::Vector<Real> Udiff = Ubg3p - Up3p;
        // for (int i=0; i<Up3p.Dim()/3; i++) {
        for (int i=0; i<5; i++) {
            std::cout << std::setprecision(10)<< sctl::sqrt<Real>(Udiff[i*3]*Udiff[i*3] + Udiff[i*3+1]*Udiff[i*3+1] + Udiff[i*3+2]*Udiff[i*3+2]) << std::endl;
        }   
    }

    sctl::Vector<Real> Ubg2p, Up2p;
    { // three-periodic
        LayerPotenOp0.SetTargetCoord(X0surf);
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XY, period_length);

        sctl::Vector<Real> sigma;
        sctl::GMRES<Real> solver(comm);
        solver(&sigma, BIO, eval_rhs(pressure_drop), gmres_tol, gmres_max_iter);
        // elem_lst0.WriteVTK("vis/Up-3p", eval_rhs(pressure_drop), comm);
        // elem_lst0.WriteVTK("vis/sigma-3p", sigma, comm);
        // std::cout << "Eval rhs in 3p" << std::endl;
        // sctl::Vector<Real> rhs_test =  eval_rhs(pressure_drop);
        // for (int i=0; i<rhs_test.Dim()/3; i++) {
        //     std::cout << rhs_test[i*3+0] << ", " << rhs_test[i*3+1] << ", " << rhs_test[i*3+2] << std::endl;
        // }

        LayerPotenOp0.SetTargetCoord(X0surf);
        sctl::Vector<Real> sigma_ubg;
        solver(&sigma_ubg, BIO, bg_flow(X0surf)*(pressure_drop/period_length), gmres_tol, gmres_max_iter);
        // elem_lst0.WriteVTK("vis/Ubg-3p", bg_flow(X0surf) * (pressure_drop/period_length), comm);

    
        { // Evaluate in interior, and write visualization
            VolumeVis<Real> vol_vis(elem_lst0, comm);
            const auto& X0 = vol_vis.GetCoord(); // set new target coordinates
            LayerPotenOp0.SetTargetCoord(X0);
            // sctl::Vector<Real> U;
            BIO(&Ubg2p, sigma_ubg);
            // vol_vis.WriteVTK("vis/U-3p-Ubg-b4", Ubg3p);
            Ubg2p -= bg_flow(X0) * (pressure_drop/period_length);
            // vol_vis.WriteVTK("vis/U-3p-Ubg", Ubg3p);
        }

        { // Evaluate in interior, and write visualization
            VolumeVis<Real> vol_vis(elem_lst0, comm);
            const auto& X0 = vol_vis.GetCoord(); // set new target coordinates
            LayerPotenOp0.SetTargetCoord(X0);
            // sctl::Vector<Real> U;
            BIO(&Up2p, sigma);
            // vol_vis.WriteVTK("vis/U-3p-b4", Up3p);
            Up2p -= eval_rhs(pressure_drop);
            // vol_vis.WriteVTK("vis/U-3p", Up3p);
        }

        // Compare diff magnitude
        std::cout << "2 peri Magnitude of difference between U from background flow and U from force densities: " << std::endl;
        sctl::Vector<Real> Udiff = Ubg2p - Up2p;
        // for (int i=0; i<Up2p.Dim()/3; i++) {
        for (int i=0; i<5; i++) {
            std::cout << std::setprecision(10)<< sctl::sqrt<Real>(Udiff[i*3]*Udiff[i*3] + Udiff[i*3+1]*Udiff[i*3+1] + Udiff[i*3+2]*Udiff[i*3+2]) << std::endl;
        }   
    }

    sctl::Vector<Real> Ubg1p, Up1p;
    { // one-periodic
        LayerPotenOp0.SetTargetCoord(X0surf);
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);

        sctl::Vector<Real> sigma;
        sctl::GMRES<Real> solver(comm);
        solver(&sigma, BIO, bg_flow(X0surf) * (pressure_drop/period_length), gmres_tol, gmres_max_iter);
        // elem_lst0.WriteVTK("vis/Ubg-1p", bg_flow(X0surf) * (pressure_drop/period_length), comm);

        sctl::Vector<Real> sigma_up;
        LayerPotenOp0.SetTargetCoord(X0surf);
        solver(&sigma_up, BIO, eval_rhs(pressure_drop), gmres_tol, gmres_max_iter);
        // elem_lst0.WriteVTK("vis/Up-1p", eval_rhs(pressure_drop), comm);
        // std::cout << "Eval rhs in 1p" << std::endl;
        // sctl::Vector<Real> rhs_test =  eval_rhs(pressure_drop);
        // for (int i=0; i<rhs_test.Dim()/3; i++) {
        //     std::cout << rhs_test[i*3+0] << ", " << rhs_test[i*3+1] << ", " << rhs_test[i*3+2] << std::endl;
        // }

        { // Evaluate in interior, and write visualization
            VolumeVis<Real> vol_vis(elem_lst0, comm);
            const auto& X0 = vol_vis.GetCoord(); // set new target coordinates
            LayerPotenOp0.SetTargetCoord(X0);
            // sctl::Vector<Real> U;
            BIO(&Ubg1p, sigma);
            // vol_vis.WriteVTK("vis/U-1p-Ubg-b4", Ubg1p);
            Ubg1p -= bg_flow(X0) * (pressure_drop/period_length);
            // vol_vis.WriteVTK("vis/U-1p-Ubg", Ubg1p);
        }

        { // Evaluate in interior, and write visualization
            VolumeVis<Real> vol_vis(elem_lst0, comm);
            const auto& X0 = vol_vis.GetCoord(); // set new target coordinates
            LayerPotenOp0.SetTargetCoord(X0);
            // sctl::Vector<Real> U;
            BIO(&Up1p, sigma_up);
            // vol_vis.WriteVTK("vis/U-1p-Up-b4", Up1p);
            Up1p -= eval_rhs(pressure_drop);
            // vol_vis.WriteVTK("vis/U-1p-Up", Up1p);
        }

        // Compare diff magnitude
        std::cout << "1 peri Magnitude of difference between U from background flow and U from force densities: " << std::endl;
        sctl::Vector<Real> Udiff = Up1p - Ubg1p;
        // for (int i=0; i<Up1p.Dim()/3; i++) {
        for (int i=0; i<5; i++) {
            std::cout << std::setprecision(10)<< sctl::sqrt<Real>(Udiff[i*3]*Udiff[i*3] + Udiff[i*3+1]*Udiff[i*3+1] + Udiff[i*3+2]*Udiff[i*3+2]) << std::endl;
        }
    }

    sctl::Vector<Real> Udiff;
    // 3 peri and 1 peri with Up
    std::cout << " Magnitude of difference between 3-peri and 1-peri, U from force densities (might be different): " << std::endl;
    Udiff = Up1p - Up3p;
    // for (int i=0; i<Up1p.Dim()/3; i++) {
    for (int i=0; i<5; i++) {
        std::cout << std::setprecision(10)<< sctl::sqrt<Real>(Udiff[i*3]*Udiff[i*3] + Udiff[i*3+1]*Udiff[i*3+1] + Udiff[i*3+2]*Udiff[i*3+2]) << std::endl;
    }

    // 3 peri and 1 peri with Ubg
    std::cout << " Magnitude of difference between 3-peri and 1-peri, U from background flow (should be the same): " << std::endl;
    Udiff = Ubg1p - Ubg3p;
    // for (int i=0; i<Ubg1p.Dim()/3; i++) {
    for (int i=0; i<5; i++) {
        std::cout << std::setprecision(10)<< sctl::sqrt<Real>(Udiff[i*3]*Udiff[i*3] + Udiff[i*3+1]*Udiff[i*3+1] + Udiff[i*3+2]*Udiff[i*3+2]) << std::endl;
    }

    // 3 peri with Urhs and 1 peri with Ubg
    std::cout << " Magnitude of difference between 3-peri from force density and 1-peri from background flow (might be different): " << std::endl;
    Udiff = Ubg1p - Up3p;
    // for (int i=0; i<Ubg1p.Dim()/3; i++) {
    for (int i=0; i<5; i++) {
        std::cout << std::setprecision(10)<< sctl::sqrt<Real>(Udiff[i*3]*Udiff[i*3] + Udiff[i*3+1]*Udiff[i*3+1] + Udiff[i*3+2]*Udiff[i*3+2]) << std::endl;
    }

    // 3 peri and 1 peri with Up
    std::cout << " Magnitude of difference between 3-peri and 1-peri, U from force densities (might be different): " << std::endl;
    Udiff = Up1p - Up3p;
    // for (int i=0; i<Up1p.Dim()/3; i++) {
    for (int i=0; i<5; i++) {
        std::cout << std::setprecision(10)<< sctl::sqrt<Real>(Udiff[i*3]*Udiff[i*3] + Udiff[i*3+1]*Udiff[i*3+1] + Udiff[i*3+2]*Udiff[i*3+2]) << std::endl;
    }

    // 2 peri and 1 peri with Ubg
    std::cout << " Magnitude of difference between 2-peri and 1-peri, U from background flow (should be the same): " << std::endl;
    Udiff = Ubg1p - Ubg2p;
    // for (int i=0; i<Ubg1p.Dim()/3; i++) {
    for (int i=0; i<5; i++) {
        std::cout << std::setprecision(10)<< sctl::sqrt<Real>(Udiff[i*3]*Udiff[i*3] + Udiff[i*3+1]*Udiff[i*3+1] + Udiff[i*3+2]*Udiff[i*3+2]) << std::endl;
    }
}

int main(int argc, char** argv) {
    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Profile::Enable(true);
        const sctl::Comm comm = sctl::Comm::World();
        test<Real>(comm);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}