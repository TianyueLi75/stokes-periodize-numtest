/*
    Time the precomputation time corresponding to certain quadrature tolerance levels.
    No MPI parallel
    Based on https://github.com/dmalhotra/stokes-periodize test/test2.cpp
*/

// Boundary integral operators
#include "stokes_bio.hpp" 

// Geometry for tests
#include "utils_geom.hpp"

// Other util functions
#include "utils_tests.cpp" 

template <class Real> void precompute_time(
    const Real tol, 
    const Real SL_scal, 
    const Real DL_scal, 
    sctl::Comm comm) {

    const Real period_length = 1;

    const sctl::Long Nelem = 4;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 28;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs, ptcls_rs, NormalOrient;
    sctl::SlenderElemList<Real> elem_lst0;
    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_straight(Nelem, ElemOrder, FourierOrder, 0.2, comm, ptcls, ptcls_rs, ptcls_Xcs, 0);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);

    sctl::Vector<Real> X0, wts;
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    Real surface_area;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        surface_area = surface_area_[0];
    }

    StokesBIO<Real> LayerPotenOp0(SL_scal, DL_scal, comm);
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length); 

    // Periodized layer potential operator
    const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> sigma_mean, sigma0;
        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_;
            elem_lst0.GetFarFieldDensity(sigma_, sigma);
            SurfaceIntegral(sigma_mean, sigma_, wts);
            sigma_mean *= (1/surface_area);
            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);
        }

        // Compute the periodic solution using sigma0
        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma0.Dim()) (*U) -= sigma0 * 0.5 * NormalOrient * DL_scal; // for double-layer on interior of channel

        // Add back sigma_mean
        AddConstVec(*U, sigma_mean);
    };

    // Arbitrary sigma
    sctl::Vector<Real> sigma(X0.Dim());
    srand48(2);
    for (sctl::Long i=0; i<X0.Dim(); i++) {
        sigma[i] = (drand48() - 0.5);
    }
    sctl::Vector<Real> Ueval;
    BIO(&Ueval, sigma);
    sctl::Profile::print(&comm);

}

int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    sctl::Comm comm = sctl::Comm::World();
    sctl::Profile::Enable(true);
    double quadr_tol = std::stod(argv[1]); // Quadrature accuracy that corresponds to different multipole order.
    double SL_scal = std::stod(argv[2]);
    double DL_scal = std::stod(argv[3]);
    precompute_time<Real>(quadr_tol, SL_scal, DL_scal, comm);    
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}