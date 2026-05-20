// Time the precomputation time corresponding to certain quadrature tolerance levels.
// No MPI parallel
// Based on code in test_3peri, which corresponds to https://github.com/dmalhotra/stokes-periodize test/test2.cpp

#include "utils.hpp"

/**
* Compute the weighted sum of vals with weights wts.
* vals: [x1,y1,z1, x2,y2,z2, ..., xN,yN,zN]
* wts:  [w1, w2, ..., wN]
* I: [mx, my, mz]
*/
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

/**
 * Add vector c0 to each point in vals.
 * vals: [x1,y1,z1, x2,y2,z2, ..., xN,yN,zN]
 * c0:   [cx,cy,cz]
 */
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
 * Build element list for a straight channel.
 */
template <class Real> sctl::SlenderElemList<Real> build_elem_lst(const sctl::Long Nelem_channel, const sctl::Long ElemOrder, const sctl::Long FourierOrder) {
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
    sctl::SlenderElemList<Real> elem_lst0(ElemOrderVec, FourierOrderVec, Xc, eps, orient);
    return elem_lst0;
}

template <class Real> void test(const Real tol, const Real SL_scal, const Real DL_scal, sctl::Comm comm) {

    const Real period_length = 1;

    const sctl::Long Nelem_channel = 4;
    const sctl::Long ElemOrder = 10;
    const sctl::Long FourierOrder = 28;

    sctl::SlenderElemList<Real> elem_lst0 = build_elem_lst<Real>(Nelem_channel, ElemOrder, FourierOrder); // geometry in the unit box [0,1]^3

    Real surface_area;
    sctl::Vector<Real> X0, wts;
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Vector<Real> X0surf = X0;
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
    LayerPotenOp0.SetTargetCoord(X0surf);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length); // TODO: are all matrices precomputed?

    // Define the boundary-integral operator: K[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> sigma_mean, sigma0;
        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_;
            elem_lst0.GetFarFieldDensity(sigma_, sigma);
            SurfaceIntegral(sigma_mean, sigma_, wts);
            sigma_mean *= (1/surface_area);
            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*DL_scal; // for double-layer on interior of channel

        AddConstVec(*U, sigma_mean);
    };

    // Arbitrary sigma
    sctl::Vector<Real> sigma(X0surf.Dim());
    srand48(2);
    for (sctl::Long i=0; i<X0surf.Dim(); i++) {
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
    double quadr_tol = std::stod(argv[1]); 
    double SL_scal = std::stod(argv[2]);
    double DL_scal = std::stod(argv[3]);
    test<Real>(quadr_tol, SL_scal, DL_scal, comm);    
  }

  sctl::Comm::MPI_Finalize();
  return 0;
}