#include "periodize.hpp"
#include "utils.hpp"
#include <fstream>

/**
 * Background flow with unit pressure gradient along X-axis.
 */
template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
  const sctl::Long N = X.Dim()/3;
  sctl::Vector<Real> U(N*3);
  for (sctl::Long i = 0; i < N; i++) {
    const auto x = X.begin() + i*3;
    U[i*3+0] = -((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
    // U[i*3+0] = -((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/25;
    U[i*3+1] = 0;
    U[i*3+2] = 0;
  }
  return U;
}

/**
 * Reference solution for checking error.
 */
template <class Real> sctl::Vector<Real> u_ref(const sctl::Vector<Real>& X) {
  const sctl::Long N = X.Dim()/3;
  sctl::Vector<Real> U(N*3);
  for (sctl::Long i = 0; i < N; i++) {
    const auto x = X.begin() + i*3;
    U[i*3+0] = 1e-2 - ((x[1]-0.4)*(x[1]-0.4) + (x[2]-0.3)*(x[2]-0.3))/4;
    U[i*3+1] = 0;
    U[i*3+2] = 0;
  }
  return U;
}

// Self convergence solution read from file.
template <class Real> sctl::Vector<Real> Read_u_ref() {
  sctl::Vector<Real> U;
  // U.Read("out/U_8_16.txt");
  U.Read("out/U_8_16_1.txt");
  return U;
}

template <class Real> void test(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Comm comm) {

  // std::cout << "size of overall comm is " << comm.Size() <<"; rank in overall comm is "<< comm.Rank() << std::endl;

  const Real SL_scal = 1.0;
  const Real DL_scal = 1.0;

  const Real tol = 1e-15;
  const Real gmres_tol = 1e-14;
  const sctl::Long ElemOrder = 10;

  const auto build_elem_lst_nbr = [](const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range){
    sctl::Vector<Real> Xc, eps, orient;
    sctl::Vector<sctl::Long> ElemOrderVec, FourierOrderVec;
    for (sctl::Long k0 = -nbr_range; k0 <= nbr_range; k0++) {
      for (sctl::Long k1 = 0; k1 <= 0; k1++) {
        for (sctl::Long k2 = 0; k2 <= 0; k2++) {
          for (sctl::Long i = 0; i < Nelem; i++) {
            ElemOrderVec.PushBack(ElemOrder);
            FourierOrderVec.PushBack(FourierOrder);
            const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
            for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
              const Real x = (i+nodes[j])/Nelem;
              Xc.PushBack(k0+x);
              // // Straight channel
              // Xc.PushBack(k1+0.4);
              // Xc.PushBack(k2+0.3);
              // eps.PushBack(0.2);

              // // wavy channel 1
              // Xc.PushBack(k1 + 0.1*cos(2*sctl::const_pi<Real>()*x)+0.5);
              // Xc.PushBack(k2+0.5);
              // eps.PushBack(0.1);

              // wavy channel 2
              Xc.PushBack(k1 + 0.3*cos(2*sctl::const_pi<Real>()*x)+0.5);
              Xc.PushBack(k2+0.5);
              eps.PushBack(0.1);

              orient.PushBack(0);
              orient.PushBack(0);
              orient.PushBack(1);
            }
          }
        }
      }
    }
    sctl::SlenderElemList<Real> elem_lst(ElemOrderVec, FourierOrderVec, Xc, eps, orient);
    return elem_lst;
  };
  const auto elem_lst0 = build_elem_lst_nbr(Nelem, ElemOrder, FourierOrder, 0); // geometry in the unit box [0,1]^3
  const auto elem_lst_nbr = build_elem_lst_nbr(Nelem, ElemOrder, FourierOrder, 1); // geometry with one set of images in each direction
  const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); // should be 3^3 = 27

  sctl::Vector<Real> X0; // target coordinates
  elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
  const auto X_proxy = Periodize<Real>::GetProxySurf(); // proxy points coordinates

  StokesBIO LayerPotenOp0(SL_scal, DL_scal,comm); // potential from elem_lst_nbr to X0
  LayerPotenOp0.AddElemList(elem_lst_nbr);
  LayerPotenOp0.SetTargetCoord(X0);
  LayerPotenOp0.SetAccuracy(tol);

  StokesBIO LayerPotenOp_proxy(SL_scal, DL_scal,comm); // potential from elem_lst0 to proxy points
  LayerPotenOp_proxy.AddElemList(elem_lst0);
  LayerPotenOp_proxy.SetTargetCoord(X_proxy);
  LayerPotenOp_proxy.SetAccuracy(tol);

  // periodized layer potential operator
  const auto BIO = [&comm,&DL_scal,&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    // std::cout <<"ID " << comm.Rank() << " in BIO." << std::endl;
    const sctl::Long N = sigma.Dim();

    sctl::Vector<Real> sigma_nbr(Nrepeat*N); // repeat sigma Nrepeat times
    for (sctl::Long k = 0; k < Nrepeat; k++) {
      for (sctl::Long i = 0; i < N; i++) {
        sigma_nbr[k*N+i] = sigma[i];
      }
    }

    U->SetZero();
    LayerPotenOp0.ComputePotential(*U, sigma_nbr);
    if (DL_scal && U->Dim() == N) (*U) -= sigma*0.5 * DL_scal; // for double-layer

    { // Add far-field
      // std::cout << "ID " << comm.Rank() << " in far eval." << std::endl;
      sctl::Vector<Real> U_proxy, U_far;
      LayerPotenOp_proxy.ComputePotential(U_proxy, sigma);
      Periodize<Real>::EvalFarField(U_far, X0, U_proxy);
      (*U) += U_far;
    }
  };

  // Solve for sigma to satisfy no-slip boundary conditions: BIO(sigma) + bg_flow = 0
  sctl::Vector<Real> sigma;
  sctl::GMRES<Real> solver(comm);
  solver(&sigma, BIO, -bg_flow(X0), gmres_tol);
  // elem_lst0.WriteVTK("vis/bgflow1",bg_flow(X0));
  // elem_lst0.WriteVTK("vis/sigma", sigma);

  { // Evaluate in interior, compute error and write visualization
    std::cout << "in eval block, sigma dim = " << sigma.Dim() << std::endl;
    const sctl::Long Nelem_trg = 8;
    const sctl::Long FourierOrder_trg = 16;
    const auto elem_lst_trg = build_elem_lst_nbr(Nelem_trg, ElemOrder, FourierOrder_trg, 0);

    VolumeVis<Real> cube(elem_lst0, comm);
    X0 = cube.GetCoord(); // set new target coordinates
    LayerPotenOp0.SetTargetCoord(X0);
    sctl::Vector<Real> U;
    std::cout <<"size of BIO is " << LayerPotenOp0.Dim(0) << ", " << LayerPotenOp0.Dim(1) << std::endl;
    BIO(&U, sigma);
    U += bg_flow(X0);

    // cube.WriteVTK("vis/U_ref", U);
    // cube.WriteVTK("vis/err", err);

    // if this is the reference parameters, write to file.
    if (write_ref) {
      std::cout << "in writing block" << std::endl;
      U.Write("out/U_8_16_1.txt");
      cube.WriteVTK("vis/U_ref_1", U);
      elem_lst_nbr.WriteVTK("vis/S-nbr_1");
      // U.Write("out/U_8_16.txt");
      // cube.WriteVTK("vis/U_ref", U);
      // elem_lst_nbr.WriteVTK("vis/S-nbr");
    } else {
      std::cout << "not write-ref, reading U_ref" << std::endl;
      sctl::Vector<double> U_ref = Read_u_ref<double>();
      // error code for self-convergence.
      double max_err = 0;
      const auto err = U - U_ref;
      for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
      std::cout<<"Max error = "<<max_err<<'\n';
    }

  }
}

int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    const sctl::Comm comm = sctl::Comm::World();
    long Nelem  = std::stol(argv[1]); // number of elements
    long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
    int write_ref = std::stoi(argv[3]); // whether the parameters are considered ''artifical true solution''.
    test<Real>(Nelem, FourierOrder, (write_ref==1), comm);
  }

  sctl::Comm::MPI_Finalize();
  // std::cout << "after finalizing MPI" << std::endl;
  return 0;
}

