#include "periodize.hpp"
#include "utils.hpp"

int main(int argc, char** argv) {
    sctl::Comm::MPI_Init(&argc, &argv);
    // using Real = double;
    using Real = sctl::QuadReal;

    {
      sctl::Comm comm = sctl::Comm::World();

      sctl::Long level_in = std::stol(argv[1]);
      sctl::Long m0_in = std::stol(argv[2]);
      
      // Set ProxySurf
      sctl::Vector<Real> X_proxy = Periodize1D<Real>::GetProxySurf(level_in, m0_in); // proxy points coordinates
      sctl::Vector<Real> U_proxy = X_proxy;

      // Set target for eval (arbitrary)
      sctl::Vector<Real> X0(3);
      X0 = 0.5;

      sctl::Vector<Real> U_far;
      Periodize1D<Real>::EvalFarField(U_far, X0, U_proxy, level_in, m0_in);

    }

    sctl::Comm::MPI_Finalize();
    return 0;
}