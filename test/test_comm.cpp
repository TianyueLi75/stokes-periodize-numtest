
#include "csbq.hpp"

int main(int argc, char** argv) {
    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;
  
    {
      sctl::Comm comm = sctl::Comm::World();
      std::cout << "Rank " << comm.Rank() << " printing." << std::endl;
      
      sctl::Vector<sctl::Long> size_loc(1);
      size_loc[0] = 7;
      sctl::Vector<sctl::Long> size_all(1);
      size_all[0] = 0;
      comm.Allreduce((sctl::Iterator<sctl::Long>) size_loc.begin(), (sctl::Iterator<sctl::Long>) size_all.begin(), 1, sctl::CommOp::SUM);
      std::cout << "rank " << comm.Rank() << " size loc = " << size_loc[0] << ", size all is " << size_all[0] << std::endl;
    
      // sctl::Vector<Real> vec_loc(7);
      // vec_loc = (Real) comm.Rank();
      // sctl::Vector<Real> vec_all(vec_loc.Dim() * comm.Size());
      // comm.Allgather((sctl::Iterator<Real>) vec_loc.begin(), size_loc[0], (sctl::Iterator<Real>) vec_all.begin(), size_all[0]);
      // if (comm.Rank() == 0) {
      //   for (int i=0; i < size_all[0]; i++) {
      //       std::cout << vec_all[i] << std::endl;
      //   }
      // }

      // if (comm.Rank() == 1) {
      //   vec_all.ReInit(0);
      // }
      // comm.Barrier();
      // comm.PartitionN(vec_all,size_loc[0]);
      // comm.Barrier();
      // if (comm.Rank() == 0) {
      //   std::cout << "rank 0" << std::endl;
      //   for (int i=0; i < vec_all.Dim(); i++) {
      //       std::cout << vec_all[i] << std::endl;
      //   }
      // }
      // comm.Barrier();
      // if (comm.Rank() == 1) {
      //   std::cout << "rank 1" << std::endl;
      //   for (int i=0; i < vec_all.Dim(); i++) {
      //       std::cout << vec_all[i] << std::endl;
      //   }
      // }
    }
  
    sctl::Comm::MPI_Finalize();
    return 0;
}

