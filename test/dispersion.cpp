// Solves an empty periodic channel problem, then save files of time evolution of Lagrangian particles following the solution flow.
// Using preconditioner.

#include "utils_vis.hpp"
#include "stokes_bio.hpp"
#include "bio_operator.hpp"
#include "utils_geom.hpp"
#include "utils_tests.cpp"


template <class Real> void trefoil_dispersion(sctl::Long Nelem_channel, sctl::Long FourierOrder, sctl::Comm comm) {
    
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real tol = 1e-14;
    const Real gmres_tol = 1e-9;
    const sctl::Long ElemOrder = 10;
    const Real period_length = 1.;
    const Real pressure_drop = -15.;
    const sctl::Long gmres_max_iter = 400;

    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient;
    sctl::Long geom_mode = 0;
    sctl::Long ptcl_ord = 1;

    std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.build_trefoil(Nelem_channel, ElemOrder, FourierOrder, comm, ptcls, ptcls_rs, ptcls_Xcs, ptcl_ord, geom_mode);
    elem_lst0 = std::get<0>(build0);
    NormalOrient = std::get<1>(build0);

    sctl::Vector<Real> X0;
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    elem_lst0.WriteVTK("vis/Trefoil_dispersion_geometry",X0,comm);
    Real surface_area;
    sctl::Vector<Real> wts;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area = sa_all[0];
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); 
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);

    MeanCorrectedStokesBIOOperator<Real> BIO(LayerPotenOp0, NormalOrient, DL_scal, comm);
    BIO.AddSurface(elem_lst0, wts, surface_area);

    // Old lambda-based BIO block retained for reference:
    // const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,&NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    //     sctl::Vector<Real> sigma_mean, sigma0;
    //     { 
    //         sctl::Vector<Real> sigma_;
    //         elem_lst0.GetFarFieldDensity(sigma_, sigma);
    //         SurfaceIntegral(sigma_mean, sigma_, wts);
    //         // MPI for total surface area
    //         sctl::Vector<Real> sa_loc = sigma_mean;
    //         sctl::Vector<Real> sa_all(3);
    //         sa_all = 0;
    //         comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
    //         comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+1, (sctl::Iterator<Real>) sa_all.begin()+1, 1, sctl::CommOp::SUM);
    //         comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+2, (sctl::Iterator<Real>) sa_all.begin()+2, 1, sctl::CommOp::SUM);
    //         sigma_mean = sa_all;
    //         sigma_mean *= (1./surface_area);
    //
    //         sigma0 = sigma;
    //         AddConstVec(sigma0, -sigma_mean);
    //     }
    //
    //     U->SetZero();
    //     LayerPotenOp0.ComputePotential(*U, sigma0);
    //     if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer
    //
    //     AddConstVec(*U, sigma_mean);
    // };

    std::string sigma_file = "out/trefoil_density_"+std::to_string(Nelem_channel)+"_"+std::to_string(FourierOrder)+"_"+std::to_string(comm.Rank())+".txt";
    sctl::Vector<Real> sigma;
    sigma.Read(sigma_file.c_str());

    if (!sigma.Dim()) {
        std::cout << "rank " << comm.Rank() << "couldn't read density, making new." << std::endl;
        // Preconditioning
        sctl::Matrix<Real> PrecondMat0, PrecondMat1;
        sctl::Long A11size;
        Real channel_radius = 0.035; // for trefoil geometry
        A11size = precond_channel(PrecondMat0, PrecondMat1, Nelem_channel, ElemOrder, FourierOrder, channel_radius, SL_scal, DL_scal, comm);
        // Apply A11inv to each panel of a vector.
        const auto AinvApply = [&PrecondMat0,&PrecondMat1,&A11size](const sctl::Vector<Real>& vec) {
            sctl::Vector<Real> AinvVec(vec.Dim());
            sctl::Long N = vec.Dim();
            sctl::Long Npanels = N / A11size; 
            for (sctl::Long i=0; i<Npanels; i++) {
                sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
                sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
                for (sctl::Long j=0; j<A11size; j++) {
                    AinvVec[i*A11size + j] = AinvVecMat(j,0);
                }
            }
            return AinvVec;
        };
        // Left diagonal preconditioner
        const auto BIO_precond = [&BIO,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
            sctl::Vector<Real> Uloc;
            BIO(&Uloc,sigma);
            (*U) = AinvApply(Uloc);
        };


        auto bg_flow = [](const sctl::Vector<Real>& X) {
            const sctl::Long N = X.Dim()/3;
            sctl::Vector<Real> U(N*3);
            for (sctl::Long i = 0; i < N; i++) {
                const auto x = X.begin() + i*3;
                U[i*3+0] = - ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
                U[i*3+1] = 0;
                U[i*3+2] = 0;
            }
            return U;
        };

        sctl::GMRES<Real> solver(comm);
        sctl::KrylovPrecond<Real> krylov;
        sctl::Vector<Real> A11invF = AinvApply(bg_flow(X0)*pressure_drop/period_length);
        solver(&sigma, BIO_precond, A11invF, gmres_tol, gmres_max_iter, false, nullptr, &krylov);

        std::cout << "rank " << comm.Rank() << "done with solve, writing sigma to file named " << sigma_file << std::endl;
        sigma.Write(sigma_file.c_str());
    } 

    { 
        PeriodicGeom<Real> trg;
        sctl::Long Nelem_trg=200;  
        const sctl::Long FourierOrder_trg = 8; // not used
        sctl::SlenderElemList<Real> elem_lst_trg;
        sctl::Vector<sctl::Long> ptcls_trg;
        sctl::Vector<Real> ptcls_Xcs_trg;
        sctl::Vector<Real> ptcls_rs_trg;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_trg = trg.build_trefoil(Nelem_trg, ElemOrder, FourierOrder_trg, comm, ptcls_trg, ptcls_rs_trg, ptcls_Xcs_trg, 1, geom_mode);
        elem_lst_trg = std::get<0>(build_trg);
        // Form targets at Ngroups cross sections, divided evenly among processes
        XsectionVis<Real> XsectVis(elem_lst_trg, comm);
        X0 = XsectVis.GetCoord();
        sctl::Vector<Real> U0;
        LayerPotenOp0.SetTargetCoord(X0);
        BIO(&U0, sigma);
        U0 -= bg_flow(X0)*pressure_drop/period_length;
        XsectVis.WriteVTK("vis/disperison/XsectionVis_t0",U0);

        Real T = 100000.;
        sctl::Long Nt = 10000;
        Real dt = T / Nt; 
        for (sctl::Long tind = 1; tind <= Nt; tind++) {
            sctl::Vector<Real> U;
            if (tind > 1) {
                LayerPotenOp0.SetTargetCoord(X0);
                BIO(&U, sigma);
                U -= bg_flow(X0)*pressure_drop/period_length;
            } else {
                U = U0;
            }
            for (sctl::Long xind=0; xind<X0.Dim()/3; xind++) {
                const Real current_x = X0[xind*3+0]+dt*U[xind*3+0];
                const Real current_y = X0[xind*3+1]+dt*U[xind*3+1];
                const Real current_z = X0[xind*3+2]+dt*U[xind*3+2];
                auto [is_in_trefoil, xc, yc, zc] = trg.in_trefoil(current_x,current_y,current_z);
                if (is_in_trefoil) {
                    if (current_x > 1+1e-5 || current_x < -1e-5) { 
                        X0[xind*3+0] = current_x - std::floor(current_x);
                        X0[xind*3+1] = current_y;
                        X0[xind*3+2] = current_z;
                    } else {
                        X0[xind*3+0] = current_x;
                        X0[xind*3+1] = current_y;
                        X0[xind*3+2] = current_z;
                    }
                } else {
                    // if new point would be out of trefoil, negate normal vector for repulsive force
                    Real nx = xc - current_x;
                    Real ny = yc - current_y;
                    Real nz = zc - current_z;
                    Real n2 = nx*nx + ny*ny + nz*nz;
                    Real udotn = U[xind*3+0] * nx + U[xind*3+1] * ny + U[xind*3+2] * nz;
                    udotn = udotn / n2; 
                    // buffer: u_new = u - 1.01*u_n
                    Real U_x = U[xind*3+0] - 1.01*udotn*nx;
                    Real U_y = U[xind*3+1] - 1.01*udotn*ny;
                    Real U_z = U[xind*3+2] - 1.01*udotn*nz;
                    Real new_x = X0[xind*3+0]+dt*U_x;
                    Real new_y = X0[xind*3+1]+dt*U_y;
                    Real new_z = X0[xind*3+2]+dt*U_z;
                    auto [is_in_trefoil2, xc2, yc2, zc2] = trg.in_trefoil(new_x, new_y, new_z);
                    if (is_in_trefoil2) { // double check new position back in trefoil, otherwise stays put.
                        if (new_x > 1+1e-5 || new_x < -1e-5) { 
                            X0[xind*3+0] = new_x - std::floor(new_x);
                            X0[xind*3+1] = new_y;
                            X0[xind*3+2] = new_z;
                        } else {
                            X0[xind*3+0] = new_x;
                            X0[xind*3+1] = new_y;
                            X0[xind*3+2] = new_z;
                        }
                    }
                }
            }
            XsectVis.SetCoord(X0);
            if (tind % 1000 == 0) {
                XsectVis.WriteVTK("vis/dispersion/XsectionVis_t"+std::to_string(tind),U);
            }
        }
    }
}

int main(int argc, char** argv) {
  sctl::Comm::MPI_Init(&argc, &argv);
  using Real = double;

  {
    sctl::Comm comm = sctl::Comm::World();
    long Nelem_channel = std::stol(argv[1]); // number of elements
    long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes

    trefoil_dispersion<Real>(Nelem_channel, FourierOrder, comm);

  }

  sctl::Comm::MPI_Finalize();
  return 0;
}

