/**
 Porous media flow in a periodic geometry
 May 15: 
    First try one ptcl inside straight channel.
*/

#include "periodize.hpp"
#include <fstream>
using namespace sctl;

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

// Stokes combined field operator (S+D)
template <sctl::Long SL_scal> struct Stokes3D_CF_ {
    static const std::string& Name() {
        // Name determines what quadrature tables to use.
        // Single-layer quadrature tables also works for combined fields.
        static const std::string name = "Stokes3D-FxU";
        return name;
    }
  
    static constexpr sctl::Integer FLOPS() {
        return 50;
    }
  
    template <class Real> static constexpr Real uKerScaleFactor() {
        return 1 / (8 * sctl::const_pi<Real>());
    }
  
    template <sctl::Integer digits, class VecType>
    static void uKerMatrix(VecType (&u)[3][3], const VecType (&r)[3], const VecType (&n)[3], const void* ctx_ptr) {
        using Real = typename VecType::ScalarType;
        const auto SL_scal_ = VecType((Real)SL_scal);
        const auto r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        const auto rinv = sctl::approx_rsqrt<digits>(r2, r2 > VecType::Zero()); // Compute inverse square root
        const auto rinv2 = rinv * rinv;
        const auto rinv3 = rinv2 * rinv;
        const auto rinv5 = rinv3 * rinv2;
        const auto rdotn = r[0] * n[0] + r[1] * n[1] + r[2] * n[2];
        const auto rdotn_rinv5_6 = VecType((Real)6) * rdotn * rinv5;
        for (sctl::Integer i = 0; i < 3; i++) {
            for (sctl::Integer j = 0; j < 3; j++) {
                const auto ri_rj = r[i] * r[j];
                const auto ker_dl = ri_rj * rdotn_rinv5_6; // Double-layer kernel
                const auto ker_sl = (i == j ? rinv + ri_rj * rinv3 : ri_rj * rinv3); // Single-layer kernel
                u[i][j] = ker_dl + ker_sl * SL_scal_; // Combine kernels
                // u[i][j] = ker_sl; // only SL
            }
        }
    }
};
using Stokes3D_CF = sctl::GenericKernel<Stokes3D_CF_<1>>;
  
/**
* Visualize volume inside SlenderElemList.
*/
template <class Real> class VolumeVis {
        static constexpr sctl::Integer COORD_DIM = 3;
        static constexpr sctl::Integer s_order = 20;
        static constexpr sctl::Integer t_order = 60;
        static constexpr sctl::Integer r_order = 12;
    public:

        VolumeVis() = default;

        /**
        * @brief Construct a new VolumeVis object.
        *
        * @param elem_lst the geometry.
        * @param comm MPI communicator.
        */
        VolumeVis(const sctl::SlenderElemList<Real>& elem_lst, const sctl::Comm& comm = sctl::Comm::Self()) : comm_(comm) {
            Nelem = elem_lst.Size();
            sctl::Vector<Real> s_param, sin_theta, cos_theta;
            for (sctl::Long i = 0; i < s_order; i++) {
                const Real t = i/(Real)(s_order-1);
                s_param.PushBack(t);
            }
            for (sctl::Long i = 0; i < t_order; i++) {
                const Real t = i/(Real)t_order;
                sin_theta.PushBack(sctl::sin<Real>(2*sctl::const_pi<Real>()*t));
                cos_theta.PushBack(sctl::cos<Real>(2*sctl::const_pi<Real>()*t));
            }
            for (sctl::Long elem_idx = 0; elem_idx < Nelem; elem_idx++) {
                const Real t_order_inv = 1/(Real)t_order;
                const Real r_order_inv = (1-1e-6)/(Real)(r_order-1);
                sctl::Vector<Real> X_, Xc(COORD_DIM);
                elem_lst.GetGeom(&X_, nullptr, nullptr, nullptr, nullptr, s_param, sin_theta, cos_theta, elem_idx);
                for (sctl::Long i = 0; i < s_order; i++) {
                    Xc = 0;
                    for (sctl::Long j = 0; j < t_order; j++) {
                        for (sctl::Long l = 0; l < COORD_DIM; l++) {
                            Xc[l] += X_[(i*t_order+j)*COORD_DIM+l] * t_order_inv;
                        }
                    }
                    for (sctl::Long j = 0; j < t_order; j++) {
                        for (sctl::Long k = 0; k < r_order; k++) {
                            for (sctl::Long l = 0; l < COORD_DIM; l++) {
                                coord.PushBack((X_[(i*t_order+j)*COORD_DIM+l]-Xc[l])*k*r_order_inv + Xc[l]);
                            }
                        }
                    }
                }
            }
        }
  
        /**
        * @brief Get the coordinates of the discretization points.
        *
        * @return const Vector<Real>& Vector containing the coordinates.
        */
        const sctl::Vector<Real>& GetCoord() const {
            return coord;
        }
  
        /**
        * @brief Write the volume to a VTK file.
        *
        * @param fname File name.
        * @param F Data associated with the discretization points.
        */
        void WriteVTK(const std::string& fname, const sctl::Vector<Real>& F) const {
            sctl::VTUData vtu_data;
            GetVTUData(vtu_data, F);
            vtu_data.WriteVTK(fname, comm_);
        }
    
  
        /**
        * @brief Get VTU data.
        *
        * @param vtu_data VTU data object.
        * @param F Data associated with the discretization points.
        */
        void GetVTUData(sctl::VTUData& vtu_data, const sctl::Vector<Real>& F) const {
            for (const auto& x : coord) vtu_data.coord.PushBack((float)x);
            for (const auto& x :     F) vtu_data.value.PushBack((float)x);
            for (sctl::Long l = 0; l < Nelem; l++) {
                const sctl::Long offset = l * s_order*t_order*r_order;
                for (sctl::Long i = 0; i < s_order-1; i++) {
                    for (sctl::Long j = 0; j < t_order; j++) {
                        for (sctl::Long k = 0; k < r_order-1; k++) {
                            auto idx = [this,&offset](sctl::Long i, sctl::Long j, sctl::Long k) {
                                return offset+(i*t_order+(j%t_order))*r_order+k;
                            };
                            vtu_data.connect.PushBack(idx(i+0,j+0,k+0));
                            vtu_data.connect.PushBack(idx(i+0,j+0,k+1));
                            vtu_data.connect.PushBack(idx(i+0,j+1,k+1));
                            vtu_data.connect.PushBack(idx(i+0,j+1,k+0));
                            vtu_data.connect.PushBack(idx(i+1,j+0,k+0));
                            vtu_data.connect.PushBack(idx(i+1,j+0,k+1));
                            vtu_data.connect.PushBack(idx(i+1,j+1,k+1));
                            vtu_data.connect.PushBack(idx(i+1,j+1,k+0));
                            vtu_data.offset.PushBack(vtu_data.connect.Dim());;
                            vtu_data.types.PushBack(12);
                        }
                    }
                }
            }
        }
  
    private:
  
      sctl::Comm comm_;
      sctl::Long Nelem;
      sctl::Vector<Real> coord;
};

template <class Real> void InitGeom(Vector<Real>& X, Vector<Real>& R, Vector<Real>& OrientVec, Vector<Long>& ElemOrder, Vector<Long>& FourierOrder, Vector<Real>& panel_len, Vector<Real>& Mr_lst, Vector<Long>& cnt, Vector<Long>& dsp, const Long Nobj, const Real loop_rad, const Geom& geom_type) {
    srand48(2);
    const Long ElemOrder0 = 10;
    // Set ElemOrder, FourierOrder, panel_len, Mr_lst, cnt, dsp (Uniform discretization)
    const Long Npanel = 16;
    const Long FourierOrder0 = 16;

    cnt.ReInit(Nobj); cnt = Npanel;
    dsp.ReInit(Nobj); dsp = 0;
    omp_par::scan(cnt.begin(), dsp.begin(), Nobj);

    ElemOrder.ReInit(Nobj * Npanel);
    FourierOrder.ReInit(Nobj * Npanel);
    panel_len.ReInit(Nobj * Npanel);

    ElemOrder = ElemOrder0;
    FourierOrder = FourierOrder0;
    panel_len = 1/(Real)Npanel;

    auto loop_geom = [&loop_rad](Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta){
        x = loop_rad * cos<Real>(theta);
        y = loop_rad * sin<Real>(theta);
        z = 0;
        ex = 0;
        ey = 0;
        ez = 1;
        r = 0.025;
    };
    auto bacteria_geom = [&loop_rad](Real& x, Real& y, Real& z, Real& ex, Real& ey, Real& ez, Real& r, const Real theta){
        Real t = theta/const_pi<Real>()-1; // -1:1
        Real aspect = const_pi<Real>()*3/2+1;

        Real L = aspect-1+const_pi<Real>()/2;
        Real scal = loop_rad/(1+L-const_pi<Real>()/2) * 0.7;
        if (L*(1+t) < const_pi<Real>()/2) z = scal * (-cos<Real>(L*(1+t)) - L+const_pi<Real>()/2);
        else if (L*(1-t) < const_pi<Real>()/2) z = scal * (cos<Real>(L*(1-t)) + L-const_pi<Real>()/2);
        else z = scal * L * t;

        y = 0;
        x = 0;
        ex = 1;
        ey = 0;
        ez = 0;

        if (L*(1+t) < const_pi<Real>()/2) r = scal * sin<Real>(L*(1+t));
        else if (L*(1-t) < const_pi<Real>()/2) r = scal * sin<Real>(L*(1-t));
        else r = scal;
    };

    X.ReInit(0);
    R.ReInit(0);
    Mr_lst.ReInit(0);
    for (Long i = 0; i < Nobj; i++) {
        Real X0, Y0, Z0;
        { // Set offsets X0, Y0, Z0
            const Long N = (Long)ceil((double)pow<Real>((Real)Nobj,1/(Real)3));
            X0 = (i/pow<0>(N))%N;
            Y0 = (i/pow<1>(N))%N;
            Z0 = (i/pow<2>(N))%N;

            if (geom_type == Geom::Loop) {
            if (Nobj>2) Z0 += drand48()*0.5;
            } else if (geom_type == Geom::Bacteria) {
            if (Nobj>2) {
                X0 = X0*0.8 + drand48()*0.4;
                Y0 = Y0*0.8 + drand48()*0.4;
                Z0 = Z0*1.6 + drand48()*0.4;
            } else {
                X0 = X0*0.2;
            }
            } else {
                SCTL_ASSERT(false); // not implemented
            }
        }

        Real s_dsp = 0;
        // TODO: change Init Geom to handle 1D neighbors; add wall class?
        // for (sctl::Long k0 = -nbr_range; k0 <= nbr_range; k0++) {
        //     for (sctl::Long k1 = 0; k1 <= 0; k1++) {
        //         for (sctl::Long k2 = 0; k2 <= 0; k2++) {
        //             for (sctl::Long i = 0; i < Nelem; i++) {
        //                 ElemOrderVec.PushBack(ElemOrder);
        //                 FourierOrderVec.PushBack(FourierOrder);
        //                 const sctl::Vector<Real>& nodes = sctl::SlenderElemList<Real>::CenterlineNodes(ElemOrderVec[i]);
        //                 for (sctl::Long j = 0; j < ElemOrderVec[i]; j++) {
        //                     const Real x = (i+nodes[j])/Nelem;
        //                     //TODO: particle Xc and orientation.
        //                 }
        //             }
        //         }
        //     }
        // } 
        for (Long j = 0; j < cnt[i]; j++) { // Set X, OrientVec, R
            const Long ElemOrder_ = ElemOrder[dsp[i]+j];
            const auto& nds = SlenderElemList<Real>::CenterlineNodes(ElemOrder_);
            for (Long k = 0; k < ElemOrder_; k++) {
                Real x, y, z, ex, ey, ez, r;
                Real s = s_dsp + nds[k]*panel_len[dsp[i]+j];
                if (geom_type == Geom::Loop) {
                    loop_geom(x, y, z, ex, ey, ez, r, 2*const_pi<Real>()*s);
                } else if (geom_type == Geom::Bacteria) {
                    bacteria_geom(x, y, z, ex, ey, ez, r, 2*const_pi<Real>()*s);
                } else {
                    SCTL_ASSERT(false); // not implemented
                }
                X.PushBack(x+X0);
                X.PushBack(y+Y0);
                X.PushBack(z+Z0);
                OrientVec.PushBack(ex);
                OrientVec.PushBack(ey);
                OrientVec.PushBack(ez);
                R.PushBack(r);
            }
            s_dsp += panel_len[dsp[i]+j];
        }

        for (Long j = 0; j < COORD_DIM; j++) { // Set Mr
            for (Long k = 0; k < COORD_DIM; k++) {
                Mr_lst.PushBack(j==k?1:0);
            }
        }
    }
}

// TODO: copy over InitElemList
  

template <class Real> void test(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref) {
    using KerType = Stokes3D_CF;
    //using KerType = sctl::Stokes3D_FxU; // unstable if gmres_tol is too small
  
    const Real tol = 1e-14;
    const Real gmres_tol = 1e-14;
    const sctl::Long ElemOrder = 10;
  
    const sctl::Comm comm = sctl::Comm::Self();
    const KerType stokes_ker;
  
    const auto build_elem_lst_nbr_wall = [](const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range){
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
                            // straight channel
                            Xc.PushBack(k1+0.4);
                            Xc.PushBack(k2+0.3);
                            eps.PushBack(0.2);
                            // wavy channel
                            // Xc.PushBack(k1 + 0.1*cos(2*sctl::const_pi<Real>()*x)+0.5);
                            // Xc.PushBack(k2+0.5);
                            // eps.PushBack(0.1);
            
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
    const auto elem_lst0 = build_elem_lst_nbr_wall(Nelem, ElemOrder, FourierOrder, 0); // geometry in the unit box [0,1]^3
    const auto elem_lst_nbr = build_elem_lst_nbr_wall(Nelem, ElemOrder, FourierOrder, 1); // geometry with one set of images in each direction
    const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); // should be 3
  
    // TODO
    long Nobj = 1;
    double loop_rad = 0.2;
    string geom_type = "bacteria";
    Vector<double> X_ptcl, R_ptcl, OrientVec_ptcl, ElemOrder_ptcl, FourierOrder_ptcl, panel_len_ptcl, Mr_lst_ptcl, obj_elem_cnt_ptcl, obj_elem_dsp_ptcl;
    InitGeom(X, R, OrientVec, ElemOrder, FourierOrder, panel_len, Mr_lst, obj_elem_cnt, obj_elem_dsp, Nobj, loop_rad, geom_type);
    Vector<double> loc_elem_cnt_ptcl, loc_elem_dsp_ptcl, elem_lst_ptcl;
    InitElemList(loc_elem_cnt, loc_elem_dsp, elem_lst, ElemOrder_ptcl, FourierOrder_ptcl, X_ptcl, R_ptcl, OrientVec_ptcl, comm);
    Vector<double> Xc_ptcl;
    GetXc(Xc_ptcl, elem_lst_ptcl, obj_elem_cnt_ptcl, obj_elem_dsp_ptcl, comm);
    const auto build_elem_lst_nbr_ptcl = [](const sctl::Long Nelem, const sctl::Long ElemOrder, const sctl::Long FourierOrder, const sctl::Integer nbr_range){
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
                            //TODO: particle Xc and orientation.
                        }
                    }
                }
            }
        } 
        sctl::SlenderElemList<Real> elem_lst(ElemOrderVec, FourierOrderVec, Xc, eps, orient);
        return elem_lst;
    };
    const auto elem_lst0_ptcl = build_elem_lst_nbr_ptcl(Nelem, ElemOrder, FourierOrder, 0); // geometry in the unit box [0,1]^3
    const auto elem_lst_nbr_ptcl = build_elem_lst_nbr_ptcl(Nelem, ElemOrder, FourierOrder, 1); // geometry with one set of images in each direction

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    const auto X_proxy = Periodize<Real>::GetProxySurf(); // proxy points coordinates
  
    sctl::BoundaryIntegralOp<Real, KerType> LayerPotenOp0(stokes_ker); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst_nbr);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
  
    sctl::BoundaryIntegralOp<Real, KerType> LayerPotenOp_proxy(stokes_ker); // potential from elem_lst0 to proxy points
    LayerPotenOp_proxy.AddElemList(elem_lst0);
    LayerPotenOp_proxy.SetTargetCoord(X_proxy);
    LayerPotenOp_proxy.SetAccuracy(tol);
  
    // periodized layer potential operator
    const auto BIO = [&LayerPotenOp0,&LayerPotenOp_proxy,&X0,&Nrepeat](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
      const sctl::Long N = sigma.Dim();
  
      // std::cout << "Object input dim = " << LayerPotenOp0.Dim(1) << ", output = " << LayerPotenOp0.Dim(0) <<"; Proxy input dim = " << LayerPotenOp_proxy.Dim(1) << ", output = " << LayerPotenOp_proxy.Dim(0) << std::endl;
  
      sctl::Vector<Real> sigma_nbr(Nrepeat*N); // repeat sigma Nrepeat times
      for (sctl::Long k = 0; k < Nrepeat; k++) {
        for (sctl::Long i = 0; i < N; i++) {
          sigma_nbr[k*N+i] = sigma[i];
        }
      }
  
      U->SetZero();
      LayerPotenOp0.ComputePotential(*U, sigma_nbr);
      if (U->Dim() == N && std::is_same<KerType,Stokes3D_CF>::value) (*U) -= sigma*0.5; // for double-layer
  
      { // Add far-field
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
    // std::cout << "done solving for sigma" << std::endl;
    // elem_lst0.WriteVTK("vis/bgflow1",bg_flow(X0));
    // elem_lst0.WriteVTK("vis/sigma", sigma);
  
    { // Evaluate in interior, compute error and write visualization
      // NEW: use densest grid on surface to get same target points for comparison.
      const sctl::Long Nelem_trg = 8;
      const sctl::Long FourierOrder_trg = 16;
      const auto elem_lst_trg = build_elem_lst_nbr(Nelem_trg, ElemOrder, FourierOrder_trg, 0); 
  
      VolumeVis<Real> cube(elem_lst_trg, comm);
      // VolumeVis<Real> cube(elem_lst0, comm);
      X0 = cube.GetCoord(); // set new target coordinates
      LayerPotenOp0.SetTargetCoord(X0);
      sctl::Vector<Real> U;
      BIO(&U, sigma);
      U += bg_flow(X0);
  
      // cube.WriteVTK("vis/U_ref", U);
      // cube.WriteVTK("vis/err", err);
  
      // if this is the reference parameters, write to file.
      if (write_ref) {
        // Open a file for writing
        std::ofstream outFile("out/U_"+std::to_string(Nelem_trg)+"_"+std::to_string(FourierOrder_trg)+".txt");
  
        // Check if the file opened successfully
        if (!outFile) {
            std::cerr << "Error: Could not open U.txt for writing." << std::endl;
        } else {
          // Write the contents of the vector to the file
          for (Real value : U) {
            outFile << value << std::endl; // each value on a new line
          }
          // Close the file
          outFile.close();
        }
        cube.WriteVTK("vis/U_ref", U);
        elem_lst_nbr.WriteVTK("vis/S-nbr");
      } else {
        // std::cout << "not write-ref, reading U_ref" << std::endl;
        sctl::Vector<double> U_ref = Read_u_ref<double>();
        // error code for self-convergence.
        double max_err = 0;
        const auto err = U - U_ref;
        for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
        std::cout<<"Max error = "<<max_err<<'\n';
      }
  
      // // error code when straight pipe with exact solutions.
      // Real max_err = 0;
      // const auto err = U - u_ref(X0);
      // for (const auto e : err) max_err = std::max<Real>(max_err, sctl::fabs(e));
      // std::cout<<"Max error = "<<max_err<<'\n';
  
    }
}
  
