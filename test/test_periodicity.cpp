#include "periodize.hpp"
#include "utils.hpp"
#include "planeNaive.hpp"

// Test script for calculation and timing of 1, 2, and 3 periodic problems with background pressure flow. 

/**
 * Background flow with unit pressure gradient along X-axis.
 */
template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
    const Real dpdx = -1;
    const sctl::Long N = X.Dim()/3;
    sctl::Vector<Real> U(N*3);
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        U[i*3+0] = dpdx * ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
        U[i*3+1] = 0;
        U[i*3+2] = 0;
    }
    return U;
}

template <class Real> sctl::Vector<Real> bg_pres_flow(const sctl::Vector<Real>& X) {
    const sctl::Long N = X.Dim()/3;
    sctl::Vector<Real> U(N*3);
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        U[i*3+0] = - 0.5 * ((x[2]-0.5)*(x[2]-0.5)); // 2-periodic flow between plates.
        U[i*3+1] = 0.;
        // U[i*3+0] = -0.5*((x[2]-0.01)*(0.99-x[2]))/0.2401;
        // U[i*3+1] = -0.8*((x[2]-0.01)*(0.99-x[2]))/0.2401; // devide by z dir max when no scaling in x,y dir.
        U[i*3+2] = 0.;
    }
    return U;
}

// Uniform background flow in x direction.
template <class Real> sctl::Vector<Real> bg_unif_flow(const sctl::Vector<Real>& X) {
    sctl::Vector<Real> U = X;
    const sctl::Long N = X.Dim() /3;
    // U = 1.; // background flow diagonal to avoid planes of unaffected flows between periods.
    for (sctl::Long i = 0; i < N; i++) {
        U[i*3+0] = 1.; 
        U[i*3+1] = 0.; 
        U[i*3+2] = 0.; 
    }
    return U;
}

template <class Real> sctl::Vector<sctl::Vector<Real>> get_rot_mat(const sctl::Vector<Real> Xc) {
    sctl::Vector<Real> center;
    center = {0.5,0.5,0.5};
    sctl::Vector<Real> r1 = Xc - center;
    // std::cout << "Xc - center = r1 = " << r1[0] << ", " << r1[1] << ", " << r1[2] << std::endl;
    Real r1norm = r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2];
    sctl::Vector<Real> r2, r3;
    if (r1norm > 1e-5) {
        // std::cout << "center not at (0.5,0.5,0.5)." << std::endl;
        r2 = {r1[1], -r1[0], 0.};
        Real r2norm = r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2];
        r2 = r2 / sctl::sqrt<Real>(r2norm);
        r1 = r1 / sctl::sqrt<Real>(r1norm);
        // std::cout << "new r1 = "<< r1[0] << ", " << r1[1] << ", " << r1[2] << std::endl;
        r3 = { \
            r1[1]*r2[2] - r1[2]*r2[1], \
            -r1[0]*r2[2] + r1[2]*r2[0], \
            r1[0]*r2[1] - r1[1]*r2[0]
        };
        r3 = -r3;
        // std::cout << "r2 = "<< r2[0] << ", " << r2[1] << ", " << r2[2] << std::endl;
        // std::cout << "r3 = "<< r3[0] << ", " << r3[1] << ", " << r3[2] << std::endl;
    } else {
        r1 = {1.,0.,0.};
        r2 = {0.,1.,0.};
        r3 = {0.,0.,1.};
    }
    
    sctl::Vector<sctl::Vector<Real>> R;
    // R = {r1,r2,r3}; // NOTE: R = [ -r1T- ; -r2T- ; -r3T- ], actually the COB from standard to new basis.
    R = {r1, r3, r2}; // Same order as x-y-z.
    return R;
}

template <class Real> Real get_rot_mat_direction(const sctl::Vector<Real> Ftot, sctl::Vector<sctl::Vector<Real>>* R) {
    Real utilde = Ftot[0]*Ftot[0] + Ftot[1]*Ftot[1] + Ftot[2]*Ftot[2]; // technically radius * translational_velocity.
    sctl::Vector<Real> r1, r2, r3; // unit vectors for rotation matrix
    r1 = {-Ftot[1], Ftot[0], 0.}; // normal to Ftot, but rotated cw instead of ccw.
    Real r1norm = r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2];
    r2 = Ftot / sctl::sqrt<Real>(utilde);
    r1 = r1 / sctl::sqrt<Real>(r1norm);
    // std::cout << "new r1 = "<< r1[0] << ", " << r1[1] << ", " << r1[2] << std::endl;
    r3 = { \
        r1[1]*r2[2] - r1[2]*r2[1], \
        -r1[0]*r2[2] + r1[2]*r2[0], \
        r1[0]*r2[1] - r1[1]*r2[0]
    };
    r3 = -r3;
    // std::cout << "r2 = "<< r2[0] << ", " << r2[1] << ", " << r2[2] <<"); norm of Ftot = " << utilde << std::endl;
    // std::cout << "r3 = "<< r3[0] << ", " << r3[1] << ", " << r3[2] << std::endl;
    
    // R = {r1,r2,r3}; // NOTE: R = [ -r1T- ; -r2T- ; -r3T- ], actually the COB from standard to new basis.
    (*R) = {r1, r3, r2}; // Same order as x-y-z.
    
    return utilde;
}

template <class Real> sctl::Vector<Real> vslip(const sctl::Vector<Real> Xtrg, const sctl::Vector<Real> Xc, const Real r) {
    sctl::Long Ntrg = Xtrg.Dim()/3;
    sctl::Vector<Real> Utrg(Xtrg.Dim());
    sctl::Vector<sctl::Vector<Real>> R = get_rot_mat(Xc);
    auto COB = [&R](sctl::Vector<Real> v, bool RT) {
        sctl::Vector<Real> vr(3);
        if (RT) { // If using R transposed
            vr[0] = R[0][0] * v[0] + R[1][0] * v[1] + R[2][0] * v[2];
            vr[1] = R[0][1] * v[0] + R[1][1] * v[1] + R[2][1] * v[2];
            vr[2] = R[0][2] * v[0] + R[1][2] * v[1] + R[2][2] * v[2];
        } else {
            vr[0] = R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2];
            vr[1] = R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2];
            vr[2] = R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2];
        }
        return vr;
    };
    for (sctl::Long i=0; i < Ntrg; i++) {
        sctl::Vector<Real> Xtrg_here(3,(sctl::Iterator<Real>)Xtrg.begin()+i*3,true);
        sctl::Vector<Real> XtoXc_here = Xtrg_here - Xc;
        sctl::Vector<Real> Xtrg_rot = COB(XtoXc_here, false);
        // std::cout << "Xtrg before rot: (" << Xtrg_here[0]<<","<<Xtrg_here[1]<<","<<Xtrg_here[2] <<"), vector from center: (" << XtoXc_here[0]<<","<<XtoXc_here[1]<<","<<XtoXc_here[2] << "); after rot = ("<< Xtrg_rot[0]<<","<<Xtrg_rot[1]<<","<<Xtrg_rot[2]<<")." <<std::endl;
        Real phi = sctl::atan2<Real>(Xtrg_rot[1],Xtrg_rot[0]);
        Real theta = sctl::acos<Real>((Xtrg_rot[2])/r);
        // std::cout << "angles in body frame: phi = " << phi << ", theta = " << theta << std::endl;
        sctl::Vector<Real> vslip_here;
        vslip_here = { \
            - sctl::sin<Real>(theta) * sctl::cos<Real>(theta) * sctl::cos<Real>(phi), \
            - sctl::sin<Real>(theta) * sctl::cos<Real>(theta) * sctl::sin<Real>(phi), \
            sctl::sin<Real>(theta) * sctl::sin<Real>(theta) 
        };
        sctl::Vector<Real> Utrg_here = COB(vslip_here,true);
        // std::cout << "vslip in body frame: (" << vslip_here[0]<<","<<vslip_here[1]<<","<<vslip_here[2] <<"), in lab frame = ("<< Utrg_here[0]<<","<<Utrg_here[1]<<","<<Utrg_here[2]<<")." <<std::endl;
        Utrg[i*3+0] = Utrg_here[0];
        Utrg[i*3+1] = Utrg_here[1];
        Utrg[i*3+2] = Utrg_here[2];
    }
    return Utrg;
}

template <class Real> sctl::Vector<Real> vslip_direction(const sctl::Vector<Real> Xtrg, const sctl::Vector<Real> Xc, const Real r, const sctl::Vector<Real> Ftot) {
    sctl::Long Ntrg = Xtrg.Dim()/3;
    sctl::Vector<Real> Utrg(Xtrg.Dim());
    sctl::Vector<sctl::Vector<Real>> R;
    // std::cout << "Ftot in vslip direction: " << Ftot[0] << ", " << Ftot[1] << ", " << Ftot[2] << std::endl;
    Real utilde = get_rot_mat_direction(-Ftot, &R); 
    auto COB = [&R](sctl::Vector<Real> v, bool RT) {
        sctl::Vector<Real> vr(3);
        if (RT) { // If using R transposed
            vr[0] = R[0][0] * v[0] + R[1][0] * v[1] + R[2][0] * v[2];
            vr[1] = R[0][1] * v[0] + R[1][1] * v[1] + R[2][1] * v[2];
            vr[2] = R[0][2] * v[0] + R[1][2] * v[1] + R[2][2] * v[2];
        } else {
            vr[0] = R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2];
            vr[1] = R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2];
            vr[2] = R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2];
        }
        return vr;
    };
    for (sctl::Long i=0; i < Ntrg; i++) {
        sctl::Vector<Real> Xtrg_here(3,(sctl::Iterator<Real>)Xtrg.begin()+i*3,true);
        sctl::Vector<Real> XtoXc_here = Xtrg_here - Xc;
        sctl::Vector<Real> Xtrg_rot = COB(XtoXc_here, false);
        // std::cout << "Xtrg before rot: (" << Xtrg_here[0]<<","<<Xtrg_here[1]<<","<<Xtrg_here[2] <<"), vector from center: (" << XtoXc_here[0]<<","<<XtoXc_here[1]<<","<<XtoXc_here[2] << "); after rot = ("<< Xtrg_rot[0]<<","<<Xtrg_rot[1]<<","<<Xtrg_rot[2]<<")." <<std::endl;
        Real phi = sctl::atan2<Real>(Xtrg_rot[1],Xtrg_rot[0]);
        Real theta = sctl::acos<Real>((Xtrg_rot[2])/r);
        // std::cout << "angles in body frame: phi = " << phi << ", theta = " << theta << std::endl;
        sctl::Vector<Real> vslip_here;
        vslip_here = { \
            - sctl::sin<Real>(theta) * sctl::cos<Real>(theta) * sctl::cos<Real>(phi), \
            - sctl::sin<Real>(theta) * sctl::cos<Real>(theta) * sctl::sin<Real>(phi), \
            sctl::sin<Real>(theta) * sctl::sin<Real>(theta) 
        };
        vslip_here = sctl::sqrt<Real>(utilde) / r * vslip_here; // rescaled by required translational velocity
        sctl::Vector<Real> Utrg_here = COB(vslip_here,true);
        // std::cout << "vslip in body frame: (" << vslip_here[0]<<","<<vslip_here[1]<<","<<vslip_here[2] <<"), in lab frame = ("<< Utrg_here[0]<<","<<Utrg_here[1]<<","<<Utrg_here[2]<<")." <<std::endl;
        Utrg[i*3+0] = Utrg_here[0];
        Utrg[i*3+1] = Utrg_here[1];
        Utrg[i*3+2] = Utrg_here[2];
    }
    return Utrg;
}


template <class Real> sctl::Vector<Real> total_vslip(const sctl::Vector<Real> X0, const sctl::Long ptcl_gridsize, const sctl::Long Nptcl, const sctl::Vector<Real> ptcls_Xcs, const sctl::Vector<Real> ptcls_rs) {
    sctl::Vector<Real> Uslip(X0.Dim());
    // Set up for total force calculation
    sctl::Vector<Real> Ftot(3);
    Ftot = 0.;
    sctl::Vector<Real> center;
    center = {0.5,0.5,0.5};
    sctl::Vector<Real> r1, r2;
    /////////////////////////
    for (int ptcl_ind = 0; ptcl_ind < Nptcl-1; ptcl_ind++) {
        sctl::Vector<Real> Xtrg_here(ptcl_gridsize, (sctl::Iterator<Real>) X0.begin()+ptcl_ind * ptcl_gridsize,true);
        sctl::Vector<Real> Xc_here(3, (sctl::Iterator<Real>) ptcls_Xcs.begin()+ptcl_ind * 3,true);
        Real r_here = ptcls_rs[ptcl_ind];
        sctl::Vector<Real> Uslip_here = vslip(Xtrg_here, Xc_here, r_here);
        for (int i=0; i<ptcl_gridsize; i++) {
            Uslip[ptcl_ind*ptcl_gridsize + i] = Uslip_here[i];
        }
        // Sum up total force by Stokes drag law from this particle going U=1 velocity in angular direction
        r1 = Xc_here - center;
        Real r1norm = r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2];
        if (r1norm > 1e-5) {
            r2 = {r1[1], -r1[0], 0.};
            Real r2norm = r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2];
            r2 = r2 / sctl::sqrt<Real>(r2norm);
        } else {
            r2 = {0.,0.,1.}; // TODO: check this.
        }
        // std::cout << "ptcl ind = " << ptcl_ind << ", Ftot r2 = (" << r2[0] << ", " << r2[1] << ", " << r2[2] << ")." << std::endl;
        Ftot += r_here * r2;
    }
    // Set vslip on last sphere such that total force is zero in a periodic box.
    sctl::Vector<Real> Xtrg_here(ptcl_gridsize, (sctl::Iterator<Real>) X0.begin()+(Nptcl-1) * ptcl_gridsize,true);
    sctl::Vector<Real> Xc_here(3, (sctl::Iterator<Real>) ptcls_Xcs.begin()+(Nptcl-1) * 3,true);
    Real r_here = ptcls_rs[Nptcl-1];
    sctl::Vector<Real> Uslip_here = vslip_direction(Xtrg_here, Xc_here, r_here, Ftot);
    for (int i=0; i<ptcl_gridsize; i++) {
        Uslip[(Nptcl-1)*ptcl_gridsize + i] = Uslip_here[i];
    }
    // // CHECK that total F adds up to (0,0,0).
    // sctl::Vector<sctl::Vector<Real>> Rtemp;
    // Real utilde = get_rot_mat_direction(-Ftot, &Rtemp); 
    // std::cout << "F before last particle: (" << Ftot[0] << ", " << Ftot[1] << ", " << Ftot[2] <<");" <<std::endl;
    // std::cout << "r2 from inside get_rot_mat is (" << Rtemp[2][0] << ", " << Rtemp[2][1] << ", " << Rtemp[2][2] << ")" << std::endl;
    // Ftot = Ftot + sctl::sqrt<Real>(utilde)*Rtemp[2];
    // std::cout << "F after adding Uslip of last particle: (" << Ftot[0] << ", " << Ftot[1] << ", " << Ftot[2] <<");" <<std::endl;
    
    return Uslip;
}

template <class Real> void plot_setup(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Integer peri_mode, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, const Real gmres_tol, const Real tol) {

    // // Combine single-layer and double-layer kernels in these proportions
    // const Real SL_scal = 1.0;
    // const Real DL_scal = 1.0;

    const sctl::Long ElemOrder = 10;
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0;
    sctl::Vector<Real> NormalOrient;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else { 
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 0, 1, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    }
    // const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); 
    Nptcl = ptcls_rs.Dim(); 
    // std::cout << "periodic mode is " << peri_mode << ", Nrepeat is " << Nrepeat << std::endl;

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Long ptcl_gridsize = Nelem * ElemOrder * FourierOrder * 3; 
    sctl::Long Nptcl_slip = elem_lst0.Size() / Nelem; 
    sctl::Vector<Real> ptcls_Xcs_slip(Nptcl_slip * 3, (sctl::Iterator<Real>)ptcls_Xcs.begin() + comm.Rank()*Nptcl_slip*3, true);
    sctl::Vector<Real> ptcls_rs_slip(Nptcl_slip,  (sctl::Iterator<Real>)ptcls_rs.begin()+comm.Rank()*Nptcl_slip, true);
    sctl::Vector<Real> Uslip = total_vslip(X0, ptcl_gridsize, Nptcl_slip, ptcls_Xcs_slip, ptcls_rs_slip);
    if (write_ref) {
        elem_lst0.WriteVTK("vis/"+std::to_string(Nptcl)+"spheres",Uslip,comm);
    }  
}

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

template <class Real> void test(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Integer peri_mode, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, const Real gmres_tol, const Real tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real pressure_drop = -1.0;
    const Real period_length = 1;

    const sctl::Long ElemOrder = 10;
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else { 
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 0, peri_mode, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    }
    Nptcl = ptcls_rs.Dim(); 

    // sctl::Vector<Real> X0; // target coordinates
    // elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    sctl::Vector<Real> X0_ptcl; // target coordinates
    elem_lst0.GetNodeCoord(&X0_ptcl, nullptr, nullptr);
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
        // surface_area = surface_area_[0];
    }

    // sctl::Long ptcl_gridsize = Nelem * ElemOrder * FourierOrder * 3; 
    // sctl::Long Nptcl_slip = elem_lst0.Size() / Nelem; 
    // sctl::Vector<Real> ptcls_Xcs_slip(Nptcl_slip * 3, (sctl::Iterator<Real>)ptcls_Xcs.begin() + comm.Rank()*Nptcl_slip*3, true);
    // sctl::Vector<Real> ptcls_rs_slip(Nptcl_slip,  (sctl::Iterator<Real>)ptcls_rs.begin()+comm.Rank()*Nptcl_slip, true);
    // sctl::Vector<Real> Uslip = total_vslip(X0_ptcl, ptcl_gridsize, Nptcl_slip, ptcls_Xcs_slip, ptcls_rs_slip);
    // if (write_ref) {
    //     elem_lst0.WriteVTK("vis/"+std::to_string(Nptcl)+"spheres_streamline_vslip/"+std::to_string(Nptcl)+"spheres",Uslip,comm);
    // }  

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0,"1");

    sctl::Vector<Real> X0;
    StokesBIO LayerPotenOp1(0., DL_scal, comm); 
    
    // NEW: Add plane elements
    sctl::Long gl_order = 40;
    sctl::Long Nelem_x = 2;
    sctl::Long Nelem_y = 2;
    sctl::PlaneIntegral<Real> plane(gl_order, Nelem_x, Nelem_y, 0.01);
    sctl::Vector<Real> X0_wall;
    plane.GetNodeCoord(&X0_wall, nullptr, nullptr);
    // LayerPotenOp0.AddElemList(plane,"2",false,true); 
    // DEBUGGING: Separate Op for plane and particle.
    LayerPotenOp1.AddElemList(plane);
    LayerPotenOp1.SetTargetCoord(X0_wall);
    LayerPotenOp1.SetAccuracy(tol);
    LayerPotenOp1.SetPeriodicity(sctl::Periodicity::XY, 1.0);

    // Add plane nodes as target nodes after on-particle nodes
    if (peri_mode==2) {
        X0.ReInit(X0_ptcl.Dim() + X0_wall.Dim());
        for (int j=0; j<X0_ptcl.Dim(); j++) {
            X0[j] = X0_ptcl[j];
        }
        for (int j=0; j<X0_wall.Dim(); j++) {
            X0[j+X0_ptcl.Dim()] = X0_wall[j];
        }
    } else {
        X0 = X0_ptcl;
    }
    
    // // Add plane normal orient as well 
    // sctl::Vector<Real> NormalOrient_(NormalOrient.Dim() + X0_wall.Dim());
    // NormalOrient_ = -1.; // Normal orient = -1 (-sign below) means all normals point into fluid (exterior problem)
    // NormalOrient_.Swap(NormalOrient);
    
    // LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    if (peri_mode==1) {
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);
    } else if (peri_mode==3) {
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length);
    } else if (peri_mode==2) {
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XY, period_length);
    } else {
        SCTL_ASSERT(false);
    }

    // // periodized layer potential operator
    // const auto BIO = [DL_scal,&LayerPotenOp0,NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    //     U->SetZero();
    //     LayerPotenOp0.ComputePotential(*U, sigma);
    //     if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma*0.5*NormalOrient * DL_scal; // for double-layer
    // };
    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,&NormalOrient, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> sigma_mean, sigma0;
        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_;
            elem_lst0.GetFarFieldDensity(sigma_, sigma);
            SurfaceIntegral(sigma_mean, sigma_, wts);
            // MPI
            sctl::Vector<Real> sa_loc(1);
            sa_loc[0] = sigma_mean[0];
            sctl::Vector<Real> sa_all(1);
            sa_all[0] = 0;
            comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
            // sigma_mean *= (1/surface_area);
            sigma_mean *= (1/sa_all[0]);

            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    const auto BIO_solve = [&LayerPotenOp0,&LayerPotenOp1,&X0_ptcl,&X0_wall,&tol](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        U->SetZero();
        sctl::Vector<Real> sig_ptcl(X0_ptcl.Dim(), (sctl::Iterator<Real>)sigma.begin(),true);
        sctl::Vector<Real> sig_wall(X0_wall.Dim(), (sctl::Iterator<Real>)sigma.begin()+X0_ptcl.Dim(),true);

        sctl::Vector<Real> Uptcl_self,Uptcl_wall,Uwall_self,Uwall_ptcl;
        LayerPotenOp0.SetTargetCoord(X0_ptcl);
        LayerPotenOp0.ComputePotential(Uptcl_self, sig_ptcl);
        Uptcl_self += 0.5*sig_ptcl; // exterior problem to particle -> PV + 1/2*density
        LayerPotenOp0.SetTargetCoord(X0_wall);
        LayerPotenOp0.ComputePotential(Uptcl_wall, sig_ptcl);

        LayerPotenOp1.SetTargetCoord(X0_wall);
        LayerPotenOp1.ComputePotential(Uwall_self, sig_wall);
        Uwall_self += 0.5*sig_wall; // exterior problem to walls -> PV + 1/2*density
        LayerPotenOp1.SetTargetCoord(X0_ptcl);
        LayerPotenOp1.ComputePotential(Uwall_ptcl, sig_wall);

        for (sctl::Long ind=0; ind<sigma.Dim(); ind++) {
            if (ind < X0_ptcl.Dim()) { // target on particle
                (*U)[ind] = Uptcl_self[ind] + Uwall_ptcl[ind];
            } else { // target on wall
                sctl::Long ind_2 = ind - X0_ptcl.Dim();
                (*U)[ind] = Uptcl_wall[ind_2] + Uwall_self[ind_2];
            }
        }

    };

    const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
        sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
        AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

        sctl::Vector<Real> U0;
        LayerPotenOp0.ComputeSL(U0, force_density);
        return U0;
    };

    // first gmres to remove timing for matrix loading, and set Krylov preconditioner.
    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    if (peri_mode == 2) {
        solver(&sigma,BIO_solve, -bg_pres_flow(X0), gmres_tol);
    } else if (peri_mode == 3) {
        LayerPotenOp0.SetTargetCoord(X0);
        // solver(&sigma,BIO, Uslip, gmres_tol);
        solver(&sigma, BIO, eval_rhs(pressure_drop), gmres_tol);
    } else {
        LayerPotenOp0.SetTargetCoord(X0);
        solver(&sigma, BIO, bg_flow(X0) * (pressure_drop/period_length), gmres_tol);
    }  

    if (write_ref) { 
        if (peri_mode==2) {
            sctl::Vector<Real> ptcl_dens(X0_ptcl.Dim(), (sctl::Iterator<Real>) sigma.begin(), true);
            sctl::Vector<Real> wall_dens(X0_wall.Dim(), (sctl::Iterator<Real>) sigma.begin()+ptcl_dens.Dim(), true);
            elem_lst0.WriteVTK("vis/2peri_"+std::to_string(Nptcl)+"spheres",ptcl_dens,comm);
            plane.WriteVTK("vis/2peri_plane_density",wall_dens,comm);
        } else {
            elem_lst0.WriteVTK("vis/"+std::to_string(peri_mode)+"peri_"+std::to_string(Nptcl)+"spheres",sigma,comm);
        }
        sctl::Long Ntrg_side = 5;
        Real gap = 1./(Ntrg_side+1);
        X0.ReInit(Ntrg_side * Ntrg_side * 3);
        // X symmetry
        sctl::Vector<Real> X1(Ntrg_side * Ntrg_side * 3);
        for (sctl::Long yind=0; yind<Ntrg_side; yind++) {
            for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
                sctl::Long Ntrg_nodeind = yind*Ntrg_side + zind;
                X0[Ntrg_nodeind * 3 + 0] = 0.;
                X0[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
                X0[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
                X1[Ntrg_nodeind * 3 + 0] = 1.;
                X1[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
                X1[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
            }
        }
        const auto BIO_eval = [&LayerPotenOp0,&LayerPotenOp1,&X0_ptcl,&X0_wall](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
            U->SetZero();
            sctl::Vector<Real> sig_ptcl(X0_ptcl.Dim(), (sctl::Iterator<Real>)sigma.begin(),true);
            sctl::Vector<Real> sig_wall(X0_wall.Dim(), (sctl::Iterator<Real>)sigma.begin()+X0_ptcl.Dim(),true);
            sctl::Vector<Real> Uptcl,Uwall;
            LayerPotenOp0.ComputePotential(Uptcl, sig_ptcl);
            LayerPotenOp1.ComputePotential(Uwall, sig_wall);
            (*U) = Uptcl + Uwall;
        };
        sctl::Vector<Real> UX0(X0.Dim());
        LayerPotenOp0.SetTargetCoord(X0);
        if (peri_mode == 2) {
            LayerPotenOp1.SetTargetCoord(X0);
            BIO_eval(&UX0, sigma); 
        } else {
            BIO(&UX0,sigma);
        }
        sctl::Vector<Real> UX1(X1.Dim());
        LayerPotenOp0.SetTargetCoord(X1);
        if (peri_mode == 2) {
            LayerPotenOp1.SetTargetCoord(X1);
            BIO_eval(&UX1, sigma); 
        } else {
            BIO(&UX1,sigma);
        }
        if (peri_mode==1) {
            UX0 -= bg_flow(X0) * (pressure_drop/period_length);
            UX1 -= bg_flow(X1) * (pressure_drop/period_length);
        } else if (peri_mode==3) {
            UX0 -= eval_rhs(pressure_drop);
            UX1 -= eval_rhs(pressure_drop);
        }
        std::cout << "============ X periodicity =================" << std::endl;
        sctl::Vector<Real> UdiffX = UX0-UX1;
        for (int i=0; i<UdiffX.Dim()/3; i++) {
            std::cout << std::setprecision(8) << UdiffX[i*3+0]<< ", " << UdiffX[i*3+1]<< ", " << UdiffX[i*3+2]<< ". " << std::endl;
        }

        if (peri_mode > 1) {
            // Y symmetry
            sctl::Vector<Real> Y0(Ntrg_side * Ntrg_side * 3);
            sctl::Vector<Real> Y1(Ntrg_side * Ntrg_side * 3);
            for (sctl::Long xind=0; xind<Ntrg_side; xind++) {
                for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
                    sctl::Long Ntrg_nodeind = xind*Ntrg_side + zind;
                    Y0[Ntrg_nodeind * 3 + 0] = (xind+1)*gap;
                    Y0[Ntrg_nodeind * 3 + 1] = 0.0;
                    Y0[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
                    Y1[Ntrg_nodeind * 3 + 0] = (xind+1)*gap;
                    Y1[Ntrg_nodeind * 3 + 1] = 1.;
                    Y1[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
                }
            }
            sctl::Vector<Real> UY0(Y0.Dim());
            LayerPotenOp0.SetTargetCoord(Y0);
            if (peri_mode == 2) {
                LayerPotenOp1.SetTargetCoord(Y0);
                BIO_eval(&UY0, sigma); 
            } else {
                BIO(&UY0,sigma);
            }
            sctl::Vector<Real> UY1(Y1.Dim());
            LayerPotenOp0.SetTargetCoord(Y1);
            if (peri_mode == 2) {
                LayerPotenOp1.SetTargetCoord(Y1);
                BIO_eval(&UY1, sigma); 
            } else {
                BIO(&UY1,sigma);
            }
            std::cout << "============ Y periodicity =================" << std::endl;
            sctl::Vector<Real> UdiffY = UY0-UY1;
            for (int i=0; i<UdiffY.Dim()/3; i++) {
                std::cout << std::setprecision(8) << UdiffY[i*3+0]<< ", " << UdiffY[i*3+1]<< ", " << UdiffY[i*3+2]<< ". " << std::endl;
            }
        }
        if (peri_mode > 2) {
            // Z symmetry
            sctl::Vector<Real> Z0(Ntrg_side * Ntrg_side * 3);
            sctl::Vector<Real> Z1(Ntrg_side * Ntrg_side * 3);
            for (sctl::Long xind=0; xind<Ntrg_side; xind++) {
                for (sctl::Long yind=0; yind<Ntrg_side; yind++) {
                    sctl::Long Ntrg_nodeind = xind*Ntrg_side + yind;
                    Z0[Ntrg_nodeind * 3 + 0] = (xind+1)*gap;
                    Z0[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
                    Z0[Ntrg_nodeind * 3 + 2] = 0.0;
                    Z1[Ntrg_nodeind * 3 + 0] = (xind+1)*gap;
                    Z1[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
                    Z1[Ntrg_nodeind * 3 + 2] = 1.0;
                }
            }
            sctl::Vector<Real> UZ0(Z0.Dim());
            LayerPotenOp0.SetTargetCoord(Z0);
            BIO(&UZ0,sigma);
            sctl::Vector<Real> UZ1(Z1.Dim());
            LayerPotenOp0.SetTargetCoord(Z1);
            BIO(&UZ1,sigma);
            std::cout << "============ Z periodicity =================" << std::endl;
            sctl::Vector<Real> UdiffZ = UZ0-UZ1;
            for (int i=0; i<UdiffZ.Dim()/3; i++) {
                std::cout << std::setprecision(8) << UdiffZ[i*3+0]<< ", " << UdiffZ[i*3+1]<< ", " << UdiffZ[i*3+2]<< ". " << std::endl;
            }
        }

        // Streamlines
        PeriodicGeom<Real> trg;    
        CubeVolumeVisShifted<Real> vol_vis(15, 0.92, comm);
        // X0 = vol_vis.GetCoord();
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);
        std::cout << "number of target points: " << X0.Dim() << std::endl;

        sctl::Vector<Real> Utrg;
        LayerPotenOp0.SetTargetCoord(X0);
        if (peri_mode == 2) {
            LayerPotenOp1.SetTargetCoord(X0);
            Utrg.ReInit(X0.Dim());
            BIO_eval(&Utrg, sigma); 
            Utrg += bg_pres_flow(X0); 
        } else if (peri_mode == 3) {
            BIO(&Utrg,sigma);
        } else {
            BIO(&Utrg,sigma);
            Utrg += bg_unif_flow(X0); 
        }
        
        sctl::Vector<Real> U_vis(X0_all.Dim());
        U_vis = 0.;
        sctl::Long X1_ptr = 0;
        for (sctl::Long i=0; i<X0_all.Dim()/3; i++) {
            if (filtered_inds[i] == 0) {
                U_vis[i*3] = Utrg[X1_ptr*3];
                U_vis[i*3+1] = Utrg[X1_ptr*3+1];
                U_vis[i*3+2] = Utrg[X1_ptr*3+2];
                X1_ptr += 1;
            }
        }
        // vol_vis.WriteVTK("vis/"+std::to_string(Nptcl)+"spheres_streamline_vslip/"+std::to_string(Nptcl)+"streamlines", U_vis); 
        vol_vis.WriteVTK("vis/"+std::to_string(peri_mode)+"peri_"+std::to_string(Nptcl)+"streamlines", U_vis); 
    } 
}

template <class Real> void test1peri(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, const Real gmres_tol, const Real tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real pressure_drop = -1.0;
    const Real period_length = 1;

    const sctl::Long ElemOrder = 10;
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else { 
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 0, 1, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    }
    Nptcl = ptcls_rs.Dim(); 

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    Real surface_area;
    sctl::Vector<Real> wts;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        surface_area = surface_area_[0];
    }

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, period_length);
    

    // // periodized layer potential operator
    // const auto BIO = [DL_scal,&LayerPotenOp0,NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    //     U->SetZero();
    //     LayerPotenOp0.ComputePotential(*U, sigma);
    //     if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma*0.5*NormalOrient * DL_scal; // for double-layer
    // };

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,&NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> sigma_mean, sigma0;
        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_;
            elem_lst0.GetFarFieldDensity(sigma_, sigma);
            SurfaceIntegral(sigma_mean, sigma_, wts);
            sigma_mean *= (1/surface_area);

            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);

            // DEBUG: check that sigma-sigma_mean has surface integral = 0:
            sctl::Vector<Real> sigma1 = sigma_;
            AddConstVec(sigma1, -sigma_mean);
            sctl::Vector<Real> sigma_test_;
            SurfaceIntegral(sigma_test_, sigma1, wts);
            std::cout << "Surface integral of sigma - sigma bar = " << sigma_test_[0] << ", "<< sigma_test_[1] << ", " << sigma_test_[2] << ". "<< std::endl;
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    LayerPotenOp0.SetTargetCoord(X0);
    solver(&sigma, BIO, bg_pres_flow(X0) * (pressure_drop/period_length), gmres_tol); 

    // const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
    //   sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
    //   AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

    //   sctl::Vector<Real> U0;
    //   LayerPotenOp0.ComputeSL(U0, force_density);
    //   return U0;
    // };
    // sctl::Vector<Real> sigma_up;
    // solver(&sigma_up, BIO, eval_rhs(pressure_drop), gmres_tol); 

    sctl::Long Ntrg_side = 5;
    Real gap = 1./(Ntrg_side+1);
    X0.ReInit(Ntrg_side * Ntrg_side * 3);
    // X symmetry
    sctl::Vector<Real> X1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long yind=0; yind<Ntrg_side; yind++) {
        for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
            sctl::Long Ntrg_nodeind = yind*Ntrg_side + zind;
            X0[Ntrg_nodeind * 3 + 0] = 0.;
            X0[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
            X0[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
            X1[Ntrg_nodeind * 3 + 0] = 1.;
            X1[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
            X1[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
        }
    }
    sctl::Vector<Real> UX0(X0.Dim());
    LayerPotenOp0.SetTargetCoord(X0);
    BIO(&UX0,sigma);
    sctl::Vector<Real> UX1(X1.Dim());
    LayerPotenOp0.SetTargetCoord(X1);
    BIO(&UX1,sigma);
    UX0 -= bg_flow(X0) * (pressure_drop/period_length);
    UX1 -= bg_flow(X1) * (pressure_drop/period_length);
    std::cout << "============ X periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffX = UX0-UX1;
    for (int i=0; i<UdiffX.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffX[i*3+0]<< ", " << UdiffX[i*3+1]<< ", " << UdiffX[i*3+2]<< ". " << std::endl;
    }
}

template <class Real> void test2peri(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, const Real gmres_tol, const Real tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real pressure_drop = -10.0;
    const Real period_length = 1;

    const sctl::Long ElemOrder = 10;
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else { 
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 0, 2, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    }
    Nptcl = ptcls_rs.Dim(); 

    sctl::Vector<Real> X0_ptcl; // target coordinates
    elem_lst0.GetNodeCoord(&X0_ptcl, nullptr, nullptr);
    Real surface_area;
    sctl::Vector<Real> wts;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        // surface_area = surface_area_[0];
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area = sa_all[0];
    }
    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0,"1");

    sctl::Vector<Real> X0;
    
    // NEW: Add plane elements
    // sctl::Long gl_order = 49;
    // sctl::Long Nelem_x = 2;
    // sctl::Long Nelem_y = 2;
    sctl::Long gl_order = 20;
    sctl::Long Nelem_x = 1;
    sctl::Long Nelem_y = 1;
    sctl::PlaneIntegral<Real> plane(gl_order, Nelem_x, Nelem_y, 0.01);
    sctl::Vector<Real> X0_wall;
    plane.GetNodeCoord(&X0_wall, nullptr, nullptr);
    LayerPotenOp0.AddElemList(plane,"2",false,true); 

    // Extra SL op for plane with singularity subtraction.
    StokesBIO LayerPotenOp2(1.0, 0.0, comm);
    LayerPotenOp2.AddElemList(plane);
    LayerPotenOp2.SetAccuracy(tol);
    LayerPotenOp2.SetPeriodicity(sctl::Periodicity::XY, period_length);
    // LayerPotenOp2.SetPeriodicity(sctl::Periodicity::XYZ, period_length);

    Real surface_area_wall;
    sctl::Vector<Real> wts_wall;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        plane.GetFarFieldNodes(X, Xn, wts_wall, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts_wall*0+1, wts_wall);
        // surface_area_wall = surface_area_[0];
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area_wall = sa_all[0];
        std::cout << "DEBUG wall surface area, computed to be " << surface_area_wall << std::endl;
    }

    X0.ReInit(X0_ptcl.Dim() + X0_wall.Dim());
    for (int j=0; j<X0_ptcl.Dim(); j++) {
        X0[j] = X0_ptcl[j];
    }
    for (int j=0; j<X0_wall.Dim(); j++) {
        X0[j+X0_ptcl.Dim()] = X0_wall[j];
    }
    
    // Add plane normal orient as well 
    sctl::Vector<Real> NormalOrient_(NormalOrient.Dim() + X0_wall.Dim());
    NormalOrient_ = -1.; // Normal orient = -1 (-sign below) means all normals point into fluid (exterior problem)
    NormalOrient_.Swap(NormalOrient);
    
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XY, period_length);
    // LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length);


    // Real eps = 1e-5;

    // const auto reg_sl_func = [&eps, &wts_wall](const sctl::Vector<Real> Xself, const sctl::Vector<Real> sigma) {
    //     SCTL_ASSERT(Xself.Dim() == sigma.Dim());
    //     sctl::Long N = Xself.Dim() / 3;
    //     sctl::Vector<Real> SL_eps(N*3);
    //     Real eps2 = eps*eps;
    //     for (sctl::Long i=0; i<N; i++) { // trg idx
    //         sctl::Vector<Real> xtrg(3, (sctl::Iterator<Real>) Xself.begin() + i*3, false);
    //         for (sctl::Long j=0; j<N; j++) { // src idx
    //             sctl::Vector<Real> xsrc(3, (sctl::Iterator<Real>) Xself.begin() + j*3, false);
    //             sctl::Vector<Real> r = xtrg - xsrc;
    //             Real r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
    //             Real sqrt_r2e2 = sctl::sqrt<Real>(r2 + eps2);
    //             Real inv_r2e2 = 1./sqrt_r2e2;
    //             Real inv3_r2e2 = inv_r2e2*inv_r2e2*inv_r2e2;
    //             for (sctl::Long k1 = 0; k1 < 3; k1++) { // trg dim
    //                 for (sctl::Long k2 = 0; k2 < 3; k2++) { // src dim
    //                     SL_eps[i*3 + k1] += (k1==k2 ? ( (r2+2.*eps2)*inv3_r2e2)*sigma[j*3+k2]*wts_wall[j] : 0.) + r[k1]*r[k2]*inv3_r2e2*sigma[j*3+k2]*wts_wall[j];
    //                 }
    //             }
    //         }
    //     }       
    //     return SL_eps;
    // };

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&wts_wall,&surface_area_wall,&plane,&LayerPotenOp0,&LayerPotenOp2,&DL_scal,&X0_ptcl,&X0_wall,&NormalOrient, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        
        sctl::Vector<Real> ptcl_dens(X0_ptcl.Dim(), (sctl::Iterator<Real>) sigma.begin(), true);
        sctl::Vector<Real> wall_dens(X0_wall.Dim(), (sctl::Iterator<Real>) sigma.begin()+ptcl_dens.Dim(), true);
            
        sctl::Vector<Real> sigma_mean, sigma0;

        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_mean_ptcl, sigma_mean_wall;
            sctl::Vector<Real> sigma_ptcl_;
            elem_lst0.GetFarFieldDensity(sigma_ptcl_, ptcl_dens);
            SurfaceIntegral(sigma_mean_ptcl, sigma_ptcl_, wts);
            sctl::Vector<Real> sigma_wall_ = wall_dens;
            SurfaceIntegral(sigma_mean_wall, sigma_wall_, wts_wall);
            // MPI
            sctl::Vector<Real> sa_loc(6);
            for (int i=0; i<3; i++) {
                sa_loc[i] = sigma_mean_ptcl[i];
            }
            for (int i=0; i<3; i++) {
                sa_loc[i+3] = sigma_mean_wall[i];
            }
            sctl::Vector<Real> sa_all(6);
            sa_all = 0.;
            for (int i=0; i<3; i++) {
                comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+i, (sctl::Iterator<Real>) sa_all.begin()+i, 1, sctl::CommOp::SUM);
            }
            for (int i=0; i<3; i++) {
                comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+i+3, (sctl::Iterator<Real>) sa_all.begin()+i+3, 1, sctl::CommOp::SUM);
            }
            // sigma_mean = sigma_mean_ptcl + sigma_mean_wall; 
            sigma_mean.ReInit(sigma_mean_ptcl.Dim());
            for (int i=0; i<3; i++) {
                sigma_mean[i] = sa_all[i] + sa_all[i+3]; // adding total sigma in {x,y,z} coord separately
            }
            sigma_mean *= (1/(surface_area + surface_area_wall));
            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);

            // // DEBUG: check that sigma-sigma_mean has surface integral = 0:
            // AddConstVec(sigma_ptcl_, -sigma_mean); // note: ptcl_dens created to be copy not pointer, so this shouldn't change sigma values.
            // AddConstVec(sigma_wall_, -sigma_mean); 
            // sctl::Vector<Real> sigma_test_ptcl, sigma_test_wall;
            // SurfaceIntegral(sigma_test_ptcl, sigma_ptcl_, wts);
            // SurfaceIntegral(sigma_test_wall, sigma_wall_, wts_wall);
            // std::cout << "size of sigma_test wall is: " << sigma_test_wall.Dim() <<  "; Surface integral of sigma - sigma bar = " << sigma_test_wall[0] + sigma_test_ptcl[0] << ", "<< sigma_test_wall[1] + sigma_test_ptcl[1] << ", " << sigma_test_wall[2] + sigma_test_ptcl[2] << ". "<< std::endl;
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma0.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        // single layer on walls
        {
            sctl::Vector<Real> wall_dens_0(X0_wall.Dim(), (sctl::Iterator<Real>) sigma0.begin()+ptcl_dens.Dim(), true);
            // First compute SL[walls](ptcl)
            sctl::Vector<Real> U_sl_1(X0_ptcl.Dim());
            LayerPotenOp2.SetTargetCoord(X0_ptcl);
            LayerPotenOp2.ComputePotential(U_sl_1, wall_dens_0);
            for (int i=0; i<X0_ptcl.Dim(); i++) {
                (*U)[i] += U_sl_1[i];
            }
            sctl::Vector<Real> reg_sl_U;
            LayerPotenOp2.SetTargetCoord(X0_wall);
            LayerPotenOp2.ComputePotential(reg_sl_U, wall_dens_0);
            for (int i=X0_ptcl.Dim(); i<U->Dim(); i++) {
                (*U)[i] += reg_sl_U[i-ptcl_dens.Dim()];
            }
            // // Singularity subtraction for self-to-self on top AND bottom planes
            // sctl::Long Nentries_plane = X0_wall.Dim()/2; 
            // sctl::Long Nentries_ptcl = ptcl_dens.Dim(); 
            // sctl::Vector<Real> U2(3); 
            // // Top plane first.
            // for (int i=0; i<X0_wall.Dim()/6; i++) {
            //     sctl::Vector<Real> sigma_i(3, (sctl::Iterator<Real>) wall_dens_0.begin() + i*3, true);
            //     sctl::Vector<Real> x_i(3, (sctl::Iterator<Real>) X0_wall.begin() + i*3, true);
            //     sctl::Vector<Real> sigma_all_topsub = wall_dens_0;
            //     sctl::Vector<Real> sigma_top(Nentries_plane, (sctl::Iterator<Real>) sigma_all_topsub.begin(), false); // only sigma on first plane subtracted by sigma_i.
            //     AddConstVec(sigma_top, -sigma_i);
            //     // sctl::Vector<Real> sigma_top(Nentries_plane, (sctl::Iterator<Real>) wall_dens_0.begin(), true); // only sigma on first plane subtracted by sigma_i.
            //     // AddConstVec(sigma_top, -sigma_i);
            //     // sctl::Vector<Real> sigma_all_topsub = wall_dens_0;
            //     // for (int j=0; j<Nentries_plane; j++) {
            //     //     sigma_all_topsub[j] = sigma_top[j];
            //     // }
            //     LayerPotenOp2.SetTargetCoord(x_i);
            //     LayerPotenOp2.ComputePotential(U2, sigma_all_topsub); // SL1[sigma_1 - sigma_1(x_i)](x_i) + SL2[sigma_2](x_i), x_i on top plane, smooth
            //     // Allegedly: int_R2 G(x,y) dy = 0... 
            //     // Otherwise add to U2 then put into U
            //     for (int k=0; k<3; k++) {
            //         (*U)[i*3+k+Nentries_ptcl] += U2[k];
            //     }
            // }
            // // Bottom plane
            // for (int i=0; i<X0_wall.Dim()/6; i++) {
            //     sctl::Vector<Real> sigma_i(3, (sctl::Iterator<Real>) wall_dens_0.begin() + Nentries_plane + i*3, true);
            //     sctl::Vector<Real> x_i(3, (sctl::Iterator<Real>) X0_wall.begin() + Nentries_plane + i*3, true);
            //     sctl::Vector<Real> sigma_all_botsub = wall_dens_0;
            //     sctl::Vector<Real> sigma_bot(Nentries_plane, (sctl::Iterator<Real>) sigma_all_botsub.begin() + Nentries_plane, false); // only sigma on second plane subtracted by sigma_i.
            //     AddConstVec(sigma_bot, -sigma_i);
            //     // sctl::Vector<Real> sigma_bot(Nentries_plane, (sctl::Iterator<Real>) wall_dens_0.begin() + Nentries_plane, true); // only sigma on second plane subtracted by sigma_i.
            //     // AddConstVec(sigma_bot, -sigma_i);
            //     // sctl::Vector<Real> sigma_all_botsub = wall_dens_0;
            //     // for (int j=0; j<Nentries_plane; j++) {
            //     //     sigma_all_botsub[j + Nentries_plane] = sigma_bot[j];
            //     // }
            //     LayerPotenOp2.SetTargetCoord(x_i);
            //     LayerPotenOp2.ComputePotential(U2, sigma_all_botsub); // SL2[sigma_2 - sigma_2(x_i)](x_i) + SL1[sigma_1](x_i), x_i on bottom plane, smooth
            //     // Allegedly: int_R2 G(x,y) dy = 0... 
            //     // Otherwise add to U2 then put into U
            //     for (int k=0; k<3; k++) {
            //         (*U)[i*3+k+Nentries_plane+Nentries_ptcl] += U2[k];
            //     }
            // }
        }

        AddConstVec(*U, sigma_mean);
    };


    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO_eval = [&wts,&surface_area,&elem_lst0,&wts_wall,&surface_area_wall,&plane,&LayerPotenOp0,&LayerPotenOp2,&DL_scal,&X0_ptcl,&X0_wall,&NormalOrient, &comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        
        sctl::Vector<Real> ptcl_dens(X0_ptcl.Dim(), (sctl::Iterator<Real>) sigma.begin(), true);
        sctl::Vector<Real> wall_dens(X0_wall.Dim(), (sctl::Iterator<Real>) sigma.begin()+ptcl_dens.Dim(), true);
            
        sctl::Vector<Real> sigma_mean, sigma0;

        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_mean_ptcl, sigma_mean_wall;
            sctl::Vector<Real> sigma_ptcl_;
            elem_lst0.GetFarFieldDensity(sigma_ptcl_, ptcl_dens);
            SurfaceIntegral(sigma_mean_ptcl, sigma_ptcl_, wts);
            sctl::Vector<Real> sigma_wall_ = wall_dens;
            SurfaceIntegral(sigma_mean_wall, sigma_wall_, wts_wall);
            // MPI
            sctl::Vector<Real> sa_loc(6);
            for (int i=0; i<3; i++) {
                sa_loc[i] = sigma_mean_ptcl[i];
            }
            for (int i=0; i<3; i++) {
                sa_loc[i+3] = sigma_mean_wall[i];
            }
            sctl::Vector<Real> sa_all(6);
            sa_all = 0.;
            for (int i=0; i<3; i++) {
                comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+i, (sctl::Iterator<Real>) sa_all.begin()+i, 1, sctl::CommOp::SUM);
            }
            for (int i=0; i<3; i++) {
                comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+i+3, (sctl::Iterator<Real>) sa_all.begin()+i+3, 1, sctl::CommOp::SUM);
            }
            // sigma_mean = sigma_mean_ptcl + sigma_mean_wall; 
            sigma_mean.ReInit(sigma_mean_ptcl.Dim());
            for (int i=0; i<3; i++) {
                sigma_mean[i] = sa_all[i] + sa_all[i+3]; // adding total sigma in {x,y,z} coord separately
            }
            sigma_mean *= (1/(surface_area + surface_area_wall));
            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);

        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        sctl::Vector<Real> wall_dens_0(X0_wall.Dim(), (sctl::Iterator<Real>) sigma0.begin()+ptcl_dens.Dim(), true);
        sctl::Vector<Real> U2(U->Dim());
        LayerPotenOp2.ComputePotential(U2, wall_dens_0);
        (*U) += U2;

        AddConstVec(*U, sigma_mean);
    };

    // const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
    //     sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
    //     AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

    //     sctl::Vector<Real> U0;
    //     LayerPotenOp0.ComputeSL(U0, force_density);
    //     return U0;
    // };

    // /*
    // Getting singular values
    sctl::Long Nsrc = X0_ptcl.Dim()/3 + X0_wall.Dim()/3;
    sctl::Vector<Real> sigma_eye(Nsrc*3);
    sctl::Vector<sctl::Vector<Real>> LPOvecvec(Nsrc*3);
    for (sctl::Long i=0; i<Nsrc; i++) {
        for (sctl::Long k=0; k<3; k++) {
            sigma_eye.SetZero();
            sigma_eye[i*3+k] = 1.;
            // std::cout << "Node number is = " << i << ", dimension = " << k << std::endl;
            BIO(LPOvecvec.begin()+i*3+k, sigma_eye);
        }
    }
    // SVD
    sctl::Matrix<Real> LPOmat(Nsrc*3,Nsrc*3);
    for (long i=0; i < Nsrc*3; i++) {
        for (long j = 0; j < Nsrc*3; j++) {
            LPOmat(j,i) = LPOvecvec[i][j];
        }
    }      
    sctl::Matrix<Real> Usvd_p, VT_p, S_p;
    sctl::Matrix<Real> LPOforSVD = sctl::Matrix<Real>(LPOmat);
    LPOforSVD.SVD(Usvd_p, S_p, VT_p);
    std::cout << "debug by printing the last 10 matrix singular values:" << std::endl;
    std::cout << "shape of S_p is " << S_p.Dim(0) << ", " << S_p.Dim(1) << std::endl;
    for (long i=S_p.Dim(0)-10; i<S_p.Dim(0); i++) {
        std::cout << S_p(i,i) << std::endl;
    }
    // */

    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    solver(&sigma,BIO, bg_pres_flow(X0) * (pressure_drop/period_length), gmres_tol);
    // elem_lst0.WriteVTK("vis/sigma_ptcl2p", sigma, comm);

    sctl::Long Ntrg_side = 5;
    // Real gap = 1./(Ntrg_side+1);
    Real gap = 1./(Ntrg_side+5); // TRY: taking targets farther from plates
    X0.ReInit(Ntrg_side * Ntrg_side * 3);
    // X symmetry
    sctl::Vector<Real> X1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long yind=0; yind<Ntrg_side; yind++) {
        for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
            sctl::Long Ntrg_nodeind = yind*Ntrg_side + zind;
            X0[Ntrg_nodeind * 3 + 0] = 0.001;
            // X0[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
            // X0[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
            X0[Ntrg_nodeind * 3 + 1] = (yind+3)*gap; // Shift to start further from the plates.
            X0[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
            X1[Ntrg_nodeind * 3 + 0] = 0.999;
            X1[Ntrg_nodeind * 3 + 1] = (yind+3)*gap;
            X1[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
            
        }
        // std::cout << "debug: check y adn z values: at yind = " << yind << ", y value is " << (yind+3)*gap << std::endl;
    }
    sctl::Vector<Real> UX0(X0.Dim());
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp2.SetTargetCoord(X0);
    BIO_eval(&UX0,sigma);
    sctl::Vector<Real> UX1(X1.Dim());
    LayerPotenOp0.SetTargetCoord(X1);
    LayerPotenOp2.SetTargetCoord(X1);
    BIO_eval(&UX1,sigma);
    UX0 -= bg_pres_flow(X0) * (pressure_drop/period_length);
    UX1 -= bg_pres_flow(X1) * (pressure_drop/period_length);
    std::cout << "============ X periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffX = UX0-UX1;
    for (int i=0; i<UdiffX.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffX[i*3+0]<< ", " << UdiffX[i*3+1]<< ", " << UdiffX[i*3+2]<< ". " << std::endl;
    }

    // Y symmetry
    sctl::Vector<Real> Y0(Ntrg_side * Ntrg_side * 3);
    sctl::Vector<Real> Y1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long xind=0; xind<Ntrg_side; xind++) {
        for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
            sctl::Long Ntrg_nodeind = xind*Ntrg_side + zind;
            Y0[Ntrg_nodeind * 3 + 0] = (xind+3)*gap; // TRY: Changed here as well.
            Y0[Ntrg_nodeind * 3 + 1] = 0.001;
            Y0[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
            Y1[Ntrg_nodeind * 3 + 0] = (xind+3)*gap;
            Y1[Ntrg_nodeind * 3 + 1] = 0.999;
            Y1[Ntrg_nodeind * 3 + 2] = (zind+3)*gap;
        }
    }
    sctl::Vector<Real> UY0(Y0.Dim());
    LayerPotenOp0.SetTargetCoord(Y0);
    LayerPotenOp2.SetTargetCoord(Y0);
    BIO_eval(&UY0, sigma); 
    sctl::Vector<Real> UY1(Y1.Dim());
    LayerPotenOp0.SetTargetCoord(Y1);
    LayerPotenOp2.SetTargetCoord(Y1);
    BIO_eval(&UY1, sigma); 
    std::cout << "============ Y periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffY = UY0-UY1;
    for (int i=0; i<UdiffY.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffY[i*3+0]<< ", " << UdiffY[i*3+1]<< ", " << UdiffY[i*3+2]<< ". " << std::endl;
    }

    {
        PeriodicGeom<Real> trg;    
        CubeVolumeVisShifted<Real> vol_vis(10, 0.9, comm);
        // X0 = vol_vis.GetCoord();
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, 0);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);

        LayerPotenOp0.SetTargetCoord(X0);
        LayerPotenOp2.SetTargetCoord(X0);
        sctl::Vector<Real> U(X0.Dim());
        BIO_eval(&U, sigma);
        U -= bg_pres_flow(X0) * (pressure_drop/period_length);
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
        vol_vis.WriteVTK("vis/plane-U-2p-Ubg", U_vis);
    }
}

template <class Real> void test2peri_noplanes(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, const Real gmres_tol, const Real tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real pressure_drop = -1.0;
    const Real period_length = 1;

    const sctl::Long ElemOrder = 10;
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else { 
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 0, 2, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    }
    Nptcl = ptcls_rs.Dim(); 

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    Real surface_area;
    sctl::Vector<Real> wts;
    { // get wts and surface area
        sctl::Vector<Real> X, Xn, dist_far, surface_area_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        elem_lst0.GetFarFieldNodes(X, Xn, wts, dist_far, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts*0+1, wts);
        // surface_area = surface_area_[0];
        sctl::Vector<Real> sa_loc(1);
        sa_loc[0] = surface_area_[0];
        sctl::Vector<Real> sa_all(1);
        sa_all[0] = 0;
        comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
        surface_area = sa_all[0];
    }
    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0,"1");
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XY, period_length);

    // // periodized layer potential operator
    // const auto BIO = [DL_scal,&LayerPotenOp0,NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    //     U->SetZero();
    //     LayerPotenOp0.ComputePotential(*U, sigma);
    //     if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma*0.5*NormalOrient * DL_scal; // for double-layer
    // };
    
    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,&NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> sigma_mean, sigma0;
        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_;
            elem_lst0.GetFarFieldDensity(sigma_, sigma);
            SurfaceIntegral(sigma_mean, sigma_, wts);
            sigma_mean *= (1/surface_area);

            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);

            // // DEBUG: check that sigma-sigma_mean has surface integral = 0:
            // sctl::Vector<Real> sigma1 = sigma_;
            // AddConstVec(sigma1, -sigma_mean);
            // sctl::Vector<Real> sigma_test_;
            // SurfaceIntegral(sigma_test_, sigma1, wts);
            // std::cout << "Surface integral of sigma - sigma bar = " << sigma_test_[0] << ", "  << sigma_test_[1] << ", "  << sigma_test_[2] << std::endl;
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    // const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
    //     sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
    //     AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

    //     sctl::Vector<Real> U0;
    //     LayerPotenOp0.ComputeSL(U0, force_density);
    //     return U0;
    // };

    // first gmres to remove timing for matrix loading, and set Krylov preconditioner.
    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
    solver(&sigma,BIO, bg_pres_flow(X0) * (pressure_drop/period_length), gmres_tol);

    sctl::Long Ntrg_side = 5;
    Real gap = 1./(Ntrg_side+1);
    X0.ReInit(Ntrg_side * Ntrg_side * 3);
    // X symmetry
    sctl::Vector<Real> X1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long yind=0; yind<Ntrg_side; yind++) {
        for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
            sctl::Long Ntrg_nodeind = yind*Ntrg_side + zind;
            X0[Ntrg_nodeind * 3 + 0] = 0.;
            X0[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
            X0[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
            X1[Ntrg_nodeind * 3 + 0] = 1.;
            X1[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
            X1[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
        }
    }
    sctl::Vector<Real> UX0(X0.Dim());
    LayerPotenOp0.SetTargetCoord(X0);
    BIO(&UX0,sigma);
    sctl::Vector<Real> UX1(X1.Dim());
    LayerPotenOp0.SetTargetCoord(X1);
    BIO(&UX1,sigma);
    UX0 -= bg_pres_flow(X0) * (pressure_drop/period_length);
    UX1 -= bg_pres_flow(X1) * (pressure_drop/period_length);
    std::cout << "============ X periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffX = UX0-UX1;
    for (int i=0; i<UdiffX.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffX[i*3+0]<< ", " << UdiffX[i*3+1]<< ", " << UdiffX[i*3+2]<< ". " << std::endl;
    }

    // Y symmetry
    sctl::Vector<Real> Y0(Ntrg_side * Ntrg_side * 3);
    sctl::Vector<Real> Y1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long xind=0; xind<Ntrg_side; xind++) {
        for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
            sctl::Long Ntrg_nodeind = xind*Ntrg_side + zind;
            Y0[Ntrg_nodeind * 3 + 0] = (xind+1)*gap;
            Y0[Ntrg_nodeind * 3 + 1] = 0.0;
            Y0[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
            Y1[Ntrg_nodeind * 3 + 0] = (xind+1)*gap;
            Y1[Ntrg_nodeind * 3 + 1] = 1.;
            Y1[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
        }
    }
    sctl::Vector<Real> UY0(Y0.Dim());
    LayerPotenOp0.SetTargetCoord(Y0);
    BIO(&UY0, sigma); 
    sctl::Vector<Real> UY1(Y1.Dim());
    LayerPotenOp0.SetTargetCoord(Y1);
    BIO(&UY1, sigma); 
    std::cout << "============ Y periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffY = UY0-UY1;
    for (int i=0; i<UdiffY.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffY[i*3+0]<< ", " << UdiffY[i*3+1]<< ", " << UdiffY[i*3+2]<< ". " << std::endl;
    }
}

template <class Real> void test3peri(sctl::Long Nelem, sctl::Long FourierOrder, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, const Real gmres_tol, const Real tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

    const Real pressure_drop = -1.0;
    const Real period_length = 1;

    const sctl::Long ElemOrder = 10;
    
    PeriodicGeom<Real> obj;
    sctl::Vector<sctl::Long> ptcls;
    sctl::Vector<Real> ptcls_Xcs;
    sctl::Vector<Real> ptcls_rs;
    sctl::SlenderElemList<Real> elem_lst0, elem_lst_nbr;
    sctl::Vector<Real> NormalOrient;
    if (Nptcl == 1) {
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    } else { 
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build0 = obj.many_ptcls2(Nelem, ElemOrder, FourierOrder, 0, 3, comm, Nptcl, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        elem_lst0 = std::get<0>(build0);
        NormalOrient = std::get<1>(build0);
    }
    Nptcl = ptcls_rs.Dim(); 

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    Real surface_area;
    sctl::Vector<Real> wts;
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

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0,"1");
    
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, period_length);

    // // periodized layer potential operator
    // const auto BIO = [DL_scal,&LayerPotenOp0,NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
    //     U->SetZero();
    //     LayerPotenOp0.ComputePotential(*U, sigma);
    //     if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma*0.5*NormalOrient * DL_scal; // for double-layer
    // };

    // Define the boundary-integral operator: (I/2 + D + S)[sigma-sigma_mean] + sigma_mean
    const auto BIO = [&wts,&surface_area,&elem_lst0,&LayerPotenOp0,&DL_scal,&NormalOrient,&comm](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> sigma_mean, sigma0;
        { // compute sigma_mean and sigma0 = sigma - sigma_mean
            sctl::Vector<Real> sigma_;
            elem_lst0.GetFarFieldDensity(sigma_, sigma);
            SurfaceIntegral(sigma_mean, sigma_, wts);
            //MPI
            sctl::Vector<Real> sa_loc = sigma_mean;
            sctl::Vector<Real> sa_all(3);
            sa_all = 0;
            comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin(), (sctl::Iterator<Real>) sa_all.begin(), 1, sctl::CommOp::SUM);
            comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+1, (sctl::Iterator<Real>) sa_all.begin()+1, 1, sctl::CommOp::SUM);
            comm.Allreduce((sctl::Iterator<Real>) sa_loc.begin()+2, (sctl::Iterator<Real>) sa_all.begin()+2, 1, sctl::CommOp::SUM);
            sigma_mean = sa_all;
            sigma_mean *= (1/surface_area);

            sigma0 = sigma;
            AddConstVec(sigma0, -sigma_mean);

            // DEBUG: check that sigma-sigma_mean has surface integral = 0:
            sctl::Vector<Real> sigma1 = sigma_;
            AddConstVec(sigma1, -sigma_mean);
            sctl::Vector<Real> sigma_test_;
            SurfaceIntegral(sigma_test_, sigma1, wts);
            std::cout << "Surface integral of sigma - sigma bar = " << sigma_test_[0] << ", "<< sigma_test_[1] << ", " << sigma_test_[2] << ". "<< std::endl;
        
        }

        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma0);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma0*0.5*NormalOrient * DL_scal; // for double-layer

        AddConstVec(*U, sigma_mean);
    };

    const auto eval_rhs = [&LayerPotenOp0,surface_area,period_length](const Real pressure_drop) { // BIOpSL( -pressure_drop * cross_sectional_area / surface_area )
        sctl::Vector<Real> force_density(LayerPotenOp0.Dim(0)); force_density = 0;
        AddConstVec(force_density, sctl::Vector<Real>{-pressure_drop * period_length*period_length / surface_area, 0, 0});

        sctl::Vector<Real> U0;
        LayerPotenOp0.ComputeSL(U0, force_density);
        return U0;
    };

    // first gmres to remove timing for matrix loading, and set Krylov preconditioner.
    sctl::GMRES<Real> solver(comm);
    sctl::Vector<Real> sigma;
   LayerPotenOp0.SetTargetCoord(X0);
    // solver(&sigma,BIO, Uslip, gmres_tol);
    solver(&sigma, BIO, eval_rhs(pressure_drop), gmres_tol);

    sctl::Long Ntrg_side = 5;
    Real gap = 1./(Ntrg_side+1);
    X0.ReInit(Ntrg_side * Ntrg_side * 3);
    // X symmetry
    sctl::Vector<Real> X1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long yind=0; yind<Ntrg_side; yind++) {
        for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
            sctl::Long Ntrg_nodeind = yind*Ntrg_side + zind;
            X0[Ntrg_nodeind * 3 + 0] = 0.;
            X0[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
            X0[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
            X1[Ntrg_nodeind * 3 + 0] = 1.;
            X1[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
            X1[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
        }
    }
    sctl::Vector<Real> UX0(X0.Dim());
    LayerPotenOp0.SetTargetCoord(X0);
    BIO(&UX0,sigma);
    sctl::Vector<Real> UX1(X1.Dim());
    LayerPotenOp0.SetTargetCoord(X1);
    BIO(&UX1,sigma);
    UX0 -= eval_rhs(pressure_drop);
    UX1 -= eval_rhs(pressure_drop);
    std::cout << "============ X periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffX = UX0-UX1;
    for (int i=0; i<UdiffX.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffX[i*3+0]<< ", " << UdiffX[i*3+1]<< ", " << UdiffX[i*3+2]<< ". " << std::endl;
    }
    // Y symmetry
    sctl::Vector<Real> Y0(Ntrg_side * Ntrg_side * 3);
    sctl::Vector<Real> Y1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long xind=0; xind<Ntrg_side; xind++) {
        for (sctl::Long zind=0; zind<Ntrg_side; zind++) {
            sctl::Long Ntrg_nodeind = xind*Ntrg_side + zind;
            Y0[Ntrg_nodeind * 3 + 0] = (xind+1)*gap;
            Y0[Ntrg_nodeind * 3 + 1] = 0.0;
            Y0[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
            Y1[Ntrg_nodeind * 3 + 0] = (xind+1)*gap;
            Y1[Ntrg_nodeind * 3 + 1] = 1.;
            Y1[Ntrg_nodeind * 3 + 2] = (zind+1)*gap;
        }
    }
    sctl::Vector<Real> UY0(Y0.Dim());
    LayerPotenOp0.SetTargetCoord(Y0);
    BIO(&UY0,sigma);
    sctl::Vector<Real> UY1(Y1.Dim());
    LayerPotenOp0.SetTargetCoord(Y1);
    BIO(&UY1,sigma);
    std::cout << "============ Y periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffY = UY0-UY1;
    for (int i=0; i<UdiffY.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffY[i*3+0]<< ", " << UdiffY[i*3+1]<< ", " << UdiffY[i*3+2]<< ". " << std::endl;
    }
    // Z symmetry
    sctl::Vector<Real> Z0(Ntrg_side * Ntrg_side * 3);
    sctl::Vector<Real> Z1(Ntrg_side * Ntrg_side * 3);
    for (sctl::Long xind=0; xind<Ntrg_side; xind++) {
        for (sctl::Long yind=0; yind<Ntrg_side; yind++) {
            sctl::Long Ntrg_nodeind = xind*Ntrg_side + yind;
            Z0[Ntrg_nodeind * 3 + 0] = (xind+1)*gap;
            Z0[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
            Z0[Ntrg_nodeind * 3 + 2] = 0.0;
            Z1[Ntrg_nodeind * 3 + 0] = (xind+1)*gap;
            Z1[Ntrg_nodeind * 3 + 1] = (yind+1)*gap;
            Z1[Ntrg_nodeind * 3 + 2] = 1.0;
        }
    }
    sctl::Vector<Real> UZ0(Z0.Dim());
    LayerPotenOp0.SetTargetCoord(Z0);
    BIO(&UZ0,sigma);
    sctl::Vector<Real> UZ1(Z1.Dim());
    LayerPotenOp0.SetTargetCoord(Z1);
    BIO(&UZ1,sigma);
    std::cout << "============ Z periodicity =================" << std::endl;
    sctl::Vector<Real> UdiffZ = UZ0-UZ1;
    for (int i=0; i<UdiffZ.Dim()/3; i++) {
        std::cout << std::setprecision(8) << UdiffZ[i*3+0]<< ", " << UdiffZ[i*3+1]<< ", " << UdiffZ[i*3+2]<< ". " << std::endl;
}
}


int main(int argc, char** argv) {

    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Comm comm = sctl::Comm::World();
        // sctl::Profile::Enable(true);
        long Nelem_ptcl = std::stol(argv[1]); // number of elements
        long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
        int write_ref = std::stoi(argv[3]);
        int peri_mode = std::stoi(argv[4]); // what kind of periodicity does the system have; peri_mode = j for j-periodic.
        long Nptcl = std::stol(argv[5]); // number of particles inside
        long geom_mode = std::stol(argv[6]); // =0: spheres; =1: spheroids; =3: bacteria; =4: loop.
        double gmres_tol = std::stod(argv[7]);
        double tol = std::stod(argv[8]);

        // test<Real>(Nelem_ptcl, FourierOrder, (write_ref==1), peri_mode, comm, Nptcl, geom_mode, gmres_tol, tol);
        // plot_setup<Real>(Nelem_ptcl, FourierOrder, (write_ref==1), peri_mode, comm, Nptcl, geom_mode, gmres_tol, tol);
        if (peri_mode==1) {
            test1peri<Real>(Nelem_ptcl, FourierOrder, comm, Nptcl, geom_mode, gmres_tol, tol);
        } else if (peri_mode == 2) {
            test2peri<Real>(Nelem_ptcl, FourierOrder, comm, Nptcl, geom_mode, gmres_tol, tol);
        } else {
            test3peri<Real>(Nelem_ptcl, FourierOrder, comm, Nptcl, geom_mode, gmres_tol, tol);
        }
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
