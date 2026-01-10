#include "periodize.hpp"
#include "utils.hpp"
// Test script for calculation and timing of particle-only 1 and 3 periodic problems with background pressure flow. 

/**
 * Background flow with unit pressure gradient along X-axis.
 */
template <class Real> sctl::Vector<Real> bg_flow(const sctl::Vector<Real>& X) {
    const Real pdrive = 1;
    const sctl::Long N = X.Dim()/3;
    sctl::Vector<Real> U(N*3);
    for (sctl::Long i = 0; i < N; i++) {
        const auto x = X.begin() + i*3;
        U[i*3+0] = -pdrive * ((x[1]-0.5)*(x[1]-0.5) + (x[2]-0.5)*(x[2]-0.5))/4;
        U[i*3+1] = 0;
        U[i*3+2] = 0;
    }
    return U;
}

// Uniform background flow in x direction.
template <class Real> sctl::Vector<Real> bg_unif_flow(const sctl::Vector<Real>& X) {
    sctl::Vector<Real> U = X;
    // const sctl::Long N = X.Dim() /3;
    U = 1.; // background flow diagonal to avoid planes of unaffected flows between periods.
    // for (sctl::Long i = 0; i < N; i++) {
    //     U[i*3+0] = 1.; // background flow in only x direction for timing runs.
    // }
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


template <class Real> void test(sctl::Long Nelem, sctl::Long FourierOrder, bool write_ref, sctl::Integer peri_mode, sctl::Comm comm, sctl::Long Nptcl, sctl::Long geom_mode, const Real gmres_tol, const Real tol) {

    // Combine single-layer and double-layer kernels in these proportions
    const Real SL_scal = 1.0;
    const Real DL_scal = 1.0;

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
    // const sctl::Long Nrepeat = elem_lst_nbr.Size() / elem_lst0.Size(); 
    Nptcl = ptcls_rs.Dim(); 
    // std::cout << "periodic mode is " << peri_mode << ", Nrepeat is " << Nrepeat << std::endl;

    sctl::Vector<Real> X0; // target coordinates
    elem_lst0.GetNodeCoord(&X0, nullptr, nullptr);
    // if (write_ref) {
    //     elem_lst0.WriteVTK("vis/"+std::to_string(Nptcl)+"spheres",X0,comm);
    // }  
    sctl::Long ptcl_gridsize = Nelem * ElemOrder * FourierOrder * 3; 
    sctl::Long Nptcl_slip = elem_lst0.Size() / Nelem; 
    sctl::Vector<Real> ptcls_Xcs_slip(Nptcl_slip * 3, (sctl::Iterator<Real>)ptcls_Xcs.begin() + comm.Rank()*Nptcl_slip*3, true);
    sctl::Vector<Real> ptcls_rs_slip(Nptcl_slip,  (sctl::Iterator<Real>)ptcls_rs.begin()+comm.Rank()*Nptcl_slip, true);
    sctl::Vector<Real> Uslip = total_vslip(X0, ptcl_gridsize, Nptcl_slip, ptcls_Xcs_slip, ptcls_rs_slip);
    if (write_ref) {
        elem_lst0.WriteVTK("vis/"+std::to_string(Nptcl)+"spheres_streamline_vslip/"+std::to_string(Nptcl)+"spheres",Uslip,comm);
    }  

    StokesBIO LayerPotenOp0(SL_scal, DL_scal, comm); // potential from elem_lst_nbr to X0
    LayerPotenOp0.AddElemList(elem_lst0);
    LayerPotenOp0.SetTargetCoord(X0);
    LayerPotenOp0.SetAccuracy(tol);
    if (peri_mode==1) {
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::X, 1.0);
    } else if (peri_mode==3) {
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XYZ, 1.0);
    } else if (peri_mode==2) {
        LayerPotenOp0.SetPeriodicity(sctl::Periodicity::XY, 1.0);
    } else {
        SCTL_ASSERT(false);
    }

    // =============== PRECONDITIONING =======================================
    // Store preconditioner matrix, or make new if not present.
    std::string precond0_file = "data/precond0_ptcl_Np"+std::to_string(Nelem)+"_Nf"+std::to_string(FourierOrder)+".mat";
    std::string precond1_file = "data/precond1_ptcl_Np"+std::to_string(Nelem)+"_Nf"+std::to_string(FourierOrder)+".mat";
    sctl::Matrix<Real> PrecondMat0, PrecondMat1;
    PrecondMat0.template Read<Real>(precond0_file.c_str());

    sctl::Long A11size;

    comm.Barrier();
    if (PrecondMat0.Dim(0) || PrecondMat0.Dim(1)) {
        std::cout << " successfully read file." << std::endl;
        PrecondMat1.template Read<Real>(precond1_file.c_str());
        A11size = PrecondMat0.Dim(1);
    } else {
        std::cout << "making precond matrices" << std::endl;
        sctl::Vector<sctl::Long> ptcls_pre;
        sctl::Vector<Real> ptcls_Xcs_pre, ptcls_rs_pre;
        std::tuple<sctl::SlenderElemList<Real>,sctl::Vector<Real>> build_precond = obj.many_ptcls1(Nelem, ElemOrder, FourierOrder, 0, 1, comm.Self(), ptcls_pre, ptcls_rs_pre, ptcls_Xcs_pre, geom_mode);
        sctl::SlenderElemList<Real> elem_lst_precond = std::get<0>(build_precond);
        sctl::Vector<Real> X0_precond; // target coordinates
        elem_lst_precond.GetNodeCoord(&X0_precond, nullptr, nullptr);
        StokesBIO Precond_bio(SL_scal, DL_scal, comm.Self());
        Precond_bio.SetAccuracy(tol); // set quadrature accuracy
        Precond_bio.AddElemList(elem_lst_precond);
        Precond_bio.SetTargetCoord(X0_precond);
        const auto BIO_1ptcl = [&DL_scal,&Precond_bio](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
            U->SetZero();
            Precond_bio.ComputePotential(*U, sigma);
            (*U) += sigma * 0.5 * DL_scal;
        };
        A11size = 3*ElemOrder*FourierOrder*Nelem;
        sctl::Vector<sctl::Vector<Real>> PrecondMat(A11size);
        sctl::Vector<Real> SigmaCol_precond(A11size);
        for (sctl::Long col=0; col < A11size; col ++) {
            SigmaCol_precond = 0.;
            SigmaCol_precond[col] = 1.;
            BIO_1ptcl(PrecondMat.begin() + col,SigmaCol_precond);
        }
        sctl::Matrix<Real> A11(A11size,A11size);
        for (long col=0; col < A11size; col++) {
            for (long row = 0; row < A11size; row++) {
                A11(row,col) = PrecondMat[col][row];
            }
        }      
        sctl::Matrix<Real> Usvd, VT, S, SforInv;
        sctl::Matrix<Real> A11forSVD = sctl::Matrix<Real>(A11);
        A11forSVD.SVD(Usvd, S, VT);
        SforInv = sctl::Matrix<Real>(S);
        sctl::Matrix<Real> Sinv = SforInv.pinv(1e-16);

        PrecondMat0 = VT.Transpose();
        PrecondMat1 = Sinv * Usvd.Transpose();
        if (!comm.Rank()) {
            PrecondMat0.template Write<Real>(precond0_file.c_str());
            PrecondMat1.template Write<Real>(precond1_file.c_str());
        }
    }

    // periodized layer potential operator
    const auto BIO = [&DL_scal,&LayerPotenOp0,&X0,NormalOrient](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        U->SetZero();
        LayerPotenOp0.ComputePotential(*U, sigma);
        if (DL_scal && U->Dim() == sigma.Dim()) (*U) -= sigma*0.5*NormalOrient * DL_scal; // for double-layer
    };

    // Apply A11inv to each panel of vec.
    const auto AinvApply = [&PrecondMat0,&PrecondMat1,&A11size, &comm](const sctl::Vector<Real>& vec) {
        sctl::Long N = vec.Dim();
        sctl::Long Nptcl = N / A11size; 
        sctl::Vector<Real> AinvVec(N);
        for (sctl::Long i=0; i<Nptcl; i++) {
            // for each particle, apply A11inv.
            sctl::Matrix<Real> vecMat(A11size,1,(sctl::Iterator<Real>) vec.begin() + i*A11size,true);
            sctl::Matrix<Real> AinvVecMat = PrecondMat0 * (PrecondMat1 * vecMat);
            for (sctl::Long j=0; j<A11size; j++) {
                AinvVec[i*A11size + j] = AinvVecMat(j,0);
            }
        }
        return AinvVec;
    };

    const auto BIO_precond = [&BIO,&AinvApply](sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) {
        sctl::Vector<Real> Uloc;
        BIO(&Uloc,sigma);
        // LEFT PRECONDITIONER: u -> A11inv*u
        (*U) = AinvApply(Uloc);
    };

    // first gmres to remove timing for matrix loading, and set Krylov preconditioner.
    sctl::GMRES<Real> solver(comm);
    sctl::KrylovPrecond<Real> krylov_precond;
    sctl::Vector<Real> A11invF = AinvApply(Uslip);

    // sctl::Vector<Real> sigma_temp;
    // solver(&sigma_temp, BIO_precond, A11invF, 1e-2);
    // // solver(&sigma_temp, BIO, -bg_flow(X0), 1e0);
    // sctl::Profile::reset();

    // LayerPotenOp0.ClearSetup();
    // sctl::Profile::Tic("Setup SurfOP");
    // LayerPotenOp0.Setup();
    // sctl::Profile::Toc();
    // sctl::Profile::print(&comm);
    // sctl::Profile::reset();

    // // // Weak scaling data did not include this.
    // // sctl::Profile::Tic("Setup ProxyOP");
    // // LayerPotenOp_proxy.ClearSetup();
    // // LayerPotenOp_proxy.Setup(); 
    // // sctl::Profile::Toc();
    // // sctl::Profile::print(&comm);
    // // sctl::Profile::reset();

    // sctl::Profile::Tic("Solve without Precond");
    // sctl::Vector<Real> sigma_temp2;
    // solver(&sigma_temp2, BIO, -bg_unif_flow(X0), gmres_tol, -1, false);
    // sctl::Profile::Toc();
    // sctl::Profile::print(&comm, {"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    // sctl::Profile::reset();
    // comm.Barrier();

    sctl::Vector<Real> sigma;
    sctl::Profile::Tic("Solver: KrylovPrecond_setup");
    solver(&sigma, BIO_precond, A11invF, gmres_tol, -1, false, nullptr, &krylov_precond);
    // solver(&sigma,BIO,Uslip, gmres_tol, -1, false, nullptr, &krylov_precond);
    // solver(&sigma, BIO_precond, A11invF, gmres_tol, -1, false);
    sctl::Profile::Toc();
    sctl::Profile::print(&comm, {"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    sctl::Profile::reset();
    comm.Barrier();

    // sctl::Vector<Real> sigma1;
    // sctl::Profile::Tic("Solver1");
    // // PRECOND with Krylov
    // solver(&sigma1, BIO_precond, A11invF, gmres_tol, -1, false, nullptr, &krylov_precond);
    // // solver(&sigma1,BIO,-bg_flow(X0), gmres_tol, -1, false, nullptr, &krylov_precond);
    // sctl::Profile::Toc();
    // sctl::Profile::print(&comm, {"t_avg", "t_max", "f_avg", "f_max", "m_min", "m_avg", "m_max"});
    // sctl::Profile::reset();
    // comm.Barrier();
    if (!comm.Rank()) {
        std::cout << "------------------- DONE WITH SOLVE ======================" << std::endl;
    }

    if (write_ref) { 
        PeriodicGeom<Real> trg;    
        CubeVolumeVisShifted<Real> vol_vis(100, 0.9, comm);
        // VolumeVis<Real> vol_vis(elem_lst_trg, comm); 
        // X0 = vol_vis.GetCoord();
        sctl::Vector<Real> X0_all = vol_vis.GetCoord();
        sctl::Vector<sctl::Long> filtered_inds(X0_all.Dim()/3);
        std::tuple<sctl::Vector<Real>,sctl::Vector<sctl::Long>> trg_tuple = trg.filter_target(X0_all, ptcls, ptcls_rs, ptcls_Xcs, geom_mode);
        X0 = std::get<0>(trg_tuple);
        filtered_inds = std::get<1>(trg_tuple);
        std::cout << "number of target points: " << X0.Dim() << std::endl;

        LayerPotenOp0.SetTargetCoord(X0);
        sctl::Vector<Real> U;
        BIO(&U, sigma);

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
        vol_vis.WriteVTK("vis/"+std::to_string(Nptcl)+"spheres_streamline_vslip/"+std::to_string(Nptcl)+"streamlines", U_vis); 
    } 
}


int main(int argc, char** argv) {

    sctl::Comm::MPI_Init(&argc, &argv);
    using Real = double;

    {
        sctl::Comm comm = sctl::Comm::World();
        sctl::Profile::Enable(true);
        long Nelem_ptcl = std::stol(argv[1]); // number of elements
        long FourierOrder = std::stol(argv[2]);  // number of Fourier nodes
        int write_ref = std::stoi(argv[3]);
        int peri_mode = std::stoi(argv[4]); // what kind of periodicity does the system have; peri_mode = j for j-periodic.
        long Nptcl = std::stol(argv[5]); // number of particles inside
        long geom_mode = std::stol(argv[6]); // =0: spheres; =1: spheroids; =3: bacteria; =4: loop.
        double gmres_tol = std::stod(argv[7]);
        double tol = std::stod(argv[8]);

        test<Real>(Nelem_ptcl, FourierOrder, (write_ref==1), peri_mode, comm, Nptcl, geom_mode, gmres_tol, tol);
        // plot_setup<Real>(Nelem_ptcl, FourierOrder, (write_ref==1), peri_mode, comm, Nptcl, geom_mode, gmres_tol, tol);
    }

    sctl::Comm::MPI_Finalize();
    return 0;
}
