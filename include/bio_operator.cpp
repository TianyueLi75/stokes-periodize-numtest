#ifndef _BIO_OPERATOR_CPP_
#define _BIO_OPERATOR_CPP_

#include "bio_operator.hpp"

template <class Real, class ElemLstType> 
MeanCorrectedStokesBIOOperator<Real, ElemLstType>::MeanCorrectedStokesBIOOperator(StokesBIO<Real>& bio,
                                                                     const sctl::Vector<Real>& normal_orient,
                                                                     const Real dl_scal,
                                                                     const sctl::Comm& comm,
                                                                     const Real double_layer_sign)
    : bio_(bio),
      normal_orient_(normal_orient),
      dl_scal_(dl_scal),
      double_layer_sign_(double_layer_sign),
      comm_(comm) {}

template <class Real, class ElemLstType>
void MeanCorrectedStokesBIOOperator<Real, ElemLstType>::AddSurface(const ElemLstType& elem_lst) {
    surfaces_.PushBack(elem_lst);
}

template <class Real, class ElemLstType> 
void MeanCorrectedStokesBIOOperator<Real, ElemLstType>::ApplyDoubleLayerCorrection(sctl::Vector<Real>& U,
                                                                      const sctl::Vector<Real>& sigma) const {
    if (!dl_scal_ || U.Dim() != sigma.Dim() || normal_orient_.Dim() != sigma.Dim()) {
        return;
    }
    U -= sigma * 0.5 * normal_orient_ * dl_scal_ * double_layer_sign_;
}

template <class Real, class ElemLstType> 
void MeanCorrectedStokesBIOOperator<Real, ElemLstType>::AccumulateSurfaceMean(sctl::Vector<Real>& sigma_mean,
                                                                 Real& total_surface_area,
                                                                 const sctl::Vector<Real>& sigma) const {
    sigma_mean.ReInit(3);
    sigma_mean = 0;
    total_surface_area = 0;

    for (const auto& surface : surfaces_) {
        // Get surface area
        sctl::Vector<Real> X, Xn, dist_far_, surface_area_, wts_;
        sctl::Vector<sctl::Long> element_wise_node_cnt;
        surface.GetFarFieldNodes(X, Xn, wts_, dist_far_, element_wise_node_cnt, 1);
        SurfaceIntegral(surface_area_, wts_*0+1, wts_);

        // TODO: delete X, Xn, dist_far.
        sctl::Vector<Real> sigma_far;
        surface.GetFarFieldDensity(sigma_far, sigma);

        // Get integral of density
        sctl::Vector<Real> mean_local;
        SurfaceIntegral(mean_local, sigma_far, wts_);

        sctl::Vector<Real> mean_global(3);
        sctl::Vector<Real> surface_area_global(1);
        mean_global = 0;
        surface_area_global = 0;
        for (sctl::Long i = 0; i < 3; i++) {
            sctl::Vector<Real> loc(1);
            loc[0] = mean_local[i];
            sctl::Vector<Real> all(1);
            all[0] = 0;
            comm_.Allreduce((sctl::Iterator<Real>) loc.begin(),
                            (sctl::Iterator<Real>) all.begin(),
                            1,
                            sctl::CommOp::SUM);
            mean_global[i] = all[0];
        }
        // collect surface area
        sctl::Vector<Real> loc(1);
        loc[0] = surface_area_[0];
        sctl::Vector<Real> all(1);
        all[0] = 0.;
        comm_.Allreduce((sctl::Iterator<Real>) loc.begin(),
                        (sctl::Iterator<Real>) all.begin(),
                        1,
                        sctl::CommOp::SUM);
        surface_area_global[0] = all[0];

        sigma_mean += mean_global;
        total_surface_area += surface_area_global[0];
    }

    if (total_surface_area != 0) {
        sigma_mean *= (1. / total_surface_area);
    }
}

template <class Real, class ElemLstType> 
void MeanCorrectedStokesBIOOperator<Real, ElemLstType>::Apply(sctl::Vector<Real>& U,
                                                 const sctl::Vector<Real>& sigma) const {
    U.SetZero();

    sctl::Vector<Real> sigma_eval = sigma;
    sctl::Vector<Real> sigma_mean;

    if (surfaces_.Dim()) {
        Real total_surface_area = 0;
        AccumulateSurfaceMean(sigma_mean, total_surface_area, sigma);
        sigma_eval = sigma;
        AddConstVec(sigma_eval, -sigma_mean);
    }

    bio_.ComputePotential(U, sigma_eval);
    ApplyDoubleLayerCorrection(U, sigma_eval);

    if (surfaces_.Dim()) {
        AddConstVec(U, sigma_mean);
    }
}

template <class Real, class ElemLstType> 
PreconditionedMeanCorrectedStokesBIOOperator<Real, ElemLstType>::PreconditionedMeanCorrectedStokesBIOOperator(
    const MeanCorrectedStokesBIOOperator<Real, ElemLstType>& bio,
    const sctl::Matrix<Real>& precond_mat0,
    const sctl::Matrix<Real>& precond_mat1,
    const sctl::Long A11size)
    : bio_(bio),
      precond_mat0_(precond_mat0),
      precond_mat1_(precond_mat1),
      A11size_(A11size) {}

template <class Real, class ElemLstType> 
void PreconditionedMeanCorrectedStokesBIOOperator<Real, ElemLstType>::Apply(sctl::Vector<Real>& U,
                                                               const sctl::Vector<Real>& sigma) const {
    sctl::Vector<Real> Uloc;
    bio_.Apply(Uloc, sigma);
    U = ApplyAinv(Uloc);
}

template <class Real, class ElemLstType> 
sctl::Vector<Real> PreconditionedMeanCorrectedStokesBIOOperator<Real, ElemLstType>::ApplyAinv(
    const sctl::Vector<Real>& vec) const {
    sctl::Vector<Real> AinvVec(vec.Dim());
    const sctl::Long N = vec.Dim();
    const sctl::Long Npanels = N / A11size_;

    for (sctl::Long i = 0; i < Npanels; i++) {
        sctl::Matrix<Real> vecMat(A11size_, 1, (sctl::Iterator<Real>) vec.begin() + i * A11size_, true);
        sctl::Matrix<Real> AinvVecMat = precond_mat0_ * (precond_mat1_ * vecMat);
        for (sctl::Long j = 0; j < A11size_; j++) {
            AinvVec[i * A11size_ + j] = AinvVecMat(j, 0);
        }
    }

    return AinvVec;
}

#endif
