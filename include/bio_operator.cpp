#ifndef _BIO_OPERATOR_CPP_
#define _BIO_OPERATOR_CPP_

#include "bio_operator.hpp"

template <class Real>
MeanCorrectedStokesBIOOperator<Real>::MeanCorrectedStokesBIOOperator(StokesBIO<Real>& bio,
                                                                     const sctl::Vector<Real>& normal_orient,
                                                                     const Real dl_scal,
                                                                     const sctl::Comm& comm,
                                                                     const Real double_layer_sign)
    : bio_(bio),
      normal_orient_(normal_orient),
      dl_scal_(dl_scal),
      double_layer_sign_(double_layer_sign),
      comm_(comm) {}

template <class Real>
void MeanCorrectedStokesBIOOperator<Real>::AddSurface(const DensityGetter& density_getter,
                                                      const sctl::Vector<Real>& weights,
                                                      const Real surface_area) {
    surfaces_.push_back({density_getter, weights, surface_area});
}

template <class Real>
template <class ElemLstType>
void MeanCorrectedStokesBIOOperator<Real>::AddSurface(const ElemLstType& elem_lst,
                                                      const sctl::Vector<Real>& weights,
                                                      const Real surface_area) {
    AddSurface([&elem_lst](sctl::Vector<Real>& sigma_far, const sctl::Vector<Real>& sigma) {
                   elem_lst.GetFarFieldDensity(sigma_far, sigma);
               },
               weights,
               surface_area);
}

template <class Real>
void MeanCorrectedStokesBIOOperator<Real>::ApplyDoubleLayerCorrection(sctl::Vector<Real>& U,
                                                                      const sctl::Vector<Real>& sigma) const {
    if (!dl_scal_ || U.Dim() != sigma.Dim() || normal_orient_.Dim() != sigma.Dim()) {
        return;
    }
    U -= sigma * 0.5 * normal_orient_ * dl_scal_ * double_layer_sign_;
}

template <class Real>
void MeanCorrectedStokesBIOOperator<Real>::AccumulateSurfaceMean(sctl::Vector<Real>& sigma_mean,
                                                                 Real& total_surface_area,
                                                                 const sctl::Vector<Real>& sigma) const {
    sigma_mean.ReInit(3);
    sigma_mean = 0;
    total_surface_area = 0;

    for (const auto& surface : surfaces_) {
        sctl::Vector<Real> sigma_far;
        surface.density_getter(sigma_far, sigma);

        sctl::Vector<Real> mean_local;
        SurfaceIntegral(mean_local, sigma_far, surface.weights);

        sctl::Vector<Real> mean_global(3);
        mean_global = 0;
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

        sigma_mean += mean_global;
        total_surface_area += surface.surface_area;
    }

    if (total_surface_area != 0) {
        sigma_mean *= (1. / total_surface_area);
    }
}

template <class Real>
void MeanCorrectedStokesBIOOperator<Real>::Apply(sctl::Vector<Real>& U,
                                                 const sctl::Vector<Real>& sigma) const {
    U.SetZero();

    sctl::Vector<Real> sigma_eval = sigma;
    sctl::Vector<Real> sigma_mean;

    if (!surfaces_.empty()) {
        Real total_surface_area = 0;
        AccumulateSurfaceMean(sigma_mean, total_surface_area, sigma);
        sigma_eval = sigma;
        AddConstVec(sigma_eval, -sigma_mean);
    }

    bio_.ComputePotential(U, sigma_eval);
    ApplyDoubleLayerCorrection(U, sigma_eval);

    if (!surfaces_.empty()) {
        AddConstVec(U, sigma_mean);
    }
}

template <class Real>
PreconditionedMeanCorrectedStokesBIOOperator<Real>::PreconditionedMeanCorrectedStokesBIOOperator(
    const MeanCorrectedStokesBIOOperator<Real>& bio,
    const sctl::Matrix<Real>& precond_mat0,
    const sctl::Matrix<Real>& precond_mat1,
    const sctl::Long A11size)
    : bio_(bio),
      precond_mat0_(precond_mat0),
      precond_mat1_(precond_mat1),
      A11size_(A11size) {}

template <class Real>
void PreconditionedMeanCorrectedStokesBIOOperator<Real>::Apply(sctl::Vector<Real>& U,
                                                               const sctl::Vector<Real>& sigma) const {
    sctl::Vector<Real> Uloc;
    bio_.Apply(Uloc, sigma);
    U = ApplyAinv(Uloc);
}

template <class Real>
sctl::Vector<Real> PreconditionedMeanCorrectedStokesBIOOperator<Real>::ApplyAinv(
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
