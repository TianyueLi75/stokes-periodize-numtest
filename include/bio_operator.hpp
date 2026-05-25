#ifndef _BIO_OPERATOR_HPP_
#define _BIO_OPERATOR_HPP_

#include "stokes_bio.hpp"

#include <functional>
#include <vector>

template <class Real> void SurfaceIntegral(sctl::Vector<Real>& I,
                                            const sctl::Vector<Real>& vals,
                                            const sctl::Vector<Real>& wts) {
    const sctl::Long dof = vals.Dim() / wts.Dim();
    SCTL_ASSERT(vals.Dim() == wts.Dim() * dof);
    if (I.Dim() != dof) I.ReInit(dof);
    I = 0;
    for (sctl::Long i = 0; i < wts.Dim(); i++) {
        for (sctl::Long j = 0; j < dof; j++) {
            I[j] += vals[i * dof + j] * wts[i];
        }
    }
}

template <class Real> void AddConstVec(sctl::Vector<Real>& vals,
                                        const sctl::Vector<Real>& c0) {
    const sctl::Long dof = c0.Dim();
    const sctl::Long N = vals.Dim() / dof;
    SCTL_ASSERT(vals.Dim() == N * dof);
    for (sctl::Long i = 0; i < N; i++) {
        for (sctl::Long j = 0; j < dof; j++) {
            vals[i * dof + j] += c0[j];
        }
    }
}

template <class Real>
class MeanCorrectedStokesBIOOperator {
  public:
    using DensityGetter = std::function<void(sctl::Vector<Real>&, const sctl::Vector<Real>&)>;

    struct SurfaceSpec {
        DensityGetter density_getter;
        sctl::Vector<Real> weights;
        Real surface_area;
    };

    MeanCorrectedStokesBIOOperator(StokesBIO<Real>& bio,
                                   const sctl::Vector<Real>& normal_orient,
                                   const Real dl_scal,
                                   const sctl::Comm& comm,
                                   const Real double_layer_sign = -1.0);

    void AddSurface(const DensityGetter& density_getter,
                    const sctl::Vector<Real>& weights,
                    const Real surface_area);

    template <class ElemLstType>
    void AddSurface(const ElemLstType& elem_lst,
                    const sctl::Vector<Real>& weights,
                    const Real surface_area);

    void Apply(sctl::Vector<Real>& U, const sctl::Vector<Real>& sigma) const;

    void operator()(sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) const {
        Apply(*U, sigma);
    }

  private:
    void ApplyDoubleLayerCorrection(sctl::Vector<Real>& U,
                                    const sctl::Vector<Real>& sigma) const;

    void AccumulateSurfaceMean(sctl::Vector<Real>& sigma_mean,
                               Real& total_surface_area,
                               const sctl::Vector<Real>& sigma) const;

    StokesBIO<Real>& bio_;
    sctl::Vector<Real> normal_orient_;
    Real dl_scal_;
    Real double_layer_sign_;
    sctl::Comm comm_;
    std::vector<SurfaceSpec> surfaces_;
};

template <class Real>
class PreconditionedMeanCorrectedStokesBIOOperator {
  public:
    PreconditionedMeanCorrectedStokesBIOOperator(const MeanCorrectedStokesBIOOperator<Real>& bio,
                                                 const sctl::Matrix<Real>& precond_mat0,
                                                 const sctl::Matrix<Real>& precond_mat1,
                                                 const sctl::Long A11size);

    void Apply(sctl::Vector<Real>& U, const sctl::Vector<Real>& sigma) const;

    void operator()(sctl::Vector<Real>* U, const sctl::Vector<Real>& sigma) const {
        Apply(*U, sigma);
    }

    sctl::Vector<Real> ApplyAinv(const sctl::Vector<Real>& vec) const;

  private:
    const MeanCorrectedStokesBIOOperator<Real>& bio_;
    sctl::Matrix<Real> precond_mat0_;
    sctl::Matrix<Real> precond_mat1_;
    sctl::Long A11size_;
};

#include "bio_operator.cpp"

#endif
