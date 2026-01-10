#include <sctl.hpp>
#include "planeNaive.hpp"

namespace sctl {
    template <class Real> const std::pair<Vector<Real>,Vector<Real>>& PlaneIntegral<Real>::LegendreQuad_plane(Integer ORDER) {
        constexpr Integer max_order = 50;
        auto compute_nds_wts = [max_order]() {
            Vector<std::pair<Vector<Real>,Vector<Real>>> nds_wts(max_order);
            for (Integer order = 1; order < max_order; order++) {
                auto& x_ = nds_wts[order].first;
                auto& w_ = nds_wts[order].second;
                LegQuadRule<Real>::ComputeNdsWts(&x_, &w_, order);
            }
            return nds_wts;
        };
        static const auto nds_wts = compute_nds_wts();

        SCTL_ASSERT(ORDER < max_order);
        return nds_wts[ORDER];
    }

    template <class Real> PlaneIntegral<Real>::PlaneIntegral(const Long order, const Long Nelem_x, const Long Nelem_y, const Real z_offset) {
        // constructor
        order_ = order;
        z_offset_ = z_offset;
        auto gl_grid = LegendreQuad_plane(order_);
        gl_nodes_ = gl_grid.first;
        gl_wts_ = gl_grid.second;
        Nelem_x_ = Nelem_x;
        Nelem_y_ = Nelem_y;
        Long g = gl_nodes_.Dim();

        Long Nnodes = Nelem_x * g * Nelem_y * g * 2;
        Xsrc_.ReInit(Nnodes*KDIM);
        Xsrc_n_.ReInit(Nnodes*KDIM);
        Xsrc_n_.SetZero();
        Xwts_.ReInit(Nnodes);
        Real elem_len_x = 1./Nelem_x; // Hardcoded size of unit box being 1.
        Real elem_len_y = 1./Nelem_y;
        // Top plane
        for (Long xind=0; xind < Nelem_x; xind ++) {
            for (Long yind=0; yind < Nelem_y; yind ++) {
                for (Long node_xind = 0; node_xind < g; node_xind ++) {
                    for (Long node_yind = 0; node_yind < g; node_yind ++) {
                        Long cur_node_ind = xind * Nelem_y * g * g + \
                                            yind * g * g + node_xind * g + node_yind;
                        Real x_corner = xind * elem_len_x;
                        Real y_corner = yind * elem_len_y;
                        Xwts_[cur_node_ind] = gl_wts_[node_xind] * gl_wts_[node_yind] * elem_len_x * elem_len_y;
                        Xsrc_[KDIM*cur_node_ind + 0] = x_corner + elem_len_x * gl_nodes_[node_xind];
                        Xsrc_[KDIM*cur_node_ind + 1] = y_corner + elem_len_y * gl_nodes_[node_yind];
                        Xsrc_[KDIM*cur_node_ind + 2] = 1. - z_offset_; // top plane first
                        Xsrc_n_[KDIM*cur_node_ind + 2] = -1.; // inward normal points down
                    }
                }
            }
        }
        // Bottom plane 
        for (Long xind=0; xind < Nelem_x; xind ++) {
            for (Long yind=0; yind < Nelem_y; yind ++) {
                for (Long node_xind = 0; node_xind < gl_nodes_.Dim(); node_xind ++) {
                    for (Long node_yind = 0; node_yind < gl_nodes_.Dim(); node_yind ++) {
                        Long cur_node_ind = xind * Nelem_y * g * g + \
                                            yind * g * g + node_xind * g + node_yind + \
                                            Nelem_x * g * Nelem_y * g;
                        Real x_corner = xind * elem_len_x;
                        Real y_corner = yind * elem_len_y;
                        Xwts_[cur_node_ind] = gl_wts_[node_xind] * gl_wts_[node_yind] * elem_len_x * elem_len_y;
                        Xsrc_[KDIM*cur_node_ind + 0] = x_corner + elem_len_x * gl_nodes_[node_xind];
                        Xsrc_[KDIM*cur_node_ind + 1] = y_corner + elem_len_y * gl_nodes_[node_yind];
                        Xsrc_[KDIM*cur_node_ind + 2] = z_offset_; 
                        Xsrc_n_[KDIM*cur_node_ind + 2] = 1.; // inward normal points up
                    }
                }
            }
        }
    }


    template <class Real> void PlaneIntegral<Real>::GetNodeCoord(Vector<Real>* X, Vector<Real>* Xn, Vector<Long>* element_wise_node_cnt) const {
        (*X) = Xsrc_;
        if (Xn) {
            (*Xn) = Xsrc_n_;
        }
        if (element_wise_node_cnt) { // expecting a scan or a vector with number of elements per node?
            element_wise_node_cnt->ReInit(Nelem_x_*Nelem_y_*2);
            (*element_wise_node_cnt) = gl_nodes_.Dim() * gl_nodes_.Dim(); 
        }
    }

    template <class Real> void PlaneIntegral<Real>:: GetFarFieldNodes(Vector<Real>& X, Vector<Real>& Xn, Vector<Real>& wts, Vector<Real>& dist_far, Vector<Long>& element_wise_node_cnt, const Real tol) const {
        // first no change to gl_grid for far
        X = Xsrc_;
        Xn = Xsrc_n_;
        wts = Xwts_;
        if (element_wise_node_cnt.Dim() != Nelem_x_*Nelem_y_*2) {
            element_wise_node_cnt.ReInit(Nelem_x_ * Nelem_y_*2);
        }
        element_wise_node_cnt = gl_nodes_.Dim() * gl_nodes_.Dim();
        
        // dist_far by gauss rule, from slender_element.cpp
        if (dist_far.Dim() != Xwts_.Dim()) {
            dist_far.ReInit(Xwts_.Dim()); // TODO: one far-eval threshold for each node on object? Or one for each element?
        }

        dist_far.SetZero();
    }

    template <class Real> Vector<Real> PlaneIntegral<Real>::glNodes(const Integer Order) const {
        return gl_nodes_;
    }

    template <class Real> Long PlaneIntegral<Real>::Size() const {
        return Nelem_x_ * Nelem_y_ * 2;
    }

    template <class Real> void PlaneIntegral<Real>::GetVTUData(sctl::VTUData& vtu_data, const sctl::Vector<Real>& F) const {
        for (const auto& x : Xsrc_) vtu_data.coord.PushBack((float)x);
        for (const auto& x :     F) vtu_data.value.PushBack((float)x);
        // sctl::Long N = Nelem_x_*order_;
        // for (sctl::Long j = 0; j < N-1; j++) {
        //     for (sctl::Long k = 0; k < N-1; k++) {
        //         auto idx = [this,N](sctl::Long j, sctl::Long k) {
        //             return j*N+k;
        //         };
        //         vtu_data.connect.PushBack(idx(j+0,k+0));
        //         vtu_data.connect.PushBack(idx(j+0,k+1));
        //         vtu_data.connect.PushBack(idx(j+1,k+1));
        //         vtu_data.connect.PushBack(idx(j+1,k+0));
        //         vtu_data.connect.PushBack(idx(j+0,k+0));
        //         vtu_data.connect.PushBack(idx(j+0,k+1));
        //         vtu_data.connect.PushBack(idx(j+1,k+1));
        //         vtu_data.connect.PushBack(idx(j+1,k+0));
        //         vtu_data.offset.PushBack(vtu_data.connect.Dim());;
        //         vtu_data.types.PushBack(12);
        //     }
        // }
    }

    template <class Real> void PlaneIntegral<Real>::WriteVTK(const std::string& fname, const sctl::Vector<Real>& F, sctl::Comm comm) const {
        sctl::VTUData vtu_data;
        GetVTUData(vtu_data, F);
        vtu_data.WriteVTK(fname, comm);
    }

}
