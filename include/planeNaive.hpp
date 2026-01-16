#ifndef _PLANENAIVE_HPP_
#define _PLANENAIVE_HPP_

#include <sctl.hpp>

namespace sctl {
    // Element List that includes BOTH top and bottom planes with z_offset from z=0 and z=1 (symmetric for now).
    template <class Real> class PlaneIntegral : public ElementListBase<Real> { 
            static constexpr Long KDIM = 3;

        public: 

            PlaneIntegral() {}

            /**
             * @brief Construct a new Plane Integral object with panel-based Chebyshev discretization
             *        Assumptions: 
             *                  unit box: X,Y = [0,1] x [0,1]
             *                  uniform order on each panel
             *                  same discretization for top and bottom
             *                  normal vector points toward each other -- fluid domain in between so exterior problem on each plane.
             * 
             * @param order : of Chebyshev on each panel, in each variable (x and y)
             * @param Nelem_x : Number of panels in x 
             * @param Nelem_y : Number of panels in y
             * @param z_offset : to avoid touching the periodic boundary z=0,1, offset by z_offset symmetrically.
             */
            PlaneIntegral(const Long order, const Long Nelem_x, const Long Nelem_y, const Real z_offset);

            virtual ~PlaneIntegral() {}

            /**
             * @brief Number of panels on BOTH top and bottom
             * 
             * @return Long 
             */
            Long Size() const override;

            /**
             * @brief Number of nodes on EITHER top and bottom
             *          TODO: update to Order_x and Order_y, and change from one value to a vectro of values.
             * 
             * @return Long 
             */
            Long Order_const() const;

            /**
             * @brief Populate destination vectors with list of nodes, normals, and number of nodes on each panel, top plane first, dimensions fast node slow.
             * 
             * @param X 
             * @param Xn 
             * @param element_wise_node_cnt 
             */
            void GetNodeCoord(Vector<Real>* X, Vector<Real>* Xn, Vector<Long>* element_wise_node_cnt) const override;

            /**
             * @brief Get the Far Field ndoes; Currently just returns all nodes (i.e. no near eval)
             * 
             * @param X 
             * @param Xn 
             * @param wts 
             * @param dist_far 
             * @param element_wise_node_cnt 
             * @param tol 
             */
            void GetFarFieldNodes(Vector<Real>& X, Vector<Real>& Xn, Vector<Real>& wts, Vector<Real>& dist_far, Vector<Long>& element_wise_node_cnt, const Real tol) const override;

            Vector<Real> glNodes(const Integer Order) const;

            void GetVTUData(sctl::VTUData& vtu_data, const sctl::Vector<Real>& F) const;

            void WriteVTK(const std::string& fname, const sctl::Vector<Real>& F, sctl::Comm comm) const;

            template <class Kernel> static void SelfInterac(sctl::Vector<sctl::Matrix<Real>>& M_lst, const Kernel& ker, Real tol, bool trg_dot_prod, const sctl::ElementListBase<Real>* self);

        private: 
            static const std::pair<Vector<Real>,Vector<Real>>& LegendreQuad_plane(Integer ORDER);
            static void reg_sl(sctl::Matrix<Real>& SL_eps, const sctl::Vector<Real> Xsrc, const Real eps);
            static void subtr_sl(sctl::Matrix<Real>& SL_subtr, const sctl::Vector<Real> Xsrc);

            Long order_, Nelem_x_, Nelem_y_; // order of gl grids on each element in both directions; number of elements in the x and y directions on EACH plane.
            Real z_offset_; // assume two symmetric flat planes, so only specify z offset only.
            Vector<Real> gl_nodes_, gl_wts_; // 1D nodes (and weights) of gaussian order order_
            Vector<Real> Xsrc_, Xsrc_n_, Xwts_; // collection of all nodes on square (both planes)
            // Vector<Real> dens_; // of dimension KDIM * N_nodes * N_nodes
    };
}

#include <planeNaive.cpp>

#endif