#ifndef CONVECTIVE_FLUX_RECONSTRUCTOR_SECOND_ORDER_HLLC_HPP
#define CONVECTIVE_FLUX_RECONSTRUCTOR_SECOND_ORDER_HLLC_HPP

#include "flow/convective_flux_reconstructors/ConvectiveFluxReconstructor.hpp"

class ConvectiveFluxReconstructorSecondOrderHLLC: public ConvectiveFluxReconstructor
{
    public:
        ConvectiveFluxReconstructorSecondOrderHLLC(
            const std::string& object_name,
            const tbox::Dimension& dim,
            const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
            const int& num_eqn,
            const boost::shared_ptr<FlowModel>& flow_model,
            const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db);
        
        ~ConvectiveFluxReconstructorSecondOrderHLLC() {}
        
        /*
         * Print all characteristics of the convective flux reconstruction class.
         */
        void
        printClassData(std::ostream& os) const;
        
        /*
         * Put the characteristics of the convective flux reconstruction class
         * into the restart database.
         */
        void
        putToRestart(
            const boost::shared_ptr<tbox::Database>& restart_db) const;
        
        /*
         * Compute the convective flux and source due to splitting of convective term on a patch.
         */
        void computeConvectiveFluxAndSourceOnPatch(
            hier::Patch& patch,
            const boost::shared_ptr<pdat::SideVariable<double> >& variable_convective_flux,
            const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
            const boost::shared_ptr<hier::VariableContext>& data_context,
            const double time,
            const double dt,
            const int RK_step_number);
        
    private:
        std::vector<EQN_FORM::TYPE> d_eqn_form;
        bool d_has_advective_eqn_form;
        
        /*
         * boost::shared_ptr to the Riemann solver object.
         */
        boost::shared_ptr<FlowModelRiemannSolver> d_riemann_solver;

        //void Limiter(const double v_ll, const double v_l, const double v_r, const double v_rrdouble & v_min, double & v_plus)
        
};

#endif /* CONVECTIVE_FLUX_RECONSTRUCTOR_SECOND_ORDER_HLLC_HPP */
