#include "flow/nonconservative_diffusive_flux_divergence_operators/sixth_order/NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "util/basic_geometry/WallTreatment.hpp"
#include <map>

NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& nonconservative_diffusive_flux_divergence_operator_db):
        NonconservativeDiffusiveFluxDivergenceOperator(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model,
            nonconservative_diffusive_flux_divergence_operator_db)
{
    d_num_diff_ghosts = hier::IntVector::getOne(d_dim)*3;
}


/*
 * Print all characteristics of the non-conservative diffusive flux divergence operator class.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::printClassData(
    std::ostream& os) const
{
    os << "\nPrint NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder object..."
       << std::endl;
    
    os << std::endl;
    
    os << "NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder: this = "
       << (NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Put the characteristics of the non-conservative diffusive flux divergence operator class into
 * the restart database.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putString("d_nonconservative_diffusive_flux_divergence_operator", "SIXTH_ORDER");
}


/*
 * Compute the non-conservative diffusive flux divergence on a patch.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeNonconservativeDiffusiveFluxDivergenceOnPatch(
    hier::Patch& patch,
    const boost::shared_ptr<pdat::CellVariable<double> >& variable_diffusive_flux_divergence,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(variable_diffusive_flux_divergence);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the cell data of diffusive flux divergence.
    boost::shared_ptr<pdat::CellData<double> > diffusive_flux_divergence(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(variable_diffusive_flux_divergence, data_context)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(diffusive_flux_divergence);
    TBOX_ASSERT(diffusive_flux_divergence->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Initialize the data of diffusive flux to zero.
    diffusive_flux_divergence->fillAll(double(0));

    /*
     * Register the patch and derived cell variables in the flow model and compute the corresponding cell data.
     */

    d_flow_model->registerPatchWithDataContext(patch, data_context);

    /*
     * Compute the cell status.
     */

    boost::shared_ptr<pdat::CellData<double> > cell_status
            = d_flow_model->getGlobalCellStatus();


    d_flow_model->registerDiffusiveFlux(d_num_diff_ghosts);

    d_flow_model->computeGlobalDerivedCellData();

    boost::shared_ptr<pdat::CellData<double> > velocity =
            d_flow_model->getGlobalCellData("VELOCITY");

    boost::shared_ptr<pdat::CellData<double> > temperature =
            d_flow_model->getGlobalCellData("TEMPERATURE");


    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimension and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        

        
        /*
         * Delcare containers for computing flux derivatives in different directions.
         */
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_x;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_x;
        
        std::vector<std::vector<int> > var_component_idx_x;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_x;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xx;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_x;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_x_computed;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xx_computed;
        
        /*
         * Compute the derivatives for diffusive flux in x-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        // Get the diffusivities for the terms in x-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        
        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);
        
        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x);
        
        // Compute the second derivatives of variables in x-direction.
        computeSecondDerivativesInX(
            patch,
            var_derivative_xx,
            derivative_xx_computed,
            var_data_x,
            var_component_idx_x);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_xx[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudx = var_derivative_x[ei][vi]->getPointer(0);
                double* d2udxdx = var_derivative_xx[ei][vi]->getPointer(0);
                double* dmudx = diffusivities_derivative_x[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0; i++)
                {
                    // Compute the linear indices.
                    const int idx_nghost = i;
                    
                    const int idx_diffusivity = i + num_subghosts_0_diffusivity;
                    
                    const int idx_diff = i + num_diff_ghosts_0;
                    
                    nabla_F[idx_nghost] += dt*(dmudx[idx_diff]*dudx[idx_diff] +
                        mu[idx_diffusivity]*d2udxdx[idx_diff]);
                }
            }
        }
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xx.clear();
        
//        /*
//         * Unregister the patch and data of all registered derived cell variables in the flow model.
//         */
//
//        d_flow_model->unregisterPatch();
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];

        
        /*
         * Delcare containers for computing flux derivatives in different directions.
         */
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_y;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_y;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        
        std::vector<std::vector<int> > var_derivative_component_idx_x;
        std::vector<std::vector<int> > var_derivative_component_idx_y;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_y;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xx;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xy;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_yx;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_yy;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_y;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_x_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_y_computed;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xx_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xy_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_yx_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_yy_computed;

        /*
         * Build ghost cell map
         */
        std::vector<std::vector<std::array<int,3> > > ghost_cell_maps;
        ghost_cell_maps.resize(2);
        buildGhostCellMap2D(cell_status, ghost_cell_maps);
        
        /*
         * (1) Compute the derivatives for diffusive flux in x-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell2D(velocity, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        
        // Get the diffusivities for the terms in x-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION, true);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell2D(velocity, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell2D(temperature, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        
        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);

        
        // Compute the second derivatives of variables in x-direction.
        computeSecondDerivativesInX(
            patch,
            var_derivative_xx,
            derivative_xx_computed,
            var_data_x,
            var_component_idx_x);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_xx[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudx = var_derivative_x[ei][vi]->getPointer(0);
                double* d2udxdx = var_derivative_xx[ei][vi]->getPointer(0);
                double* dmudx = diffusivities_derivative_x[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_nghost = i +
                            j*(interior_dim_0);
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        const int idx_diff = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        nabla_F[idx_nghost] += dt*(dmudx[idx_diff]*dudx[idx_diff] +
                            mu[idx_diffusivity]*d2udxdx[idx_diff]);
                    }
                }
            }
        }
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xx.clear();
        
        // Get the variables for the terms in y-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);

        /*
        * mirror ghost cell data for computing the direvative in the x-direction.
        */
        mirrorGhostCell2D(velocity, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        
        // Get the diffusivities for the terms in y-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION, true);//todo need recompute???
        
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y);


        /*
        * mirror ghost cell data for computing the direvative in the y-direction.
        */
        mirrorGhostCell2D(velocity, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell2D(temperature, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);

        //Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            patch,
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);

        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_y.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_y[ei].resize(var_derivative_y[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_y[ei].size()); vi++)
            {
                var_derivative_component_idx_y[ei][vi] = 0;
            }
        }

        {
            /*
             * Mirror the derivatives in x-direction to compute du/dydx in the x direction.
             */
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_y;
            du_y.reserve(2); //du/dy, dv/dy
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
                du_y.push_back(var_derivative_y[ei][vi]);
            mirrorGhostCellDerivative2D(du_y, ghost_cell_maps, DIRECTION::X_DIRECTION);
        }
        
        computeFirstDerivativesInX(
            patch,
            var_derivative_xy,
            derivative_xy_computed,
            var_derivative_y,
            var_derivative_component_idx_y);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_xy[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudy = var_derivative_y[ei][vi]->getPointer(0);
                double* d2udxdy = var_derivative_xy[ei][vi]->getPointer(0);
                double* dmudx = diffusivities_derivative_x[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_nghost = i +
                            j*(interior_dim_0);
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        const int idx_diff = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        nabla_F[idx_nghost] += dt*(dmudx[idx_diff]*dudy[idx_diff] +
                            mu[idx_diffusivity]*d2udxdy[idx_diff]);
                    }
                }
            }
        }
        
        var_data_y.clear();
        
        diffusivities_data_y.clear();
        
        var_component_idx_y.clear();
        
        diffusivities_component_idx_y.clear();
        
        var_derivative_component_idx_y.clear();
        
        var_derivative_y.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xy.clear();

        
        /*
         * (2) Compute the derivatives for diffusive flux in y-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in y-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell2D(velocity, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);


        // Get the diffusivities for the terms in x-direction in the diffusive flux in y-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION, true);

        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
             patch,
             diffusivities_derivative_y,
             derivative_y_computed,
             diffusivities_data_x,
             diffusivities_component_idx_x);



        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell2D(velocity, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell2D(temperature, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);

        
        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);

        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_x.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_x[ei].resize(var_derivative_x[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_x[ei].size()); vi++)
            {
                var_derivative_component_idx_x[ei][vi] = 0;
            }
        }
        {
            /*
             * Mirror the derivatives in y-direction to compute flux du/dxdy in the y direction.
             */
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_x;
            du_x.reserve(2); //du/dx, dv/dx
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
                du_x.push_back(var_derivative_x[ei][vi]);
            mirrorGhostCellDerivative2D(du_x, ghost_cell_maps, DIRECTION::Y_DIRECTION);
        }
        
        computeFirstDerivativesInY(
            patch,
            var_derivative_yx,
            derivative_yx_computed,
            var_derivative_x,
            var_derivative_component_idx_x);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_yx[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_y[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudx = var_derivative_x[ei][vi]->getPointer(0);
                double* d2udydx = var_derivative_yx[ei][vi]->getPointer(0);
                double* dmudy = diffusivities_derivative_y[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_nghost = i +
                            j*(interior_dim_0);
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        const int idx_diff = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        nabla_F[idx_nghost] += dt*(dmudy[idx_diff]*dudx[idx_diff] +
                            mu[idx_diffusivity]*d2udydx[idx_diff]);
                    }
                }
            }
        }
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_y.clear();
        var_derivative_yx.clear();
        
        // Get the variables for the terms in y-direction in the diffusive flux in y-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);

        /*
        * mirror ghost cell data for computing the direvative in the x-direction.
        */
        mirrorGhostCell2D(velocity, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        
        // Get the diffusivities for the terms in y-direction in the diffusive flux in y-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION, true);
        
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);


        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            patch,
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y);

        /*
        * mirror ghost cell data for computing the direvative in the y-direction.
        */
        mirrorGhostCell2D(velocity, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell2D(temperature, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        
        //Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            patch,
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);

        // Compute the second derivatives of variables in y-direction.
        computeSecondDerivativesInY(
            patch,
            var_derivative_yy,
            derivative_yy_computed,
            var_data_y,
            var_component_idx_y);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_yy[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudy = var_derivative_y[ei][vi]->getPointer(0);
                double* d2udydy = var_derivative_yy[ei][vi]->getPointer(0);
                double* dmudy = diffusivities_derivative_y[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                
                for (int j = 0; j < interior_dim_1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_nghost = i +
                            j*(interior_dim_0);
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        const int idx_diff = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        nabla_F[idx_nghost] += dt*(dmudy[idx_diff]*dudy[idx_diff] +
                            mu[idx_diffusivity]*d2udydy[idx_diff]);
                    }
                }
            }
        }
        
        var_data_y.clear();
        
        diffusivities_data_y.clear();
        
        var_component_idx_y.clear();
        
        diffusivities_component_idx_y.clear();
        
        var_derivative_y.clear();
        diffusivities_derivative_y.clear();
        var_derivative_yy.clear();

        
//        /*
//         * Unregister the patch and data of all registered derived cell variables in the flow model.
//         */
//
//        d_flow_model->unregisterPatch();
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        /*
         * Delcare containers for computing flux derivatives in different directions.
         */
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_y;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_z;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_y;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_z;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        std::vector<std::vector<int> > var_component_idx_z;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        std::vector<std::vector<int> > diffusivities_component_idx_z;
        
        std::vector<std::vector<int> > var_derivative_component_idx_x;
        std::vector<std::vector<int> > var_derivative_component_idx_y;
        std::vector<std::vector<int> > var_derivative_component_idx_z;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_y;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_z;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xx;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xy;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_xz;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_yx;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_yy;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_yz;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_zx;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_zy;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_derivative_zz;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_y;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_derivative_z;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_x_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_y_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_z_computed;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xx_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xy_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_xz_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_yx_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_yy_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_yz_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_zx_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_zy_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_zz_computed;

        /*
         * Build ghost cell map
         */
        std::vector<std::vector<std::array<int,4> > > ghost_cell_maps;
        ghost_cell_maps.resize(3);
        buildGhostCellMap3D(cell_status, ghost_cell_maps);
        
        /*
         * (1) Compute the derivatives for diffusive flux in x-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);

        /*
        * mirror ghost cell data for computing the derivative of diffusivities in the x-direction.
        */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        
        // Get the diffusivities for the terms in x-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION,true);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x);

        /*
       * mirror ghost cell data for computing the derivative in the x-direction.
       */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell3D(temperature, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        
        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);
        
        // Compute the second derivatives of variables in x-direction.
        computeSecondDerivativesInX(
            patch,
            var_derivative_xx,
            derivative_xx_computed,
            var_data_x,
            var_component_idx_x);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_xx[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudx = var_derivative_x[ei][vi]->getPointer(0);
                double* d2udxdx = var_derivative_xx[ei][vi]->getPointer(0);
                double* dmudx = diffusivities_derivative_x[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_diff = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            nabla_F[idx_nghost] += dt*(dmudx[idx_diff]*dudx[idx_diff] +
                                mu[idx_diffusivity]*d2udxdx[idx_diff]);
                        }
                    }
                }
            }
        }
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xx.clear();
        
        // Get the variables for the terms in y-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);

        /*
        * mirror ghost cell data for computing the derivative in the x-direction.
        */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        
        // Get the diffusivities for the terms in y-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION, true);
        
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y);

        /*
         * mirror ghost cell data for computing the derivative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell3D(temperature, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);


        //Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            patch,
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);
        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_y.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_y[ei].resize(var_derivative_y[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_y[ei].size()); vi++)
            {
                var_derivative_component_idx_y[ei][vi] = 0;
            }
        }

        {
            /*
             * Mirror the derivatives in x-direction to compute ddu/dydx in the x direction.
             */
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_y;
            du_y.reserve(2); //du/dy, dv/dy
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
                du_y.push_back(var_derivative_y[ei][vi]);
            mirrorGhostCellDerivative3D(du_y, ghost_cell_maps, DIRECTION::X_DIRECTION);
        }
        
        computeFirstDerivativesInX(
            patch,
            var_derivative_xy,
            derivative_xy_computed,
            var_derivative_y,
            var_derivative_component_idx_y);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_xy[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudy = var_derivative_y[ei][vi]->getPointer(0);
                double* d2udxdy = var_derivative_xy[ei][vi]->getPointer(0);
                double* dmudx = diffusivities_derivative_x[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_diff = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            nabla_F[idx_nghost] += dt*(dmudx[idx_diff]*dudy[idx_diff] +
                                mu[idx_diffusivity]*d2udxdy[idx_diff]);
                        }
                    }
                }
            }
        }
        
        var_data_y.clear();
        
        diffusivities_data_y.clear();
        
        var_component_idx_y.clear();
        
        diffusivities_component_idx_y.clear();
        
        var_derivative_component_idx_y.clear();
        
        var_derivative_y.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xy.clear();
        
        // Get the variables for the terms in z-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::X_DIRECTION,
            DIRECTION::Z_DIRECTION);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);

        // Get the diffusivities for the terms in z-direction in the diffusive flux in x-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::X_DIRECTION,
            DIRECTION::Z_DIRECTION, true);
        
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in x-direction.
        computeFirstDerivativesInX(
            patch,
            diffusivities_derivative_x,
            derivative_x_computed,
            diffusivities_data_z,
            diffusivities_component_idx_z);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell3D(temperature, ghost_cell_maps, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);


        //Compute the first derivatives of variables in z-direction.
        computeFirstDerivativesInZ(
            patch,
            var_derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z);
        

        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_z.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_z[ei].resize(var_derivative_z[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_z[ei].size()); vi++)
            {
                var_derivative_component_idx_z[ei][vi] = 0;
            }
        }

        {
            /*
             * Mirror the derivatives in x-direction to compute ddu/dzdx in the x direction.
             */
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_z;
            du_z.reserve(2); //du/dz, dw/dz
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
                du_z.push_back(var_derivative_z[ei][vi]);
            mirrorGhostCellDerivative3D(du_z, ghost_cell_maps, DIRECTION::X_DIRECTION);
        }
        
        computeFirstDerivativesInX(
            patch,
            var_derivative_xz,
            derivative_xz_computed,
            var_derivative_z,
            var_derivative_component_idx_z);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_z[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_xz[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudz = var_derivative_z[ei][vi]->getPointer(0);
                double* d2udxdz = var_derivative_xz[ei][vi]->getPointer(0);
                double* dmudx = diffusivities_derivative_x[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_diff = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            nabla_F[idx_nghost] += dt*(dmudx[idx_diff]*dudz[idx_diff] +
                                mu[idx_diffusivity]*d2udxdz[idx_diff]);
                        }
                    }
                }
            }
        }
        
        var_data_z.clear();
        
        diffusivities_data_z.clear();
        
        var_component_idx_z.clear();
        
        diffusivities_component_idx_z.clear();
        
        var_derivative_component_idx_z.clear();
        
        var_derivative_z.clear();
        diffusivities_derivative_x.clear();
        var_derivative_xz.clear();
        
        /*
         * (2) Compute the derivatives for diffusive flux in y-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in y-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);
        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);

        // Get the diffusivities for the terms in x-direction in the diffusive flux in y-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION, true);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
                patch,
                diffusivities_derivative_y,
                derivative_y_computed,
                diffusivities_data_x,
                diffusivities_component_idx_x);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell3D(temperature, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);

        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);
        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_x.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_x[ei].resize(var_derivative_x[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_x[ei].size()); vi++)
            {
                var_derivative_component_idx_x[ei][vi] = 0;
            }
        }

        {
            /*
             * Mirror the derivatives in y-direction to compute ddu/dxdy in the y direction.
             */
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_x;
            du_x.reserve(2); //du/dy, dv/dy
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
                du_x.push_back(var_derivative_x[ei][vi]);
            mirrorGhostCellDerivative3D(du_x, ghost_cell_maps, DIRECTION::Y_DIRECTION);
        }
        
        computeFirstDerivativesInY(
            patch,
            var_derivative_yx,
            derivative_yx_computed,
            var_derivative_x,
            var_derivative_component_idx_x);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_yx[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_y[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudx = var_derivative_x[ei][vi]->getPointer(0);
                double* d2udydx = var_derivative_yx[ei][vi]->getPointer(0);
                double* dmudy = diffusivities_derivative_y[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_diff = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            nabla_F[idx_nghost] += dt*(dmudy[idx_diff]*dudx[idx_diff] +
                                mu[idx_diffusivity]*d2udydx[idx_diff]);
                        }
                    }
                }
            }
        }
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_y.clear();
        var_derivative_yx.clear();

        // Get the variables for the terms in y-direction in the diffusive flux in y-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        
        // Get the diffusivities for the terms in y-direction in the diffusive flux in y-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION, true);
        
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
                patch,
                diffusivities_derivative_y,
                derivative_y_computed,
                diffusivities_data_y,
                diffusivities_component_idx_y);

        /*
        * mirror ghost cell data for computing the direvative in the x-direction.
        */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell3D(temperature, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);

        //Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            patch,
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);
        
        // Compute the second derivatives of variables in y-direction.
        computeSecondDerivativesInY(
            patch,
            var_derivative_yy,
            derivative_yy_computed,
            var_data_y,
            var_component_idx_y);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_yy[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudy = var_derivative_y[ei][vi]->getPointer(0);
                double* d2udydy = var_derivative_yy[ei][vi]->getPointer(0);
                double* dmudy = diffusivities_derivative_y[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_diff = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            nabla_F[idx_nghost] += dt*(dmudy[idx_diff]*dudy[idx_diff] +
                                mu[idx_diffusivity]*d2udydy[idx_diff]);
                        }
                    }
                }
            }
        }
        
        var_data_y.clear();
        
        diffusivities_data_y.clear();
        
        var_component_idx_y.clear();
        
        diffusivities_component_idx_y.clear();
        
        var_derivative_y.clear();
        diffusivities_derivative_y.clear();
        var_derivative_yy.clear();
        
        // Get the variables for the terms in z-direction in the diffusive flux in y-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Z_DIRECTION);

        /*
        * mirror ghost cell data for computing the direvative in the x-direction.
        */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);

        // Get the diffusivities for the terms in z-direction in the diffusive flux in y-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Z_DIRECTION,true);
        
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in y-direction.
        computeFirstDerivativesInY(
            patch,
            diffusivities_derivative_y,
            derivative_y_computed,
            diffusivities_data_z,
            diffusivities_component_idx_z);


        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell3D(temperature, ghost_cell_maps, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);


        //Compute the first derivatives of variables in z-direction.
        computeFirstDerivativesInZ(
            patch,
            var_derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z);
        

        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_z.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_z[ei].resize(var_derivative_z[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_z[ei].size()); vi++)
            {
                var_derivative_component_idx_z[ei][vi] = 0;
            }
        }

        {
            /*
             * Mirror the derivatives in y-direction to compute ddu/dzdy in the y direction.
             */
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_z;
            du_z.reserve(2); //dv/dz, dw/dz
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
                du_z.push_back(var_derivative_z[ei][vi]);
            mirrorGhostCellDerivative3D(du_z, ghost_cell_maps, DIRECTION::Y_DIRECTION);
        }
        
        computeFirstDerivativesInY(
            patch,
            var_derivative_yz,
            derivative_yz_computed,
            var_derivative_z,
            var_derivative_component_idx_z);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_z[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_yz[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_y[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudz = var_derivative_z[ei][vi]->getPointer(0);
                double* d2udydz = var_derivative_yz[ei][vi]->getPointer(0);
                double* dmudy = diffusivities_derivative_y[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_diff = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            nabla_F[idx_nghost] += dt*(dmudy[idx_diff]*dudz[idx_diff] +
                                mu[idx_diffusivity]*d2udydz[idx_diff]);
                        }
                    }
                }
            }
        }
        
        var_data_z.clear();
        
        diffusivities_data_z.clear();
        
        var_component_idx_z.clear();
        
        diffusivities_component_idx_z.clear();
        
        var_derivative_component_idx_z.clear();
        
        var_derivative_z.clear();
        diffusivities_derivative_y.clear();
        var_derivative_yz.clear();
        
        /*
         * (3) Compute the derivatives for diffusive flux in z-direction.
         */
        
        // Get the variables for the terms in x-direction in the diffusive flux in z-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::Z_DIRECTION,
            DIRECTION::X_DIRECTION);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);

        // Get the diffusivities for the terms in x-direction in the diffusive flux in z-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Z_DIRECTION,
            DIRECTION::X_DIRECTION, true);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in z-direction.
        computeFirstDerivativesInZ(
            patch,
            diffusivities_derivative_z,
            derivative_z_computed,
            diffusivities_data_x,
            diffusivities_component_idx_x);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell3D(temperature, ghost_cell_maps, DIRECTION::X_DIRECTION, WALL_NO_SLIP);

        //Compute the first derivatives of variables in x-direction.
        computeFirstDerivativesInX(
            patch,
            var_derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);

        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_x.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_x[ei].resize(var_derivative_x[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_x[ei].size()); vi++)
            {
                var_derivative_component_idx_x[ei][vi] = 0;
            }
        }

        {
            /*
             * Mirror the derivatives in z-direction to compute ddu/dxdz in the z direction.
             */
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_x;
            du_x.reserve(2); //du/dx, dw/dx
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
                du_x.push_back(var_derivative_x[ei][vi]);
            mirrorGhostCellDerivative3D(du_x, ghost_cell_maps, DIRECTION::Z_DIRECTION);
        }
        
        computeFirstDerivativesInZ(
            patch,
            var_derivative_zx,
            derivative_zx_computed,
            var_derivative_x,
            var_derivative_component_idx_x);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_x[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_zx[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_z[ei].size()) ==
                        static_cast<int>(var_component_idx_x[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudx = var_derivative_x[ei][vi]->getPointer(0);
                double* d2udzdx = var_derivative_zx[ei][vi]->getPointer(0);
                double* dmudz = diffusivities_derivative_z[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_diff = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            nabla_F[idx_nghost] += dt*(dmudz[idx_diff]*dudx[idx_diff] +
                                mu[idx_diffusivity]*d2udzdx[idx_diff]);
                        }
                    }
                }
            }
        }
        
        var_data_x.clear();
        
        diffusivities_data_x.clear();
        
        var_component_idx_x.clear();
        
        diffusivities_component_idx_x.clear();
        
        var_derivative_component_idx_x.clear();
        
        var_derivative_x.clear();
        diffusivities_derivative_z.clear();
        var_derivative_zx.clear();
        
        // Get the variables for the terms in y-direction in the diffusive flux in z-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Y_DIRECTION);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);

        // Get the diffusivities for the terms in y-direction in the diffusive flux in z-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Y_DIRECTION, true);
        
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in z-direction.
        computeFirstDerivativesInZ(
            patch,
            diffusivities_derivative_z,
            derivative_z_computed,
            diffusivities_data_y,
            diffusivities_component_idx_y);


        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell3D(temperature, ghost_cell_maps, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);


        //Compute the first derivatives of variables in y-direction.
        computeFirstDerivativesInY(
            patch,
            var_derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);
        
        // Compute the mixed derivatives of variables.
        var_derivative_component_idx_y.resize(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            var_derivative_component_idx_y[ei].resize(var_derivative_y[ei].size());
            for (int vi = 0; vi < static_cast<int>(var_derivative_component_idx_y[ei].size()); vi++)
            {
                var_derivative_component_idx_y[ei][vi] = 0;
            }
        }

        {
            /*
             * Mirror the derivatives in z-direction to compute ddu/dydz in the z direction.
             */
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_y;
            du_y.reserve(2); //dv/dy, dw/dy
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
                du_y.push_back(var_derivative_y[ei][vi]);
            mirrorGhostCellDerivative3D(du_y, ghost_cell_maps, DIRECTION::Z_DIRECTION);
        }
        
        computeFirstDerivativesInZ(
            patch,
            var_derivative_zy,
            derivative_zy_computed,
            var_derivative_y,
            var_derivative_component_idx_y);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_y[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_zy[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_z[ei].size()) ==
                        static_cast<int>(var_component_idx_y[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudy = var_derivative_y[ei][vi]->getPointer(0);
                double* d2udzdy = var_derivative_zy[ei][vi]->getPointer(0);
                double* dmudz = diffusivities_derivative_z[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_y[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_diff = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            nabla_F[idx_nghost] += dt*(dmudz[idx_diff]*dudy[idx_diff] +
                                mu[idx_diffusivity]*d2udzdy[idx_diff]);
                        }
                    }
                }
            }
        }
        
        var_data_y.clear();
        
        diffusivities_data_y.clear();
        
        var_component_idx_y.clear();
        
        diffusivities_component_idx_y.clear();
        
        var_derivative_component_idx_y.clear();
        
        var_derivative_y.clear();
        diffusivities_derivative_z.clear();
        var_derivative_zy.clear();
        
        // Get the variables for the terms in z-direction in the diffusive flux in z-direction.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Z_DIRECTION);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);


        // Get the diffusivities for the terms in z-direction in the diffusive flux in z-direction.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Z_DIRECTION, true);
        
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);

        // Compute the first derivatives of diffusivities in z-direction.
        computeFirstDerivativesInZ(
            patch,
            diffusivities_derivative_z,
            derivative_z_computed,
            diffusivities_data_z,
            diffusivities_component_idx_z);

        /*
        * mirror ghost cell data for computing the direvative in the x-direction.
        */
        mirrorGhostCell3D(velocity, ghost_cell_maps, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell3D(temperature, ghost_cell_maps, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);

        //Compute the first derivatives of variables in z-direction.
        computeFirstDerivativesInZ(
            patch,
            var_derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z);
        
        // Compute the second derivatives of variables in z-direction.
        computeSecondDerivativesInZ(
            patch,
            var_derivative_zz,
            derivative_zz_computed,
            var_data_z,
            var_component_idx_z);
        
        /*
         * Add the derivatives to the divergence of diffusive flux.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_z[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(var_derivative_zz[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_derivative_z[ei].size()) ==
                        static_cast<int>(var_component_idx_z[ei].size()));
            
            double* nabla_F = diffusive_flux_divergence->getPointer(ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivatives.
                double* dudz = var_derivative_z[ei][vi]->getPointer(0);
                double* d2udzdz = var_derivative_zz[ei][vi]->getPointer(0);
                double* dmudz = diffusivities_derivative_z[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width and ghost box dimensions of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostCellWidth();
                
                hier::IntVector subghostcell_dims_diffusivity =
                    diffusivities_data_z[ei][vi]->getGhostBox().numberCells();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                const int num_subghosts_1_diffusivity = num_subghosts_diffusivity[1];
                const int num_subghosts_2_diffusivity = num_subghosts_diffusivity[2];
                const int subghostcell_dim_0_diffusivity = subghostcell_dims_diffusivity[0];
                const int subghostcell_dim_1_diffusivity = subghostcell_dims_diffusivity[1];
                
                for (int k = 0; k < interior_dim_2; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_nghost = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_diff = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            nabla_F[idx_nghost] += dt*(dmudz[idx_diff]*dudz[idx_diff] +
                                mu[idx_diffusivity]*d2udzdz[idx_diff]);
                        }
                    }
                }
            }
        }
        
        var_data_z.clear();
        
        diffusivities_data_z.clear();
        
        var_component_idx_z.clear();
        
        diffusivities_component_idx_z.clear();
        
        var_derivative_z.clear();
        diffusivities_derivative_z.clear();
        var_derivative_zz.clear();
        
//        /*
//         * Unregister the patch and data of all registered derived cell variables in the flow model.
//         */
//
//        d_flow_model->unregisterPatch();
    }
    /*
     * Unregister the patch and data of all registered derived cell variables in the flow model.
     */

    d_flow_model->unregisterPatch();
}


/*
 * Compute the first derivatives in the x-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeFirstDerivativesInX(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_x,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_x_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_x[ei].size()) ==
                    static_cast<int>(data_component_idx_x[ei].size()));
    }
#endif
    
    derivative_x.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_0 = dx[0];
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudx = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_diff_ghosts_0;
                        
                        const int idx_data_LLL = i - 3 + num_subghosts_0_data;
                        const int idx_data_LL  = i - 2 + num_subghosts_0_data;
                        const int idx_data_L   = i - 1 + num_subghosts_0_data;
                        const int idx_data_R   = i + 1 + num_subghosts_0_data;
                        const int idx_data_RR  = i + 2 + num_subghosts_0_data;
                        const int idx_data_RRR = i + 3 + num_subghosts_0_data;
                        
                        dudx[idx] = (double(3)/double(4)*(u[idx_data_R] - u[idx_data_L]) +
                                     double(-3)/double(20)*(u[idx_data_RR] - u[idx_data_LL]) +
                                     double(1)/double(60)*(u[idx_data_RRR] - u[idx_data_LLL]))/
                                        dx_0;
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudx = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = -3; j < interior_dim_1 + 3; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            dudx[idx] = (double(3)/double(4)*(u[idx_data_R] - u[idx_data_L]) +
                                         double(-3)/double(20)*(u[idx_data_RR] - u[idx_data_LL]) +
                                         double(1)/double(60)*(u[idx_data_RRR] - u[idx_data_LLL]))/
                                            dx_0;
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudx = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = -3; k < interior_dim_2 + 3; k++)
                    {
                        for (int j = -3; j < interior_dim_1 + 3; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                dudx[idx] = (double(3)/double(4)*(u[idx_data_R] - u[idx_data_L]) +
                                             double(-3)/double(20)*(u[idx_data_RR] - u[idx_data_LL]) +
                                             double(1)/double(60)*(u[idx_data_RRR] - u[idx_data_LLL]))/
                                                dx_0;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the first derivatives in the y-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeFirstDerivativesInY(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_y,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_y_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_y[ei].size()) ==
                    static_cast<int>(data_component_idx_y[ei].size()));
    }
#endif
    
    derivative_y.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_1 = dx[1];
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeFirstDerivativesInY()\n"
            << "There isn't y-direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudy = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = -3; i < interior_dim_0 + 3; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_BBB = (i + num_subghosts_0_data) +
                                (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BB = (i + num_subghosts_0_data) +
                                (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_B = (i + num_subghosts_0_data) +
                                (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_T = (i + num_subghosts_0_data) +
                                (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TT = (i + num_subghosts_0_data) +
                                (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTT = (i + num_subghosts_0_data) +
                                (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            dudy[idx] = (double(3)/double(4)*(u[idx_data_T] - u[idx_data_B]) +
                                         double(-3)/double(20)*(u[idx_data_TT] - u[idx_data_BB]) +
                                         double(1)/double(60)*(u[idx_data_TTT] - u[idx_data_BBB]))/
                                            dx_1;
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudy = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = -3; k < interior_dim_2 + 3; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = -3; i < interior_dim_0 + 3; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_T = (i + num_subghosts_0_data) +
                                    (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TT = (i + num_subghosts_0_data) +
                                    (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTT = (i + num_subghosts_0_data) +
                                    (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                dudy[idx] = (double(3)/double(4)*(u[idx_data_T] - u[idx_data_B]) +
                                             double(-3)/double(20)*(u[idx_data_TT] - u[idx_data_BB]) +
                                             double(1)/double(60)*(u[idx_data_TTT] - u[idx_data_BBB]))/
                                                dx_1;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the first derivatives in the z-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeFirstDerivativesInZ(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_z,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_z_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_z[ei].size()) ==
                    static_cast<int>(data_component_idx_z[ei].size()));
    }
#endif
    
    derivative_z.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_2 = dx[2];
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_z[ei].reserve(static_cast<int>(data_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_z[ei][vi];
                
                if (derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))
                    == derivative_z_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_z[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* dudz = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_z[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_z[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = -3; j < interior_dim_1 + 3; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = -3; i < interior_dim_0 + 3; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_F = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                dudz[idx] = (double(3)/double(4)*(u[idx_data_F] - u[idx_data_B]) +
                                             double(-3)/double(20)*(u[idx_data_FF] - u[idx_data_BB]) +
                                             double(1)/double(60)*(u[idx_data_FFF] - u[idx_data_BBB]))/
                                                dx_2;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_z_computed.insert(derivative_pair);
                }
                
                derivative_z[ei].push_back(
                    derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the second derivatives in the x-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeSecondDerivativesInX(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_x,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_x_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_x[ei].size()) ==
                    static_cast<int>(data_component_idx_x[ei].size()));
    }
#endif
    
    derivative_x.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dx_sq = dx[0]*dx[0];
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udx2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx = i + num_diff_ghosts_0;
                        
                        const int idx_data_LLL = i - 3 + num_subghosts_0_data;
                        const int idx_data_LL  = i - 2 + num_subghosts_0_data;
                        const int idx_data_L   = i - 1 + num_subghosts_0_data;
                        const int idx_data     = i     + num_subghosts_0_data;
                        const int idx_data_R   = i + 1 + num_subghosts_0_data;
                        const int idx_data_RR  = i + 2 + num_subghosts_0_data;
                        const int idx_data_RRR = i + 3 + num_subghosts_0_data;
                        
                        d2udx2[idx] = (double(-49)/double(18)*u[idx_data] +
                                       double(3)/double(2)*(u[idx_data_L] + u[idx_data_R]) +
                                       double(-3)/double(20)*(u[idx_data_LL] + u[idx_data_RR]) +
                                       double(1)/double(90)*(u[idx_data_LLL] + u[idx_data_RRR]))/
                                        dx_sq;
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udx2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            d2udx2[idx] = (double(-49)/double(18)*u[idx_data] +
                                           double(3)/double(2)*(u[idx_data_L] + u[idx_data_R]) +
                                           double(-3)/double(20)*(u[idx_data_LL] + u[idx_data_RR]) +
                                           double(1)/double(90)*(u[idx_data_LLL] + u[idx_data_RRR]))/
                                            dx_sq;
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_x[ei].reserve(static_cast<int>(data_x[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_x[ei][vi];
                
                if (derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))
                    == derivative_x_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_x[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udx2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_x[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_x[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_LLL = (i - 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_LL = (i - 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_L = (i - 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_R = (i + 1 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RR = (i + 2 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_RRR = (i + 3 + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                d2udx2[idx] = (double(-49)/double(18)*u[idx_data] +
                                               double(3)/double(2)*(u[idx_data_L] + u[idx_data_R]) +
                                               double(-3)/double(20)*(u[idx_data_LL] + u[idx_data_RR]) +
                                               double(1)/double(90)*(u[idx_data_LLL] + u[idx_data_RRR]))/
                                                dx_sq;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_x_computed.insert(derivative_pair);
                }
                
                derivative_x[ei].push_back(
                    derivative_x_computed.find(data_x[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the second derivatives in the y-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeSecondDerivativesInY(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_y,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_y_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_y[ei].size()) ==
                    static_cast<int>(data_component_idx_y[ei].size()));
    }
#endif
    
    derivative_y.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dy_sq = dx[1]*dx[1];
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeSecondDerivativesInY()\n"
            << "There isn't y-direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udy2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                            
                            const int idx_data_BBB = (i + num_subghosts_0_data) +
                                (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_BB = (i + num_subghosts_0_data) +
                                (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_B = (i + num_subghosts_0_data) +
                                (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data = (i + num_subghosts_0_data) +
                                (j + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_T = (i + num_subghosts_0_data) +
                                (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TT = (i + num_subghosts_0_data) +
                                (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            const int idx_data_TTT = (i + num_subghosts_0_data) +
                                (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data;
                            
                            d2udy2[idx] = (double(-49)/double(18)*u[idx_data] +
                                           double(3)/double(2)*(u[idx_data_B] + u[idx_data_T]) +
                                           double(-3)/double(20)*(u[idx_data_BB] + u[idx_data_TT]) +
                                           double(1)/double(90)*(u[idx_data_BBB] + u[idx_data_TTT]))/
                                            dy_sq;
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_y[ei].reserve(static_cast<int>(data_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_y[ei][vi];
                
                if (derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))
                    == derivative_y_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_y[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udy2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_y[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_y[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j - 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j - 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j - 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_T = (i + num_subghosts_0_data) +
                                    (j + 1 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TT = (i + num_subghosts_0_data) +
                                    (j + 2 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_TTT = (i + num_subghosts_0_data) +
                                    (j + 3 + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                d2udy2[idx] = (double(-49)/double(18)*u[idx_data] +
                                               double(3)/double(2)*(u[idx_data_B] + u[idx_data_T]) +
                                               double(-3)/double(20)*(u[idx_data_BB] + u[idx_data_TT]) +
                                               double(1)/double(90)*(u[idx_data_BBB] + u[idx_data_TTT]))/
                                                dy_sq;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_y_computed.insert(derivative_pair);
                }
                
                derivative_y[ei].push_back(
                    derivative_y_computed.find(data_y[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}


/*
 * Compute the second derivatives in the z-direction.
 */
void
NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::computeSecondDerivativesInZ(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_z,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_z_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei ++)
    {
        TBOX_ASSERT(static_cast<int>(data_z[ei].size()) ==
                    static_cast<int>(data_component_idx_z[ei].size()));
    }
#endif
    
    derivative_z.resize(d_num_eqn);
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    
    const double dz_sq = dx[2]*dx[2];
    
    if (d_dim == tbox::Dimension(1))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeSecondDerivativesInZ()\n"
            << "There isn't z-direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": NonconservativeDiffusiveFluxDivergenceOperatorSixthOrder::"
            << "computeSecondDerivativesInZ()\n"
            << "There isn't z-direction for 2D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the dimensions and number of ghost cells.
         */
        
        const int interior_dim_0 = interior_dims[0];
        const int interior_dim_1 = interior_dims[1];
        const int interior_dim_2 = interior_dims[2];
        
        const int num_diff_ghosts_0 = d_num_diff_ghosts[0];
        const int num_diff_ghosts_1 = d_num_diff_ghosts[1];
        const int num_diff_ghosts_2 = d_num_diff_ghosts[2];
        
        const int diff_ghostcell_dim_0 = diff_ghostcell_dims[0];
        const int diff_ghostcell_dim_1 = diff_ghostcell_dims[1];
        
        for (int ei = 0; ei < d_num_eqn; ei ++)
        {
            derivative_z[ei].reserve(static_cast<int>(data_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int u_idx = data_component_idx_z[ei][vi];
                
                if (derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))
                    == derivative_z_computed.end())
                {
                    // Get the pointer to variable for derivative.
                    double* u = data_z[ei][vi]->getPointer(u_idx);
                    
                    // Declare container to store the derivative.
                    boost::shared_ptr<pdat::CellData<double> > derivative(
                        new pdat::CellData<double>(
                            interior_box, 1, d_num_diff_ghosts));
                    
                    // Get the pointer to the derivative.
                    double* d2udz2 = derivative->getPointer(0);
                    
                    /*
                     * Get the sub-ghost cell width and ghost box dimensions of the variable.
                     */
                    
                    hier::IntVector num_subghosts_data =
                        data_z[ei][vi]->getGhostCellWidth();
                    
                    hier::IntVector subghostcell_dims_data =
                        data_z[ei][vi]->getGhostBox().numberCells();
                    
                    const int num_subghosts_0_data = num_subghosts_data[0];
                    const int num_subghosts_1_data = num_subghosts_data[1];
                    const int num_subghosts_2_data = num_subghosts_data[2];
                    const int subghostcell_dim_0_data = subghostcell_dims_data[0];
                    const int subghostcell_dim_1_data = subghostcell_dims_data[1];
                    
                    for (int k = 0; k < interior_dim_2; k++)
                    {
                        for (int j = 0; j < interior_dim_1; j++)
                        {
#ifdef HAMERS_ENABLE_SIMD
                            #pragma omp simd
#endif
                            for (int i = 0; i < interior_dim_0; i++)
                            {
                                // Compute the linear indices.
                                const int idx = (i + num_diff_ghosts_0) +
                                    (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                    (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                        diff_ghostcell_dim_1;
                                
                                const int idx_data_BBB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_BB = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_B = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k - 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_F = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 1 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 2 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                const int idx_data_FFF = (i + num_subghosts_0_data) +
                                    (j + num_subghosts_1_data)*subghostcell_dim_0_data +
                                    (k + 3 + num_subghosts_2_data)*subghostcell_dim_0_data*
                                        subghostcell_dim_1_data;
                                
                                d2udz2[idx] = (double(-49)/double(18)*u[idx_data] +
                                               double(3)/double(2)*(u[idx_data_B] + u[idx_data_F]) +
                                               double(-3)/double(20)*(u[idx_data_BB] + u[idx_data_FF]) +
                                               double(1)/double(90)*(u[idx_data_BBB] + u[idx_data_FFF]))/
                                                dz_sq;
                            }
                        }
                    }
                    
                    std::pair<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_pair(
                        u,
                        derivative);
                    
                    derivative_z_computed.insert(derivative_pair);
                }
                
                derivative_z[ei].push_back(
                    derivative_z_computed.find(data_z[ei][vi]->getPointer(u_idx))->
                        second);
            }
        }
    }
}
