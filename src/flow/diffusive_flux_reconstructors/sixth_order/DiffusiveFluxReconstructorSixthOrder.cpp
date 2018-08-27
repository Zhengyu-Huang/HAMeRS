#include "flow/diffusive_flux_reconstructors/sixth_order/DiffusiveFluxReconstructorSixthOrder.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "util/basic_geometry/WallTreatment.hpp"
#include <map>

DiffusiveFluxReconstructorSixthOrder::DiffusiveFluxReconstructorSixthOrder(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& diffusive_flux_reconstructor_db):
        DiffusiveFluxReconstructor(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model,
            diffusive_flux_reconstructor_db)
{
    d_num_diff_ghosts = hier::IntVector::getOne(d_dim)*6;
}


/*
 * Print all characteristics of the diffusive flux reconstruction class.
 */
void
DiffusiveFluxReconstructorSixthOrder::printClassData(
    std::ostream& os) const
{
    os << "\nPrint DiffusiveFluxReconstructorSixthOrder object..."
       << std::endl;
    
    os << std::endl;
    
    os << "DiffusiveFluxReconstructorSixthOrder: this = "
       << (DiffusiveFluxReconstructorSixthOrder *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Put the characteristics of the diffusive flux reconstruction class
 * into the restart database.
 */
void
DiffusiveFluxReconstructorSixthOrder::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putString("d_diffusive_flux_reconstructor", "SIXTH_ORDER");
}


/*
 * Compute the diffusive flux on a patch.
 */
void
DiffusiveFluxReconstructorSixthOrder::computeDiffusiveFluxOnPatch(
    hier::Patch& patch,
    const boost::shared_ptr<pdat::SideVariable<double> >& variable_diffusive_flux,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(variable_diffusive_flux);
#endif
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // diffusive ghost cells.
    hier::Box diff_ghost_box = interior_box;
    diff_ghost_box.grow(d_num_diff_ghosts);
    const hier::IntVector diff_ghostcell_dims = diff_ghost_box.numberCells();
    
    // Get the side data of diffusive flux.
    boost::shared_ptr<pdat::SideData<double> > diffusive_flux(
        BOOST_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(variable_diffusive_flux, data_context)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(diffusive_flux);
    TBOX_ASSERT(diffusive_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Initialize the data of diffusive flux to zero.
    diffusive_flux->fillAll(double(0));

    /*
     * Compute the cell status.
     */
    d_flow_model->registerPatchWithDataContext(patch, data_context);

    boost::shared_ptr<pdat::CellData<double> > cell_status
            = d_flow_model->getGlobalCellStatus();

    double* cell_status_data = cell_status->getPointer(0);
    const hier::IntVector num_ghosts_cell_status = cell_status->getGhostCellWidth();
    const hier::IntVector ghostcell_dims_cell_status = cell_status->getGhostBox().numberCells();

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
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_x;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_x;
        
        std::vector<std::vector<int> > var_component_idx_x;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > derivative_x;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_x_computed;
        
        /*
         * (1) Compute the flux in the x-direction.
         */
        
        // Get the variables for the derivatives in the diffusive flux.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        // Get the diffusivities in the diffusive flux.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        
        /*
         * Compute the derivatives in x-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInX(
            patch,
            derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);
        
        /*
         * Reconstruct the flux in x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            double* F_face_x = diffusive_flux->getPointer(0, ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivative_x[ei][vi]->getPointer(0);
                
                /*
                 * Get the sub-ghost cell width of the diffusivity.
                 */
                
                hier::IntVector num_subghosts_diffusivity =
                    diffusivities_data_x[ei][vi]->getGhostCellWidth();
                
                const int num_subghosts_0_diffusivity = num_subghosts_diffusivity[0];
                
#ifdef HAMERS_ENABLE_SIMD
                #pragma omp simd
#endif
                for (int i = 0; i < interior_dim_0 + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_x = i;
                    const int idx_diffusivity = i + num_subghosts_0_diffusivity;
                    const int idx_node_LLL = i - 3 + num_diff_ghosts_0;
                    const int idx_node_LL  = i - 2 + num_diff_ghosts_0;
                    const int idx_node_L   = i - 1 + num_diff_ghosts_0;
                    const int idx_node_R   = i + num_diff_ghosts_0;
                    const int idx_node_RR  = i + 1 + num_diff_ghosts_0;
                    const int idx_node_RRR = i + 2 + num_diff_ghosts_0;
                    
                    F_face_x[idx_face_x] += dt*mu[idx_diffusivity]*(
                        double(37)/double(60)*(dudx[idx_node_L] + dudx[idx_node_R]) +
                        double(-2)/double(15)*(dudx[idx_node_LL] + dudx[idx_node_RR]) +
                        double(1)/double(60)*(dudx[idx_node_LLL] + dudx[idx_node_RRR]));
                }
            }
        }
        
        var_data_x.clear();
        diffusivities_data_x.clear();
        var_component_idx_x.clear();
        diffusivities_component_idx_x.clear();
        derivative_x.clear();

        
    } // if (d_dim == tbox::Dimension(1))
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
         * Delcare containers for computing fluxes in different directions.
         */
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > var_data_y;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > diffusivities_data_y;
        
        std::vector<std::vector<int> > var_component_idx_x;
        std::vector<std::vector<int> > var_component_idx_y;
        
        std::vector<std::vector<int> > diffusivities_component_idx_x;
        std::vector<std::vector<int> > diffusivities_component_idx_y;
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > derivative_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > derivative_y;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_x_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_y_computed;
        
        /*
         * (1) Compute the flux in the x-direction.
         */
        
        // Get the variables for the derivatives in the diffusive flux.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);


        /*
         * mirror ghost cell data for computing flux in the x-direction.
         */
        mirrorGhostCell(velocity, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);

        // Get the diffusivities in the diffusive flux.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION, true);
        
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */

        mirrorGhostCell(velocity, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);

        /*
         * Compute the derivatives in x-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInX(
            patch,
            derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);


        /*
        * mirror ghost cell data for computing the direvative in the y-direction.
        */

        mirrorGhostCell(velocity, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInY(
            patch,
            derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);

        /*
         * Reconstruct the flux in x-direction.
         */

        {
           /*
            * Mirror the derivatives in x-direction to compute flux in the x direction.
            */
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_x;
            du_x.reserve(5); //du/dx, dv/dx, dT/dx, du/dy, dv/dy,
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
                du_x.push_back(derivative_x[d_num_eqn - 1][vi]);
            for(int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
                du_x.push_back(derivative_y[d_num_eqn - 1][vi]);
            mirrorGhostCellDerivative(du_x, cell_status, DIRECTION::X_DIRECTION);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            double* F_face_x = diffusive_flux->getPointer(0, ei);

            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++) {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];

                // Get the pointer to diffusivity.
                double *mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);

                // Get the pointer to derivative.
                double *dudx = derivative_x[ei][vi]->getPointer(0);

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

                for (int j = 0; j < interior_dim_1; j++) {
#ifdef HAMERS_ENABLE_SIMD
#pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0 + 1; i++) {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                                               j * (interior_dim_0 + 1);

                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                                    (j + num_subghosts_1_diffusivity) * subghostcell_dim_0_diffusivity;

                        const int idx_node_LLL = (i - 3 + num_diff_ghosts_0) +
                                                 (j + num_diff_ghosts_1) * diff_ghostcell_dim_0;

                        const int idx_node_LL = (i - 2 + num_diff_ghosts_0) +
                                                (j + num_diff_ghosts_1) * diff_ghostcell_dim_0;

                        const int idx_node_L = (i - 1 + num_diff_ghosts_0) +
                                               (j + num_diff_ghosts_1) * diff_ghostcell_dim_0;

                        const int idx_node_R = (i + num_diff_ghosts_0) +
                                               (j + num_diff_ghosts_1) * diff_ghostcell_dim_0;

                        const int idx_node_RR = (i + 1 + num_diff_ghosts_0) +
                                                (j + num_diff_ghosts_1) * diff_ghostcell_dim_0;

                        const int idx_node_RRR = (i + 2 + num_diff_ghosts_0) +
                                                 (j + num_diff_ghosts_1) * diff_ghostcell_dim_0;

                        F_face_x[idx_face_x] += dt * (
                                double(37) / double(60) *
                                (mu[idx_node_L] * dudx[idx_node_L] + mu[idx_node_R] * dudx[idx_node_R]) +
                                double(-2) / double(15) *
                                (mu[idx_node_LL] * dudx[idx_node_LL] + mu[idx_node_RR] * dudx[idx_node_RR]) +
                                double(1) / double(60) *
                                (mu[idx_node_LLL] * dudx[idx_node_LLL] + mu[idx_node_RRR] * dudx[idx_node_RRR]));

                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivative_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivative_y[ei][vi]->getPointer(0);
                
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
                    for (int i = 0; i < interior_dim_0 + 1; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                            j*(interior_dim_0 + 1);
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        const int idx_node_LLL = (i - 3 + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_LL  = i - 2 + num_diff_ghosts_0 +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_L   = i - 1 + num_diff_ghosts_0 +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_R   = i + num_diff_ghosts_0 +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_RR  = i + 1 + num_diff_ghosts_0 +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_RRR = i + 2 + num_diff_ghosts_0 +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;

                        F_face_x[idx_face_x] += dt*(
                            double(37)/double(60)*(mu[idx_node_L]*dudy[idx_node_L] + mu[idx_node_R]*dudy[idx_node_R]) +
                            double(-2)/double(15)*(mu[idx_node_LL]*dudy[idx_node_LL] + mu[idx_node_RR]*dudy[idx_node_RR]) +
                            double(1)/double(60)*(mu[idx_node_LLL]*dudy[idx_node_LLL] + mu[idx_node_RRR]*dudy[idx_node_RRR]));


//                        if ((cell_status_data[idx_node_L] < 0.5 && cell_status_data[idx_node_R] > 0.5) ||
//                            (cell_status_data[idx_node_L] > 0.5 && cell_status_data[idx_node_R] < 0.5)) {
//                            if ((ei <= d_num_eqn - 2) &&(fabs(mu[idx_node_L]*dudy[idx_node_L] - mu[idx_node_R]*dudy[idx_node_R]) > 1e-8 ||
//                                                         fabs(mu[idx_node_LL]*dudy[idx_node_LL] - mu[idx_node_RR]*dudy[idx_node_RR]) > 1e-8 ||
//                                                         fabs(mu[idx_node_LLL]*dudy[idx_node_LLL] - mu[idx_node_LLL]*dudy[idx_node_RRR]) > 1e-8)) {
//                                std::cout << "ei " << ei << " vi " << vi << " i " << i << " j " << j  << std::endl;
//                                std::cout << " num_ghosts_cell_status[0] " << num_ghosts_cell_status[0]
//                                          << " num_subghosts_0_diffusivity " << num_subghosts_0_diffusivity
//                                          << " num_diff_ghosts_0 " << num_diff_ghosts_0
//                                          << " diff_ghostcell_dim_0 " << diff_ghostcell_dim_0
//                                          << " diff_ghostcell_dim_1 "<< diff_ghostcell_dims[1] <<std::endl;
//
//                                std::cout << (velocity->getPointer(1))[idx_node_R + 3*diff_ghostcell_dim_0] << " " << (velocity->getPointer(1))[idx_node_R + 2*diff_ghostcell_dim_0] << " "
//                                          << (velocity->getPointer(1))[idx_node_R + 1*diff_ghostcell_dim_0] << " " << (velocity->getPointer(1))[idx_node_R - 1*diff_ghostcell_dim_0] << " "
//                                          << (velocity->getPointer(1))[idx_node_R - 2*diff_ghostcell_dim_0] << " " << (velocity->getPointer(1))[idx_node_R - 3*diff_ghostcell_dim_0] << std::endl;
//
//                                std::cout << "dudy " << dudy[idx_node_L] << " " << dudy[idx_node_R] << " " << dudy[idx_node_LL]
//                                          << " " << dudy[idx_node_RR] << " " << dudy[idx_node_LLL] << " "
//                                          << dudy[idx_node_RRR] << std::endl;
//                                std::cout << "mu " << mu[idx_node_L] << " " << mu[idx_node_R] << " " << mu[idx_node_LL]
//                                          << " " << mu[idx_node_RR] << " " << mu[idx_node_LLL] << " "
//                                          << mu[idx_node_RRR] << std::endl;
//
//                                std::cout << "dudy -  wrong!!! ei " << ei << " vi " << vi  << std::endl;
//                                exit(1);
//                            }
//
//
//                            if ((ei == d_num_eqn - 1) &&(fabs(mu[idx_node_L]*dudy[idx_node_L] + mu[idx_node_R]*dudy[idx_node_R]) > 1e-8 ||
//                                                         fabs(mu[idx_node_LL]*dudy[idx_node_LL] + mu[idx_node_RR]*dudy[idx_node_RR]) > 1e-8 ||
//                                                         fabs(mu[idx_node_LLL]*dudy[idx_node_LLL] + mu[idx_node_RRR]*dudy[idx_node_RRR]) > 1e-8)) {
//                                std::cout << "num_ghosts_cell_status[0] " << num_ghosts_cell_status[0] << std::endl;
//                                std::cout << dudy[idx_node_L] << " " << dudy[idx_node_R] << " " << dudy[idx_node_LL]
//                                          << " " << dudy[idx_node_RR] << " " << dudy[idx_node_LLL] << " "
//                                          << dudy[idx_node_RRR] << std::endl;
//                                std::cout << "dudy +  wrong!!!" << std::endl;
//                                exit(1);
//                            }
//                        }


                    }
                }
            }
        }
        
        var_data_x.clear();
        var_data_y.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        
        derivative_x.clear();
        derivative_y.clear();
        
        /*
         * (2) Compute the flux in the y-direction.
         */
        
        // Get the variables for the derivatives in the diffusive flux.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);


        /*
         * mirror ghost cell data for computing flux in the y-direction.
         */
        mirrorGhostCell(velocity, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        
        // Get the diffusivities in the diffusive flux.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION, true);
        
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);

        /*
        * mirror ghost cell data for computing the direvative in the x-direction.
        */

        mirrorGhostCell(velocity, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        
        /*
         * Compute the derivatives in x-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInX(
            patch,
            derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);


        /*
        * mirror ghost cell data for computing the direvative in the x-direction.
        */

        mirrorGhostCell(velocity, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInY(
            patch,
            derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);

        /*
         * Reconstruct the flux in y-direction.
         */

        {
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_y;
            du_y.reserve(5); //du/dx, dv/dx, du/dy, dv/dy dT/dy
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
                du_y.push_back(derivative_x[d_num_eqn - 1][vi]);
            for(int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
                du_y.push_back(derivative_y[d_num_eqn - 1][vi]);
            mirrorGhostCellDerivative(du_y, cell_status, DIRECTION::Y_DIRECTION);
        }

        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            double* F_face_y = diffusive_flux->getPointer(1, ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivative_x[ei][vi]->getPointer(0);
                
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
                
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        const int idx_node_BBB = (i + num_diff_ghosts_0) +
                            (j - 3 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_BB = (i + num_diff_ghosts_0) +
                            (j - 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_B = (i + num_diff_ghosts_0) +
                            (j - 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_T = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_TT = (i + num_diff_ghosts_0) +
                            (j + 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_TTT = (i + num_diff_ghosts_0) +
                            (j + 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0;

                        F_face_y[idx_face_y] += dt*(
                            double(37)/double(60)*(mu[idx_node_B]*dudx[idx_node_B] +  mu[idx_node_T]*dudx[idx_node_T]) +
                            double(-2)/double(15)*(mu[idx_node_BB]*dudx[idx_node_BB] + mu[idx_node_TT]*dudx[idx_node_TT]) +
                            double(1)/double(60)*(mu[idx_node_BBB]*dudx[idx_node_BBB] + mu[idx_node_TTT]*dudx[idx_node_TTT]));


//                        if((cell_status_data[idx_node_B] < 0.5 && cell_status_data[idx_node_T] > 0.5)||
//                           (cell_status_data[idx_node_B] > 0.5 && cell_status_data[idx_node_T] < 0.5))
//                        {
//
//
//                            if((ei <= d_num_eqn - 2) &&
//                                   (fabs(mu[idx_node_B]*dudx[idx_node_B] - mu[idx_node_T]*dudx[idx_node_T]) > 1e-8 ||
//                                    fabs(mu[idx_node_BB]*dudx[idx_node_BB] - mu[idx_node_TT]*dudx[idx_node_TT]) > 1e-8 ||
//                                    fabs(mu[idx_node_BBB]*dudx[idx_node_BBB] - mu[idx_node_TTT]*dudx[idx_node_TTT]) > 1e-8)) {
//                                    std::cout << " i "  << i << " j " << j << std::endl;
//                                    std::cout << mu[idx_node_B] << " " << mu[idx_node_T] << " " << mu[idx_node_BB] << " "
//                                              << mu[idx_node_TT] << " " << mu[idx_node_BBB] << " " << mu[idx_node_TTT]
//                                              << std::endl;
//                                    std::cout << dudx[idx_node_B] << " " << dudx[idx_node_T] << " " << dudx[idx_node_BB] << " "
//                                          << dudx[idx_node_TT] << " " << dudx[idx_node_BBB] << " " << dudx[idx_node_TTT]
//                                          << std::endl;
//                                    std::cout << " ei " << ei << " vi " << vi << " dudx2 mu wrong!!!" << std::endl;
//                                    exit(1);
//                                }
//                            if((ei == d_num_eqn - 1) &&
//                                      (fabs(mu[idx_node_B]*dudx[idx_node_B] + mu[idx_node_T]*dudx[idx_node_T]) > 1e-8 ||
//                                       fabs(mu[idx_node_BB]*dudx[idx_node_BB] + mu[idx_node_TT]*dudx[idx_node_TT]) > 1e-8 ||
//                                       fabs(mu[idx_node_BBB]*dudx[idx_node_BBB] + mu[idx_node_TTT]*dudx[idx_node_TTT]) > 1e-8)) {
//                                std::cout << "num_ghosts_cell_status[0] " << num_ghosts_cell_status[0] << std::endl;
//                                std::cout << mu[idx_node_B] << " " << mu[idx_node_T] << " " << mu[idx_node_BB] << " "
//                                          << mu[idx_node_TT] << " " << mu[idx_node_BBB] << " " << mu[idx_node_TTT]
//                                          << std::endl;
//                                std::cout  << " ei " << ei << "dudx2 mu wrong!!!" << std::endl;
//                                exit(1);
//                            }
//                        }


                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivative_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivative_y[ei][vi]->getPointer(0);
                
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
                
                for (int j = 0; j < interior_dim_1 + 1; j++)
                {
#ifdef HAMERS_ENABLE_SIMD
                    #pragma omp simd
#endif
                    for (int i = 0; i < interior_dim_0; i++)
                    {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                            j*interior_dim_0;
                        
                        const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                            (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity;
                        
                        const int idx_node_BBB = (i + num_diff_ghosts_0) +
                            (j - 3 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_BB = (i + num_diff_ghosts_0) +
                            (j - 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_B = (i + num_diff_ghosts_0) +
                            (j - 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_T = (i + num_diff_ghosts_0) +
                            (j + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_TT = (i + num_diff_ghosts_0) +
                            (j + 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        const int idx_node_TTT = (i + num_diff_ghosts_0) +
                            (j + 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0;
                        
                        F_face_y[idx_face_y] += dt*(
                            double(37)/double(60)*(mu[idx_node_B]*dudy[idx_node_B] + mu[idx_node_T]*dudy[idx_node_T]) +
                            double(-2)/double(15)*(mu[idx_node_BB]*dudy[idx_node_BB] + mu[idx_node_TT]*dudy[idx_node_TT]) +
                            double(1)/double(60)*(mu[idx_node_BBB]*dudy[idx_node_BBB] + mu[idx_node_TTT]*dudy[idx_node_TTT]));


                    }
                }
            }
        }
        
        var_data_x.clear();
        var_data_y.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        
        derivative_x.clear();
        derivative_y.clear();
        
//        /*
//         * Unregister the patch and data of all registered derived cell variables in the flow model.
//         */
//
//        d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(2))
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
         * Delcare containers for computing fluxes in different directions.
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
        
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > derivative_x;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > derivative_y;
        std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > > derivative_z;
        
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_x_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_y_computed;
        std::map<double*, boost::shared_ptr<pdat::CellData<double> > > derivative_z_computed;
        
        /*
         * (1) Compute the flux in the x-direction.
         */
        
        // Get the variables for the derivatives in the diffusive flux.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::X_DIRECTION,
            DIRECTION::Z_DIRECTION);

        /*
         * mirror ghost cell data for computing flux in the x-direction.
         */
        mirrorGhostCell(velocity, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        
        // Get the diffusivities in the diffusive flux.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::X_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::X_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::X_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */

        mirrorGhostCell(velocity, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        
        /*
         * Compute the derivatives in x-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInX(
            patch,
            derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */

        mirrorGhostCell(velocity, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInY(
            patch,
            derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */

        mirrorGhostCell(velocity, cell_status, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);

        /*
         * Compute the derivatives in z-direction for diffusive flux in x-direction.
         */
        
        computeFirstDerivativesInZ(
            patch,
            derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z);
        
        /*
         * Reconstruct the flux in x-direction.
         */

        {
            /*
             * Mirror the derivatives in x-direction to compute flux in the x direction.
             */
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_x;
            du_x.reserve(8); //du/dx, dv/dx, dw/dx, dT/dx, du/dy, dv/dy, du/dz, dw/dz
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
                du_x.push_back(derivative_x[d_num_eqn - 1][vi]);
            for(int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
                du_x.push_back(derivative_y[d_num_eqn - 1][vi]);
            for(int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
                du_x.push_back(derivative_z[d_num_eqn - 1][vi]);
            mirrorGhostCellDerivative(du_x, cell_status, DIRECTION::X_DIRECTION);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            double* F_face_x = diffusive_flux->getPointer(0, ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivative_x[ei][vi]->getPointer(0);
                
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
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_node_LLL = (i - 3 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_LL = (i - 2 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_L = (i - 1 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_R = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_RR = (i + 1 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_RRR = (i + 2 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            F_face_x[idx_face_x] += dt*(
                                double(37)/double(60)*(mu[idx_node_L]*dudx[idx_node_L] + mu[idx_node_R]*dudx[idx_node_R]) +
                                double(-2)/double(15)*(mu[idx_node_LL]*dudx[idx_node_LL] + mu[idx_node_RR]*dudx[idx_node_RR]) +
                                double(1)/double(60)*(mu[idx_node_LLL]*dudx[idx_node_LLL] + mu[idx_node_RRR]*dudx[idx_node_RRR]));


                            if ((cell_status_data[idx_node_L] < 0.5 && cell_status_data[idx_node_R] > 0.5) ||
                                (cell_status_data[idx_node_L] > 0.5 && cell_status_data[idx_node_R] < 0.5)) {
                                if ((ei <= d_num_eqn - 2) &&(fabs(mu[idx_node_L]*dudx[idx_node_L] - mu[idx_node_R]*dudx[idx_node_R]) > 1e-8 ||
                                                             fabs(mu[idx_node_LL]*dudx[idx_node_LL] - mu[idx_node_RR]*dudx[idx_node_RR]) > 1e-8 ||
                                                             fabs(mu[idx_node_LLL]*dudx[idx_node_LLL] - mu[idx_node_LLL]*dudx[idx_node_RRR]) > 1e-8)) {
                                    std::cout << "num_ghosts_cell_status[0] " << num_ghosts_cell_status[0] << std::endl;
                                    std::cout << dudx[idx_node_L] << " " << dudx[idx_node_R] << " " << dudx[idx_node_LL]
                                              << " " << dudx[idx_node_RR] << " " << dudx[idx_node_LLL] << " "
                                              << dudx[idx_node_RRR] << std::endl;

                                    std::cout << mu[idx_node_L] << " " << mu[idx_node_R] << " " << mu[idx_node_LL]
                                              << " " << mu[idx_node_RR] << " " << mu[idx_node_LLL] << " "
                                              << mu[idx_node_RRR] << std::endl;


                                    std::cout << "dudx wrong!!!" << " ei " << ei << " vi " << vi  << std::endl;
                                    exit(1);
                                }


                                if ((ei == d_num_eqn - 1) &&(fabs(mu[idx_node_L]*dudx[idx_node_L] + mu[idx_node_R]*dudx[idx_node_R]) > 1e-8 ||
                                                             fabs(mu[idx_node_LL]*dudx[idx_node_LL] + mu[idx_node_RR]*dudx[idx_node_RR]) > 1e-8 ||
                                                             fabs(mu[idx_node_LLL]*dudx[idx_node_LLL] + mu[idx_node_RRR]*dudx[idx_node_RRR]) > 1e-8)) {
                                    std::cout << "num_ghosts_cell_status[0] " << num_ghosts_cell_status[0] << std::endl;
                                    std::cout << "dudx " << dudx[idx_node_L] << " " << dudx[idx_node_R] << " " << dudx[idx_node_LL]
                                              << " " << dudx[idx_node_RR] << " " << dudx[idx_node_LLL] << " "
                                              << dudx[idx_node_RRR] << std::endl;

                                    std::cout << "mu " << mu[idx_node_L] << " " << mu[idx_node_R] << " " << mu[idx_node_LL]
                                              << " " << mu[idx_node_RR] << " " << mu[idx_node_LLL] << " "
                                              << mu[idx_node_RRR] << std::endl;

                                    std::cout << mu[idx_node_L]*dudx[idx_node_L] + mu[idx_node_R]*dudx[idx_node_R]  << " "
                                              << mu[idx_node_LL]*dudx[idx_node_LL] + mu[idx_node_RR]*dudx[idx_node_RR] << " "
                                              << mu[idx_node_LLL]*dudx[idx_node_LLL] + mu[idx_node_RRR]*dudx[idx_node_RRR] << std::endl;

                                    std::cout << "dudx wrong!!!" << " ei " << ei << " vi " << vi  << std::endl;
                                    exit(1);
                                }
                            }


                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivative_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivative_y[ei][vi]->getPointer(0);
                
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
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_node_LLL = (i - 3 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_LL = (i - 2 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_L = (i - 1 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_R = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_RR = (i + 1 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_RRR = (i + 2 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            F_face_x[idx_face_x] += dt*(
                                double(37)/double(60)*(mu[idx_node_L]*dudy[idx_node_L] + mu[idx_node_R]*dudy[idx_node_R]) +
                                double(-2)/double(15)*(mu[idx_node_LL]*dudy[idx_node_LL] + mu[idx_node_RR]*dudy[idx_node_RR]) +
                                double(1)/double(60)*(mu[idx_node_LLL]*dudy[idx_node_LLL] + mu[idx_node_RRR]*dudy[idx_node_RRR]));
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivative_z[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudz = derivative_z[ei][vi]->getPointer(0);
                
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
                        for (int i = 0; i < interior_dim_0 + 1; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_x = i +
                                j*(interior_dim_0 + 1) +
                                k*(interior_dim_0 + 1)*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_node_LLL = (i - 3 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_LL = (i - 2 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_L = (i - 1 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_R = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_RR = (i + 1 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_RRR = (i + 2 + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            F_face_x[idx_face_x] += dt*(
                                double(37)/double(60)*(mu[idx_node_L]*dudz[idx_node_L] + mu[idx_node_R]*dudz[idx_node_R]) +
                                double(-2)/double(15)*(mu[idx_node_LL]*dudz[idx_node_LL] + mu[idx_node_RR]*dudz[idx_node_RR]) +
                                double(1)/double(60)*(mu[idx_node_LLL]*dudz[idx_node_LLL] + mu[idx_node_RRR]*dudz[idx_node_RRR]));
                        }
                    }
                }
            }
        }
        
        var_data_x.clear();
        var_data_y.clear();
        var_data_z.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        diffusivities_data_z.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        var_component_idx_z.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        diffusivities_component_idx_z.clear();
        
        derivative_x.clear();
        derivative_y.clear();
        derivative_z.clear();
        
        /*
         * (2) Compute the flux in the y-direction.
         */
        
        // Get the variables for the derivatives in the diffusive flux.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Z_DIRECTION);

        /*
         * mirror ghost cell data for computing flux in the y-direction. need to recompute diffusivities
         */
        mirrorGhostCell(velocity, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        
        // Get the diffusivities in the diffusive flux.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Y_DIRECTION,
            DIRECTION::X_DIRECTION, true);
        
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::Y_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);

        /*
        * mirror ghost cell data for computing the direvative in the x-direction.
        */

        mirrorGhostCell(velocity, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        
        /*
         * Compute the derivatives in x-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInX(
            patch,
            derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);

        /*
        * mirror ghost cell data for computing the direvative in the x-direction.
        */

        mirrorGhostCell(velocity, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInY(
            patch,
            derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);

        /*
       * mirror ghost cell data for computing the direvative in the x-direction.
       */

        mirrorGhostCell(velocity, cell_status, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);
        
        /*
         * Compute the derivatives in z-direction for diffusive flux in y-direction.
         */
        
        computeFirstDerivativesInZ(
            patch,
            derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z);
        
        /*
         * Reconstruct the flux in y-direction.
         */

        {
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_y;
            du_y.reserve(8); //du/dx, dv/dx, du/dy, dv/dy, dw/dy, dT/dy, dv/dz, dw/dz
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
                du_y.push_back(derivative_x[d_num_eqn - 1][vi]);
            for(int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
                du_y.push_back(derivative_y[d_num_eqn - 1][vi]);
            for(int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
                du_y.push_back(derivative_z[d_num_eqn - 1][vi]);
            mirrorGhostCellDerivative(du_y, cell_status, DIRECTION::Y_DIRECTION);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            double* F_face_y = diffusive_flux->getPointer(1, ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivative_x[ei][vi]->getPointer(0);
                
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
                    for (int j = 0; j < interior_dim_1 + 1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_node_BBB = (i + num_diff_ghosts_0) +
                                (j - 3 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_BB = (i + num_diff_ghosts_0) +
                                (j - 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_B = (i + num_diff_ghosts_0) +
                                (j - 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_T = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_TT = (i + num_diff_ghosts_0) +
                                (j + 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_TTT = (i + num_diff_ghosts_0) +
                                (j + 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            F_face_y[idx_face_y] += dt*(
                                double(37)/double(60)*(mu[idx_node_B]*dudx[idx_node_B] + mu[idx_node_T]*dudx[idx_node_T]) +
                                double(-2)/double(15)*(mu[idx_node_BB]*dudx[idx_node_BB] + mu[idx_node_TT]*dudx[idx_node_TT]) +
                                double(1)/double(60)*(mu[idx_node_BBB]*dudx[idx_node_BBB] + mu[idx_node_TTT]*dudx[idx_node_TTT]));
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivative_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivative_y[ei][vi]->getPointer(0);
                
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
                    for (int j = 0; j < interior_dim_1 + 1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_node_BBB = (i + num_diff_ghosts_0) +
                                (j - 3 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_BB = (i + num_diff_ghosts_0) +
                                (j - 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_B = (i + num_diff_ghosts_0) +
                                (j - 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_T = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_TT = (i + num_diff_ghosts_0) +
                                (j + 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_TTT = (i + num_diff_ghosts_0) +
                                (j + 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            F_face_y[idx_face_y] += dt*(
                                double(37)/double(60)*(mu[idx_node_B]*dudy[idx_node_B] + mu[idx_node_T]*dudy[idx_node_T]) +
                                double(-2)/double(15)*(mu[idx_node_BB]*dudy[idx_node_BB] + mu[idx_node_TT]*dudy[idx_node_TT]) +
                                double(1)/double(60)*(mu[idx_node_BBB]*dudy[idx_node_BBB] + mu[idx_node_TTT]*dudy[idx_node_TTT]));



                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivative_z[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudz = derivative_z[ei][vi]->getPointer(0);
                
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
                    for (int j = 0; j < interior_dim_1 + 1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_y = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*(interior_dim_1 + 1);
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_node_BBB = (i + num_diff_ghosts_0) +
                                (j - 3 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_BB = (i + num_diff_ghosts_0) +
                                (j - 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_B = (i + num_diff_ghosts_0) +
                                (j - 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_T = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_TT = (i + num_diff_ghosts_0) +
                                (j + 1 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_TTT = (i + num_diff_ghosts_0) +
                                (j + 2 + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            F_face_y[idx_face_y] += dt*(
                                double(37)/double(60)*(mu[idx_node_B]*dudz[idx_node_B] + mu[idx_node_T]*dudz[idx_node_T]) +
                                double(-2)/double(15)*(mu[idx_node_BB]*dudz[idx_node_BB] + mu[idx_node_TT]*dudz[idx_node_TT]) +
                                double(1)/double(60)*(mu[idx_node_BBB]*dudz[idx_node_BBB] + mu[idx_node_TTT]*dudz[idx_node_TTT]));
                        }
                    }
                }
            }
        }
        
        var_data_x.clear();
        var_data_y.clear();
        var_data_z.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        diffusivities_data_z.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        var_component_idx_z.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        diffusivities_component_idx_z.clear();
        
        derivative_x.clear();
        derivative_y.clear();
        derivative_z.clear();
        
        /*
         * (3) Compute the flux in the z-direction.
         */
        
        // Get the variables for the derivatives in the diffusive flux.
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_x,
            var_component_idx_x,
            DIRECTION::Z_DIRECTION,
            DIRECTION::X_DIRECTION);
        
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_y,
            var_component_idx_y,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        d_flow_model->getDiffusiveFluxVariablesForDerivative(
            var_data_z,
            var_component_idx_z,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Z_DIRECTION);

        /*
         * mirror ghost cell data for computing flux in the z-direction.
         */
        mirrorGhostCell(velocity, cell_status, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);
        
        // Get the diffusivities in the diffusive flux.
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_x,
            diffusivities_component_idx_x,
            DIRECTION::Z_DIRECTION,
            DIRECTION::X_DIRECTION, true);
        
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_y,
            diffusivities_component_idx_y,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Y_DIRECTION);
        
        d_flow_model->getDiffusiveFluxDiffusivities(
            diffusivities_data_z,
            diffusivities_component_idx_z,
            DIRECTION::Z_DIRECTION,
            DIRECTION::Z_DIRECTION);
        
        TBOX_ASSERT(static_cast<int>(var_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(var_component_idx_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_data_z.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_x.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_y.size()) == d_num_eqn);
        TBOX_ASSERT(static_cast<int>(diffusivities_component_idx_z.size()) == d_num_eqn);

        /*
        * mirror ghost cell data for computing the direvative in the x-direction.
        */

        mirrorGhostCell(velocity, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::X_DIRECTION, WALL_NO_SLIP);
        
        /*
         * Compute the derivatives in x-direction for diffusive flux in z-direction.
         */
        
        computeFirstDerivativesInX(
            patch,
            derivative_x,
            derivative_x_computed,
            var_data_x,
            var_component_idx_x);
        
        /*
         * Compute the derivatives in y-direction for diffusive flux in z-direction.
         */


        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */

        mirrorGhostCell(velocity, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::Y_DIRECTION, WALL_NO_SLIP);


        computeFirstDerivativesInY(
            patch,
            derivative_y,
            derivative_y_computed,
            var_data_y,
            var_component_idx_y);

        /*
         * mirror ghost cell data for computing the direvative in the x-direction.
         */

        mirrorGhostCell(velocity, cell_status, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);
        mirrorGhostCell(temperature, cell_status, DIRECTION::Z_DIRECTION, WALL_NO_SLIP);

        /*
         * Compute the derivatives in z-direction for diffusive flux in z-direction.
         */
        
        computeFirstDerivativesInZ(
            patch,
            derivative_z,
            derivative_z_computed,
            var_data_z,
            var_component_idx_z);
        
        /*
         * Reconstruct the flux in z-direction.
         */

        /*
         * Reconstruct the flux in x-direction.
         */

        {
            /*
             * Mirror the derivatives in z-direction to compute flux in the z direction.
             */
            std::vector<boost::shared_ptr<pdat::CellData<double> > > du_z;
            du_z.reserve(8); //du/dx, dw/dx, dv/dy, dw/dy, du/dz, dv/dz, dw/dz, dT/dz
            int ei = d_num_eqn - 1;
            for(int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
                du_z.push_back(derivative_x[d_num_eqn - 1][vi]);
            for(int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
                du_z.push_back(derivative_y[d_num_eqn - 1][vi]);
            for(int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
                du_z.push_back(derivative_z[d_num_eqn - 1][vi]);
            mirrorGhostCellDerivative(du_z, cell_status, DIRECTION::Z_DIRECTION);
        }


        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            TBOX_ASSERT(static_cast<int>(derivative_x[ei].size()) ==
                        static_cast<int>(diffusivities_data_x[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_x[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_x[ei].size()));
            
            double* F_face_z = diffusive_flux->getPointer(2, ei);
            
            for (int vi = 0; vi < static_cast<int>(var_data_x[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_x[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_x[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudx = derivative_x[ei][vi]->getPointer(0);
                
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
                
                for (int k = 0; k < interior_dim_2 + 1; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_node_BBB = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k - 3 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_BB = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k - 2 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_B = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k - 1 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_F = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_FF = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + 1 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_FFF = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + 2 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            F_face_z[idx_face_z] += dt*(
                                double(37)/double(60)*(mu[idx_node_B]*dudx[idx_node_B] + mu[idx_node_F]*dudx[idx_node_F]) +
                                double(-2)/double(15)*(mu[idx_node_BB]*dudx[idx_node_BB] + mu[idx_node_FF]*dudx[idx_node_FF]) +
                                double(1)/double(60)*(mu[idx_node_BBB]*dudx[idx_node_BBB] + mu[idx_node_FFF]*dudx[idx_node_FFF]));
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivative_y[ei].size()) ==
                        static_cast<int>(diffusivities_data_y[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_y[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_y[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_y[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_y[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_y[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudy = derivative_y[ei][vi]->getPointer(0);
                
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
                
                for (int k = 0; k < interior_dim_2 + 1; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_node_BBB = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k - 3 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_BB = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k - 2 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_B = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k - 1 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_F = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_FF = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + 1 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_FFF = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + 2 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            F_face_z[idx_face_z] += dt*(
                                double(37)/double(60)*(mu[idx_node_B]*dudy[idx_node_B] + mu[idx_node_F]*dudy[idx_node_F]) +
                                double(-2)/double(15)*(mu[idx_node_BB]*dudy[idx_node_BB] + mu[idx_node_FF]*dudy[idx_node_FF]) +
                                double(1)/double(60)*(mu[idx_node_BBB]*dudy[idx_node_BBB] + mu[idx_node_FFF]*dudy[idx_node_FFF]));
                        }
                    }
                }
            }
            
            TBOX_ASSERT(static_cast<int>(derivative_z[ei].size()) ==
                        static_cast<int>(diffusivities_data_z[ei].size()));
            
            TBOX_ASSERT(static_cast<int>(diffusivities_data_z[ei].size()) ==
                        static_cast<int>(diffusivities_component_idx_z[ei].size()));
            
            for (int vi = 0; vi < static_cast<int>(var_data_z[ei].size()); vi++)
            {
                // Get the index of variable for derivative.
                const int mu_idx = diffusivities_component_idx_z[ei][vi];
                
                // Get the pointer to diffusivity.
                double* mu = diffusivities_data_z[ei][vi]->getPointer(mu_idx);
                
                // Get the pointer to derivative.
                double* dudz = derivative_z[ei][vi]->getPointer(0);
                
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
                
                for (int k = 0; k < interior_dim_2 + 1; k++)
                {
                    for (int j = 0; j < interior_dim_1; j++)
                    {
#ifdef HAMERS_ENABLE_SIMD
                        #pragma omp simd
#endif
                        for (int i = 0; i < interior_dim_0; i++)
                        {
                            // Compute the linear indices.
                            const int idx_face_z = i +
                                j*interior_dim_0 +
                                k*interior_dim_0*interior_dim_1;
                            
                            const int idx_diffusivity = (i + num_subghosts_0_diffusivity) +
                                (j + num_subghosts_1_diffusivity)*subghostcell_dim_0_diffusivity +
                                (k + num_subghosts_2_diffusivity)*subghostcell_dim_0_diffusivity*
                                    subghostcell_dim_1_diffusivity;
                            
                            const int idx_node_BBB = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k - 3 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_BB = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k - 2 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_B = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k - 1 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_F = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_FF = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + 1 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            const int idx_node_FFF = (i + num_diff_ghosts_0) +
                                (j + num_diff_ghosts_1)*diff_ghostcell_dim_0 +
                                (k + 2 + num_diff_ghosts_2)*diff_ghostcell_dim_0*
                                    diff_ghostcell_dim_1;
                            
                            F_face_z[idx_face_z] += dt*(
                                double(37)/double(60)*(mu[idx_node_B]*dudz[idx_node_B] + mu[idx_node_F]*dudz[idx_node_F]) +
                                double(-2)/double(15)*(mu[idx_node_BB]*dudz[idx_node_BB] + mu[idx_node_FF]*dudz[idx_node_FF]) +
                                double(1)/double(60)*(mu[idx_node_BBB]*dudz[idx_node_BBB] + mu[idx_node_FFF]*dudz[idx_node_FFF]));


                            if((cell_status_data[idx_node_B] < 0.5 && cell_status_data[idx_node_F] > 0.5)||
                               (cell_status_data[idx_node_B] > 0.5 && cell_status_data[idx_node_F] < 0.5))
                            {


                                if((ei <= d_num_eqn - 2) &&
                                   (fabs(mu[idx_node_B]*dudz[idx_node_B] - mu[idx_node_F]*dudz[idx_node_F]) > 1e-8 ||
                                    fabs(mu[idx_node_BB]*dudz[idx_node_BB] - mu[idx_node_FF]*dudz[idx_node_FF]) > 1e-8 ||
                                    fabs(mu[idx_node_BBB]*dudz[idx_node_BBB] - mu[idx_node_FFF]*dudz[idx_node_FFF]) > 1e-8)) {
                                    std::cout << "num_ghosts_cell_status[0] " << num_ghosts_cell_status[0] << std::endl;
                                    std::cout << mu[idx_node_B] << " " << mu[idx_node_F] << " " << mu[idx_node_BB] << " "
                                              << mu[idx_node_FF] << " " << mu[idx_node_BBB] << " " << mu[idx_node_FFF]
                                              << std::endl;

                                    std::cout << dudz[idx_node_B] << " " <<  dudz[idx_node_F] << " " <<  dudz[idx_node_BB] << " "
                                              <<  dudz[idx_node_FF] << " " <<  dudz[idx_node_BBB] << " " <<  dudz[idx_node_FFF]
                                              << std::endl;

                                    std::cout << "dudz2 mu wrong!!!" << " ei " << ei << " vi " << vi << std::endl;
                                    exit(1);
                                }
                                if((ei == d_num_eqn - 1) &&
                                   (fabs(mu[idx_node_B]*dudz[idx_node_B] + mu[idx_node_F]*dudz[idx_node_F]) > 1e-8 ||
                                    fabs(mu[idx_node_BB]*dudz[idx_node_BB] + mu[idx_node_FF]*dudz[idx_node_FF]) > 1e-8 ||
                                    fabs(mu[idx_node_BBB]*dudz[idx_node_BBB] + mu[idx_node_FFF]*dudz[idx_node_FFF]) > 1e-8)) {
                                    std::cout << "num_ghosts_cell_status[0] " << num_ghosts_cell_status[0] << std::endl;
                                    std::cout << mu[idx_node_B] << " " << mu[idx_node_F] << " " << mu[idx_node_BB] << " "
                                              << mu[idx_node_FF] << " " << mu[idx_node_BBB] << " " << mu[idx_node_FFF]
                                              << std::endl;

                                    std::cout << "vx " << (velocity->getPointer(0))[idx_node_B] << " " << (velocity->getPointer(0))[idx_node_F] << " " << (velocity->getPointer(0))[idx_node_BB] << " "
                                              << (velocity->getPointer(0))[idx_node_FF] << " " << (velocity->getPointer(0))[idx_node_BBB] << " " << (velocity->getPointer(0))[idx_node_FFF]
                                              << std::endl;

                                    std::cout << "vy " << (velocity->getPointer(1))[idx_node_B] << " " << (velocity->getPointer(1))[idx_node_F] << " " << (velocity->getPointer(1))[idx_node_BB] << " "
                                              << (velocity->getPointer(1))[idx_node_FF] << " " << (velocity->getPointer(1))[idx_node_BBB] << " " << (velocity->getPointer(1))[idx_node_FFF]
                                              << std::endl;

                                    std::cout << "vz " << (velocity->getPointer(2))[idx_node_B] << " " << (velocity->getPointer(2))[idx_node_F] << " " << (velocity->getPointer(2))[idx_node_BB] << " "
                                              << (velocity->getPointer(2))[idx_node_FF] << " " << (velocity->getPointer(2))[idx_node_BBB] << " " << (velocity->getPointer(2))[idx_node_FFF]
                                              << std::endl;

                                    std::cout << dudz[idx_node_B] << " " <<  dudz[idx_node_F] << " " <<  dudz[idx_node_BB] << " "
                                              <<  dudz[idx_node_FF] << " " <<  dudz[idx_node_BBB] << " " <<  dudz[idx_node_FFF]
                                              << std::endl;
                                    std::cout << "dudz2 mu wrong!!!" << " ei " << ei << " vi " << vi  << " in " << static_cast<int>(var_data_z[ei].size()) <<  " mu_idx " << mu_idx << std::endl;
                                    std::cout << " i " << i <<  " j "  << j  << " k " << k <<  std::endl;
                                    exit(1);
                                }
                            }

                        }
                    }
                }
            }
        }
        
        var_data_x.clear();
        var_data_y.clear();
        var_data_z.clear();
        
        diffusivities_data_x.clear();
        diffusivities_data_y.clear();
        diffusivities_data_z.clear();
        
        var_component_idx_x.clear();
        var_component_idx_y.clear();
        var_component_idx_z.clear();
        
        diffusivities_component_idx_x.clear();
        diffusivities_component_idx_y.clear();
        diffusivities_component_idx_z.clear();
        
        derivative_x.clear();
        derivative_y.clear();
        derivative_z.clear();

        
    } // if (d_dim == tbox::Dimension(3))

    /*
     * Unregister the patch and data of all registered derived cell variables in the flow model.
     */

     d_flow_model->unregisterPatch();
}


/*
 * Compute the first derivatives in the x-direction.
 */
void
DiffusiveFluxReconstructorSixthOrder::computeFirstDerivativesInX(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_x,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_x_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_x,
    const std::vector<std::vector<int> >& data_component_idx_x)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
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
                    for (int i = -3; i < interior_dim_0 + 3; i++)
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
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
                        for (int i = -3; i < interior_dim_0 + 3; i++)
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
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
                            for (int i = -3; i < interior_dim_0 + 3; i++)
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
DiffusiveFluxReconstructorSixthOrder::computeFirstDerivativesInY(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_y,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_y_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_y,
    const std::vector<std::vector<int> >& data_component_idx_y)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
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
            << ": DiffusiveFluxReconstructorSixthOrder::"
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
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
                    
                    for (int j = -3; j < interior_dim_1 + 3; j++)
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
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
DiffusiveFluxReconstructorSixthOrder::computeFirstDerivativesInZ(
    hier::Patch& patch,
    std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& derivative_z,
    std::map<double*, boost::shared_ptr<pdat::CellData<double> > >& derivative_z_computed,
    const std::vector<std::vector<boost::shared_ptr<pdat::CellData<double> > > >& data_z,
    const std::vector<std::vector<int> >& data_component_idx_z)
{
#ifdef HAMERS_DEBUG_CHECK_DEV_ASSERTIONS
    for (int ei = 0; ei < d_num_eqn; ei++)
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
            << ": DiffusiveFluxReconstructorSixthOrder::"
            << "computeFirstDerivativesInZ()\n"
            << "There isn't z-direction for 1D problem."
            << std::endl);
    }
    else if (d_dim == tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": DiffusiveFluxReconstructorSixthOrder::"
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
        
        for (int ei = 0; ei < d_num_eqn; ei++)
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
                    
                    for (int k = -3; k < interior_dim_2 + 3; k++)
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
