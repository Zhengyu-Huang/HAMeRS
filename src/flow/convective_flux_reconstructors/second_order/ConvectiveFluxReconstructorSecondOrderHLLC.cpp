#include "flow/convective_flux_reconstructors/second_order/ConvectiveFluxReconstructorSecondOrderHLLC.hpp"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "util/basic_geometry/WallTreatment.hpp"

#include <limits>

ConvectiveFluxReconstructorSecondOrderHLLC::ConvectiveFluxReconstructorSecondOrderHLLC(
    const std::string& object_name,
    const tbox::Dimension& dim,
    const boost::shared_ptr<geom::CartesianGridGeometry>& grid_geometry,
    const int& num_eqn,
    const boost::shared_ptr<FlowModel>& flow_model,
    const boost::shared_ptr<tbox::Database>& convective_flux_reconstructor_db):
        ConvectiveFluxReconstructor(
            object_name,
            dim,
            grid_geometry,
            num_eqn,
            flow_model,
            convective_flux_reconstructor_db)
{
    d_num_conv_ghosts = hier::IntVector::getOne(d_dim)*2;
    
    d_eqn_form = d_flow_model->getEquationsForm();
    d_has_advective_eqn_form = false;
    for (int ei = 0; ei < d_num_eqn; ei++)
    {
        if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
        {
            d_has_advective_eqn_form = true;
        }
    }
}

/*
 * Print all characteristics of the convective flux reconstruction class.
 */
void
ConvectiveFluxReconstructorSecondOrderHLLC::printClassData(
    std::ostream& os) const
{
    os << "\nPrint ConvectiveFluxReconstructorSecondOrderHLLC object..."
       << std::endl;
    
    os << std::endl;
    
    os << "ConvectiveFluxReconstructorSecondOrderHLLC: this = "
       << (ConvectiveFluxReconstructorSecondOrderHLLC *)this
       << std::endl;
    os << "d_object_name = "
       << d_object_name
       << std::endl;
}


/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorSecondOrderHLLC::putToRestart(
   const boost::shared_ptr<tbox::Database>& restart_db) const
{
    restart_db->putString("d_shock_capturing_scheme", "FIRST_ORDER_HLLC");
}


/*
 * Compute the convective flux and source due to splitting of convective term on a patch.
 */
void
ConvectiveFluxReconstructorSecondOrderHLLC::computeConvectiveFluxAndSourceOnPatch(
    hier::Patch& patch,
    const boost::shared_ptr<pdat::SideVariable<double> >& variable_convective_flux,
    const boost::shared_ptr<pdat::CellVariable<double> >& variable_source,
    const boost::shared_ptr<hier::VariableContext>& data_context,
    const double time,
    const double dt,
    const int RK_step_number)
{
    NULL_USE(time);
    NULL_USE(RK_step_number);
    
    d_flow_model->setupRiemannSolver();
    d_riemann_solver = d_flow_model->getFlowModelRiemannSolver();
    
    // Get the dimensions of box that covers the interior of patch.
    hier::Box interior_box = patch.getBox();
    const hier::IntVector interior_dims = interior_box.numberCells();
    
    // Get the dimensions of box that covers interior of patch plus
    // convective ghost cells.
    hier::Box conv_ghost_box = interior_box;
    conv_ghost_box.grow(d_num_conv_ghosts);
    const hier::IntVector conv_ghostcell_dims = conv_ghost_box.numberCells();
    
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
        BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
            patch.getPatchGeometry()));
    
    const double* const dx = patch_geom->getDx();
    const double* const x_lo = patch_geom->getXLower();

    /*
     * Compute the cell status.
     */
    d_flow_model->registerPatchWithDataContext(patch, data_context);
    boost::shared_ptr<pdat::CellData<double> > cell_status
            = d_flow_model->getGlobalCellStatus();

    
    // Get the side data of convective flux.
    boost::shared_ptr<pdat::SideData<double> > convective_flux(
        BOOST_CAST<pdat::SideData<double>, hier::PatchData>(
            patch.getPatchData(variable_convective_flux, data_context)));
    
    // Get the cell data of source.
    boost::shared_ptr<pdat::CellData<double> > source(
        BOOST_CAST<pdat::CellData<double>, hier::PatchData>(
            patch.getPatchData(variable_source, data_context)));
    
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(convective_flux);
    TBOX_ASSERT(convective_flux->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
    
    TBOX_ASSERT(source);
    TBOX_ASSERT(source->getGhostCellWidth() == hier::IntVector::getZero(d_dim));
#endif
    
    // Allocate temporary patch data.
    boost::shared_ptr<pdat::SideData<double> > velocity_intercell;

    if (d_has_advective_eqn_form)
    {
        velocity_intercell.reset(new pdat::SideData<double>(
            interior_box, d_dim.getValue(), hier::IntVector::getZero(d_dim)));
    }
    
    if (d_dim == tbox::Dimension(1))
    {
        /*
         * Get the pointer to the convective flux side data.
         */
        
        std::vector<double*> F_face_x;
        F_face_x.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_face_x.push_back(convective_flux->getPointer(0, ei));
        }
        
        /*
         * Get the pointers to the conservative variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
            d_flow_model->getGlobalCellDataConservativeVariables();
        
        std::vector<hier::IntVector> num_subghosts_conservative_var;
        num_subghosts_conservative_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_conservative_var;
        subghostcell_dims_conservative_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
        Q.reserve(d_num_eqn);
        
        int count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++)
        {
            int depth = conservative_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the conservative variable vector is not in the system of equations, ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                Q.push_back(conservative_variables[vi]->getPointer(di));
                num_subghosts_conservative_var.push_back(conservative_variables[vi]->getGhostCellWidth());
                subghostcell_dims_conservative_var.push_back(conservative_variables[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        /*
         * Declare temporary data containers for computing the fluxes at cell edges.
         */
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > conservative_variables_minus;
        std::vector<boost::shared_ptr<pdat::SideData<double> > > conservative_variables_plus;
        
        conservative_variables_minus.reserve(d_num_eqn);
        conservative_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            conservative_variables_minus.push_back(boost::make_shared<pdat::SideData<double> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));
            
            conservative_variables_plus.push_back(boost::make_shared<pdat::SideData<double> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));
        }
        
        /*
         * Initialize temporary data containers for computing the flux in the x-direction.
         */
        
        std::vector<double*> Q_minus;
        std::vector<double*> Q_plus;
        Q_minus.resize(d_num_eqn);
        Q_plus.resize(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            Q_minus[ei] = conservative_variables_minus[ei]->getPointer(0);
            Q_plus[ei] = conservative_variables_plus[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int i = 0; i < interior_dims[0] + 1; i++)
            {
                // Compute the linear indices.
                const int idx_face_x = i;
                const int idx_L = i - 1 + num_subghosts_conservative_var[ei][0];
                const int idx_R = i + num_subghosts_conservative_var[ei][0];
                
                Q_minus[ei][idx_face_x] = Q[ei][idx_L];
                Q_plus[ei][idx_face_x] = Q[ei][idx_R];
            }
        }
        
        /*
         * Compute flux in the x-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            d_riemann_solver->computeConvectiveFluxAndVelocityFromConservativeVariables(
                convective_flux,
                velocity_intercell,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            d_riemann_solver->computeConvectiveFluxFromConservativeVariables(
                convective_flux,
                conservative_variables_minus,
                conservative_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int i = 0; i < interior_dims[0] + 1; i++)
            {
                // Compute the linear index.
                const int idx_face_x = i;
                
                F_face_x[ei][idx_face_x] *= dt;
            }
        }

        /*
         * Compute the source.
         */

        if (d_has_advective_eqn_form)
        {
            for (int ei = 0; ei < d_num_eqn; ei ++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    double* S = source->getPointer(ei);

                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute the linear indices.
                        const int idx_cell_wghost = i + num_subghosts_conservative_var[ei][0];
                        const int idx_cell_nghost = i;
                        const int idx_face_x_L = i;
                        const int idx_face_x_R = i + 1;

                        const double& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                        const double& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];

                        S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*(u_R - u_L)/dx[0];
                    }
                }
            }
        }
        
        /*
         * Unregister the patch in the flow model.
         */
        
        //d_flow_model->unregisterPatch();
        
    } // if (d_dim == tbox::Dimension(1))
    else if (d_dim == tbox::Dimension(2))
    {
        /*
         * Get the pointers to the convective flux side data.
         */
        
        std::vector<double*> F_face_x;
        std::vector<double*> F_face_y;
        F_face_x.reserve(d_num_eqn);
        F_face_y.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_face_x.push_back(convective_flux->getPointer(0, ei));
            F_face_y.push_back(convective_flux->getPointer(1, ei));
        }


        /*
         * Register primitive variables
         */
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRIMITIVE_VARIABLES", d_num_conv_ghosts));
        d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
        d_flow_model->computeGlobalDerivedCellData();

        /*
         * Build ghost cell map
         */
        std::vector<std::vector<std::array<int,3> > > ghost_cell_maps;
        ghost_cell_maps.resize(2);
        buildGhostCellMap2D(cell_status, ghost_cell_maps);
        
        /*
         * Get the pointers to the conservative/primitive variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
            d_flow_model->getGlobalCellDataConservativeVariables();

        std::vector<boost::shared_ptr<pdat::CellData<double> > > primitive_variables =
                d_flow_model->getGlobalCellDataPrimitiveVariables();
        
        std::vector<hier::IntVector> num_subghosts_primitive_var;
        num_subghosts_primitive_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_primitive_var;
        subghostcell_dims_primitive_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
        Q.reserve(d_num_eqn);

        std::vector<double*> V;
        V.reserve(d_num_eqn);
        
        int count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the primitive variable vector is not in the system of equations, ignore it.
                if (count_eqn >= d_num_eqn)
                    break;
                
                V.push_back(primitive_variables[vi]->getPointer(di));
                Q.push_back(conservative_variables[vi]->getPointer(di));
                num_subghosts_primitive_var.push_back(primitive_variables[vi]->getGhostCellWidth());
                subghostcell_dims_primitive_var.push_back(primitive_variables[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        /*
         * Declare temporary data containers for computing the fluxes at cell edges.
         */
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_minus;
        std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_plus;

        primitive_variables_minus.reserve(d_num_eqn);
        primitive_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            primitive_variables_minus.push_back(boost::make_shared<pdat::SideData<double> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));

            primitive_variables_plus.push_back(boost::make_shared<pdat::SideData<double> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));
        }
        
        std::vector<double*> V_minus;
        std::vector<double*> V_plus;
        V_minus.resize(d_num_eqn);
        V_plus.resize(d_num_eqn);

        boost::shared_ptr<pdat::SideData<int> > bounded_flag_minus;
        boost::shared_ptr<pdat::SideData<int> > bounded_flag_plus;

        bounded_flag_minus.reset(
                new pdat::SideData<int>(interior_box, 1, hier::IntVector::getZero(d_dim)));

        bounded_flag_plus.reset(
                new pdat::SideData<int>(interior_box, 1, hier::IntVector::getZero(d_dim)));

        int* flag_minus = nullptr;
        int* flag_plus = nullptr;
        
        /*
         * mirroring ghost cell data for computing the flux in the x-direction.
         */
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {

            mirrorGhostCell2D(primitive_variables[vi], ghost_cell_maps, DIRECTION::X_DIRECTION);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++) {
            for (int j = 0; j < interior_dims[1]; j++) {
                for (int i = 0; i < interior_dims[0] + 1; i++) {
                    // Compute the linear indices.
                    const int idx_face_x = i +
                                           j * (interior_dims[0] + 1);

                    const int idx_L = (i - 1 + num_subghosts_primitive_var[ei][0]) +
                                      (j + num_subghosts_primitive_var[ei][1]) *
                                      subghostcell_dims_primitive_var[ei][0];

                    const int idx_LL = (i - 2 + num_subghosts_primitive_var[ei][0]) +
                                       (j + num_subghosts_primitive_var[ei][1]) *
                                       subghostcell_dims_primitive_var[ei][0];

                    const int idx_R = (i + num_subghosts_primitive_var[ei][0]) +
                                      (j + num_subghosts_primitive_var[ei][1]) *
                                      subghostcell_dims_primitive_var[ei][0];

                    const int idx_RR = (i + 1 + num_subghosts_primitive_var[ei][0]) +
                                       (j + num_subghosts_primitive_var[ei][1]) *
                                       subghostcell_dims_primitive_var[ei][0];
                    //limiter reconstruction
                    limiterReconstruction(V[ei][idx_LL], V[ei][idx_L], V[ei][idx_R], V[ei][idx_RR],
                                          V_minus[ei][idx_face_x], V_plus[ei][idx_face_x]);


                }
            }
        }

        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */

        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
                bounded_flag_minus,
                primitive_variables_minus,
                DIRECTION::X_DIRECTION);

        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
                bounded_flag_plus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION);

        /*
         * Use first order interpolation if interpolated side primitive variables in x-direction
         * are out of bounds.
         */

        flag_minus = bounded_flag_minus->getPointer(0);
        flag_plus = bounded_flag_plus->getPointer(0);

        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];

            for (int j = 0; j < interior_dims[1]; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
#pragma omp simd
#endif
                for (int i = 0; i < interior_dims[0] + 1; i++)
                {
                    // Compute the linear indices.
                    const int idx_midpoint_x = i + j*(interior_dims[0] + 1);

                    const int idx_cell_L = (i - 1 + num_subghosts_0_primitive_var) +
                                           (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;

                    const int idx_cell_R = (i + num_subghosts_0_primitive_var) +
                                           (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;

                    if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0)
                    {
                        //std::cout <<"negative pressure or density" << std::endl;
                        V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                        V_plus[ei][idx_midpoint_x] = V[ei][idx_cell_R];
                    }
                }
            }
        }
        
        /*
         * Compute flux in the x-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            d_riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux,
                velocity_intercell,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            d_riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int j = 0; j < interior_dims[1]; j++)
            {
                for (int i = 0; i < interior_dims[0] + 1; i++)
                {
                    // Compute the linear index.
                    const int idx_face_x = i +
                        j*(interior_dims[0] + 1);

                    F_face_x[ei][idx_face_x] *= dt;
                }
            }
        }
        
        /*
         * Initialize temporary data containers for computing the flux in the y-direction.
         */
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {

            mirrorGhostCell2D(primitive_variables[vi], ghost_cell_maps, DIRECTION::Y_DIRECTION);
        }

        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(1);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int j = 0; j < interior_dims[1] + 1; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute the linear indices.
                    const int idx_face_y = i +
                        j*interior_dims[0];
                    
                    const int idx_B = (i + num_subghosts_primitive_var[ei][0]) +
                        (j - 1 + num_subghosts_primitive_var[ei][1])*
                            subghostcell_dims_primitive_var[ei][0];

                    const int idx_BB = (i + num_subghosts_primitive_var[ei][0]) +
                                      (j - 2 + num_subghosts_primitive_var[ei][1])*
                                      subghostcell_dims_primitive_var[ei][0];
                    
                    const int idx_T = (i + num_subghosts_primitive_var[ei][0]) +
                        (j + num_subghosts_primitive_var[ei][1])*
                            subghostcell_dims_primitive_var[ei][0];

                    const int idx_TT = (i + num_subghosts_primitive_var[ei][0]) +
                                      (j + 1 + num_subghosts_primitive_var[ei][1])*
                                      subghostcell_dims_primitive_var[ei][0];
                    //limiter reconstruction
                    limiterReconstruction(V[ei][idx_BB], V[ei][idx_B], V[ei][idx_T], V[ei][idx_TT],
                                          V_minus[ei][idx_face_y], V_plus[ei][idx_face_y]);
                }
            }
        }

        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */

        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
                bounded_flag_minus,
                primitive_variables_minus,
                DIRECTION::Y_DIRECTION);

        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
                bounded_flag_plus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION);


        /*
         * Use first order interpolation if interpolated side primitive variables in y-direction
         * are out of bounds.
         */

        flag_minus = bounded_flag_minus->getPointer(1);
        flag_plus = bounded_flag_plus->getPointer(1);

        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];

            for (int j = 0; j < interior_dims[1] + 1; j++)
            {
#ifdef HAMERS_ENABLE_SIMD
#pragma omp simd
#endif
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute the linear indices.
                    const int idx_midpoint_y = i + j*interior_dims[0];

                    const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                                           (j - 1 + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;

                    const int idx_cell_T = (i + num_subghosts_0_primitive_var) +
                                           (j + num_subghosts_1_primitive_var)*subghostcell_dim_0_primitive_var;

                    if (flag_minus[idx_midpoint_y] == 0 || flag_plus[idx_midpoint_y] == 0)
                    {
                        //std::cout <<"negative pressure or density" << std::endl;
                        V_minus[ei][idx_midpoint_y] = V[ei][idx_cell_B];
                        V_plus[ei][idx_midpoint_y] = V[ei][idx_cell_T];
                    }
                }
            }
        }
        
        /*
         * Compute flux in the y-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            d_riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux,
                velocity_intercell,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            d_riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int j = 0; j < interior_dims[1] + 1; j++)
            {
                for (int i = 0; i < interior_dims[0]; i++)
                {
                    // Compute the linear index.
                    const int idx_face_y = i +
                        j*interior_dims[0];
                    
                    F_face_y[ei][idx_face_y] *= dt;
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        if (d_has_advective_eqn_form)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    double* S = source->getPointer(ei);
                    
                    for (int j = 0; j < interior_dims[1]; j++)
                    {
                        for (int i = 0; i < interior_dims[0]; i++)
                        {
                            // Compute the linear indices.
                            const int idx_cell_wghost = (i + num_subghosts_primitive_var[ei][0]) +
                                (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0];
                            
                            const int idx_cell_nghost = i + j*interior_dims[0];
                            
                            const int idx_face_x_L = i +
                                j*(interior_dims[0] + 1);
                            
                            const int idx_face_x_R = (i + 1) +
                                j*(interior_dims[0] + 1);
                            
                            const int idx_face_y_B = i +
                                j*interior_dims[0];
                            
                            const int idx_face_y_T = i +
                                (j + 1)*interior_dims[0];
                            
                            const double& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                            const double& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                            
                            const double& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                            const double& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                            
                            S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*((u_R - u_L)/dx[0] + (v_T - v_B)/dx[1]);
                        }
                    }
                }
            }
        }

        
    } // if (d_dim == tbox::Dimension(2))
    else if (d_dim == tbox::Dimension(3))
    {
        /*
         * Get the pointers to the convective flux side data.
         */
        
        std::vector<double*> F_face_x;
        std::vector<double*> F_face_y;
        std::vector<double*> F_face_z;
        F_face_x.reserve(d_num_eqn);
        F_face_y.reserve(d_num_eqn);
        F_face_z.reserve(d_num_eqn);
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            F_face_x.push_back(convective_flux->getPointer(0, ei));
            F_face_y.push_back(convective_flux->getPointer(1, ei));
            F_face_z.push_back(convective_flux->getPointer(2, ei));
        }

        /*
         * Register primitive variables
         */
        std::unordered_map<std::string, hier::IntVector> num_subghosts_of_data;
        num_subghosts_of_data.insert(std::pair<std::string, hier::IntVector>("PRIMITIVE_VARIABLES", d_num_conv_ghosts));

        d_flow_model->registerDerivedCellVariable(num_subghosts_of_data);
        d_flow_model->computeGlobalDerivedCellData();

        /*
         * Build ghost cell map
         */
        std::vector<std::vector<std::array<int,4> > > ghost_cell_maps;
        ghost_cell_maps.resize(3);
        buildGhostCellMap3D(cell_status, ghost_cell_maps);
        
        /*
         * Get the pointers to the conservative variables.
         * The numbers of ghost cells and the dimensions of the ghost cell boxes are also determined.
         */
        
        std::vector<boost::shared_ptr<pdat::CellData<double> > > conservative_variables =
            d_flow_model->getGlobalCellDataConservativeVariables();

        std::vector<boost::shared_ptr<pdat::CellData<double> > > primitive_variables =
                d_flow_model->getGlobalCellDataPrimitiveVariables();
        
        std::vector<hier::IntVector> num_subghosts_primitive_var;
        num_subghosts_primitive_var.reserve(d_num_eqn);
        
        std::vector<hier::IntVector> subghostcell_dims_primitive_var;
        subghostcell_dims_primitive_var.reserve(d_num_eqn);
        
        std::vector<double*> Q;
        Q.reserve(d_num_eqn);

        std::vector<double*> V;
        V.reserve(d_num_eqn);
        
        int count_eqn = 0;
        
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {
            int depth = primitive_variables[vi]->getDepth();
            
            for (int di = 0; di < depth; di++)
            {
                // If the last element of the conservative variable vector is not in the system of equations, ignore it.
                if (count_eqn >= d_num_eqn)
                    break;

                V.push_back(primitive_variables[vi]->getPointer(di));
                Q.push_back(conservative_variables[vi]->getPointer(di));
                num_subghosts_primitive_var.push_back(primitive_variables[vi]->getGhostCellWidth());
                subghostcell_dims_primitive_var.push_back(primitive_variables[vi]->getGhostBox().numberCells());
                
                count_eqn++;
            }
        }
        
        /*
         * Declare temporary data containers for computing the fluxes at cell edges.
         */
        
        std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_minus;
        std::vector<boost::shared_ptr<pdat::SideData<double> > > primitive_variables_plus;

        primitive_variables_minus.reserve(d_num_eqn);
        primitive_variables_plus.reserve(d_num_eqn);
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            primitive_variables_minus.push_back(boost::make_shared<pdat::SideData<double> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));
            
            primitive_variables_plus.push_back(boost::make_shared<pdat::SideData<double> >(
                interior_box, 1, hier::IntVector::getZero(d_dim)));
        }
        
        std::vector<double*> V_minus;
        std::vector<double*> V_plus;
        V_minus.resize(d_num_eqn);
        V_plus.resize(d_num_eqn);

        boost::shared_ptr<pdat::SideData<int> > bounded_flag_minus;
        boost::shared_ptr<pdat::SideData<int> > bounded_flag_plus;

        bounded_flag_minus.reset(
                new pdat::SideData<int>(interior_box, 1, hier::IntVector::getZero(d_dim)));

        bounded_flag_plus.reset(
                new pdat::SideData<int>(interior_box, 1, hier::IntVector::getZero(d_dim)));

        int* flag_minus = nullptr;
        int* flag_plus = nullptr;

        /*
         * mirroring ghost cell data for computing the flux in the x-direction.
         */
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {

            mirrorGhostCell3D(primitive_variables[vi], ghost_cell_maps, DIRECTION::X_DIRECTION);
        }


        
        /*
         * Initialize temporary data containers for computing the flux in the x-direction.
         */
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(0);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(0);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++) {
            for (int k = 0; k < interior_dims[2]; k++) {
                for (int j = 0; j < interior_dims[1]; j++) {
                    for (int i = 0; i < interior_dims[0] + 1; i++) {
                        // Compute the linear indices.
                        const int idx_face_x = i +
                                               j * (interior_dims[0] + 1) +
                                               k * (interior_dims[0] + 1) * interior_dims[1];

                        const int idx_L = (i - 1 + num_subghosts_primitive_var[ei][0]) +
                                          (j + num_subghosts_primitive_var[ei][1]) *
                                          subghostcell_dims_primitive_var[ei][0] +
                                          (k + num_subghosts_primitive_var[ei][2]) *
                                          subghostcell_dims_primitive_var[ei][0] *
                                          subghostcell_dims_primitive_var[ei][1];

                        const int idx_LL = (i - 2 + num_subghosts_primitive_var[ei][0]) +
                                          (j + num_subghosts_primitive_var[ei][1]) *
                                          subghostcell_dims_primitive_var[ei][0] +
                                          (k + num_subghosts_primitive_var[ei][2]) *
                                          subghostcell_dims_primitive_var[ei][0] *
                                          subghostcell_dims_primitive_var[ei][1];

                        const int idx_R = (i + num_subghosts_primitive_var[ei][0]) +
                                          (j + num_subghosts_primitive_var[ei][1]) *
                                          subghostcell_dims_primitive_var[ei][0] +
                                          (k + num_subghosts_primitive_var[ei][2]) *
                                          subghostcell_dims_primitive_var[ei][0] *
                                          subghostcell_dims_primitive_var[ei][1];

                        const int idx_RR = (i + 1 + num_subghosts_primitive_var[ei][0]) +
                                          (j + num_subghosts_primitive_var[ei][1]) *
                                          subghostcell_dims_primitive_var[ei][0] +
                                          (k + num_subghosts_primitive_var[ei][2]) *
                                          subghostcell_dims_primitive_var[ei][0] *
                                          subghostcell_dims_primitive_var[ei][1];
                        //limiter reconstruction
                        limiterReconstruction(V[ei][idx_LL], V[ei][idx_L], V[ei][idx_R], V[ei][idx_RR],
                                              V_minus[ei][idx_face_x], V_plus[ei][idx_face_x]);

                    }
                }
            }
        }

        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */

        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
                bounded_flag_minus,
                primitive_variables_minus,
                DIRECTION::X_DIRECTION);

        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
                bounded_flag_plus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION);

        /*
         * Use first order interpolation if interpolated side primitive variables in x-direction
         * are out of bounds.
         */

        flag_minus = bounded_flag_minus->getPointer(0);
        flag_plus = bounded_flag_plus->getPointer(0);

        for (int ei = 0; ei < d_num_eqn; ei++) {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
            const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];

            for (int k = 0; k < interior_dims[2]; k++) {
                for (int j = 0; j < interior_dims[1]; j++) {
#ifdef HAMERS_ENABLE_SIMD
#pragma omp simd
#endif
                    for (int i = 0; i < interior_dims[0] + 1; i++) {
                        // Compute the linear indices.
                        const int idx_midpoint_x =
                                i + j * (interior_dims[0] + 1) + k * (interior_dims[0] + 1) * interior_dims[1];

                        const int idx_cell_L = (i - 1 + num_subghosts_0_primitive_var) +
                                               (j + num_subghosts_1_primitive_var) * subghostcell_dim_0_primitive_var +
                                               (k + num_subghosts_2_primitive_var) * subghostcell_dim_0_primitive_var *
                                               subghostcell_dim_1_primitive_var;

                        const int idx_cell_R = (i + num_subghosts_0_primitive_var) +
                                               (j + num_subghosts_1_primitive_var) * subghostcell_dim_0_primitive_var +
                                               (k + num_subghosts_2_primitive_var) * subghostcell_dim_0_primitive_var *
                                               subghostcell_dim_1_primitive_var;

                        if (flag_minus[idx_midpoint_x] == 0 || flag_plus[idx_midpoint_x] == 0) {
                            //std::cout <<"negative pressure or density" << std::endl;
                            V_minus[ei][idx_midpoint_x] = V[ei][idx_cell_L];
                            V_plus[ei][idx_midpoint_x] = V[ei][idx_cell_R];
                        }
                    }
                }
            }
        }
        
        /*
         * Compute flux in the x-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            d_riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux,
                velocity_intercell,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            d_riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::X_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0] + 1; i++) {
                        // Compute the linear index.
                        const int idx_face_x = i +
                                               j * (interior_dims[0] + 1) +
                                               k * (interior_dims[0] + 1) * interior_dims[1];

                        F_face_x[ei][idx_face_x] *= dt;
                    }
                }
            }
        }
        
        /*
         * Initialize temporary data containers for computing the flux in the y-direction.
         */

        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {

            mirrorGhostCell3D(primitive_variables[vi], ghost_cell_maps, DIRECTION::Y_DIRECTION);
        }

        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(1);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(1);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++) {
            for (int k = 0; k < interior_dims[2]; k++) {
                for (int j = 0; j < interior_dims[1] + 1; j++) {
                    for (int i = 0; i < interior_dims[0]; i++) {
                        // Compute the linear indices.
                        const int idx_face_y = i +
                                               j * interior_dims[0] +
                                               k * interior_dims[0] * (interior_dims[1] + 1);

                        // Compute the linear indices.
                        const int idx_B = (i + num_subghosts_primitive_var[ei][0]) +
                                          (j - 1 + num_subghosts_primitive_var[ei][1]) *
                                          subghostcell_dims_primitive_var[ei][0] +
                                          (k + num_subghosts_primitive_var[ei][2]) *
                                          subghostcell_dims_primitive_var[ei][0] *
                                          subghostcell_dims_primitive_var[ei][1];

                        const int idx_BB = (i + num_subghosts_primitive_var[ei][0]) +
                                           (j - 2 + num_subghosts_primitive_var[ei][1]) *
                                           subghostcell_dims_primitive_var[ei][0] +
                                           (k + num_subghosts_primitive_var[ei][2]) *
                                           subghostcell_dims_primitive_var[ei][0] *
                                           subghostcell_dims_primitive_var[ei][1];

                        const int idx_T = (i + num_subghosts_primitive_var[ei][0]) +
                                          (j + num_subghosts_primitive_var[ei][1]) *
                                          subghostcell_dims_primitive_var[ei][0] +
                                          (k + num_subghosts_primitive_var[ei][2]) *
                                          subghostcell_dims_primitive_var[ei][0] *
                                          subghostcell_dims_primitive_var[ei][1];

                        const int idx_TT = (i + num_subghosts_primitive_var[ei][0]) +
                                          (j + 1 + num_subghosts_primitive_var[ei][1]) *
                                          subghostcell_dims_primitive_var[ei][0] +
                                          (k + num_subghosts_primitive_var[ei][2]) *
                                          subghostcell_dims_primitive_var[ei][0] *
                                          subghostcell_dims_primitive_var[ei][1];

                        //limiter reconstruction
                        limiterReconstruction(V[ei][idx_BB], V[ei][idx_B], V[ei][idx_T], V[ei][idx_TT],
                                              V_minus[ei][idx_face_y], V_plus[ei][idx_face_y]);
                    }
                }
            }
        }
        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */

        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
                bounded_flag_minus,
                primitive_variables_minus,
                DIRECTION::Y_DIRECTION);

        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
                bounded_flag_plus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION);

        /*
         * Use first order interpolation if interpolated side primitive variables in y-direction
         * are out of bounds.
         */

        flag_minus = bounded_flag_minus->getPointer(1);
        flag_plus = bounded_flag_plus->getPointer(1);

        for (int ei = 0; ei < d_num_eqn; ei++) {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
            const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];

            for (int k = 0; k < interior_dims[2]; k++) {
                for (int j = 0; j < interior_dims[1] + 1; j++) {
#ifdef HAMERS_ENABLE_SIMD
#pragma omp simd
#endif
                    for (int i = 0; i < interior_dims[0]; i++) {
                        // Compute the linear indices.
                        const int idx_midpoint_y = i + j * interior_dims[0] + k*interior_dims[0]*(interior_dims[1]+1);

                        const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                                               (j - 1 + num_subghosts_1_primitive_var) * subghostcell_dim_0_primitive_var +
                                               (k + num_subghosts_2_primitive_var) * subghostcell_dim_0_primitive_var *
                                               subghostcell_dim_1_primitive_var;

                        const int idx_cell_T = (i + num_subghosts_0_primitive_var) +
                                               (j + num_subghosts_1_primitive_var) * subghostcell_dim_0_primitive_var +
                                               (k + num_subghosts_2_primitive_var) * subghostcell_dim_0_primitive_var *
                                               subghostcell_dim_1_primitive_var;

                        if (flag_minus[idx_midpoint_y] == 0 || flag_plus[idx_midpoint_y] == 0) {
                            //std::cout <<"negative pressure or density" << std::endl;
                            V_minus[ei][idx_midpoint_y] = V[ei][idx_cell_B];
                            V_plus[ei][idx_midpoint_y] = V[ei][idx_cell_T];
                        }
                    }
                }
            }
        }

        
        /*
         * Compute flux in the y-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            d_riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux,
                velocity_intercell,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            d_riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Y_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = 0; k < interior_dims[2]; k++)
            {
                for (int j = 0; j < interior_dims[1] + 1; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute the linear index.
                        const int idx_face_y = i +
                            j*interior_dims[0] +
                            k*interior_dims[0]*(interior_dims[1] + 1);

                        F_face_y[ei][idx_face_y] *= dt;
                    }
                }
            }
        }
        
        /*
         * Initialize temporary data containers for computing the flux in the z-direction.
         */
        for (int vi = 0; vi < static_cast<int>(primitive_variables.size()); vi++)
        {

            mirrorGhostCell3D(primitive_variables[vi], ghost_cell_maps, DIRECTION::Z_DIRECTION);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            V_minus[ei] = primitive_variables_minus[ei]->getPointer(2);
            V_plus[ei] = primitive_variables_plus[ei]->getPointer(2);
        }
        
        for (int ei = 0; ei < d_num_eqn; ei++) {
            for (int k = 0; k < interior_dims[2] + 1; k++) {
                for (int j = 0; j < interior_dims[1]; j++) {
                    for (int i = 0; i < interior_dims[0]; i++) {
                        // Compute the linear indices.
                        const int idx_face_z = i +
                                               j * interior_dims[0] +
                                               k * interior_dims[0] * interior_dims[1];

                        const int idx_B = (i + num_subghosts_primitive_var[ei][0]) +
                                          (j + num_subghosts_primitive_var[ei][1]) *
                                          subghostcell_dims_primitive_var[ei][0] +
                                          (k - 1 + num_subghosts_primitive_var[ei][2]) *
                                          subghostcell_dims_primitive_var[ei][0] *
                                          subghostcell_dims_primitive_var[ei][1];

                        const int idx_BB = (i + num_subghosts_primitive_var[ei][0]) +
                                           (j + num_subghosts_primitive_var[ei][1]) *
                                           subghostcell_dims_primitive_var[ei][0] +
                                           (k - 2 + num_subghosts_primitive_var[ei][2]) *
                                           subghostcell_dims_primitive_var[ei][0] *
                                           subghostcell_dims_primitive_var[ei][1];


                        const int idx_F = (i + num_subghosts_primitive_var[ei][0]) +
                                          (j + num_subghosts_primitive_var[ei][1]) *
                                          subghostcell_dims_primitive_var[ei][0] +
                                          (k + num_subghosts_primitive_var[ei][2]) *
                                          subghostcell_dims_primitive_var[ei][0] *
                                          subghostcell_dims_primitive_var[ei][1];

                        const int idx_FF = (i + num_subghosts_primitive_var[ei][0]) +
                                           (j + num_subghosts_primitive_var[ei][1]) *
                                           subghostcell_dims_primitive_var[ei][0] +
                                           (k + 1 + num_subghosts_primitive_var[ei][2]) *
                                           subghostcell_dims_primitive_var[ei][0] *
                                           subghostcell_dims_primitive_var[ei][1];

                        //limiter reconstruction
                        limiterReconstruction(V[ei][idx_BB], V[ei][idx_B], V[ei][idx_F], V[ei][idx_FF],
                                              V_minus[ei][idx_face_z], V_plus[ei][idx_face_z]);
                    }
                }
            }
        }
        /*
         * Check whether the interpolated side primitive variables are within the bounds.
         */

        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
                bounded_flag_minus,
                primitive_variables_minus,
                DIRECTION::Z_DIRECTION);

        d_flow_model->checkGlobalSideDataPrimitiveVariablesBounded(
                bounded_flag_plus,
                primitive_variables_plus,
                DIRECTION::Z_DIRECTION);

        /*
         * Use first order interpolation if interpolated side primitive variables in y-direction
         * are out of bounds.
         */

        flag_minus = bounded_flag_minus->getPointer(2);
        flag_plus = bounded_flag_plus->getPointer(2);

        for (int ei = 0; ei < d_num_eqn; ei++) {
            const int num_subghosts_0_primitive_var = num_subghosts_primitive_var[ei][0];
            const int num_subghosts_1_primitive_var = num_subghosts_primitive_var[ei][1];
            const int num_subghosts_2_primitive_var = num_subghosts_primitive_var[ei][2];
            const int subghostcell_dim_0_primitive_var = subghostcell_dims_primitive_var[ei][0];
            const int subghostcell_dim_1_primitive_var = subghostcell_dims_primitive_var[ei][1];

            for (int k = 0; k < interior_dims[2] + 1; k++) {

                for (int j = 0; j < interior_dims[1]; j++) {
#ifdef HAMERS_ENABLE_SIMD
#pragma omp simd
#endif
                    for (int i = 0; i < interior_dims[0]; i++) {
                        // Compute the linear indices.
                        const int idx_midpoint_z = i + j * interior_dims[0] + k * interior_dims[0] * interior_dims[1];

                        const int idx_cell_B = (i + num_subghosts_0_primitive_var) +
                                               (j + num_subghosts_1_primitive_var) * subghostcell_dim_0_primitive_var +
                                               (k - 1 + num_subghosts_2_primitive_var) *
                                               subghostcell_dim_0_primitive_var *
                                               subghostcell_dim_1_primitive_var;

                        const int idx_cell_F = (i + num_subghosts_0_primitive_var) +
                                               (j + num_subghosts_1_primitive_var) * subghostcell_dim_0_primitive_var +
                                               (k + num_subghosts_2_primitive_var) * subghostcell_dim_0_primitive_var *
                                               subghostcell_dim_1_primitive_var;;

                        if (flag_minus[idx_midpoint_z] == 0 || flag_plus[idx_midpoint_z] == 0) {
                            //std::cout <<"negative pressure or density" << std::endl;
                            V_minus[ei][idx_midpoint_z] = V[ei][idx_cell_B];
                            V_plus[ei][idx_midpoint_z] = V[ei][idx_cell_F];
                        }
                    }
                }
            }
        }
        
        /*
         * Compute flux in the z-direction.
         */
        
        if (d_has_advective_eqn_form)
        {
            d_riemann_solver->computeConvectiveFluxAndVelocityFromPrimitiveVariables(
                convective_flux,
                velocity_intercell,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Z_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        else
        {
            d_riemann_solver->computeConvectiveFluxFromPrimitiveVariables(
                convective_flux,
                primitive_variables_minus,
                primitive_variables_plus,
                DIRECTION::Z_DIRECTION,
                RIEMANN_SOLVER::HLLC);
        }
        
        // Multiply flux by dt.
        for (int ei = 0; ei < d_num_eqn; ei++)
        {
            for (int k = 0; k < interior_dims[2] + 1; k++)
            {
                for (int j = 0; j < interior_dims[1]; j++)
                {
                    for (int i = 0; i < interior_dims[0]; i++)
                    {
                        // Compute the linear index.
                        const int idx_face_z = i +
                            j*interior_dims[0] +
                            k*interior_dims[0]*interior_dims[1];

                        F_face_z[ei][idx_face_z] *= dt;
                    }
                }
            }
        }
        
        /*
         * Compute the source.
         */
        
        if (d_has_advective_eqn_form)
        {
            for (int ei = 0; ei < d_num_eqn; ei++)
            {
                if (d_eqn_form[ei] == EQN_FORM::ADVECTIVE)
                {
                    double* S = source->getPointer(ei);
                    
                    for (int k = 0; k < interior_dims[2]; k++)
                    {
                        for (int j = 0; j < interior_dims[1]; j++)
                        {
                            for (int i = 0; i < interior_dims[0]; i++)
                            {
                                // Compute the linear indices. 
                                const int idx_cell_wghost = (i + num_subghosts_primitive_var[ei][0]) +
                                    (j + num_subghosts_primitive_var[ei][1])*subghostcell_dims_primitive_var[ei][0] +
                                    (k + num_subghosts_primitive_var[ei][2])*subghostcell_dims_primitive_var[ei][0]*
                                        subghostcell_dims_primitive_var[ei][1];
                                
                                const int idx_cell_nghost = i +
                                    j*interior_dims[0] +
                                    k*interior_dims[0]*interior_dims[1];
                                
                                const int idx_face_x_L = i +
                                    j*(interior_dims[0] + 1) +
                                    k*(interior_dims[0] + 1)*interior_dims[1];
                                
                                const int idx_face_x_R = (i + 1) +
                                    j*(interior_dims[0] + 1) +
                                    k*(interior_dims[0] + 1)*interior_dims[1];
                                
                                const int idx_face_y_B = i +
                                    j*interior_dims[0] +
                                    k*interior_dims[0]*(interior_dims[1] + 1);
                                
                                const int idx_face_y_T = i +
                                    (j + 1)*interior_dims[0] +
                                    k*interior_dims[0]*(interior_dims[1] + 1);
                                
                                const int idx_face_z_B = i +
                                    j*interior_dims[0] +
                                    k*interior_dims[0]*interior_dims[1];
                                
                                const int idx_face_z_F = i +
                                    j*interior_dims[0] +
                                    (k + 1)*interior_dims[0]*interior_dims[1];
                                
                                const double& u_L = velocity_intercell->getPointer(0, 0)[idx_face_x_L];
                                const double& u_R = velocity_intercell->getPointer(0, 0)[idx_face_x_R];
                                
                                const double& v_B = velocity_intercell->getPointer(1, 1)[idx_face_y_B];
                                const double& v_T = velocity_intercell->getPointer(1, 1)[idx_face_y_T];
                                
                                const double& w_B = velocity_intercell->getPointer(2, 2)[idx_face_z_B];
                                const double& w_F = velocity_intercell->getPointer(2, 2)[idx_face_z_F];
                                
                                S[idx_cell_nghost] += dt*Q[ei][idx_cell_wghost]*(
                                    (u_R - u_L)/dx[0] + (v_T - v_B)/dx[1] + (w_F - w_B)/dx[2]);
                            }
                        }
                    }
                }
            }
        }
        
        /*
         * Unregister the patch in the flow model.
         */
        

        
    } // if (d_dim == tbox::Dimension(3))
    d_flow_model->unregisterPatch();
}

/*
 * Put the characteristics of the convective flux reconstruction class
 * into the restart database.
 */
void
ConvectiveFluxReconstructorSecondOrderHLLC::limiterReconstruction
    (const double v_ll, const double v_l, const double v_r, const double v_rr,
     double & v_l_rec, double & v_r_rec)
{
    double sigma = std::numeric_limits<double>::infinity();
    if(fabs(v_l - v_r) > sigma*fabs(v_l + v_r)/2.0) {
        v_l_rec = v_l;
        v_r_rec = v_r;
        return;
    }

    double eps = 1.e-15;
    double ind_l = (v_r - v_l)/(v_l - v_ll + eps);
    double ind_r = (v_r - v_l)/(v_rr - v_r + eps);

    double phi_l = ind_l > 0. ? 2*ind_l/(ind_l + 1.): 0.;
    double phi_r = ind_r > 0. ? 2*ind_r/(ind_r + 1.): 0.;


    v_l_rec = v_l + 0.5*phi_l*(v_l - v_ll);
    v_r_rec = v_r - 0.5*phi_r*(v_rr - v_r);

}