#include "apps/Euler/EulerInitialConditions.hpp"

/*
 * Set the data on the patch interior to some initial values.
 */
void
EulerInitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
    const double data_time,
    const bool initial_time)
{
    NULL_USE(data_time);
    
    if (d_project_name != "2D advection of density wave")
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Can only initialize data for 'project_name' = '2D advection of density wave'!\n"
            << "'project_name' = '"
            << d_project_name
            << "' is given."
            << std::endl);
    }
    
    if (d_dim != tbox::Dimension(2))
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Dimension of problem should be 2!"
            << std::endl);
    }
    
    if (d_flow_model_type != FLOW_MODEL::SINGLE_SPECIES)
    {
        TBOX_ERROR(d_object_name
            << ": "
            << "Flow model should be single-species!"
            << std::endl);
    }
    
    if (initial_time)
    {
        const double* const domain_xlo = d_grid_geometry->getXLower();
        const double* const domain_xhi = d_grid_geometry->getXUpper();
        
        const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                patch.getPatchGeometry()));
        
#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(patch_geom);
#endif
        
        const double* const dx = patch_geom->getDx();
        const double* const patch_xlo = patch_geom->getXLower();
        
        // Get the dimensions of box that covers the interior of Patch.
        hier::Box patch_box = patch.getBox();
        const hier::IntVector patch_dims = patch_box.numberCells();
        
        /*
         * Initialize data for a 2D density wave advection problem.
         */
        
        boost::shared_ptr<pdat::CellData<double> > density      = conservative_variables[0];
        boost::shared_ptr<pdat::CellData<double> > momentum     = conservative_variables[1];
        boost::shared_ptr<pdat::CellData<double> > total_energy = conservative_variables[2];
        
        double* rho   = density->getPointer(0);
        double* rho_u = momentum->getPointer(0);
        double* rho_v = momentum->getPointer(1);
        double* E     = total_energy->getPointer(0);
        
        const double x_a = double(1)/double(3)*(domain_xlo[0] + domain_xhi[0]);
        const double x_b = double(2)/double(3)*(domain_xlo[0] + domain_xhi[0]);
        
        const double y_a = double(1)/double(3)*(domain_xlo[1] + domain_xhi[1]);
        const double y_b = double(2)/double(3)*(domain_xlo[1] + domain_xhi[1]);
        
        double gamma = double(7)/double(5);
        
        // Initial conditions inside the square.
        double rho_i = double(10);
        double u_i   = double(1);
        double v_i   = double(1);
        double p_i   = double(1);
        
        // Initial conditions outside the square.
        double rho_o = double(1);
        double u_o   = double(1);
        double v_o   = double(1);
        double p_o   = double(1);
        
        for (int j = 0; j < patch_dims[1]; j++)
        {
            for (int i = 0; i < patch_dims[0]; i++)
            {
                // Compute index into linear data array.
                int idx_cell = i + j*patch_dims[0];
                
                // Compute the coordinates.
                double x[2];
                x[0] = patch_xlo[0] + (double(i) + double(1)/double(2))*dx[0];
                x[1] = patch_xlo[1] + (double(j) + double(1)/double(2))*dx[1];
                
                if ((x[0] >= x_a) && (x[0] <= x_b) && (x[1] >= y_a) && (x[1] <= y_b))
                {
                    rho[idx_cell]   = rho_i;
                    rho_u[idx_cell] = rho_i*u_i;
                    rho_v[idx_cell] = rho_i*v_i;
                    E[idx_cell]     = p_i/(gamma - double(1)) + double(1)/double(2)*rho_i*
                        (u_i*u_i + v_i*v_i);
                }
                else
                {
                    rho[idx_cell]   = rho_o;
                    rho_u[idx_cell] = rho_o*u_o;
                    rho_v[idx_cell] = rho_o*v_o;
                    E[idx_cell]     = p_o/(gamma - double(1)) + double(1)/double(2)*rho_o*
                        (u_o*u_o + v_o*v_o);
                }
            }
        }
    }
}
