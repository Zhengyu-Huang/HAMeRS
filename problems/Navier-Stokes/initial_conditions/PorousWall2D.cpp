#include "apps/Navier-Stokes/NavierStokesInitialConditions.hpp"

/*
 * Set the data on the patch interior to some initial values.
 */
void
NavierStokesInitialConditions::initializeDataOnPatch(
    hier::Patch& patch,
    const std::vector<boost::shared_ptr<pdat::CellData<double> > >& conservative_variables,
    const double data_time,
    const bool initial_time) {
    NULL_USE(data_time);

    if ((d_project_name != "2D porous wall NASA-exp") &&
        (d_project_name != "2D porous wall Mars")) {
        TBOX_ERROR(d_object_name
                           << ": "
                           << "Can only initialize data for 'project_name' = '2D porous wall NASA-exp' or "
                           << "'2D porous wall Mars'!\n"
                           << "'project_name' = '"
                           << d_project_name
                           << "' is given."
                           << std::endl);
    }

    if (d_dim != tbox::Dimension(2)) {
        TBOX_ERROR(d_object_name
                           << ": "
                           << "Dimension of problem should be 2!"
                           << std::endl);
    }

    if (d_flow_model_type != FLOW_MODEL::SINGLE_SPECIES) {
        TBOX_ERROR(d_object_name
                           << ": "
                           << "Flow model should be single-species!"
                           << std::endl);
    }

    if (initial_time) {
        const boost::shared_ptr <geom::CartesianPatchGeometry> patch_geom(
                BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                        patch.getPatchGeometry()));

#ifdef HAMERS_DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(patch_geom);
#endif

        const double *const dx = patch_geom->getDx();
        const double *const patch_xlo = patch_geom->getXLower();

// Get the dimensions of box that covers the interior of Patch.
        hier::Box patch_box = patch.getBox();
        const hier::IntVector patch_dims = patch_box.numberCells();

/*
 * Initialize data for a 2D shock-vortex interaction problem.
 */

        boost::shared_ptr <pdat::CellData<double>> density = conservative_variables[0];
        boost::shared_ptr <pdat::CellData<double>> momentum = conservative_variables[1];
        boost::shared_ptr <pdat::CellData<double>> total_energy = conservative_variables[2];

        double *rho = density->getPointer(0);
        double *rho_u = momentum->getPointer(0);
        double *rho_v = momentum->getPointer(1);
        double *E = total_energy->getPointer(0);
        if (d_project_name == "2D porous wall NASA-exp") {
            const double gamma = double(7) / double(5);


// Post-shock condition.
            const double rho_post = double(1.0);
            //const double p_post = double(790.3842940685046);
            const double p_post = double(946.25);
            const double u_post = double(0.0);
            const double v_post = double(0.0);

// Pre-shock condition.
            //const double rho_pre = double(0.987347926730015);
            //const double p_pre = double(780.3842940685046);
            const double rho_pre = double(0.9947159841479526);
            const double p_pre = double(941.25);
            const double u_pre = double(0.0);
            const double v_pre = double(0.0);

//seperation between preshock and postshock
            double x0 = 0.0;

            for (int j = 0; j < patch_dims[1]; j++) {
                for (int i = 0; i < patch_dims[0]; i++) {
// Compute index into linear data array.
                    int idx_cell = i + j * patch_dims[0];

// Compute the coordinates.
                    double x[2];
                    x[0] = patch_xlo[0] + (double(i) + double(1) / double(2)) * dx[0];
                    x[1] = patch_xlo[1] + (double(j) + double(1) / double(2)) * dx[1];

                    if (x[0] < x0) {
                        rho[idx_cell] = rho_post;
                        rho_u[idx_cell] = rho_post * u_post;
                        rho_v[idx_cell] = rho_post * v_post;
                        E[idx_cell] = p_post / (gamma - double(1)) + double(1) / double(2) * rho_post *
                                                                     (u_post * u_post + v_post * v_post);
                    } else {

                        rho[idx_cell] = rho_pre;
                        rho_u[idx_cell] = rho_pre * u_pre;
                        rho_v[idx_cell] = rho_pre * v_pre;
                        E[idx_cell] = p_pre / (gamma - double(1)) + double(1) / double(2) * rho_pre *
                                                                    (u_pre * u_pre + v_pre * v_pre);
                    }

                }
            }
        }
        if (d_project_name == "2D porous wall Mars") {
            const double gamma = 1.33;


//// Post-shock condition.
//            const double rho_post = double(1.0);
//            const double p_post = double(0.7518796992481204);
//            const double u_post = double(1.0);
//            const double v_post = double(0.0);
//// Pre-shock condition.
//            const double rho_pre = double(0.3649922771465224);
//            const double p_pre = double(0.1770957021238189);
//            const double u_pre = double(0.0);
//            const double v_pre = double(0.0);


// Post-shock condition.
            const double rho_post = double(1.0);
            const double p_post = double(0.7518796992481204);
            const double u_post = double(1.2);
            const double v_post = double(0.0);
// Pre-shock condition.
            const double rho_pre = double(0.30276434298491006);
            const double p_pre = double(0.12658085012624493);
            const double u_pre = double(0.0);
            const double v_pre = double(0.0);

//// Post-shock condition.
//            const double rho_post = double(1.0);
//            const double p_post = double(0.7518796992481204);
//            const double u_post = double(0.6);
//            const double v_post = double(0.0);
//// Pre-shock condition.
//            const double rho_pre = double(0.5419203691155778);
//            const double p_pre = double(0.32599013828916323);
//            const double u_pre = double(0.0);
//            const double v_pre = double(0.0);



//seperation between preshock and postshock
            double x0 = -0.5;

            for (int j = 0; j < patch_dims[1]; j++) {
                for (int i = 0; i < patch_dims[0]; i++) {
// Compute index into linear data array.
                    int idx_cell = i + j * patch_dims[0];

// Compute the coordinates.
                    double x[2];
                    x[0] = patch_xlo[0] + (double(i) + double(1) / double(2)) * dx[0];
                    x[1] = patch_xlo[1] + (double(j) + double(1) / double(2)) * dx[1];


                    if (x[0] < x0) {
                        rho[idx_cell] = rho_post;
                        rho_u[idx_cell] = rho_post * u_post;
                        rho_v[idx_cell] = rho_post * v_post;
                        E[idx_cell] = p_post / (gamma - double(1)) + double(1) / double(2) * rho_post *
                                                                     (u_post * u_post + v_post * v_post);
                    } else {

                        rho[idx_cell] = rho_pre;
                        rho_u[idx_cell] = rho_pre * u_pre;
                        rho_v[idx_cell] = rho_pre * v_pre;
                        E[idx_cell] = p_pre / (gamma - double(1)) + double(1) / double(2) * rho_pre *
                                                                    (u_pre * u_pre + v_pre * v_pre);
                    }
                }
            }
        }
    }
}

