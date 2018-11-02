#include "util/basic_geometry/WallTreatment.hpp"

#include <algorithm>
#include <math.h>       /* fabs */
/*
 * This function decide point (x, y, z) is in the wall or not
 * return 0 if the point is in the wall
 * return 1 if the point is outside
 * the mesh size is 1/(32N).
 *
 */
int isOutsidePorousWall1(int dim, double x, double y, double z)
{
   /*
    * The porous wall is centered at x0.
    * The computational domain is at least 32 by 32 by 32
    * There are 4 holes in y and z directions, each of size 2,
    *           4 solid in y and z directions, each of size 6;
    * 3 2 6 2 6 2 6 2 3
    */

    double x0 = 0.0, r_hole = 1./32;
    if(fabs(x - x0) > r_hole) return 1;
    if(dim == 2){
        double yy = y/r_hole;
        if(yy < 3 || (yy > 5 && yy < 11) || (yy > 13 && yy < 19) || (yy > 21 && yy < 27) || yy > 29)
            return 0;
        else
            return 1;


   }
   return 1;

}



/* This is simple debugging setup
 * This function decide point (x, y, z) is in the wall or not
 * return 0 if the point is in the wall
 * return 1 if the point is outside
 * the mesh size is 1/(8N).
 *
 */
int isOutsidePorousWall2X(int dim, double x, double y, double z)
{
    /* 10 percents!!!!
     * The porous wall is centered at x0. and this is a QUASI 2D structure
     * The computational domain is at least is 100 by 100
     * There are 1 square holes in x-y directions at x=0.0,y=0.0 of size 8,
     *           2 solid in y directions, each of size 46;
     * 46 8 46 ()
     */

    double x0 = 0.0, r_hole = 4./100;
    if(fabs(x - x0) > 2*r_hole) return 1;
    if(dim == 2 || dim == 3){
        if(y < 0.5 - r_hole || (y > 0.5 + r_hole))
            return 0;
        else
            return 1;
    }
    return 1;
}

/* This is simple debugging setup
 * This function decide point (x, y, z) is in the wall or not
 * return 0 if the point is in the wall
 * return 1 if the point is outside
 * the mesh size is 1/(8N).
 *
 */
int isOutsidePorousWall2Y(int dim, double x, double y, double z)
{
    /*
     * The porous wall is centered at x0. and this is a QUASI 2D structure
     * The computational domain is at least is 8 by 8 by 8
     * There are 1 holes in y direction at y = 0.0, of size 2,
     *           2 solid in y direction, each of size 3;
     * 3 2 3
     */

    double y0 = 0.0, r_hole = 1./8;
    if(fabs(y - y0) > r_hole) return 1;
    if(dim == 2 || dim == 3){
        double xx = x/r_hole;
        if(xx < 3 || (xx > 5))
            return 0;
        else
            return 1;
    }
    return 1;
}


/* This is simple debugging setup
 * This function decide point (x, y, z) is in the wall or not
 * return 0 if the point is in the wall
 * return 1 if the point is outside
 * the mesh size is 1/(8N).
 *
 */
int isOutsidePorousWall3X(int dim, double x, double y, double z)
{
    /*
     * The porous wall is centered at x0. this is QUASI 1D structure
     * The computational domain is at least is 8 by 8 by 8
     * A solid wall at x = 0
     */
    double x0 = 0.0, r_hole = 1./8;
    if(fabs(x - x0) > r_hole) return 1;
    else return 0;
}

int isOutsidePorousWall3Y(int dim, double x, double y, double z)
{
    /*
     * The porous wall is centered at y0. this is QUASI 1D structure
     * The computational domain is at least is 8 by 8 by 8
     * A solid wall at y = 0
     */

    double y0 = 0.0, r_hole = 1./8;
    if(fabs(y - y0) > r_hole) return 1;
    else return 0;
}

int isOutsidePorousWall3Z(int dim, double x, double y, double z)
{
    /*
     * The porous wall is centered at z0. this is QUASI 1D structure
     * The computational domain is at least is 8 by 8 by 8
     * A solid wall at y = 0
     */

    double z0 = 0.0, r_hole = 1./8;
    if(fabs(z - z0) > r_hole) return 1;
    else return 0;
}

/* This is simple debugging setup
 * This function decide point (x, y, z) is in the wall or not
 * return 0 if the point is in the wall
 * return 1 if the point is outside
 * the mesh size is 1/(8N).
 *
 */
int isOutsidePorousWall4X(int dim, double x, double y, double z)
{
    /* 8.16 percents = (4./14)**2
     * The porous wall is centered at x0. and this is a 3D structure
     * The computational domain is at least is 14 by 14 by 14
     * There are 1 holes in x direction at x = 0.0, of size 4,
     *                ########
     *                ########
     *                ########
     *                ########
     *                ########
     *                ###  ###
     *                ###  ###
     *                ###  ###
     *                ###  ###
     *                ########
     *                ########
     *                ########
     *                ########
     *                ########
     * 5 4 5
     */

    double x0 = 0.0, r_hole = 1./7, hole_half_thickness = 8./100;
    if(fabs(x - x0) > hole_half_thickness) return 1;
    if(dim == 2 || dim == 3){
        if((y < 0.5 - r_hole|| y > 0.5 + r_hole) || (z < 0.5 - r_hole|| z > 0.5 + r_hole))
            return 0;
        else
            return 1;
    }
    return 1;
}

/* This is simple debugging setup
 * This function decide point (x, y, z) is in the wall or not
 * return 0 if the point is in the wall
 * return 1 if the point is outside
 * the mesh size is 1/(8N).
 *
 */
int isOutsidePorousWall4X_2(int dim, double x, double y, double z)
{
    /* 8.0 percents = (2./10)*(4./10)
     * The porous wall is centered at x0. and this is a 3D structure
     * The computational domain is at least is 14 by 14 by 14
     * There are 1 holes in x direction at x = 0.0, of size 4,
     *                ########
     *                ########
     *                ########
     *                ###  ###
     *                ###  ###
     *                ########
     *                ########
     *                ########
     * 5 4 5
     */

    double x0 = 0.0, r_hole_y = 1./10, r_hole_z = 2./10, hole_half_thickness = 8./100;
    if(fabs(x - x0) > hole_half_thickness) return 1;
    if(dim == 2 || dim == 3){
        if((y < 0.5 - r_hole_y|| y > 0.5 + r_hole_y) || (z < 0.5 - r_hole_z|| z > 0.5 + r_hole_z))
            return 0;
        else
            return 1;
    }
    return 1;
}

int isOutsidePorousWall4X_3(int dim, double x, double y, double z)
{
    /* 8.0 percents = (1./10)*(8./10)
     * The porous wall is centered at x0. and this is a 3D structure
     * The computational domain is at least is 14 by 14 by 14
     * There are 1 holes in x direction at x = 0.0, of size 4,
     *                ########
     *                ########
     *                ########
     *                ###  ###
     *                ###  ###
     *                ########
     *                ########
     *                ########
     * 5 4 5
     */

    double x0 = 0.0, r_hole_y = 1./20, r_hole_z = 4./10, hole_half_thickness = 8./100;
    if(fabs(x - x0) > hole_half_thickness) return 1;
    if(dim == 2 || dim == 3){
        if((y < 0.5 - r_hole_y|| y > 0.5 + r_hole_y) || (z < 0.5 - r_hole_z|| z > 0.5 + r_hole_z))
            return 0;
        else
            return 1;
    }
    return 1;
}


void
initializeCellStatus(hier::Patch& patch,
                     boost::shared_ptr<pdat::CellData<double> > &cell_status) {
    // Get the grid spacing.
    const boost::shared_ptr <geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                    patch.getPatchGeometry()));


    const double *const dx = patch_geom->getDx();
    const double *const x_lo = patch_geom->getXLower();

    // Get the dimensions of box that covers the interior of Patch.
    const tbox::Dimension d_dim = patch.getDim();

    hier::Box patch_box = patch.getBox();
    const hier::IntVector patch_dims = patch_box.numberCells();


    const hier::IntVector interior_dims = cell_status->getBox().numberCells();

    const hier::IntVector d_num_conv_ghosts = cell_status->getGhostCellWidth();

    double *cell_status_data = cell_status->getPointer(0);
    if (d_dim == tbox::Dimension(1)) {}
    else if (d_dim == tbox::Dimension(2)) {
        /*
         * Compute the status of cells, not including ghost cells
         */
        for (int j = 0; j < patch_dims[1]; j++) {
            for (int i = 0; i < patch_dims[0]; i++) {
                const double x = x_lo[0] + (i + 0.5) * dx[0],
                        y = x_lo[1] + (j + 0.5) * dx[1];
                const int idx = i + j * patch_dims[0];
                cell_status_data[idx] = isOutsidePorousWall2X(2, x, y);
            }
        }


    } else if (d_dim == tbox::Dimension(3)) {
        /*
         * Compute the status of cells, not including ghost cells
         */
        for (int k = 0; k < patch_dims[2]; k++) {
            for (int j = 0; j < patch_dims[1]; j++) {
                for (int i = 0; i < patch_dims[0]; i++) {
                    const double x = x_lo[0] + (i + 0.5) * dx[0],
                                 y = x_lo[1] + (j + 0.5) * dx[1],
                                 z = x_lo[2] + (k + 0.5) * dx[2];
                    const int idx = i + j * patch_dims[0] + k*patch_dims[0]*patch_dims[1];
                    cell_status_data[idx] = isOutsidePorousWall4X(3, x, y, z);
                }
            }
        }
    }
}


void
mirrorGhostCell(boost::shared_ptr<pdat::CellData<double> > &variables,
       const boost::shared_ptr<pdat::CellData<double> > &cell_status,
       const DIRECTION::TYPE d_direction,
       const WALL_TREATMENT_CONDITION d_condition) {


    const tbox::Dimension d_dim = cell_status->getDim();
    const hier::IntVector interior_dims = cell_status->getBox().numberCells();

    const hier::IntVector d_num_status_ghosts = cell_status->getGhostCellWidth();

    const hier::IntVector d_num_var_ghosts = variables->getGhostCellWidth();


    double *cell_status_data = cell_status->getPointer(0);

    const int depth = variables->getDepth();
    std::vector<double *> V;
    V.reserve(depth);
    for (int di = 0; di < depth; di++) {
        V.push_back(variables->getPointer(di));
    }


    int ghost_count;//Set to -inf to start get rid of these block that never use in current computation
    if (d_direction == DIRECTION::X_DIRECTION) {
        if (d_dim == tbox::Dimension(1)) {}
        else if (d_dim == tbox::Dimension(2)) {
            //loop y direction
            for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {
                ghost_count = std::numeric_limits<int>::min();
                //loop x direction
                for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                           (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                           (interior_dims[0] + 2 * d_num_status_ghosts[0]);

                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                        (i - 2 * ghost_count + 1) >= 0) {
                        const int idx_mirr =
                                (i - 2 * ghost_count + 1) + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);


                        if (depth == 1) {
                            V[0][idx] = V[0][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_SLIP) {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_NO_SLIP) {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                        } else if (depth == 4 && d_condition ==
                                                 WALL_SLIP) {   //Convective Flux in x direction rhou, rhouu +p, rhouv, u(E + p)
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = -V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        } else if (depth == 4 && d_condition == WALL_NO_SLIP) {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        }
                    }
                }

                ghost_count = std::numeric_limits<int>::min();
                //loop x direction inversely
                for (int i = interior_dims[0] + 2 * d_num_var_ghosts[0] - 1; i >= 0; i--) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                           (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                           (interior_dims[0] + 2 * d_num_status_ghosts[0]);
                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                        (i + 2 * ghost_count - 1) <= interior_dims[0] + 2 * d_num_var_ghosts[0] - 1) {
                        const int idx_mirr =
                                (i + 2 * ghost_count - 1) + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);


                        if (depth == 1) {
                            V[0][idx] = V[0][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_SLIP) {   //velocity u, v
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_NO_SLIP) {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                        } else if (depth == 4 && d_condition ==
                                                 WALL_SLIP) {   //Convective Flux in x direction rhou, rhouu +p, rhouv, u(E + p)
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = -V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        } else if (depth == 4 && d_condition == WALL_NO_SLIP) {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        }

                    }

                }
            }
        } else if (d_dim == tbox::Dimension(3)) {
            //loop z direction
            for (int k = 0; k < interior_dims[2] + 2 * d_num_var_ghosts[2]; k++) {
                //loop y direction
                for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {
                    ghost_count = std::numeric_limits<int>::min();
                    //loop x direction
                    for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            (i - 2 * ghost_count + 1) >= 0) {
                            const int idx_mirr =
                                    (i - 2 * ghost_count + 1) + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                    k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                    (interior_dims[1] + 2 * d_num_var_ghosts[1]);


                            if (depth == 1) {
                                V[0][idx] = V[0][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_SLIP) {
                                //Convective Flux in x direction rhou, rhouu +p, rhouv, rhouw, u(E + p)
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                                V[3][idx] = -V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            }

                        }
                    }

                    ghost_count = std::numeric_limits<int>::min();
                    //loop x direction inversely
                    for (int i = interior_dims[0] + 2 * d_num_var_ghosts[0] - 1; i >= 0; i--) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            (i + 2 * ghost_count - 1) <= interior_dims[0] + 2 * d_num_var_ghosts[0] - 1) {
                            const int idx_mirr =
                                    (i + 2 * ghost_count - 1) + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                    k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                    (interior_dims[1] + 2 * d_num_var_ghosts[1]);

                            if (depth == 1) {
                                V[0][idx] = V[0][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_SLIP) {   //velocity u, v
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_SLIP) {
                                //Convective Flux in x direction rhou, rhouu +p, rhouv, rhouw, u(E + p)
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                                V[3][idx] = -V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            }

                        }

                    }
                }
            }

        }
    }
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else if (d_direction == DIRECTION::Y_DIRECTION) {
        if (d_dim == tbox::Dimension(2)) {
            //loop x direction
            for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {

                ghost_count = std::numeric_limits<int>::min();
                //loop y direction
                for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                           (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                           (interior_dims[0] + 2 * d_num_status_ghosts[0]);
                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                        j - 2 * ghost_count + 1 >= 0) {

                        const int idx_mirr =
                                i + (j - 2 * ghost_count + 1) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);

                        if (depth == 1) {
                            V[0][idx] = V[0][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_SLIP) {
                            V[0][idx] = V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_NO_SLIP) {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                        } else if (depth == 4 && d_condition ==
                                                 WALL_SLIP) {   //Convective Flux in x direction rho v, rho u v, rho v v + p, v(E + p)
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        } else if (depth == 4 && d_condition == WALL_NO_SLIP) {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        }
                    }
                }
                ghost_count = std::numeric_limits<int>::min();
                //loop y direction inversely
                for (int j = interior_dims[1] + 2 * d_num_var_ghosts[1] - 1; j >= 0; j--) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                           (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                           (interior_dims[0] + 2 * d_num_status_ghosts[0]);
                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                        j + 2 * ghost_count - 1 <= interior_dims[1] + 2 * d_num_var_ghosts[1] - 1) {

                        const int idx_mirr =
                                i + (j + 2 * ghost_count - 1) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);

                        if (depth == 1) {
                            V[0][idx] = V[0][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_SLIP) {
                            V[0][idx] = V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_NO_SLIP) {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                        } else if (depth == 4 && d_condition ==
                                                 WALL_SLIP) {   //Convective Flux in x direction rho v, rho u v, rho v v + p, v(E + p)
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        } else if (depth == 4 && d_condition == WALL_NO_SLIP) {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        }

                    }

                }
            }
        } else if (d_dim == tbox::Dimension(3)) {
            //loop z direction
            for (int k = 0; k < interior_dims[2] + 2 * d_num_var_ghosts[2]; k++) {
                //loop x direction
                for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {
                    ghost_count = std::numeric_limits<int>::min();
                    //loop y direction
                    for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            j - 2 * ghost_count + 1 >= 0) {

                            const int idx_mirr = i +
                                                 (j - 2 * ghost_count + 1) *
                                                 (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                                 k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                                 (interior_dims[1] + 2 * d_num_var_ghosts[1]);

                            if (depth == 1) {
                                V[0][idx] = V[0][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_SLIP) {
                                V[0][idx] = V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_SLIP) {
                                //Convective Flux in x direction rho v, rho vu, rho v v + p, rho v w, v(E + p)
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = -V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            }
                        }
                    }
                    ghost_count = std::numeric_limits<int>::min();
                    //loop y direction inversely
                    for (int j = interior_dims[1] + 2 * d_num_var_ghosts[1] - 1; j >= 0; j--) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            j + 2 * ghost_count - 1 <= interior_dims[1] + 2 * d_num_var_ghosts[1] - 1) {

                            const int idx_mirr =
                                    i + (j + 2 * ghost_count - 1) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                    k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                    (interior_dims[1] + 2 * d_num_var_ghosts[1]);

                            if (depth == 1) {
                                V[0][idx] = V[0][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_SLIP) {
                                V[0][idx] = V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_SLIP) {
                                //Convective Flux in x direction rho v, rho v u, rho v v + p, rho v w, v(E + p)
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = -V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            }
                        }
                    }
                }
            }
        }
    } else if (d_direction == DIRECTION::Z_DIRECTION) {
        if (d_dim == tbox::Dimension(3)) {
            //loop y direction
            for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {
                //loop x direction
                for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {
                    //loop z direction
                    ghost_count = std::numeric_limits<int>::min();
                    for (int k = 0; k < interior_dims[2] + 2 * d_num_var_ghosts[2]; k++) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            k - 2 * ghost_count + 1 >= 0) {

                            const int idx_mirr = i +
                                                 j *
                                                 (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                                 (k - 2 * ghost_count + 1) *
                                                 (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                                 (interior_dims[1] + 2 * d_num_var_ghosts[1]);

                            if (depth == 1) {
                                V[0][idx] = V[0][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_SLIP) {
                                V[0][idx] = V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_SLIP) {
                                //Convective Flux in x direction rho w, rho wu, rho w v, rho w w + p, w(E + p)
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                                V[3][idx] = V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            }
                        }
                    }
                    ghost_count = std::numeric_limits<int>::min();
                    //loop z direction inversely
                    for (int k = interior_dims[2] + 2 * d_num_var_ghosts[2] - 1; k >= 0; k--) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            k + 2 * ghost_count - 1 <= interior_dims[2] + 2 * d_num_var_ghosts[2] - 1) {

                            const int idx_mirr =
                                    i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                    (k + 2 * ghost_count - 1) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                    (interior_dims[1] + 2 * d_num_var_ghosts[1]);

                            if (depth == 1) {
                                V[0][idx] = V[0][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_SLIP) {
                                V[0][idx] = V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                            } else if (depth == 3 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_SLIP) {
                                //Convective Flux in x direction rho v, rho v u, rho v v + p, rho v w, v(E + p)
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                                V[3][idx] = V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            } else if (depth == 5 && d_condition == WALL_NO_SLIP) {
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                            }
                        }
                    }
                }
            }
        }

    }

}

void
mirrorGhostCellDerivative(std::vector<boost::shared_ptr<pdat::CellData<double> > >du_x,
                const boost::shared_ptr<pdat::CellData<double> > &cell_status,
                const DIRECTION::TYPE d_direction) {

    int depth = du_x.size();

    const tbox::Dimension d_dim = cell_status->getDim();
    const hier::IntVector interior_dims = cell_status->getBox().numberCells();

    const hier::IntVector d_num_status_ghosts = cell_status->getGhostCellWidth();

    double *cell_status_data = cell_status->getPointer(0);

    const hier::IntVector d_num_var_ghosts = du_x[0]->getGhostCellWidth();

    std::vector<double *> V;
    V.reserve(depth);
    for (int di = 0; di < depth; di++) {
        V.push_back(du_x[di]->getPointer(0));
    }

    int ghost_count;//Set to -inf to start get rid of these block that never use in current computation
    if (d_direction == DIRECTION::X_DIRECTION) {
        if (d_dim == tbox::Dimension(1)) {}
        else if (d_dim == tbox::Dimension(2)) {
            //loop y direction
            for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {
                ghost_count = std::numeric_limits<int>::min();
                //loop x direction
                for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                           (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                           (interior_dims[0] + 2 * d_num_status_ghosts[0]);
                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                        (i - 2 * ghost_count + 1) >= 0) {
                        const int idx_mirr =
                                (i - 2 * ghost_count + 1) + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);

                        if (depth == 5) {   // du/dx, dv/dx, dT/dx,  du/dy, dv/dy
                            V[0][idx] = V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = -V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                            V[4][idx] = -V[4][idx_mirr];
                        }
                    }
                }

                ghost_count = std::numeric_limits<int>::min();
                //loop x direction inversely
                for (int i = interior_dims[0] + 2 * d_num_var_ghosts[0] - 1; i >= 0; i--) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                           (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                           (interior_dims[0] + 2 * d_num_status_ghosts[0]);
                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                        (i + 2 * ghost_count - 1) <= interior_dims[0] + 2 * d_num_var_ghosts[0] - 1) {
                        const int idx_mirr =
                                (i + 2 * ghost_count - 1) + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                        if (depth == 5) {   // du/dx, dv/dx, dT/dx,  du/dy, dv/dy
                            V[0][idx] = V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = -V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                            V[4][idx] = -V[4][idx_mirr];
                        }
                    }
                }
            }
        } else if (d_dim == tbox::Dimension(3)) {
            //loop z direction
            for (int k = 0; k < interior_dims[2] + 2 * d_num_var_ghosts[2]; k++) {
                //loop y direction
                for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {
                    ghost_count = std::numeric_limits<int>::min();
                    //loop x direction
                    for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            (i - 2 * ghost_count + 1) >= 0) {
                            const int idx_mirr =
                                    (i - 2 * ghost_count + 1) + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                    k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                    (interior_dims[1] + 2 * d_num_var_ghosts[1]);

                            if (depth == 8) { // du/dx, dv/dx, dw/x, dT/dx,  du/dy, dv/dy, du/dz, dw/dz
                                V[0][idx] = V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = -V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                                V[5][idx] = -V[5][idx_mirr];
                                V[6][idx] = -V[6][idx_mirr];
                                V[7][idx] = -V[7][idx_mirr];
                            }
                       }
                    }

                    ghost_count = std::numeric_limits<int>::min();
                    //loop x direction inversely
                    for (int i = interior_dims[0] + 2 * d_num_var_ghosts[0] - 1; i >= 0; i--) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            (i + 2 * ghost_count - 1) <= interior_dims[0] + 2 * d_num_var_ghosts[0] - 1) {
                            const int idx_mirr =
                                    (i + 2 * ghost_count - 1) + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                    k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                    (interior_dims[1] + 2 * d_num_var_ghosts[1]);

                            if (depth == 8) { // du/dx, dv/dx, dw/x, dT/dx,  du/dy, dv/dy, du/dz, dw/dz
                                V[0][idx] = V[0][idx_mirr];
                                V[1][idx] = V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = -V[3][idx_mirr];
                                V[4][idx] = -V[4][idx_mirr];
                                V[5][idx] = -V[5][idx_mirr];
                                V[6][idx] = -V[6][idx_mirr];
                                V[7][idx] = -V[7][idx_mirr];
                            }
                       }
                    }
                }
            }
        }
    }
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else if (d_direction == DIRECTION::Y_DIRECTION) {
        if (d_dim == tbox::Dimension(2)) {
            //loop x direction
            for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {

                ghost_count = std::numeric_limits<int>::min();
                //loop y direction
                for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                           (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                           (interior_dims[0] + 2 * d_num_status_ghosts[0]);
                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                        j - 2 * ghost_count + 1 >= 0) {

                        const int idx_mirr =
                                i + (j - 2 * ghost_count + 1) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);

                        if (depth == 5) {   // du/dx, dv/dx,  du/dy, dv/dy, dT/dy
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = V[3][idx_mirr];
                            V[4][idx] = -V[4][idx_mirr];
                        }
                    }
                }
                ghost_count = std::numeric_limits<int>::min();
                //loop y direction inversely
                for (int j = interior_dims[1] + 2 * d_num_var_ghosts[1] - 1; j >= 0; j--) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                           (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                           (interior_dims[0] + 2 * d_num_status_ghosts[0]);
                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                        j + 2 * ghost_count - 1 <= interior_dims[1] + 2 * d_num_var_ghosts[1] - 1) {

                        const int idx_mirr =
                                i + (j + 2 * ghost_count - 1) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);

                        if (depth == 5) {   // du/dx, dv/dx,  du/dy, dv/dy, dT/dy
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = V[3][idx_mirr];
                            V[4][idx] = -V[4][idx_mirr];
                        }
                    }

                }
            }
        } else if (d_dim == tbox::Dimension(3)) {

            //loop z direction
            for (int k = 0; k < interior_dims[2] + 2 * d_num_var_ghosts[2]; k++) {
                //loop x direction
                for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {
                    ghost_count = std::numeric_limits<int>::min();
                    //loop y direction
                    for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            j - 2 * ghost_count + 1 >= 0) {

                            const int idx_mirr = i +
                                                 (j - 2 * ghost_count + 1) *
                                                 (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                                 k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                                 (interior_dims[1] + 2 * d_num_var_ghosts[1]);

                            if (depth == 8) { // du/dx, dv/dx, du/dy, dv/dy, dw/dy, dT/dy, dv/dz, dw/dz
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = V[3][idx_mirr];
                                V[4][idx] = V[4][idx_mirr];
                                V[5][idx] = -V[5][idx_mirr];
                                V[6][idx] = -V[6][idx_mirr];
                                V[7][idx] = -V[7][idx_mirr];
                            }
                        }
                    }
                    ghost_count = std::numeric_limits<int>::min();
                    //loop y direction inversely
                    for (int j = interior_dims[1] + 2 * d_num_var_ghosts[1] - 1; j >= 0; j--) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            j + 2 * ghost_count - 1 <= interior_dims[1] + 2 * d_num_var_ghosts[1] - 1) {

                            const int idx_mirr =
                                    i + (j + 2 * ghost_count - 1) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                    k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                    (interior_dims[1] + 2 * d_num_var_ghosts[1]);


                            if (depth == 8) { // du/dx, dv/dx, du/dy, dv/dy, dw/dy, dT/dy, dv/dz, dw/dz
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = V[2][idx_mirr];
                                V[3][idx] = V[3][idx_mirr];
                                V[4][idx] = V[4][idx_mirr];
                                V[5][idx] = -V[5][idx_mirr];
                                V[6][idx] = -V[6][idx_mirr];
                                V[7][idx] = -V[7][idx_mirr];
                            }
                       }

                    }
                }
            }


        }
    } else if (d_direction == DIRECTION::Z_DIRECTION) {
        if (d_dim == tbox::Dimension(3)) {
            //loop y direction
            for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {
                //loop x direction
                for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {
                    ghost_count = std::numeric_limits<int>::min();
                    //loop z direction
                    for (int k = 0; k < interior_dims[2] + 2 * d_num_var_ghosts[2]; k++) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            k - 2 * ghost_count + 1 >= 0) {

                            const int idx_mirr = i +
                                                 j *
                                                 (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                                 (k - 2 * ghost_count + 1) *
                                                 (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                                 (interior_dims[1] + 2 * d_num_var_ghosts[1]);


                            if (depth == 8) { // du/dx, dw/dx, dv/dy, dw/dy, du/dz, dv/dz, dw/dz, dT/dz
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                                V[3][idx] = -V[3][idx_mirr];
                                V[4][idx] = V[4][idx_mirr];
                                V[5][idx] = V[5][idx_mirr];
                                V[6][idx] = V[6][idx_mirr];
                                V[7][idx] = -V[7][idx_mirr];
                            }
                        }

                    }

                    ghost_count = std::numeric_limits<int>::min();
                    //loop z direction inversely
                    for (int k = interior_dims[2] + 2 * d_num_var_ghosts[2] - 1; k >= 0; k--) {

                        const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                        k * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                        (interior_dims[1] + 2 * d_num_var_ghosts[1]);
                        const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) +
                                               (j - d_num_var_ghosts[1] + d_num_status_ghosts[1]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                               (k - d_num_var_ghosts[2] + d_num_status_ghosts[2]) *
                                               (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                               (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                        if (cell_status_data[idx_status] < 0.5) {
                            ghost_count++;
                        } else {
                            ghost_count = 0;
                        }
                        if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] &&
                            k + 2 * ghost_count - 1 <= interior_dims[2] + 2 * d_num_var_ghosts[2] - 1) {

                            const int idx_mirr =
                                    i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                                    (k + 2 * ghost_count - 1) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                                    (interior_dims[1] + 2 * d_num_var_ghosts[1]);

                            if (depth == 8) { // du/dx, dw/dx, dv/dy, dw/dy, du/dz, dv/dz, dw/dz, dT/dz
                                V[0][idx] = -V[0][idx_mirr];
                                V[1][idx] = -V[1][idx_mirr];
                                V[2][idx] = -V[2][idx_mirr];
                                V[3][idx] = -V[3][idx_mirr];
                                V[4][idx] = V[4][idx_mirr];
                                V[5][idx] = V[5][idx_mirr];
                                V[6][idx] = V[6][idx_mirr];
                                V[7][idx] = -V[7][idx_mirr];
                            }

                        }
                    }
                }
            }
        }
    }
}


void
buildGhostCellMap2D(const boost::shared_ptr<pdat::CellData<double> > &cell_status,
                  std::vector<std::vector<std::array<int,3> > >& ghost_cell_maps) {

    const hier::IntVector interior_dims = cell_status->getBox().numberCells();

    const hier::IntVector d_num_status_ghosts = cell_status->getGhostCellWidth();

    double *cell_status_data = cell_status->getPointer(0);

    std::vector<std::array<int,3> > & x_map = ghost_cell_maps[0];
    std::vector<std::array<int,3> > & y_map = ghost_cell_maps[1];


    int ghost_count;//Set to -inf to start get rid of these block that never use in current computation

    // DIRECTION::X_DIRECTION
    // loop y direction
    for (int j = -d_num_status_ghosts[1]; j < interior_dims[1] + d_num_status_ghosts[1]; j++) {
        ghost_count = std::numeric_limits<int>::min();
        //loop x direction
        for (int i = -d_num_status_ghosts[0]; i < interior_dims[0] + d_num_status_ghosts[0]; i++) {

            const int idx_status = (i + d_num_status_ghosts[0]) +
                                   (j + d_num_status_ghosts[1]) * (interior_dims[0] + 2 * d_num_status_ghosts[0]);
            if (cell_status_data[idx_status] < 0.5) {
                ghost_count++;
            } else {
                ghost_count = 0;
            }
            if (ghost_count >= 1 && ghost_count <= d_num_status_ghosts[0] &&
                (i - 2 * ghost_count + 1) >= -d_num_status_ghosts[0]) {
                //(i,j) --- > (i - 2 * ghost_count + 1 , j)
                x_map.push_back(std::array<int, 3 >{i,j, i - 2 * ghost_count + 1});
            }
        }

        ghost_count = std::numeric_limits<int>::min();
        //loop x direction inversely
        for (int i = interior_dims[0] + d_num_status_ghosts[0] - 1; i >= -d_num_status_ghosts[0]; i--) {


            const int idx_status = (i + d_num_status_ghosts[0]) +
                                   (j + d_num_status_ghosts[1]) * (interior_dims[0] + 2 * d_num_status_ghosts[0]);
            if (cell_status_data[idx_status] < 0.5) {
                ghost_count++;
            } else {
                ghost_count = 0;
            }
            if (ghost_count >= 1 && ghost_count <= d_num_status_ghosts[0] &&
                (i + 2 * ghost_count - 1) <= interior_dims[0] + d_num_status_ghosts[0] - 1) {

                // (i,j) -- > (i + 2 * ghost_count - 1 , j)
                x_map.push_back(std::array<int, 3 >{i,j, i + 2 * ghost_count - 1});

            }

        }
    }
    // DIRECTION::Y_DIRECTION
    //loop x direction
    for (int i = -d_num_status_ghosts[0]; i < interior_dims[0] + d_num_status_ghosts[0]; i++) {
        ghost_count = std::numeric_limits<int>::min();
        //loop y direction
        for (int j = -d_num_status_ghosts[1]; j < interior_dims[1] + d_num_status_ghosts[1]; j++) {

            const int idx_status = (i + d_num_status_ghosts[0]) +
                                   (j + d_num_status_ghosts[1]) * (interior_dims[0] + 2 * d_num_status_ghosts[0]);
            if (cell_status_data[idx_status] < 0.5) {
                ghost_count++;
            } else {
                ghost_count = 0;
            }
            if (ghost_count >= 1 && ghost_count <= d_num_status_ghosts[1] &&
                j - 2 * ghost_count + 1 >= -d_num_status_ghosts[1]) {
                //(i, j) --> (i, j - 2 * ghost_count + 1)
                y_map.push_back(std::array<int, 3 >{i,j, j - 2 * ghost_count + 1});

            }
        }
        ghost_count = std::numeric_limits<int>::min();
        //loop y direction inversely
        for (int j = interior_dims[1] + d_num_status_ghosts[1] - 1; j >= -d_num_status_ghosts[1]; j--) {

            const int idx_status = (i + d_num_status_ghosts[0]) +
                                   (j + d_num_status_ghosts[1]) *
                                   (interior_dims[0] + 2 * d_num_status_ghosts[0]);
            if (cell_status_data[idx_status] < 0.5) {
                ghost_count++;
            } else {
                ghost_count = 0;
            }
            if (ghost_count >= 1 && ghost_count <= d_num_status_ghosts[1] &&
                j + 2 * ghost_count - 1 <= interior_dims[1] + d_num_status_ghosts[1] - 1) {
                // (i,j) --> (i , j + 2 * ghost_count - 1)
                y_map.push_back(std::array<int, 3 >{i,j, j + 2 * ghost_count - 1});

            }

        }
    }

}


void
buildGhostCellMap3D(const boost::shared_ptr<pdat::CellData<double> > &cell_status,
                    std::vector<std::vector<std::array<int,4> > > & ghost_cell_maps) {

    const hier::IntVector interior_dims = cell_status->getBox().numberCells();

    const hier::IntVector d_num_status_ghosts = cell_status->getGhostCellWidth();

    double *cell_status_data = cell_status->getPointer(0);

    std::vector<std::array<int,4> > & x_map = ghost_cell_maps[0];
    std::vector<std::array<int,4> > & y_map = ghost_cell_maps[1];
    std::vector<std::array<int,4> > & z_map = ghost_cell_maps[2];

    int ghost_count;//Set to -inf to start get rid of these block that never use in current computation


    // DIRECTION::X_DIRECTION
    //loop z direction
    for (int k = -d_num_status_ghosts[2]; k < interior_dims[2] + d_num_status_ghosts[2]; k++) {
        //loop y direction
        for (int j = -d_num_status_ghosts[1]; j < interior_dims[1] + d_num_status_ghosts[1]; j++) {
            ghost_count = std::numeric_limits<int>::min();
            //loop x direction
            for (int i = -d_num_status_ghosts[0]; i < interior_dims[0] + d_num_status_ghosts[0]; i++) {

                const int idx_status = (i + d_num_status_ghosts[0]) +
                                       (j + d_num_status_ghosts[1]) *
                                       (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                       (k + d_num_status_ghosts[2]) *
                                       (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                       (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                if (cell_status_data[idx_status] < 0.5) {
                    ghost_count++;
                } else {
                    ghost_count = 0;
                }
                if (ghost_count >= 1 && ghost_count <= d_num_status_ghosts[0] &&
                    (i - 2 * ghost_count + 1) >= -d_num_status_ghosts[0]) {
                    //(i,j,k) --> (i - 2 * ghost_count + 1, j, k)
                    x_map.push_back(std::array<int, 4 >{i,j,k, i - 2 * ghost_count + 1});



                }
            }

            ghost_count = std::numeric_limits<int>::min();
            //loop x direction inversely
            for (int i = interior_dims[0] + d_num_status_ghosts[0] - 1; i >= -d_num_status_ghosts[0]; i--) {


                const int idx_status = (i + d_num_status_ghosts[0]) +
                                       (j + d_num_status_ghosts[1]) *
                                       (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                       (k + d_num_status_ghosts[2]) *
                                       (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                       (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                if (cell_status_data[idx_status] < 0.5) {
                    ghost_count++;
                } else {
                    ghost_count = 0;
                }
                if (ghost_count >= 1 && ghost_count <= d_num_status_ghosts[0] &&
                    (i + 2 * ghost_count - 1) <= interior_dims[0] + d_num_status_ghosts[0] - 1) {
                    //(i,j,k) --> (i + 2 * ghost_count - 1, j, k)
                    x_map.push_back(std::array<int, 4 >{i,j,k, i + 2 * ghost_count - 1});


                }

            }
        }
    }
    // DIRECTION::Y_DIRECTION
    //loop z direction
    for (int k = -d_num_status_ghosts[2]; k < interior_dims[2] + d_num_status_ghosts[2]; k++) {
        //loop x direction
        for (int i = -d_num_status_ghosts[0]; i < interior_dims[0] + d_num_status_ghosts[0]; i++) {
            ghost_count = std::numeric_limits<int>::min();
            //loop y direction
            for (int j = -d_num_status_ghosts[1]; j < interior_dims[1] + d_num_status_ghosts[1]; j++) {

                const int idx_status = (i + d_num_status_ghosts[0]) +
                                       (j + d_num_status_ghosts[1]) *
                                       (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                       (k + d_num_status_ghosts[2]) *
                                       (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                       (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                if (cell_status_data[idx_status] < 0.5) {
                    ghost_count++;
                } else {
                    ghost_count = 0;
                }
                if (ghost_count >= 1 && ghost_count <= d_num_status_ghosts[1] &&
                    j - 2 * ghost_count + 1 >= -d_num_status_ghosts[1]) {
                    //(i,j,k) --> (i, j - 2 * ghost_count + 1, k)
                    y_map.push_back(std::array<int, 4 >{i,j,k, j - 2 * ghost_count + 1});

                }
            }
            ghost_count = std::numeric_limits<int>::min();
            //loop y direction inversely
            for (int j = interior_dims[1] + d_num_status_ghosts[1] - 1; j >= -d_num_status_ghosts[1]; j--) {

                const int idx_status = (i + d_num_status_ghosts[0]) +
                                       (j + d_num_status_ghosts[1]) *
                                       (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                       (k + d_num_status_ghosts[2]) *
                                       (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                       (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                if (cell_status_data[idx_status] < 0.5) {
                    ghost_count++;
                } else {
                    ghost_count = 0;
                }
                if (ghost_count >= 1 && ghost_count <= d_num_status_ghosts[1] &&
                    j + 2 * ghost_count - 1 <= interior_dims[1] + d_num_status_ghosts[1] - 1) {
                    //(i,j,k) --> (i, j + 2 * ghost_count - 1, k)
                    y_map.push_back(std::array<int, 4 >{i,j,k, j + 2 * ghost_count - 1});

                }

            }
        }
    }

    // DIRECTION::Z_DIRECTION
    //loop y direction
    for (int j = -d_num_status_ghosts[1]; j < interior_dims[1] + d_num_status_ghosts[1]; j++) {
    //loop x direction
        for (int i = -d_num_status_ghosts[0]; i < interior_dims[0] + d_num_status_ghosts[0]; i++) {
            //loop z direction
            ghost_count = std::numeric_limits<int>::min();
            for (int k = -d_num_status_ghosts[2]; k < interior_dims[2] + d_num_status_ghosts[2]; k++) {

                const int idx_status = (i + d_num_status_ghosts[0]) +
                                       (j + d_num_status_ghosts[1]) * (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                       (k + d_num_status_ghosts[2]) * (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                       (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                if (cell_status_data[idx_status] < 0.5) {
                    ghost_count++;
                } else {
                    ghost_count = 0;
                }
                if (ghost_count >= 1 && ghost_count <= d_num_status_ghosts[2] &&
                    k - 2 * ghost_count + 1 >= -d_num_status_ghosts[2]) {
                    //(i,j,k) --> (i + j, k - 2 * ghost_count + 1)
                    z_map.push_back(std::array<int, 4 >{i,j,k, k - 2 * ghost_count + 1});

                }
            }
            ghost_count = std::numeric_limits<int>::min();
            //loop z direction inversely
            for (int k = interior_dims[2] + d_num_status_ghosts[2] - 1; k >= -d_num_status_ghosts[2]; k--) {

                const int idx_status = (i + d_num_status_ghosts[0]) +
                                       (j + d_num_status_ghosts[1]) * (interior_dims[0] + 2 * d_num_status_ghosts[0]) +
                                       (k + d_num_status_ghosts[2]) * (interior_dims[0] + 2 * d_num_status_ghosts[0]) *
                                       (interior_dims[1] + 2 * d_num_status_ghosts[1]);
                if (cell_status_data[idx_status] < 0.5) {
                    ghost_count++;
                } else {
                    ghost_count = 0;
                }
                if (ghost_count >= 1 && ghost_count <= d_num_status_ghosts[2] &&
                    k + 2 * ghost_count - 1 <= interior_dims[2] + d_num_status_ghosts[2] - 1) {
                    //(i,j,k) -->(i,j, k + 2 * ghost_count - 1)
                    z_map.push_back(std::array<int, 4 >{i,j,k, k + 2 * ghost_count - 1});
                }

            }
        }
    }
}




void
mirrorGhostCell2D(boost::shared_ptr<pdat::CellData<double> > &variables,
                const std::vector<std::vector<std::array<int,3> > > & ghost_cell_maps,
                const DIRECTION::TYPE d_direction,
                const WALL_TREATMENT_CONDITION d_condition) {

    const hier::IntVector interior_dims = variables->getBox().numberCells();

    const hier::IntVector d_num_var_ghosts = variables->getGhostCellWidth();

    const int depth = variables->getDepth();
    std::vector<double *> V;
    V.reserve(depth);
    for (int di = 0; di < depth; di++) {
        V.push_back(variables->getPointer(di));
    }
    if (d_direction == DIRECTION::X_DIRECTION) {
        const std::vector <std::array<int, 3> > &x_map = ghost_cell_maps[0];
        int i, j, i_mirr, idx, idx_mirr;
        for (std::vector < std::array < int, 3> > ::const_iterator it = x_map.begin(); it != x_map.end();
        ++it)
        {
            i = (*it)[0];
            j = (*it)[1];
            i_mirr = (*it)[2];

            if (i < - d_num_var_ghosts[0] || i > interior_dims[0] + d_num_var_ghosts[0] - 1
                || j < - d_num_var_ghosts[1] || j > interior_dims[1] + d_num_var_ghosts[1] - 1
                || i_mirr < - d_num_var_ghosts[0] || i_mirr > interior_dims[0] + d_num_var_ghosts[0] - 1 )
                continue;

            idx = (i + d_num_var_ghosts[0]) + (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
            idx_mirr = (i_mirr + d_num_var_ghosts[0])+ (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);

            if (depth == 1) {
                V[0][idx] = V[0][idx_mirr];
            } else if (depth == 2 && d_condition == WALL_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
            } else if (depth == 2 && d_condition == WALL_NO_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
            } else if (depth == 4 && d_condition == WALL_SLIP) {
                //Convective Flux in x direction rhou, rhouu +p, rhouv, u(E + p)
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
                V[2][idx] = -V[2][idx_mirr];
                V[3][idx] = -V[3][idx_mirr];
            } else if (depth == 4 && d_condition == WALL_NO_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
                V[3][idx] = -V[3][idx_mirr];
            }
        }
    } else if (d_direction == DIRECTION::Y_DIRECTION) {
        const std::vector <std::array<int, 3>> &y_map = ghost_cell_maps[1];
        int i, j, j_mirr, idx, idx_mirr;
        for (std::vector < std::array < int, 3 > > ::const_iterator it = y_map.begin(); it != y_map.end();
        ++it)
        {
            i = (*it)[0];
            j = (*it)[1];
            j_mirr = (*it)[2];
            if (i < - d_num_var_ghosts[0] || i > interior_dims[0] + d_num_var_ghosts[0] - 1
                || j < - d_num_var_ghosts[1] || j > interior_dims[1] + d_num_var_ghosts[1] - 1
                || j_mirr < - d_num_var_ghosts[1] || j_mirr > interior_dims[1] + d_num_var_ghosts[1] - 1 )
                continue;
            idx = (i + d_num_var_ghosts[0]) + (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
            idx_mirr = (i + d_num_var_ghosts[0])+ (j_mirr + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
            if (depth == 1) {
                V[0][idx] = V[0][idx_mirr];
            } else if (depth == 2 && d_condition == WALL_SLIP) {
                V[0][idx] = V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
            } else if (depth == 2 && d_condition == WALL_NO_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
            } else if (depth == 4 && d_condition == WALL_SLIP) {
                //Convective Flux in x direction rho v, rho u v, rho v v + p, v(E + p)
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
                V[3][idx] = -V[3][idx_mirr];
            } else if (depth == 4 && d_condition == WALL_NO_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
                V[3][idx] = -V[3][idx_mirr];
            }
        }
    }
}


void
mirrorGhostCell3D(boost::shared_ptr<pdat::CellData<double> > &variables,
                const std::vector<std::vector<std::array<int,4> > >& ghost_cell_maps,
                const DIRECTION::TYPE d_direction,
                const WALL_TREATMENT_CONDITION d_condition) {

    const hier::IntVector interior_dims = variables->getBox().numberCells();

    const hier::IntVector d_num_var_ghosts = variables->getGhostCellWidth();

    const int depth = variables->getDepth();
    std::vector<double *> V;
    V.reserve(depth);
    for (int di = 0; di < depth; di++) {
        V.push_back(variables->getPointer(di));
    }

    if (d_direction == DIRECTION::X_DIRECTION) {
        const std::vector <std::array<int, 4>> &x_map = ghost_cell_maps[0];
        int i, j, k, i_mirr, idx, idx_mirr;
        for (std::vector < std::array < int, 4 > > ::const_iterator it = x_map.begin(); it != x_map.end();
        ++it)
        {
            i = (*it)[0];
            j = (*it)[1];
            k = (*it)[2];
            i_mirr = (*it)[3];
            if (i < -d_num_var_ghosts[0] || i > interior_dims[0] + d_num_var_ghosts[0] - 1
                || j < -d_num_var_ghosts[1] || j > interior_dims[1] + d_num_var_ghosts[1] - 1
                || k < -d_num_var_ghosts[2] || k > interior_dims[2] + d_num_var_ghosts[2] - 1
                || i_mirr < -d_num_var_ghosts[0] || i_mirr > interior_dims[0] + d_num_var_ghosts[0] - 1)
                continue;
            idx = (i + d_num_var_ghosts[0]) + (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                  (k + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                  (interior_dims[1] + 2 * d_num_var_ghosts[1]);
            idx_mirr = (i_mirr + d_num_var_ghosts[0]) +
                       (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                       (k + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                       (interior_dims[1] + 2 * d_num_var_ghosts[1]);

            if (depth == 1) {
                V[0][idx] = V[0][idx_mirr];
            } else if (depth == 3 && d_condition == WALL_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
            } else if (depth == 3 && d_condition == WALL_NO_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
                V[2][idx] = -V[2][idx_mirr];
            } else if (depth == 5 && d_condition == WALL_SLIP) {
                //Convective Flux in x direction rhou, rhouu +p, rhouv, rhouw, u(E + p)
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
                V[2][idx] = -V[2][idx_mirr];
                V[3][idx] = -V[3][idx_mirr];
                V[4][idx] = -V[4][idx_mirr];
            } else if (depth == 5 && d_condition == WALL_NO_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
                V[3][idx] = V[3][idx_mirr];
                V[4][idx] = -V[4][idx_mirr];
            }


        }
    } else if (d_direction == DIRECTION::Y_DIRECTION) {
        const std::vector <std::array<int, 4>> &y_map = ghost_cell_maps[1];
        int i, j, k, j_mirr, idx, idx_mirr;
        for (std::vector < std::array < int, 4 > > ::const_iterator it = y_map.begin(); it != y_map.end();
        ++it)
        {
            i = (*it)[0];
            j = (*it)[1];
            k = (*it)[2];
            j_mirr = (*it)[3];
            if (i < -d_num_var_ghosts[0] || i > interior_dims[0] + d_num_var_ghosts[0] - 1
                || j < -d_num_var_ghosts[1] || j > interior_dims[1] + d_num_var_ghosts[1] - 1
                || k < -d_num_var_ghosts[2] || k > interior_dims[2] + d_num_var_ghosts[2] - 1
                || j_mirr < -d_num_var_ghosts[1] || j_mirr > interior_dims[1] + d_num_var_ghosts[1] - 1)
                continue;
            idx = (i + d_num_var_ghosts[0]) + (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                  (k + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                  (interior_dims[1] + 2 * d_num_var_ghosts[1]);
            idx_mirr = (i + d_num_var_ghosts[0]) +
                       (j_mirr + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                       (k + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                       (interior_dims[1] + 2 * d_num_var_ghosts[1]);

            if (depth == 1) {
                V[0][idx] = V[0][idx_mirr];
            } else if (depth == 3 && d_condition == WALL_SLIP) {
                V[0][idx] = V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
            } else if (depth == 3 && d_condition == WALL_NO_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
                V[2][idx] = -V[2][idx_mirr];
            } else if (depth == 5 && d_condition == WALL_SLIP) {
                //Convective Flux in x direction rho v, rho vu, rho v v + p, rho v w, v(E + p)
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
                V[3][idx] = -V[3][idx_mirr];
                V[4][idx] = -V[4][idx_mirr];
            } else if (depth == 5 && d_condition == WALL_NO_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
                V[3][idx] = V[3][idx_mirr];
                V[4][idx] = -V[4][idx_mirr];
            }
        }
    } else if (d_direction == DIRECTION::Z_DIRECTION) {
        const std::vector <std::array<int, 4>> &z_map = ghost_cell_maps[2];
        int i, j, k, k_mirr, idx, idx_mirr;
        for (std::vector < std::array < int, 4 > > ::const_iterator it = z_map.begin(); it != z_map.end();
        ++it)
        {
            i = (*it)[0];
            j = (*it)[1];
            k = (*it)[2];
            k_mirr = (*it)[3];
            if (i < -d_num_var_ghosts[0] || i > interior_dims[0] + d_num_var_ghosts[0] - 1
                || j < -d_num_var_ghosts[1] || j > interior_dims[1] + d_num_var_ghosts[1] - 1
                || k < -d_num_var_ghosts[2] || k > interior_dims[2] + d_num_var_ghosts[2] - 1
                || k_mirr < -d_num_var_ghosts[2] || k_mirr > interior_dims[2] + d_num_var_ghosts[2] - 1)
                continue;
            idx = (i + d_num_var_ghosts[0]) + (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                  (k + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                  (interior_dims[1] + 2 * d_num_var_ghosts[1]);
            idx_mirr = (i + d_num_var_ghosts[0]) +
                       (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                       (k_mirr + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                       (interior_dims[1] + 2 * d_num_var_ghosts[1]);

            if (depth == 1) {
                V[0][idx] = V[0][idx_mirr];
            } else if (depth == 3 && d_condition == WALL_SLIP) {
                V[0][idx] = V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
                V[2][idx] = -V[2][idx_mirr];
            } else if (depth == 3 && d_condition == WALL_NO_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
                V[2][idx] = -V[2][idx_mirr];
            } else if (depth == 5 && d_condition == WALL_SLIP) {
                //Convective Flux in x direction rho w, rho wu, rho w v, rho w w + p, w(E + p)
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
                V[2][idx] = -V[2][idx_mirr];
                V[3][idx] = V[3][idx_mirr];
                V[4][idx] = -V[4][idx_mirr];
            } else if (depth == 5 && d_condition == WALL_NO_SLIP) {
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
                V[3][idx] = V[3][idx_mirr];
                V[4][idx] = -V[4][idx_mirr];
            }

        }
    }
}



void mirrorGhostCellDerivative2D(std::vector<boost::shared_ptr<pdat::CellData<double> > >du_x,
                                 const std::vector<std::vector<std::array<int,3> > >& ghost_cell_maps,
                                 const DIRECTION::TYPE d_direction) {

    const hier::IntVector interior_dims = du_x[0]->getBox().numberCells();

    const hier::IntVector d_num_status_ghosts = du_x[0]->getGhostCellWidth();

    const hier::IntVector d_num_var_ghosts = du_x[0]->getGhostCellWidth();

    const int depth = du_x.size();
    std::vector<double *> V;
    V.reserve(depth);
    for (int di = 0; di < depth; di++) {
        V.push_back(du_x[di]->getPointer(0));
    }

    if (d_direction == DIRECTION::X_DIRECTION) {
        const std::vector <std::array<int, 3>> &x_map = ghost_cell_maps[0];
        int i, j, i_mirr, idx, idx_mirr;
        for (std::vector < std::array < int, 3 > > ::const_iterator it = x_map.begin(); it != x_map.end();
        ++it)
        {
            i = (*it)[0];
            j = (*it)[1];
            i_mirr = (*it)[2];
            if (i < -d_num_var_ghosts[0] || i > interior_dims[0] + d_num_var_ghosts[0] - 1
                || j < -d_num_var_ghosts[1] || j > interior_dims[1] + d_num_var_ghosts[1] - 1
                || i_mirr < -d_num_var_ghosts[0] || i_mirr > interior_dims[0] + d_num_var_ghosts[0] - 1)
                continue;

            idx = (i + d_num_var_ghosts[0]) + (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
            idx_mirr = (i_mirr + d_num_var_ghosts[0]) +
                       (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
            if (depth == 2) {   //du/dy, dv/dy
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
            }else if (depth == 5) {   // du/dx, dv/dx, dT/dx,  du/dy, dv/dy
                V[0][idx] = V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
                V[2][idx] = -V[2][idx_mirr];
                V[3][idx] = -V[3][idx_mirr];
                V[4][idx] = -V[4][idx_mirr];
            }


        }
    } else if (d_direction == DIRECTION::Y_DIRECTION) {
        const std::vector <std::array<int, 3>> &y_map = ghost_cell_maps[1];
        int i, j, j_mirr, idx, idx_mirr;
        for (std::vector < std::array < int, 3 > > ::const_iterator it = y_map.begin(); it != y_map.end();
        ++it)
        {
            i = (*it)[0];
            j = (*it)[1];
            j_mirr = (*it)[2];
            if (i < -d_num_var_ghosts[0] || i > interior_dims[0] + d_num_var_ghosts[0] - 1
                || j < -d_num_var_ghosts[1] || j > interior_dims[1] + d_num_var_ghosts[1] - 1
                || j_mirr < -d_num_var_ghosts[1] || j_mirr > interior_dims[1] + d_num_var_ghosts[1] - 1)
                continue;
            idx = (i + d_num_var_ghosts[0]) + (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
            idx_mirr = (i + d_num_var_ghosts[0]) +
                       (j_mirr + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
            if (depth == 2) {   //du/dx, dv/dx
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
            }else if (depth == 5) {   // du/dx, dv/dx,  du/dy, dv/dy, dT/dy
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
                V[3][idx] = V[3][idx_mirr];
                V[4][idx] = -V[4][idx_mirr];
            }
        }
    }
}




void mirrorGhostCellDerivative3D(std::vector<boost::shared_ptr<pdat::CellData<double> > >du_x,
                                 const std::vector<std::vector<std::array<int,4> > >& ghost_cell_maps,
                                 const DIRECTION::TYPE d_direction) {

    const hier::IntVector interior_dims = du_x[0]->getBox().numberCells();

    const hier::IntVector d_num_status_ghosts = du_x[0]->getGhostCellWidth();

    const hier::IntVector d_num_var_ghosts = du_x[0]->getGhostCellWidth();

    const int depth = du_x.size();
    std::vector<double *> V;
    V.reserve(depth);
    for (int di = 0; di < depth; di++) {
        V.push_back(du_x[di]->getPointer(0));
    }

    if (d_direction == DIRECTION::X_DIRECTION) {
        const std::vector <std::array<int, 4>> &x_map = ghost_cell_maps[0];
        int i, j, k, i_mirr, idx, idx_mirr;
        for (std::vector < std::array < int, 4 > > ::const_iterator it = x_map.begin(); it != x_map.end();
        ++it)
        {
            i = (*it)[0];
            j = (*it)[1];
            k = (*it)[2];
            i_mirr = (*it)[3];
            if (i < -d_num_var_ghosts[0] || i > interior_dims[0] + d_num_var_ghosts[0] - 1
                || j < -d_num_var_ghosts[1] || j > interior_dims[1] + d_num_var_ghosts[1] - 1
                || k < -d_num_var_ghosts[2] || k > interior_dims[2] + d_num_var_ghosts[2] - 1
                || i_mirr < -d_num_var_ghosts[0] || i_mirr > interior_dims[0] + d_num_var_ghosts[0] - 1)
                continue;
            idx = (i + d_num_var_ghosts[0]) + (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                  (k + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                  (interior_dims[1] + 2 * d_num_var_ghosts[1]);
            idx_mirr = (i_mirr + d_num_var_ghosts[0]) +
                       (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                       (k + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                       (interior_dims[1] + 2 * d_num_var_ghosts[1]);
            if (depth == 2) {   //du/dy, dv/dy or du/dz, dw/dz
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
            } else if (depth == 8) { // du/dx, dv/dx, dw/x, dT/dx,  du/dy, dv/dy, du/dz, dw/dz
                V[0][idx] = V[0][idx_mirr];
                V[1][idx] = V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
                V[3][idx] = -V[3][idx_mirr];
                V[4][idx] = -V[4][idx_mirr];
                V[5][idx] = -V[5][idx_mirr];
                V[6][idx] = -V[6][idx_mirr];
                V[7][idx] = -V[7][idx_mirr];
            }


        }
    } else if (d_direction == DIRECTION::Y_DIRECTION) {
        const std::vector <std::array<int, 4>> &y_map = ghost_cell_maps[1];
        int i, j, k, j_mirr, idx, idx_mirr;
        for (std::vector < std::array < int, 4 > > ::const_iterator it = y_map.begin(); it != y_map.end();
        ++it)
        {
            i = (*it)[0];
            j = (*it)[1];
            k = (*it)[2];
            j_mirr = (*it)[3];
            if (i < -d_num_var_ghosts[0] || i > interior_dims[0] + d_num_var_ghosts[0] - 1
                || j < -d_num_var_ghosts[1] || j > interior_dims[1] + d_num_var_ghosts[1] - 1
                || k < -d_num_var_ghosts[2] || k > interior_dims[2] + d_num_var_ghosts[2] - 1
                || j_mirr < -d_num_var_ghosts[1] || j_mirr > interior_dims[1] + d_num_var_ghosts[1] - 1)
                continue;
            idx = (i + d_num_var_ghosts[0]) + (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                  (k + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                  (interior_dims[1] + 2 * d_num_var_ghosts[1]);
            idx_mirr = (i + d_num_var_ghosts[0]) +
                       (j_mirr + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                       (k + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                       (interior_dims[1] + 2 * d_num_var_ghosts[1]);

            if (depth == 2) {   //du/dx, dv/dx or dv/dz, dw/dz
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
            } else if (depth == 8) { // du/dx, dv/dx, du/dy, dv/dy, dw/dy, dT/dy, dv/dz, dw/dz
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
                V[2][idx] = V[2][idx_mirr];
                V[3][idx] = V[3][idx_mirr];
                V[4][idx] = V[4][idx_mirr];
                V[5][idx] = -V[5][idx_mirr];
                V[6][idx] = -V[6][idx_mirr];
                V[7][idx] = -V[7][idx_mirr];
            }
        }
    } else if (d_direction == DIRECTION::Z_DIRECTION) {
        const std::vector <std::array<int, 4>> &z_map = ghost_cell_maps[2];
        int i, j, k, k_mirr, idx, idx_mirr;
        for (std::vector < std::array < int, 4 > > ::const_iterator it = z_map.begin(); it != z_map.end();
        ++it)
        {
            i = (*it)[0];
            j = (*it)[1];
            k = (*it)[2];
            k_mirr = (*it)[3];
            if (i < -d_num_var_ghosts[0] || i > interior_dims[0] + d_num_var_ghosts[0] - 1
                || j < -d_num_var_ghosts[1] || j > interior_dims[1] + d_num_var_ghosts[1] - 1
                || k < -d_num_var_ghosts[2] || k > interior_dims[2] + d_num_var_ghosts[2] - 1
                || k_mirr < -d_num_var_ghosts[2] || k_mirr > interior_dims[2] + d_num_var_ghosts[2] - 1)
                continue;
            idx = (i + d_num_var_ghosts[0]) + (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                  (k + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                  (interior_dims[1] + 2 * d_num_var_ghosts[1]);
            idx_mirr = (i + d_num_var_ghosts[0]) +
                       (j + d_num_var_ghosts[1]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) +
                       (k_mirr + d_num_var_ghosts[2]) * (interior_dims[0] + 2 * d_num_var_ghosts[0]) *
                       (interior_dims[1] + 2 * d_num_var_ghosts[1]);

            if (depth == 2) {   //du/dx, dw/dx or dv/dy, dw/dy
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
            } else if (depth == 8) { // du/dx, dw/dx, dv/dy, dw/dy, du/dz, dv/dz, dw/dz, dT/dz
                V[0][idx] = -V[0][idx_mirr];
                V[1][idx] = -V[1][idx_mirr];
                V[2][idx] = -V[2][idx_mirr];
                V[3][idx] = -V[3][idx_mirr];
                V[4][idx] = V[4][idx_mirr];
                V[5][idx] = V[5][idx_mirr];
                V[6][idx] = V[6][idx_mirr];
                V[7][idx] = -V[7][idx_mirr];
            }
        }
    }
}
