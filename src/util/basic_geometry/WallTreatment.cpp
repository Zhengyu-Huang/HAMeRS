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
int isOutsidePorousWall2(int dim, double x, double y, double z)
{
    /*
     * The porous wall is centered at x0.
     * The computational domain is at least is 8 by 8 by 8
     * There are 1 holes in y and z directions, of size 2,
     *           2 solid in y and z directions, each of size 3;
     * 3 2 3
     */

    double x0 = 0.0, r_hole = 1./8;
    if(fabs(x - x0) > r_hole) return 1;
    if(dim == 2){
        double yy = y/r_hole;
        if(yy < 3 || (yy > 5))
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
int isOutsidePorousWall3(int dim, double x, double y, double z)
{
    /*
     * The porous wall is centered at x0.
     * The computational domain is at least is 8 by 8 by 8
     * A solid wall at x = 0
     */
    double x0 = 0.0, r_hole = 1./8;
    if(fabs(x - x0) > r_hole) return 1;
    else return 0;
}


void
initializeCellStatus(hier::Patch& patch,
                     boost::shared_ptr<pdat::CellData<double> > &cell_status)
{
    // Get the grid spacing.
    const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
            BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
                    patch.getPatchGeometry()));


    const double* const dx = patch_geom->getDx();
    const double* const x_lo = patch_geom->getXLower();

    // Get the dimensions of box that covers the interior of Patch.
    const tbox::Dimension d_dim = patch.getDim();

    hier::Box patch_box = patch.getBox();
    const hier::IntVector patch_dims = patch_box.numberCells();




    const hier::IntVector  interior_dims = cell_status->getBox().numberCells();

    const hier::IntVector d_num_conv_ghosts = cell_status->getGhostCellWidth();

    double* cell_status_data = cell_status->getPointer(0);
    if (d_dim == tbox::Dimension(1)) {}
    else if (d_dim == tbox::Dimension(2)) {
/*
 * Compute the status of cells, including ghost cells
 */

            for (int j = 0; j < patch_dims[1]; j++) {
                for (int i = 0; i < patch_dims[0]; i++) {
                    const double x = x_lo[0] + (i + 0.5) * dx[0],
                                 y = x_lo[1] + (j + 0.5) * dx[1];
                    const int idx = i + j * patch_dims[0];
                    cell_status_data[idx] = isOutsidePorousWall2(2, x, y);
                }
            }


    } else if (d_dim == tbox::Dimension(3)) {}
}


void
mirrorGhostCell(boost::shared_ptr<pdat::CellData<double> > &variables,
       const boost::shared_ptr<pdat::CellData<double> > &cell_status,
       const DIRECTION::TYPE d_direction,
       const WALL_TREATMENT_CONDITION d_condition)
{

    const tbox::Dimension d_dim = cell_status->getDim();
    const hier::IntVector  interior_dims = cell_status->getBox().numberCells();

    const hier::IntVector d_num_status_ghosts = cell_status->getGhostCellWidth();

    const hier::IntVector d_num_var_ghosts = variables->getGhostCellWidth();



    double* cell_status_data = cell_status->getPointer(0);

    const int depth = variables->getDepth();
    std::vector<double*> V;
    V.reserve(depth);
    for (int di = 0; di < depth; di++) {
        V.push_back(variables->getPointer(di));
    }

    //std::cout << " mirrorGhostCell depth is " << depth << std::endl;



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
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) + (j - d_num_var_ghosts[1] + d_num_status_ghosts[1])* (interior_dims[0] + 2 * d_num_status_ghosts[0]);
                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] && (i - 2 * ghost_count + 1) >=0 ) {
                        const int idx_mirr =
                                (i - 2 * ghost_count + 1) + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);

                        if(depth == 1)
                        {
                            V[0][idx] = V[0][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_SLIP)
                        {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_NO_SLIP)
                        {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                            if(fabs(V[0][idx_mirr]) + fabs(V[1][idx_mirr]) > 1e-8) {
                                std::cout << "X1 i is " << i << " j is " << j << " ghost is " << ghost_count
                                          << "i_m " << i - 2 * ghost_count + 1 << " j_m " << j
                                          << " u " << V[0][idx_mirr] << " v " << V[1][idx_mirr] << std::endl;
                                exit(1);
                            }
                        } else if (depth == 4 && d_condition == WALL_SLIP)
                        {   //Convective Flux in x direction rhou, rhouu +p, rhouv, u(E + p)
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = -V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        } else if (depth == 4 && d_condition == WALL_NO_SLIP)
                        {
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
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) + (j - d_num_var_ghosts[1] + d_num_status_ghosts[1])* (interior_dims[0] + 2 * d_num_status_ghosts[0]);
                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] && (i + 2 * ghost_count - 1) <= interior_dims[0] + 2 * d_num_var_ghosts[0] - 1) {
                        const int idx_mirr =
                                (i + 2 * ghost_count - 1) + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                        if(depth == 1) {
                            V[0][idx] = V[0][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_SLIP)
                        {   //velocity u, v
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                        } else if (depth == 2  && d_condition == WALL_NO_SLIP) {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                            if(fabs(V[0][idx_mirr]) + fabs(V[1][idx_mirr]) > 1e-8) {
                                std::cout << "X2 i is " << i << " j is " << j << " ghost is " << ghost_count
                                          << " u " << V[0][idx_mirr] << " v " << V[1][idx_mirr] << std::endl;
                                exit(1);
                            }
                        } else if (depth == 4 && d_condition == WALL_SLIP)
                        {   //Convective Flux in x direction rhou, rhouu +p, rhouv, u(E + p)
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = -V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        } else if (depth == 4 && d_condition == WALL_NO_SLIP)
                        {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        }

                    }

                }
            }
        } else if (d_dim == tbox::Dimension(3)) {}
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    else if (d_direction == DIRECTION::Y_DIRECTION) {
        if (d_dim == tbox::Dimension(1)) {}
        else if (d_dim == tbox::Dimension(2)) {
            //loop x direction
            for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {

                ghost_count = std::numeric_limits<int>::min();
                //loop y direction
                for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) + (j - d_num_var_ghosts[1] + d_num_status_ghosts[1])* (interior_dims[0] + 2 * d_num_status_ghosts[0]);
                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] && j - 2 * ghost_count + 1 >=0) {

                        const int idx_mirr =
                                i + (j - 2 * ghost_count + 1) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                        if(depth == 1) {
                            V[0][idx] = V[0][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_SLIP) {
                            V[0][idx] = V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                        }
                        else if (depth == 2 && d_condition == WALL_NO_SLIP) {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                            if(fabs(V[0][idx_mirr]) + fabs(V[1][idx_mirr]) > 1e-8)
                                std::cout << "Y1 i is " << i << " j is " << j << " ghost is " << ghost_count << std::endl;
                        }else if (depth == 4 && d_condition == WALL_SLIP)
                        {   //Convective Flux in x direction rho v, rho u v, rho v v + p, v(E + p)
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        } else if (depth == 4 && d_condition == WALL_NO_SLIP)
                        {
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
                    const int idx_status = (i - d_num_var_ghosts[0] + d_num_status_ghosts[0]) + (j - d_num_var_ghosts[1] + d_num_status_ghosts[1])* (interior_dims[0] + 2 * d_num_status_ghosts[0]);
                    if (cell_status_data[idx_status] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_var_ghosts[d_direction] && j + 2 * ghost_count - 1 <= interior_dims[1] + 2 * d_num_var_ghosts[1] - 1) {

                        const int idx_mirr =
                                i + (j + 2 * ghost_count - 1) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
                        if(depth == 1) {
                            V[0][idx] = V[0][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_SLIP)
                        {
                            V[0][idx] = V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                        } else if (depth == 2 && d_condition == WALL_NO_SLIP)
                        {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                            if(fabs(V[0][idx_mirr]) + fabs(V[1][idx_mirr]) > 1e-8)
                                std::cout << "Y2 i is " << i << " j is " << j << " ghost is " << ghost_count << std::endl;
                        } else if (depth == 4 && d_condition == WALL_SLIP)
                        {   //Convective Flux in x direction rho v, rho u v, rho v v + p, v(E + p)
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = -V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        } else if (depth == 4 && d_condition == WALL_NO_SLIP)
                        {
                            V[0][idx] = -V[0][idx_mirr];
                            V[1][idx] = V[1][idx_mirr];
                            V[2][idx] = V[2][idx_mirr];
                            V[3][idx] = -V[3][idx_mirr];
                        }

                    }

                }
            }
        } else if (d_dim == tbox::Dimension(3)) {}
    }


//    std::cout << " dim: " << interior_dims[0]    << " " << interior_dims[1] << " ghost " <<
//                             d_num_var_ghosts[0] << " " << d_num_var_ghosts[1] << std::endl;
//    for (int i = 0; i < interior_dims[0] + 2 * d_num_var_ghosts[0]; i++) {
//        //loop x direction
//        for (int j = 0; j < interior_dims[1] + 2 * d_num_var_ghosts[1]; j++) {
//            //loop y direction
//            int idx = i + j * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
//            int idx_ = i + (interior_dims[1] + 2 * d_num_var_ghosts[1]- j - 1) * (interior_dims[0] + 2 * d_num_var_ghosts[0]);
//            for(int d = 0; d < depth; d++) {
//                if (fabs(V[0][idx] - V[0][idx_]) > 1e-10) {
//                    std::cout << "error (i,j) " << i <<" " << j << " d " << d << " " << (interior_dims[1] + 2 * d_num_var_ghosts[1]- j - 1)
//                              << " d_drection "<< d_direction << " " << " value " << V[0][idx]  << " " << V[0][idx_]<< std::endl;
//                    exit(1);
//                }
//            }
//        }
//    }



}



void
populateGhostCellsHelper(std::vector<boost::shared_ptr<pdat::CellData<double > > > &conservative_variables,
                         const boost::shared_ptr<pdat::CellData<double> > &cell_status,
                         const WALL_TREATMENT_CONDITION d_condition) {
    const tbox::Dimension d_dim = cell_status->getDim();
    const hier::IntVector interior_dims = cell_status->getBox().numberCells();
    const hier::IntVector d_num_conv_ghosts = cell_status->getGhostCellWidth();
    double *cell_status_data = cell_status->getPointer(0);

    int d_num_eqn = d_dim.getValue() + 2;
    std::vector<double *> Q;
    Q.reserve(d_num_eqn);
    for(int vi = 0; vi < static_cast<int>(conservative_variables.size()); vi++) {
        int depth = conservative_variables[vi]->getDepth();
        for (int di = 0; di<depth; di++)
            Q.push_back(conservative_variables[vi]->getPointer(di));
    }


    int ghost_count;//Set to -inf to start get rid of these block that never use in current computation
    if (d_dim == tbox::Dimension(1)) {
        int num_cells = (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
        double *weight = new double[num_cells];
    } else if (d_dim == tbox::Dimension(2)) {
        int num_cells = (interior_dims[0] + 2 * d_num_conv_ghosts[0]) *
                        (interior_dims[1] + 2 * d_num_conv_ghosts[1]);
        double *weight = new double[num_cells];
        for (
                int wi = 0;
                wi < num_cells;
                wi++)
            weight[wi] = -1.;

/*
 * d_direction == DIRECTION::X_DIRECTION)
 */

//loop y direction
        for (int j = 0; j < interior_dims[1] + 2 * d_num_conv_ghosts[1]; j++) {
            ghost_count = std::numeric_limits<int>::min();
//loop x direction
            for (int i = 0; i < interior_dims[0] + 2 * d_num_conv_ghosts[0]; i++) {
                const int idx = i + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
                if (cell_status_data[idx] < 0.5) {
                    ghost_count++;
                } else {
                    ghost_count = 0;
                }
                if (ghost_count >= 1 && ghost_count <= d_num_conv_ghosts[0]) {
                    const int idx_mirr =
                            (i - 2 * ghost_count + 1) + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);

// constant population
// slip
// no-slip
                    if (weight[idx] < 0.) {
                        weight[idx] = 0.0;
                        for (
                                int di = 0;
                                di < d_num_eqn;
                                di++)
                            Q[di][idx] = 0.0;
                    }

                    weight[idx] += 1.0;

                    if (d_condition == WALL_NO_SLIP) {
                        Q[0][idx] += Q[0][idx_mirr];
                        Q[1][idx] += -Q[1][idx_mirr];
                        Q[2][idx] += -Q[2][idx_mirr];
                        Q[3][idx] += Q[3][idx_mirr];
                    }

                    if (d_condition == WALL_NO_SLIP) {
                        Q[0][idx] += Q[0][idx_mirr];
                        Q[1][idx] += -Q[1][idx_mirr];
                        Q[2][idx] += Q[2][idx_mirr];
                        Q[3][idx] += Q[3][idx_mirr];
                    }

                }
            }

            ghost_count = std::numeric_limits<int>::min();
//loop x direction inversely
            for (
                    int i = interior_dims[0] + 2 * d_num_conv_ghosts[0] - 1;
                    i >= 0; i--) {

                const int idx = i + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
                if (cell_status_data[idx] < 0.5) {
                    ghost_count++;
                } else {
                    ghost_count = 0;
                }
                if (ghost_count >= 1 && ghost_count <= d_num_conv_ghosts[0]) {
                    const int idx_mirr =
                            (i + 2 * ghost_count - 1) + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);


                    if (weight[idx] < 0.) {
                        weight[idx] = 0.0;
                        for (
                                int di = 0;
                                di < d_num_eqn;
                                di++)
                            Q[di][idx] = 0.0;
                    }

                    weight[idx] += 1.0;

                    if (d_condition == WALL_NO_SLIP) {
                        Q[0][idx] += Q[0][idx_mirr];
                        Q[1][idx] += -Q[1][idx_mirr];
                        Q[2][idx] += -Q[2][idx_mirr];
                        Q[3][idx] += Q[3][idx_mirr];
                    }

                    if (d_condition == WALL_NO_SLIP) {
                        Q[0][idx] += Q[0][idx_mirr];
                        Q[1][idx] += -Q[1][idx_mirr];
                        Q[2][idx] += Q[2][idx_mirr];
                        Q[3][idx] += Q[3][idx_mirr];
                    }

                }

            }
        }


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//(d_direction == DIRECTION::Y_DIRECTION) {
//loop x direction
        for (
                int i = 0;
                i < interior_dims[0] + 2 * d_num_conv_ghosts[0]; i++) {

            ghost_count = std::numeric_limits<int>::min();
//loop y direction
            for (
                    int j = 0;
                    j < interior_dims[1] + 2 * d_num_conv_ghosts[1]; j++) {

                const int idx = i + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
                if (cell_status_data[idx] < 0.5) {
                    ghost_count++;
                } else {
                    ghost_count = 0;
                }
                if (ghost_count >= 1 && ghost_count <= d_num_conv_ghosts[1]) {
                    const int idx_mirr =
                            i + (j - 2 * ghost_count + 1) * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);


                    if (weight[idx] < 0.) {
                        weight[idx] = 0.0;
                        for (
                                int di = 0;
                                di < d_num_eqn;
                                di++)
                            Q[di][idx] = 0.0;
                    }

                    weight[idx] += 1.0;

                    if (d_condition == WALL_NO_SLIP) {
                        Q[0][idx] += Q[0][idx_mirr];
                        Q[1][idx] += -Q[1][idx_mirr];
                        Q[2][idx] += -Q[2][idx_mirr];
                        Q[3][idx] += Q[3][idx_mirr];
                    }

                    if (d_condition == WALL_NO_SLIP) {
                        Q[0][idx] += Q[0][idx_mirr];
                        Q[1][idx] += Q[1][idx_mirr];
                        Q[2][idx] += -Q[2][idx_mirr];
                        Q[3][idx] += Q[3][idx_mirr];
                    }


                }
            }
            ghost_count = std::numeric_limits<int>::min();
//loop y direction inversely
            for (
                    int j = interior_dims[1] + 2 * d_num_conv_ghosts[1] - 1;
                    j >= 0; j--) {

                const int idx = i + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
                if (cell_status_data[idx] < 0.5) {
                    ghost_count++;
                } else {
                    ghost_count = 0;
                }
                if (ghost_count >= 1 && ghost_count <= d_num_conv_ghosts[1]) {
                    const int idx_mirr =
                            i + (j + 2 * ghost_count - 1) * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);


                    if (weight[idx] < 0.) {
                        weight[idx] = 0.0;
                        for (
                                int di = 0;
                                di < d_num_eqn;
                                di++)
                            Q[di][idx] = 0.0;
                    }

                    weight[idx] += 1.0;

                    if (d_condition == WALL_NO_SLIP) {
                        Q[0][idx] += Q[0][idx_mirr];
                        Q[1][idx] += -Q[1][idx_mirr];
                        Q[2][idx] += -Q[2][idx_mirr];
                        Q[3][idx] += Q[3][idx_mirr];
                    }

                    if (d_condition == WALL_NO_SLIP) {
                        Q[0][idx] += Q[0][idx_mirr];
                        Q[1][idx] += Q[1][idx_mirr];
                        Q[2][idx] += -Q[2][idx_mirr];
                        Q[3][idx] += Q[3][idx_mirr];
                    }

                }

            }
        }
        delete[]
                weight;
    } else if (d_dim == tbox::Dimension(3)) {
        int num_cells = (interior_dims[0] + 2 * d_num_conv_ghosts[0]) *
                        (interior_dims[1] + 2 * d_num_conv_ghosts[1]) *
                        (interior_dims[2] + 2 * d_num_conv_ghosts[2]);
        double *weight = new double[num_cells];
    }


}

