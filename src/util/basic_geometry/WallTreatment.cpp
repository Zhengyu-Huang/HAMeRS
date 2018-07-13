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
    * The computational domain is 32 by 32 by 32
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


void
computeCellStatus(boost::shared_ptr<pdat::CellData<double> > &cell_status,
                  const double* x_lo,  const double* dx)
{
    const tbox::Dimension d_dim = cell_status->getDim();
    const hier::IntVector  interior_dims = cell_status->getBox().numberCells();

    const hier::IntVector d_num_conv_ghosts = cell_status->getGhostCellWidth();

    double* cell_status_data = cell_status->getPointer(0);
    if (d_dim == tbox::Dimension(1)) {}
    else if (d_dim == tbox::Dimension(2)) {
/*
 * Compute the status of cells, including ghost cells
 */

        for (int i = 0; i < interior_dims[0] + 2 * d_num_conv_ghosts[0]; i++)
            for (int j = 0; j < interior_dims[1] + 2 * d_num_conv_ghosts[1]; j++) {
                const double x = x_lo[0] + (i + 0.5 - d_num_conv_ghosts[0]) * dx[0],
                             y = x_lo[1] + (j + 0.5 - d_num_conv_ghosts[1]) * dx[1];
                const int idx = i + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
                //std::cout << "(i,j) " << i << " " << j << " interior dim " << interior_dims
                //          << " ghosts " << d_num_conv_ghosts << " idx " << idx << std::endl;
                cell_status_data[idx] = isOutsidePorousWall1(2, x, y);
            }


    } else if (d_dim == tbox::Dimension(3)) {}
}


void
mirrorGhostCell(boost::shared_ptr<pdat::CellData<double> > &variables,
       const boost::shared_ptr<pdat::CellData<double> > &cell_status,
       const DIRECTION::TYPE d_direction)
{
    const tbox::Dimension d_dim = cell_status->getDim();
    const hier::IntVector  interior_dims = cell_status->getBox().numberCells();

    const hier::IntVector d_num_conv_ghosts = cell_status->getGhostCellWidth();

    double* cell_status_data = cell_status->getPointer(0);

    const int depth = variables->getDepth();
    std::vector<double*> V;
    V.reserve(depth);
    for (int di = 0; di < depth; di++) {
        V.push_back(variables->getPointer(di));
    }



    int ghost_count;//Set to -inf to start get rid of these block that never use in current computation
    if (d_direction == DIRECTION::X_DIRECTION) {
        if (d_dim == tbox::Dimension(1)) {}
        else if (d_dim == tbox::Dimension(2)) {
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
                    if (ghost_count >= 1 && ghost_count <= d_num_conv_ghosts[d_direction]) {
                        const int idx_mirr =
                                (i - 2 * ghost_count + 1) + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);


                        for (int di = 0; di < depth; di++)
                            V[di][idx] = V[di][idx_mirr];
                        if (depth > 1)
                            V[d_direction][idx] = -V[d_direction][idx_mirr];

                    }
                }

                ghost_count = std::numeric_limits<int>::min();
                //loop x direction inversely
                for (int i = interior_dims[0] + 2 * d_num_conv_ghosts[0] - 1; i >= 0; i--) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
                    if (cell_status_data[idx] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_conv_ghosts[d_direction]) {
                        const int idx_mirr =
                                (i + 2 * ghost_count - 1) + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
                        for (int di = 0; di < depth; di++)
                            V[di][idx] = V[di][idx_mirr];
                        if (depth > 1)
                            V[d_direction][idx] = -V[d_direction][idx_mirr];

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
            for (int i = 0; i < interior_dims[0] + 2 * d_num_conv_ghosts[0]; i++) {

                ghost_count = std::numeric_limits<int>::min();
                //loop y direction
                for (int j = 0; j < interior_dims[1] + 2 * d_num_conv_ghosts[1]; j++) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
                    if (cell_status_data[idx] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_conv_ghosts[d_direction]) {
                        const int idx_mirr =
                                i + (j - 2 * ghost_count + 1) * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
                        for (int di = 0; di < depth; di++)
                            V[di][idx] = V[di][idx_mirr];
                        if (depth > 1)
                            V[d_direction][idx] = -V[d_direction][idx_mirr];
                    }
                }
                ghost_count = std::numeric_limits<int>::min();
                //loop y direction inversely
                for (int j = interior_dims[1] + 2 * d_num_conv_ghosts[1] - 1; j >= 0; j--) {

                    const int idx = i + j * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
                    if (cell_status_data[idx] < 0.5) {
                        ghost_count++;
                    } else {
                        ghost_count = 0;
                    }
                    if (ghost_count >= 1 && ghost_count <= d_num_conv_ghosts[d_direction]) {
                        const int idx_mirr =
                                i + (j + 2 * ghost_count - 1) * (interior_dims[0] + 2 * d_num_conv_ghosts[0]);
                        for (int di = 0; di < depth; di++)
                            V[di][idx] = V[di][idx_mirr];
                        if (depth > 1)
                            V[d_direction][idx] = -V[d_direction][idx_mirr];

                    }

                }
            }
        } else if (d_dim == tbox::Dimension(3)) {}
    }
}





