#include "util/basic_geometry/PorousWall.hpp"

#include <algorithm>
#include <math.h>       /* fabs */
/*
 * This function decide point (x, y, z) is in the wall or not
 * return 0 if the point is in the wall
 * return 1 if the point is outside
 * the mesh size is 1/(32N).
 *
 */
int isOutside(int dim, double x, double y, double z)
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
computeCellStatus(boost::shared_ptr<pdat::CellData<int> > &cell_status,
                  const double* x_lo,  const double* dx)
{
    const tbox::Dimension d_dim = cell_status->getDim();
    const hier::IntVector  interior_dims = cell_status->getBox().numberCells();

    const hier::IntVector d_num_conv_ghosts = cell_status->getGhostCellWidth();

    int* cell_status_data = cell_status->getPointer(0);
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
                cell_status_data[idx] = isOutside(2, x, y);
            }


    } else if (d_dim == tbox::Dimension(3)) {}
}





