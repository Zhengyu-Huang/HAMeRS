#ifndef WALL_TREATMENT_CONDITIONS_HPP
#define WALL_TREATMENT_CONDITIONS_HPP

#include "HAMeRS_config.hpp"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include "util/Directions.hpp"

#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace SAMRAI;


/*
* WALL TREATMENT CONDITIONS
*   SYMMETRY    - s_g = s_i
*   NO_SLIP_VEL - v_g = -v_i
*   SLIP_VEL    - v_g n = -v_i n
*/
enum WALL_TREATMENT_CONDITION {
    WALL_NO_SLIP = 0,
    WALL_SLIP = 1
};

int isOutsidePorousWall1(int dim, double x, double y, double z = 0.0);
int isOutsidePorousWall2X(int dim, double x, double y, double z = 0.0);
int isOutsidePorousWall2Y(int dim, double x, double y, double z = 0.0);

int isOutsidePorousWall3X(int dim, double x, double y, double z = 0.0);
int isOutsidePorousWall3Y(int dim, double x, double y, double z = 0.0);
int isOutsidePorousWall3Z(int dim, double x, double y, double z);

int isOutsidePorousWall4X(int dim, double x, double y, double z);
int isOutsidePorousWall4X_2(int dim, double x, double y, double z);
int isOutsidePorousWall4X_3(int dim, double x, double y, double z);
void
        initializeCellStatus(hier::Patch& patch,
                             boost::shared_ptr <pdat::CellData<double>> &cell_status);

void
        mirrorGhostCell(boost::shared_ptr <pdat::CellData<double>> &variables,
                        const boost::shared_ptr <pdat::CellData<double>> &cell_status,
                        const DIRECTION::TYPE d_direction,
                        const WALL_TREATMENT_CONDITION d_condition = WALL_SLIP);
void
        buildGhostCellMap2D(const boost::shared_ptr<pdat::CellData<double> > &cell_status,
                            std::vector<std::vector<std::array<int,3> > >& ghost_cell_maps);


void
        buildGhostCellMap3D(const boost::shared_ptr<pdat::CellData<double> > &cell_status,
                            std::vector<std::vector<std::array<int,4> > >& ghost_cell_maps);


void
        mirrorGhostCellDerivative(std::vector<boost::shared_ptr<pdat::CellData<double> > >  du_x,
                                  const boost::shared_ptr<pdat::CellData<double> > &cell_status,
                                  const DIRECTION::TYPE d_direction);


void
        mirrorGhostCell2D(boost::shared_ptr<pdat::CellData<double> > &variables,
                          const std::vector<std::vector<std::array<int,3> > > & ghost_cell_maps,
                          const DIRECTION::TYPE d_direction,
                          const WALL_TREATMENT_CONDITION d_condition = WALL_SLIP);

void
        mirrorGhostCell3D(boost::shared_ptr<pdat::CellData<double> > &variables,
                          const std::vector<std::vector<std::array<int,4> > > & ghost_cell_maps,
                          const DIRECTION::TYPE d_direction,
                          const WALL_TREATMENT_CONDITION d_condition = WALL_SLIP);


void
        mirrorGhostCellDerivative2D(std::vector<boost::shared_ptr<pdat::CellData<double> > >du_x,
                                    const std::vector<std::vector<std::array<int,3> > >& ghost_cell_maps,
                                    const DIRECTION::TYPE d_direction);


void
        mirrorGhostCellDerivative3D(std::vector<boost::shared_ptr<pdat::CellData<double> > >du_x,
                                    const std::vector<std::vector<std::array<int,4> > >& ghost_cell_maps,
                                    const DIRECTION::TYPE d_direction);
#endif /* WALL_TREATMENT_CONDITIONS_HPP */
