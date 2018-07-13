#ifndef WALL_TREATMENT_CONDITIONS_HPP
#define WALL_TREATMENT_CONDITIONS_HPP

#include "HAMeRS_config.hpp"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"


#include "util/Directions.hpp"

#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace SAMRAI;

int
isOutsidePorousWall1(int dim, double x, double y, double z = 0.0);



void
computeCellStatus(boost::shared_ptr<pdat::CellData<double> > &cell_status,
                  const double* x_lo,  const double* dx);

void
mirrorGhostCell(boost::shared_ptr<pdat::CellData<double> > &variables,
                const boost::shared_ptr<pdat::CellData<double> > &cell_status,
                const DIRECTION::TYPE d_direction);
#endif /* WALL_TREATMENT_CONDITIONS_HPP */
