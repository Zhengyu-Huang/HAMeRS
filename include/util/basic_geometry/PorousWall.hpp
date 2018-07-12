#ifndef POROUS_WALL_CONDITIONS_HPP
#define POROUS_WALL_CONDITIONS_HPP

#include "HAMeRS_config.hpp"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Dimension.h"


#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace SAMRAI;

int
isOutside(int dim, double x, double y, double z = 0.0);



void
computeCellStatus(boost::shared_ptr<pdat::CellData<double> > &cell_status,
                  const double* x_lo,  const double* dx);
#endif /* POROUS_WALL_CONDITIONS_HPP */
