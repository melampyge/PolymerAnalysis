
/* main module to start a particular analysis */

/////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <map>
#include <set>
#include <tuple>
#include <stack>
#include <algorithm>

#include "gen_velocity_grid.hpp"
#include "AnalyseNeighbours.hpp"
#include "AnalysePairCorr.hpp"
#include "AnalyseStaticStruct.hpp"
#include "AnalyseLatticeOrder.hpp"
#include "AnalyseAreaDiff.hpp"
#include "AnalyseDensityHistogram.hpp"
#include "AnalyseInterScattering.hpp"
#include "AnalyseSpatialVelocityCorr.hpp"
#include "AnalyseVorticity.hpp"
#include "AnalyseVelocityStructure.hpp"
#include "AnalyseEnergySpectrum.hpp"
#include "AnalyseGyrationTensor.hpp"

#define pi M_PI

/////////////////////////////////////////////////////////////////////////////////

