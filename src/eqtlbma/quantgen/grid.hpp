/** \file grid.hpp
 *
 *  `Grid' is a class 
 *  Copyright (C) 2013 Timothee Flutre
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef QUANTGEN_GRID_HPP
#define QUANTGEN_GRID_HPP

#include <cstdlib>

#include <vector>
#include <string>
#include <iostream>

#include "quantgen/utils_io.hpp"

namespace quantgen {

  class Grid {
  public:
    vector<double> phi2s;
    vector<double> oma2s;
    vector<double> phi2s_fix;
    vector<double> oma2s_fix;
    vector<double> phi2s_maxh;
    vector<double> oma2s_maxh;
    Grid();
    Grid(const string & gridFile, const bool & makeFixMaxh, const int & verbose);
    size_t size (void) const { return phi2s.size(); }
  };

} // namespace quantgen

#endif // QUANTGEN_GRID_HPP
