/** \file data_loader.hpp
 *
 *  `data_loader' gathers functions useful to load gen-etics/-omics data
 *  Copyright (C) 2011-2014 Timothée Flutre
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

#ifndef QUANTGEN_DATA_LOADER_HPP
#define QUANTGEN_DATA_LOADER_HPP

#include <map>

#include "utils/utils_io.hpp"

namespace quantgen {
  
  std::map<std::string,std::string> loadTwoColumnFile(
    const std::string & file,
    const int & verbose);

} // namespace quantgen

#endif // QUANTGEN_DATA_LOADER_HPP
