/****************************************************************************
* Copyright (c) 2025, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/

#ifndef Polygon_geom_tools_included
#define Polygon_geom_tools_included

#include <cmath>
#include <cassert>

struct Polygon_geom_data
{
  double area_ = 0.;
  double moment_r_ = 0.; // only relevant for axisymmetric geometries
};

template <typename CoordTab, typename IndexFunctor>
inline Polygon_geom_data compute_polygon_geom(const CoordTab& coord,
                                              int geom_dimension,
                                              int nb_vertices,
                                              const IndexFunctor& index_of,
                                              bool compute_moment = false)
{
  Polygon_geom_data data;
  if (nb_vertices < 3)
    return data;

  if (geom_dimension == 2)
    {
      double sum_cross = 0.;
      double sum_r_cross = 0.;
      for (int i = 0; i < nb_vertices; ++i)
        {
          const int next = (i + 1 == nb_vertices) ? 0 : (i + 1);
          const auto si = index_of(i);
          const auto sj = index_of(next);

          const double ri = coord(si, 0);
          const double zi = coord(si, 1);
          const double rj = coord(sj, 0);
          const double zj = coord(sj, 1);

          const double cross = ri * zj - rj * zi;
          sum_cross += cross;
          if (compute_moment)
            sum_r_cross += (ri + rj) * cross;
        }

      data.area_ = 0.5 * std::fabs(sum_cross);
      if (compute_moment)
        data.moment_r_ = sum_r_cross / 6.0;
    }
  else
    {
      assert(geom_dimension == 3);
      const auto s0 = index_of(0);
      const double x0 = coord(s0, 0);
      const double y0 = coord(s0, 1);
      const double z0 = coord(s0, 2);

      double nx = 0.;
      double ny = 0.;
      double nz = 0.;
      for (int i = 1; i < nb_vertices - 1; ++i)
        {
          const auto s1 = index_of(i);
          const auto s2 = index_of(i + 1);

          const double ax = coord(s1, 0) - x0;
          const double ay = coord(s1, 1) - y0;
          const double az = coord(s1, 2) - z0;

          const double bx = coord(s2, 0) - x0;
          const double by = coord(s2, 1) - y0;
          const double bz = coord(s2, 2) - z0;

          nx += ay * bz - az * by;
          ny += az * bx - ax * bz;
          nz += ax * by - ay * bx;
        }

      const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
      data.area_ = 0.5 * norm;
    }

  return data;
}

#endif /* Polygon_geom_tools_included */
