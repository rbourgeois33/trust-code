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

#ifndef Ecrire_CGNS_helper_included
#define Ecrire_CGNS_helper_included

#include <TRUST_2_CGNS.h>
#include <Option_CGNS.h>
#include <TRUSTTab.h>
#include <cgns++.h>
#include <sstream>
#include <iomanip>
//#include <regex>

#ifdef HAS_CGNS

#pragma GCC diagnostic push
#if __GNUC__ > 5 || __clang_major__ > 10
#pragma GCC diagnostic ignored "-Wsuggest-override"
#endif
#ifdef MPI_
#include <pcgnslib.h>
#else
#include <cgnslib.h>
#endif
#pragma GCC diagnostic pop

#define CGNS_STR_SIZE 32
#define CGNS_DOUBLE_TYPE Option_CGNS::SINGLE_PRECISION>0?CGNS_ENUMV(RealSingle):CGNS_ENUMV(RealDouble)

enum class TYPE_ECRITURE_CGNS { SEQ , PAR_IN, PAR_OVER };
enum class TYPE_LINK_CGNS { GRID , SOLUTION, FINAL_LINK };
enum class TYPE_MODE_CGNS { WRITE, MODIFY , READ };
enum class TYPE_RUN_CGNS { SEQ, PAR };

inline void TRUST_CGNS_ERROR()
{
#ifdef MPI_
  Process::is_sequential() ? cg_error_exit() : cgp_error_exit();
#else
  Process::is_sequential() ? cg_error_exit() : Process::exit(); /* OpenMP ?? */
#endif
}

struct Ecrire_CGNS_helper
{
  template<TYPE_RUN_CGNS _TYPE_, TYPE_MODE_CGNS _MODE_ = TYPE_MODE_CGNS::WRITE>
  inline void cgns_open_file(const std::string&, int&, const bool print = true);

  template<TYPE_RUN_CGNS _TYPE_>
  inline void cgns_close_file(const std::string&, const int, const bool print = true);

  template<TYPE_ECRITURE_CGNS _TYPE_>
  inline void cgns_write_zone_grid_coord(const int, const int, const int, const char*, const cgsize_t*, int&,
                                         const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, int&, int&, int&);

  template<TYPE_ECRITURE_CGNS _TYPE_>
  inline std::enable_if_t< _TYPE_ != TYPE_ECRITURE_CGNS::SEQ, void>
  cgns_write_grid_coord_data(const int, const int, const int, const int, const int, const int, const int, const cgsize_t, const cgsize_t,
                             const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);

  template<TYPE_ECRITURE_CGNS _TYPE_>
  inline void cgns_sol_write(const int, const int, const int, const int, const double, const std::vector<int>&,
                             const std::string&, std::string&, std::string&, std::string&,
                             bool&, bool&, bool&, int&, int&, int&);

  template<TYPE_ECRITURE_CGNS _TYPE_>
  inline std::enable_if_t< _TYPE_ != TYPE_ECRITURE_CGNS::SEQ, void>
  cgns_field_write(const int, const int, const int, const int, const std::vector<int>&, const std::string&,
                   const int, const int, const int, const char*, int&, int&, int&);

  template<TYPE_ECRITURE_CGNS _TYPE_>
  inline std::enable_if_t<_TYPE_ == TYPE_ECRITURE_CGNS::SEQ, void>
  cgns_field_write_data(const int, const int, const int, const std::vector<int>&, const std::string&,
                        const int, const int, const int, const int,
                        const char * , const DoubleTab& , int& , int& , int&);

  template<TYPE_ECRITURE_CGNS _TYPE_>
  inline std::enable_if_t<_TYPE_ != TYPE_ECRITURE_CGNS::SEQ, void>
  cgns_field_write_data(const int, const int, const int, const std::vector<int>&, const std::string&,
                        const int, const int, const int, const int, const int, const int, const int,
                        const cgsize_t, const cgsize_t, const DoubleTab&);

  template<TYPE_ECRITURE_CGNS _TYPE_>
  inline void cgns_write_iters(const bool, const int, const int , const int, const int, const std::vector<int>&,
                               const std::string&, const std::string&, const std::string&, const std::string&, const std::vector<double>&);

  std::string convert_double_to_string(const double t)
  {
    // On fait comme dans les latas !
    char str_temps[100] = "0.0";
    if (t >= 0.)
      snprintf(str_temps, 100, "%.10f", t);

    return std::string(str_temps);

// Sinon je laisse cette belle methode au cas ou elle pourrait servir un jour !
//    /* KEEP GOOD PRECISION */
//    std::stringstream ss;
//    ss << std::fixed << std::setprecision(15) << t;
//
//    std::string ret_str = ss.str();
//
//    /* REMOVE USELESS 0s */
//    ret_str = std::regex_replace(ret_str, std::regex("([0-9]*\\.[0-9]*?[1-9])0+$"), "$1"); // remove useless 0s
//    ret_str = std::regex_replace(ret_str, std::regex("\\.0+$"), ""); // remove the final .0 if integer
//
//    return ret_str;
  }
};

#include <Ecrire_CGNS_helper.tpp>

#endif /* HAS_CGNS */

#endif /* Ecrire_CGNS_helper_included */
