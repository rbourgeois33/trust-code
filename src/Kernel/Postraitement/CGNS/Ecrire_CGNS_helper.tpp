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

#ifndef Ecrire_CGNS_helper_tpp_included
#define Ecrire_CGNS_helper_tpp_included

#include <TRUSTTrav.h>

template<TYPE_RUN_CGNS _TYPE_, TYPE_MODE_CGNS _MODE_>
inline void Ecrire_CGNS_helper::cgns_open_file(const std::string& fn, int& fileId, const bool is_print)
{
  constexpr bool is_SEQ = (_TYPE_ == TYPE_RUN_CGNS::SEQ), is_WRITE = (_MODE_ == TYPE_MODE_CGNS::WRITE), is_MODIFY = (_MODE_ == TYPE_MODE_CGNS::MODIFY);
  if (is_SEQ)
    {
      if (cg_open(fn.c_str(), is_WRITE ? CG_MODE_WRITE : ( is_MODIFY ? CG_MODE_MODIFY : CG_MODE_READ), &fileId) != CG_OK)
        {
          Cerr << "Error Ecrire_CGNS_helper::cgns_open_file : cg_open !" << finl;
          TRUST_CGNS_ERROR();
        }

      if (is_print)
        Cerr << "**** CGNS file " << fn << " opened !" << finl;
    }
  else
    {
#ifdef MPI_
      if (cgp_open(fn.c_str(), is_WRITE ? CG_MODE_WRITE : ( is_MODIFY ? CG_MODE_MODIFY : CG_MODE_READ), &fileId) != CG_OK)
        {
          Cerr << "Error Ecrire_CGNS_helper::cgns_open_file : cgp_open !" << finl;
          TRUST_CGNS_ERROR();
        }

      if (is_print)
        Cerr << "**** Parallel CGNS file " << fn << " opened !" << finl;
#else
      Cerr << "Parallel CGNS files need MPI installed ... " << finl;
      TRUST_CGNS_ERROR();
#endif
    }
}

template<TYPE_RUN_CGNS _TYPE_>
inline void Ecrire_CGNS_helper::cgns_close_file(const std::string& fn, const int fileId, const bool is_print)
{
  constexpr bool is_SEQ = (_TYPE_ == TYPE_RUN_CGNS::SEQ);
  if (is_SEQ)
    {
      if (cg_close(fileId) != CG_OK)
        {
          Cerr << "Error Ecrire_CGNS_helper::cgns_close_file : cg_close !" << finl;
          TRUST_CGNS_ERROR();
        }

      if (is_print)
        Cerr << "**** CGNS file " << fn << " closed !" << finl;
    }
  else
    {
#ifdef MPI_
      if (cgp_close(fileId) != CG_OK)
        {
          Cerr << "Error Ecrire_CGNS_helper::cgns_close_file : cgp_close !" << finl;
          TRUST_CGNS_ERROR();
        }

      if (is_print)
        Cerr << "**** Parallel CGNS file " << fn << " closed !" << finl;
#else
      Cerr << "Parallel CGNS files need MPI installed ... " << finl;
      TRUST_CGNS_ERROR();
#endif
    }
}

template<TYPE_ECRITURE_CGNS _TYPE_>
inline void Ecrire_CGNS_helper::cgns_write_zone_grid_coord(const int icelldim, const int fileId, const int baseId, const char *zonename, const cgsize_t *isize, int& zoneId,
                                                           const std::vector<double>& xCoords, const std::vector<double>& yCoords, const std::vector<double>& zCoords,
                                                           int& coordsIdx, int& coordsIdy, int& coordsIdz)
{
  constexpr bool is_SEQ = (_TYPE_ == TYPE_ECRITURE_CGNS::SEQ);

  if (cg_zone_write(fileId, baseId, zonename, isize, CGNS_ENUMV(Unstructured), &zoneId) != CG_OK)
    Cerr << "Error Ecrire_CGNS_helper::cgns_write_zone_grid_coord : cg_zone_write !" << finl, TRUST_CGNS_ERROR();

  if (is_SEQ)
    {
      if (cg_coord_write(fileId, baseId, zoneId, CGNS_DOUBLE_TYPE, "CoordinateX", xCoords.data(), &coordsIdx) != CG_OK)
        Cerr << "Error Ecrire_CGNS_helper::cgns_write_zone_grid_coord : cg_coord_write - X !" << finl, TRUST_CGNS_ERROR();

      if (cg_coord_write(fileId, baseId, zoneId, CGNS_DOUBLE_TYPE, "CoordinateY", yCoords.data(), &coordsIdy) != CG_OK)
        Cerr << "Error Ecrire_CGNS_helper::cgns_write_zone_grid_coord : cg_coord_write - Y !" << finl, TRUST_CGNS_ERROR();

      if (Objet_U::dimension > 2)
        if (cg_coord_write(fileId, baseId, zoneId, CGNS_DOUBLE_TYPE, "CoordinateZ", zCoords.data(), &coordsIdz) != CG_OK)
          Cerr << "Error Ecrire_CGNS_helper::cgns_write_zone_grid_coord : cg_coord_write - Z !" << finl, TRUST_CGNS_ERROR();
    }
  else
    {
#ifdef MPI_
      int gridId = -123;
      if (cg_grid_write(fileId, baseId, zoneId, "GridCoordinates", &gridId) != CG_OK)
        Cerr << "Error Ecrire_CGNS_helper::cgns_write_zone_grid_coord : cg_grid_write !" << finl, TRUST_CGNS_ERROR();

      if (cgp_coord_write(fileId, baseId, zoneId, CGNS_DOUBLE_TYPE, "CoordinateX", &coordsIdx) != CG_OK)
        Cerr << "Error Ecrire_CGNS_helper::cgns_write_zone_grid_coord : cgp_coord_write - X !" << finl, TRUST_CGNS_ERROR();

      if (cgp_coord_write(fileId, baseId, zoneId, CGNS_DOUBLE_TYPE, "CoordinateY", &coordsIdy) != CG_OK)
        Cerr << "Error Ecrire_CGNS_helper::cgns_write_zone_grid_coord : cgp_coord_write - Y !" << finl, TRUST_CGNS_ERROR();

      if (Objet_U::dimension > 2)
        if (cgp_coord_write(fileId, baseId, zoneId, CGNS_DOUBLE_TYPE, "CoordinateZ", &coordsIdz) != CG_OK)
          Cerr << "Error Ecrire_CGNS_helper::cgns_write_zone_grid_coord : cgp_coord_write - Z !" << finl, TRUST_CGNS_ERROR();
#endif
    }
}

template<TYPE_ECRITURE_CGNS _TYPE_>
inline std::enable_if_t<_TYPE_ != TYPE_ECRITURE_CGNS::SEQ, void>
Ecrire_CGNS_helper::cgns_write_grid_coord_data(const int icelldim, const int fileId, const int baseId, const int zoneId,
                                               const int coordsIdx, const int coordsIdy, const int coordsIdz, const cgsize_t min, const cgsize_t max,
                                               const std::vector<double>& xCoords, const std::vector<double>& yCoords, const std::vector<double>& zCoords)
{
#ifdef MPI_
  if (cgp_coord_write_data(fileId, baseId, zoneId, coordsIdx, &min, &max, xCoords.data()) != CG_OK)
    Cerr << "Error Ecrire_CGNS_helper::cgns_write_grid_coord_data : cgp_coord_write_data - X !" << finl, TRUST_CGNS_ERROR();

  if (cgp_coord_write_data(fileId, baseId, zoneId, coordsIdy, &min, &max, yCoords.data()) != CG_OK)
    Cerr << "Error Ecrire_CGNS_helper::cgns_write_grid_coord_data : cgp_coord_write_data - Y !" << finl, TRUST_CGNS_ERROR();

  if (icelldim > 2)
    if (cgp_coord_write_data(fileId, baseId, zoneId, coordsIdz, &min, &max, zCoords.data()) != CG_OK)
      Cerr << "Error Ecrire_CGNS_helper::cgns_write_grid_coord_data : cgp_coord_write_data - Z !" << finl, TRUST_CGNS_ERROR();
#endif
}

template<TYPE_ECRITURE_CGNS _TYPE_>
inline void Ecrire_CGNS_helper::cgns_sol_write(const int nb_zones_to_write, const int fileId, const int baseId, const int ind,
                                               const double temps, const std::vector<int>& zoneId, const std::string& LOC,
                                               std::string& solname_som, std::string& solname_elem, std::string& solname_faces,
                                               bool& solname_som_written, bool& solname_elem_written, bool& solname_faces_written,
                                               int& flowId_som, int& flowId_elem, int& flowId_faces)
{
  // uen fois par dt !!
  constexpr bool is_SEQ = (_TYPE_ == TYPE_ECRITURE_CGNS::SEQ), is_PAR_OVER = (_TYPE_ == TYPE_ECRITURE_CGNS::PAR_OVER);

  if (!solname_som_written && LOC == "SOM")
    {
      std::string solname = "FlowSolution" + convert_double_to_string(temps) + "_" + LOC;
      solname.resize(CGNS_STR_SIZE, ' ');
      solname_som += solname;

      // on boucle seulement sur les procs qui n'ont pas des nb_elem 0
      for (int ii = 0; ii != nb_zones_to_write; ii++)
        {
          if (cg_sol_write(fileId, baseId, zoneId[is_PAR_OVER ? ii : ind], solname.c_str(), CGNS_ENUMV(Vertex), &flowId_som) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_sol_write : cg_sol_write -- SOM !" << finl, TRUST_CGNS_ERROR();

          if (!is_SEQ)
            if (cg_goto(fileId, baseId, "Zone_t", zoneId[is_PAR_OVER ? ii : ind], "FlowSolution_t", flowId_som, "end") != CG_OK)
              Cerr << "Error Ecrire_CGNS_helper::cgns_sol_write : cg_goto -- SOM !" << finl, TRUST_CGNS_ERROR();

        }

      solname_som_written = true;
    }

  if (!solname_elem_written && LOC == "ELEM")
    {
      std::string solname = "FlowSolution" + convert_double_to_string(temps) + "_" + LOC;
      solname.resize(CGNS_STR_SIZE, ' ');
      solname_elem += solname;

      // on boucle seulement sur les procs qui n'ont pas des nb_elem 0
      for (int ii = 0; ii != nb_zones_to_write; ii++)
        {
          if (cg_sol_write(fileId, baseId, zoneId[is_PAR_OVER ? ii : ind], solname.c_str(), CGNS_ENUMV(CellCenter), &flowId_elem) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_sol_write : cg_sol_write  -- ELEM !" << finl, TRUST_CGNS_ERROR();

          if (!is_SEQ)
            if (cg_goto(fileId, baseId, "Zone_t", zoneId[is_PAR_OVER ? ii : ind], "FlowSolution_t", flowId_elem, "end") != CG_OK)
              Cerr << "Error Ecrire_CGNS_helper::cgns_sol_write : cg_goto  -- ELEM !" << finl, TRUST_CGNS_ERROR();

        }

      solname_elem_written = true;
    }

  if (!solname_faces_written && LOC == "FACES")
    {
      std::string solname = "FlowSolution" + convert_double_to_string(temps) + "_" + LOC;
      solname.resize(CGNS_STR_SIZE, ' ');
      solname_faces += solname;

      // on boucle seulement sur les procs qui n'ont pas des nb_elem 0
      for (int ii = 0; ii != nb_zones_to_write; ii++)
        {
          if (cg_sol_write(fileId, baseId, zoneId[is_PAR_OVER ? ii : ind], solname.c_str(), CGNS_ENUMV(CellCenter), &flowId_faces) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_sol_write : cg_sol_write  -- FACES !" << finl, TRUST_CGNS_ERROR();

          if (!is_SEQ)
            if (cg_goto(fileId, baseId, "Zone_t", zoneId[is_PAR_OVER ? ii : ind], "FlowSolution_t", flowId_faces, "end") != CG_OK)
              Cerr << "Error Ecrire_CGNS_helper::cgns_sol_write : cg_goto  -- FACES !" << finl, TRUST_CGNS_ERROR();

        }

      solname_faces_written = true;
    }
}

template<TYPE_ECRITURE_CGNS _TYPE_>
inline std::enable_if_t<_TYPE_ != TYPE_ECRITURE_CGNS::SEQ, void>
Ecrire_CGNS_helper::cgns_field_write(const int nb_zones_to_write, const int fileId, const int baseId, const int ind, const std::vector<int>& zoneId, const std::string& LOC,
                                     const int flowId_som, const int flowId_elem, const int flowId_faces,
                                     const char * id_champ, int& fieldId_som, int& fieldId_elem, int& fieldId_faces)
{
#ifdef MPI_
  constexpr bool is_PAR_OVER = (_TYPE_ == TYPE_ECRITURE_CGNS::PAR_OVER);
  for (int ii = 0; ii != nb_zones_to_write; ii++)
    {
      if (LOC == "SOM")
        {
          if (cgp_field_write(fileId, baseId, zoneId[is_PAR_OVER ? ii : ind], flowId_som, CGNS_DOUBLE_TYPE, id_champ, &fieldId_som) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_field_write : cgp_field_write  -- SOM !" << finl, TRUST_CGNS_ERROR();
        }
      else if (LOC == "ELEM")
        {
          if (cgp_field_write(fileId, baseId, zoneId[is_PAR_OVER ? ii : ind], flowId_elem, CGNS_DOUBLE_TYPE, id_champ, &fieldId_elem) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_field_write : cgp_field_write  -- ELEM !" << finl, TRUST_CGNS_ERROR();
        }
      else if (LOC == "FACES")
        {
          if (cgp_field_write(fileId, baseId, zoneId[is_PAR_OVER ? ii : ind], flowId_faces, CGNS_DOUBLE_TYPE, id_champ, &fieldId_faces) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_field_write : cgp_field_write  -- FACES !" << finl, TRUST_CGNS_ERROR();
        }
      else
        throw std::runtime_error("Ecrire_CGNS_helper::cgns_field_write => Unsupported LOC : " + LOC);
    }
#endif
}

template<TYPE_ECRITURE_CGNS _TYPE_>
inline std::enable_if_t<_TYPE_ == TYPE_ECRITURE_CGNS::SEQ, void>
Ecrire_CGNS_helper::cgns_field_write_data(const int fileId, const int baseId, const int ind, const std::vector<int>& zoneId,
                                          const std::string& LOC, const int flowId_som, const int flowId_elem, const int flowId_faces, const int comp,
                                          const char * id_champ, const DoubleTab& valeurs, int& fieldId_som, int& fieldId_elem, int& fieldId_faces)
{
  if (valeurs.dimension(1) == 1) /* No stride ! */
    {
      if (LOC == "SOM")
        {
          if (cg_field_write(fileId, baseId, zoneId[ind], flowId_som, CGNS_DOUBLE_TYPE, id_champ, valeurs.addr(), &fieldId_som) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_field_seq : cg_field_write  -- SOM !" << finl, TRUST_CGNS_ERROR();
        }
      else if (LOC == "ELEM")
        {
          if (cg_field_write(fileId, baseId, zoneId[ind], flowId_elem, CGNS_DOUBLE_TYPE, id_champ, valeurs.addr(), &fieldId_elem) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_field_seq : cg_field_write  -- ELEM !" << finl, TRUST_CGNS_ERROR();
        }
      else if (LOC == "FACES")
        {
          if (cg_field_write(fileId, baseId, zoneId[ind], flowId_faces, CGNS_DOUBLE_TYPE, id_champ, valeurs.addr(), &fieldId_faces) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_field_seq : cg_field_write  -- FACES !" << finl, TRUST_CGNS_ERROR();
        }
      else
        throw std::runtime_error("Ecrire_CGNS_helper::cgns_field_write_data => Unsupported LOC : " + LOC);
    }
  else
    {
      const int nb = valeurs.dimension(0);
      DoubleTrav field_cgns(nb);
      for (int i = 0; i < nb; i++)
        field_cgns(i) = valeurs(i, comp);

      if (LOC == "SOM")
        {
          if (cg_field_write(fileId, baseId, zoneId[ind], flowId_som, CGNS_DOUBLE_TYPE, id_champ, field_cgns.addr(), &fieldId_som) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_field_seq : cg_field_write  -- SOM !" << finl, TRUST_CGNS_ERROR();
        }
      else if (LOC == "ELEM")
        {
          if (cg_field_write(fileId, baseId, zoneId[ind], flowId_elem, CGNS_DOUBLE_TYPE, id_champ, field_cgns.addr(), &fieldId_elem) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_field_seq : cg_field_write  -- ELEM !" << finl, TRUST_CGNS_ERROR();
        }
      else if (LOC == "FACES")
        {
          if (cg_field_write(fileId, baseId, zoneId[ind], flowId_faces, CGNS_DOUBLE_TYPE, id_champ, field_cgns.addr(), &fieldId_faces) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_field_seq : cg_field_write  -- FACES !" << finl, TRUST_CGNS_ERROR();
        }
      else
        throw std::runtime_error("Ecrire_CGNS_helper::cgns_field_write_data => Unsupported LOC : " + LOC);
    }
}

template<TYPE_ECRITURE_CGNS _TYPE_>
inline std::enable_if_t<_TYPE_ != TYPE_ECRITURE_CGNS::SEQ, void>
Ecrire_CGNS_helper::cgns_field_write_data(const int fileId, const int baseId, const int ind, const std::vector<int>& zoneId, const std::string& LOC,
                                          const int flowId_som, const int flowId_elem, const int flowId_faces,
                                          const int fieldId_som, const int fieldId_elem, const int fieldId_faces,
                                          const int comp, const cgsize_t min, const cgsize_t max, const DoubleTab& valeurs)
{
#ifdef MPI_
  if (valeurs.dimension(1) == 1) /* No stride ! */
    {
      if (LOC == "SOM")
        {
          if (cgp_field_write_data(fileId, baseId, zoneId[ind], flowId_som, fieldId_som, &min, &max, valeurs.addr()) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_field_write_data : cgp_field_write_data  -- SOM !" << finl, TRUST_CGNS_ERROR();
        }
      else if (LOC == "ELEM")
        {
          if (cgp_field_write_data(fileId, baseId, zoneId[ind], flowId_elem, fieldId_elem, &min, &max, valeurs.addr()) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_field_write_data : cgp_field_write_data  -- ELEM !" << finl,  TRUST_CGNS_ERROR();
        }
      else if (LOC == "FACES")
        {
          if (cgp_field_write_data(fileId, baseId, zoneId[ind], flowId_faces, fieldId_faces, &min, &max, valeurs.addr()) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_field_write_data : cgp_field_write_data  -- FACES !" << finl,  TRUST_CGNS_ERROR();
        }
      else
        throw std::runtime_error("Ecrire_CGNS_helper::cgns_field_write_data => Unsupported LOC : " + LOC);
    }
  else
    {
      const int nb = valeurs.dimension(0);
      DoubleTrav field_cgns(nb);
      for (int i = 0; i < nb; i++)
        field_cgns(i) = valeurs(i, comp);

      if (LOC == "SOM")
        {
          if (cgp_field_write_data(fileId, baseId, zoneId[ind], flowId_som, fieldId_som, &min, &max, field_cgns.addr()) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_field_write_data : cgp_field_write_data  -- SOM !" << finl,  TRUST_CGNS_ERROR();
        }
      else if (LOC == "ELEM")
        {
          if (cgp_field_write_data(fileId, baseId, zoneId[ind], flowId_elem, fieldId_elem, &min, &max, field_cgns.addr()) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_field_write_data : cgp_field_write_data  -- ELEM !" << finl,  TRUST_CGNS_ERROR();
        }
      else if (LOC == "FACES")
        {
          if (cgp_field_write_data(fileId, baseId, zoneId[ind], flowId_faces, fieldId_faces, &min, &max, field_cgns.addr()) != CG_OK)
            Cerr << "Error Ecrire_CGNS_helper::cgns_field_write_data : cgp_field_write_data  -- FACES !" << finl,  TRUST_CGNS_ERROR();
        }
      else
        throw std::runtime_error("Ecrire_CGNS_helper::cgns_field_write_data => Unsupported LOC : " + LOC);
    }
#endif
}

template<TYPE_ECRITURE_CGNS _TYPE_>
inline void Ecrire_CGNS_helper::cgns_write_iters(const bool has_field, const int nb_zones_to_write, const int fileId, const int baseId, const int ind, const std::vector<int>& zoneId,
                                                 const std::string& LOC, const std::string& solname_som, const std::string& solname_elem, const std::string& solname_faces,const std::vector<double>& time_post)
{
  const int nsteps = static_cast<int>(time_post.size());
  constexpr bool is_PAR_OVER = (_TYPE_ == TYPE_ECRITURE_CGNS::PAR_OVER);

  /* create BaseIterativeData */
  if (cg_biter_write(fileId, baseId, "TimeIterValues", nsteps) != CG_OK)
    Cerr << "Error Ecrire_CGNS_helper::cgns_write_iters : cg_biter_write !" << finl, TRUST_CGNS_ERROR();

  /* go to BaseIterativeData level and write time values */
  if (cg_goto(fileId, baseId, "BaseIterativeData_t", 1, "end") != CG_OK)
    Cerr << "Error Ecrire_CGNS_helper::cgns_write_iters : cg_goto !" << finl, TRUST_CGNS_ERROR();

  cgsize_t nuse = static_cast<cgsize_t>(nsteps);
  if (cg_array_write("TimeValues", CGNS_DOUBLE_TYPE, 1, &nuse, time_post.data()) != CG_OK)
    Cerr << "Error Ecrire_CGNS_helper::cgns_write_iters : cg_array_write !" << finl, TRUST_CGNS_ERROR();

  cgsize_t idata[2];
  idata[0] = CGNS_STR_SIZE;
  idata[1] = nsteps;

  for (int ii = 0; ii != nb_zones_to_write; ii++)
    if (zoneId[ind] != -123)
      {
        /* create ZoneIterativeData */
        if (cg_ziter_write(fileId, baseId, zoneId[is_PAR_OVER ? ii : ind], "ZoneIterativeData") != CG_OK)
          Cerr << "Error Ecrire_CGNS_helper::cgns_write_iters : cg_ziter_write !" << finl, cg_error_exit();

        if (cg_goto(fileId, baseId, "Zone_t", zoneId[is_PAR_OVER ? ii : ind], "ZoneIterativeData_t", 1, "end") != CG_OK)
          Cerr << "Error Ecrire_CGNS_helper::cgns_write_iters : cg_goto !" << finl, TRUST_CGNS_ERROR();

        if (has_field)
          {
            if (LOC == "SOM")
              {
                if (cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, idata, solname_som.c_str()) != CG_OK)
                  Cerr << "Error Ecrire_CGNS_helper::cgns_write_iters : cg_array_write  -- SOM !" << finl, TRUST_CGNS_ERROR();
              }
            else if (LOC == "ELEM")
              {
                if (cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, idata, solname_elem.c_str()) != CG_OK)
                  Cerr << "Error Ecrire_CGNS_helper::cgns_write_iters : cg_array_write  -- ELEM !" << finl, TRUST_CGNS_ERROR();
              }
            else if (LOC == "FACES")
              {
                if (cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, idata, solname_faces.c_str()) != CG_OK)
                  Cerr << "Error Ecrire_CGNS_helper::cgns_write_iters : cg_array_write  -- FACES !" << finl, TRUST_CGNS_ERROR();
              }
            else
              throw std::runtime_error("Ecrire_CGNS_helper::cgns_write_iters => Unsupported LOC : " + LOC);
          }
      }

  if (cg_simulation_type_write(fileId, baseId, CGNS_ENUMV(TimeAccurate)) != CG_OK)
    Cerr << "Error Ecrire_CGNS_helper::cgns_write_iters : cg_simulation_type_write !" << finl, TRUST_CGNS_ERROR();
}

#endif /* Ecrire_CGNS_helper_tpp_included */
