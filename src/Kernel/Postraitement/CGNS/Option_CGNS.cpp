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

#include <Option_CGNS.h>
#include <Motcle.h>
#include <Param.h>

Implemente_instanciable(Option_CGNS, "Option_CGNS", Interprete);
// XD Option_CGNS interprete Option_CGNS 1 Class for CGNS options.

bool Option_CGNS::SINGLE_PRECISION = false; /* NOT BY DEFAULT */
bool Option_CGNS::PARALLEL_OVER_ZONE = false; /* NOT BY DEFAULT */
bool Option_CGNS::USE_LINKS = false; /* NOT BY DEFAULT */
bool Option_CGNS::FILE_PER_COMM_GROUP = false; /* NOT BY DEFAULT */
bool Option_CGNS::SINGLE_SAFE_FILE = false; /* NOT BY DEFAULT */
int Option_CGNS::CLOSE_EVERY_N = -1; /* -1 BY DEFAULT => never opened/closed */
int Option_CGNS::FLUSH_EVERY_N = 1; /* 1 BY DEFAULT => flush each dt post */

Sortie& Option_CGNS::printOn(Sortie& os) const { return Interprete::printOn(os); }
Entree& Option_CGNS::readOn(Entree& is) { return Interprete::readOn(is); }

Entree& Option_CGNS::interpreter(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter_non_std("SINGLE_PRECISION", (this)); // XD_ADD_P rien If used, data will be written with a single_precision format inside the CGNS file (it concerns both mesh coordinates and field values).
  param.ajouter_non_std("PARALLEL_OVER_ZONE", (this)); // XD_ADD_P rien If used, data will be written in separate zones (ie: one zone per processor). This is not so performant but easier to read later ...
  param.ajouter_non_std("USE_LINKS", (this)); // XD_ADD_P rien If used, data will be written in separate files; one file for mesh, and then one file for solution time. Links will be used.
  param.ajouter_non_std("FILE_PER_COMM_GROUP", (this)); // XD_ADD_P rien If used, data will be written (at each comm group) in separate files; one file for mesh, and then one file for solution time. Links will be used.
  param.ajouter_non_std("SINGLE_SAFE_FILE", (this)); // XD_ADD_P rien If used, data will be written in a single file that will be opened and closed at each dt post so that file can be visualized in live. Safer if simulation stops, the file can be used.
  param.ajouter("CLOSE_EVERY_N", &CLOSE_EVERY_N); // XD_ADD_P entier Used to fix the opening/closing frequency when the SINGLE_SAFE_FILE option is used.
  param.ajouter("FLUSH_EVERY_N", &FLUSH_EVERY_N); // XD_ADD_P entier Used to fix the flush-to-disc frequency when the SINGLE_SAFE_FILE option is used.
  param.lire_avec_accolades_depuis(is);

  if (PARALLEL_OVER_ZONE && (USE_LINKS || FILE_PER_COMM_GROUP))
    {
      Cerr << "Error in Option_CGNS :" << finl;
      Cerr << "       - You can not activate the option 'PARALLEL_OVER_ZONE' with 'USE_LINKS' and/or 'FILE_PER_COMM_GROUP' !!!" << finl;
      Process::exit();
    }

  if (SINGLE_SAFE_FILE && (USE_LINKS || FILE_PER_COMM_GROUP))
    {
      Cerr << "Error in Option_CGNS :" << finl;
      Cerr << "       - You can not activate the option 'SINGLE_SAFE_FILE' with 'USE_LINKS' and/or 'FILE_PER_COMM_GROUP' !!!" << finl;
      Process::exit();
    }

  return is;
}

int Option_CGNS::lire_motcle_non_standard(const Motcle& mot_cle, Entree& is)
{
  int retval = 1;

  if (mot_cle == "SINGLE_PRECISION")
    {
      Cerr << mot_cle << " => CGNS data will be written in a single precision format ..." << finl;
      SINGLE_PRECISION = true;
    }
  else if (mot_cle == "PARALLEL_OVER_ZONE")
    {
      Cerr << mot_cle << " => CGNS data will be written in separate zones ..." << finl;
      PARALLEL_OVER_ZONE = true;
    }
  else if (mot_cle == "USE_LINKS")
    {
      Cerr << mot_cle << " => CGNS data will be written in separate files (mesh, solution ...)" << finl;
      USE_LINKS = true;
    }
  else if (mot_cle == "FILE_PER_COMM_GROUP")
    {
      Cerr << mot_cle << " => A CGNS file will be written in each COMM group ..." << finl;
      FILE_PER_COMM_GROUP = true;
      Cerr << "USE_LINKS => CGNS data will be written in separate files (mesh, solution ...)" << finl;
      USE_LINKS = true;
    }
  else if (mot_cle == "SINGLE_SAFE_FILE")
    {
      Cerr << mot_cle << " => A single smart CGNS file will be written ..." << finl;
      SINGLE_SAFE_FILE = true;
    }
  else
    retval = -1;

  return retval;
}
