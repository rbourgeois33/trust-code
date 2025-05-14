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

#include <Extraire_surface.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <NettoieNoeuds.h>
#include <Parser_U.h>
#include <Domaine.h>
#include <Domaine_VF.h>
#include <Param.h>
#include <ParserView.h>
#include <Extraire_plan.h>

Implemente_instanciable(Extraire_surface,"Extraire_surface",Interprete_geometrique_base);
// XD extraire_surface interprete extraire_surface -3 This keyword extracts a surface mesh named domain_name (this domain should have been declared before) from the mesh of the pb_name problem. The surface mesh is defined by one or two conditions. The first condition is about elements with Condition_elements. For example: Condition_elements x*x+y*y+z*z<1 NL2 Will define a surface mesh with external faces of the mesh elements inside the sphere of radius 1 located at (0,0,0). The second condition Condition_faces is useful to give a restriction.NL2 By default, the faces from the boundaries are not added to the surface mesh excepted if option avec_les_bords is given (all the boundaries are added), or if the option avec_certains_bords is used to add only some boundaries.

Sortie& Extraire_surface::printOn(Sortie& os) const { return Interprete::printOn(os); }

Entree& Extraire_surface::readOn(Entree& is) { return Interprete::readOn(is); }

Entree& Extraire_surface::interpreter_(Entree& is)
{
  Nom nom_pb;
  Nom nom_domaine_surfacique;
  Nom expr_elements("1"),expr_faces("1");
  bool avec_les_bords = false;
  Noms noms_des_bords, groupes_faces;
  Param param(que_suis_je());
  param.ajouter("domaine",&nom_domaine_surfacique,Param::REQUIRED); // XD_ADD_P ref_domaine Domain in which faces are saved
  param.ajouter("probleme",&nom_pb,Param::REQUIRED); // XD_ADD_P ref_Pb_base Problem from which faces should be extracted
  param.ajouter("condition_elements",&expr_elements); // XD_ADD_P chaine condition on center of elements
  param.ajouter("condition_faces",&expr_faces); // XD_ADD_P chaine not_set
  param.ajouter_flag("avec_les_bords",&avec_les_bords); // XD_ADD_P rien not_set
  param.ajouter("avec_certains_bords",&noms_des_bords); // XD_ADD_P listchaine not_set
  param.ajouter("groupes_faces",&groupes_faces); // XD_ADD_P listchaine not_set
  param.lire_avec_accolades_depuis(is);

  associer_domaine(nom_domaine_surfacique);
  Domaine& domaine_surfacique=domaine();

  if (domaine_surfacique.nb_som_tot()!=0)
    {
      Cerr << "Error!" << finl;
      Cerr <<"The domain " << domaine_surfacique.le_nom() << " can't be used by the Extraire_surface keyword." <<finl;
      Cerr <<"The domain for Extraire_surface keyword should be empty and created just before." << finl;
      exit();
    }

  // on recupere le pb
  if(! sub_type(Probleme_base, objet(nom_pb)))
    {
      Cerr << nom_pb << " is of type " << objet(nom_pb).que_suis_je() << finl;
      Cerr << "and not of type Probleme_base" << finl;
      exit();
    }
  Probleme_base& pb=ref_cast(Probleme_base, objet(nom_pb));
  const Domaine_VF& domaine_vf=ref_cast(Domaine_VF,pb.domaine_dis());
  const Domaine& domaine_volumique = domaine_vf.domaine();

  extraire_surface(domaine_surfacique,domaine_volumique,nom_domaine_surfacique,domaine_vf,expr_elements,expr_faces,avec_les_bords,noms_des_bords, groupes_faces);

  return is;
}

// Extraction d'une ou plusieurs frontieres du domaine volumique selon certaines conditions

void Extraire_surface::extraire_surface(Domaine& domaine_surfacique,const Domaine& domaine_volumique, const Nom& nom_domaine_surfacique, const Domaine_VF& domaine_vf, const Nom& expr_elements,const Nom& expr_faces, bool avec_les_bords, const Noms& noms_des_bords, const Noms& groupes_faces)
{
  extraire_surface_without_cleaning(domaine_surfacique,domaine_volumique,nom_domaine_surfacique,domaine_vf,expr_elements,expr_faces,avec_les_bords,noms_des_bords, groupes_faces);
  NettoieNoeuds::nettoie(domaine_surfacique);
}
void Extraire_surface::extraire_surface_without_cleaning(Domaine& domaine_surfacique,const Domaine& domaine_volumique, const Nom& nom_domaine_surfacique, const Domaine_VF& domaine_vf, const Nom& expr_elements,const Nom& expr_faces, bool avec_les_bords, const Noms& noms_des_bords, const Noms& groupes_faces)
{
  domaine_surfacique.nommer(nom_domaine_surfacique);
  Parser_U condition_elements,condition_faces;
  condition_elements.setNbVar(dimension);
  condition_elements.addVar("x");
  condition_elements.addVar("y");
  if (dimension==3)
    condition_elements.addVar("z");
  condition_faces.setNbVar(dimension);
  condition_faces.addVar("x");
  condition_faces.addVar("y");
  if (dimension==3)
    condition_faces.addVar("z");

  condition_faces.setString(expr_faces);
  condition_faces.parseString();
  condition_elements.setString(expr_elements);
  condition_elements.parseString();

  // Copie des sommets
  domaine_surfacique.les_sommets()=domaine_volumique.les_sommets();
  const Nom& type_elem=domaine_vf.domaine().type_elem()->que_suis_je();

  if (dimension==3)
    if (type_elem==Motcle("Tetraedre"))
      domaine_surfacique.typer("Triangle");
    else if (type_elem==Motcle("Segment"))
      domaine_surfacique.typer("Point");
    else if (type_elem==Motcle("Segment_axi"))
      domaine_surfacique.typer("Point");
    else if ((type_elem==Motcle("Hexaedre"))|| (type_elem==Motcle("Hexaedre_VEF")))
      domaine_surfacique.typer("Quadrangle");
    else if (type_elem==Motcle("Polyedre"))
      domaine_surfacique.typer("Polygone");
    else
      {
        Cerr<<"WARNING "<<type_elem<< " not coded, use Quadrangle" <<finl;
        domaine_surfacique.typer("Quadrangle");
        //      exit();
      }
  else
    domaine_surfacique.typer("segment");

  int dim = Objet_U::dimension;
  int nb_elem=domaine_vf.nb_elem();
  IntTrav tab_marq_elem;
  domaine_vf.domaine().creer_tableau_elements(tab_marq_elem);

  ParserView parser_condition_elements(condition_elements);
  parser_condition_elements.parseString();
  CDoubleTabView xp = domaine_vf.xp().view_ro();
  IntArrView marq_elem = static_cast<ArrOfInt&>(tab_marq_elem).view_wo();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_elem, KOKKOS_LAMBDA(const int elem)
  {
    int threadId = parser_condition_elements.acquire();
    parser_condition_elements.setVar(0,xp(elem,0),threadId);
    parser_condition_elements.setVar(1,xp(elem,1),threadId);
    if (dim==3)
      parser_condition_elements.setVar(2,xp(elem,2),threadId);
    double res = parser_condition_elements.eval(threadId);
    marq_elem(elem) = std::fabs(res)>1e-5 ? 1 : 0;
    parser_condition_elements.release(threadId);
  });
  end_gpu_timer(__KERNEL_NAME__);
  tab_marq_elem.echange_espace_virtuel();

  int nb_faces=domaine_vf.nb_faces();

  ArrOfInt tab_marq(nb_faces);
  // on marque les joints
  int nbjoints=domaine_vf.nb_joints();
  if (nbjoints>0)
    {
      ToDo_Kokkos("critical");
    }
  for(int njoint=0; njoint<nbjoints; njoint++)
    {
      const Joint& joint_temp = domaine_vf.joint(njoint);
      int pe_voisin=joint_temp.PEvoisin();
      if (pe_voisin<me())
        {
          const IntTab& indices_faces_joint = joint_temp.joint_item(JOINT_ITEM::FACE).renum_items_communs();
          const int nbfaces = indices_faces_joint.dimension(0);
          for (int j = 0; j < nbfaces; j++)
            {
              int face_de_joint = indices_faces_joint(j, 1);
              tab_marq[face_de_joint] = -1;
            }
        }
    }
  int nb_t=0;

  // on flage les face sde bords qui nous interesse
  ArrOfInt tab_face_bord_int(nb_faces);
  if (avec_les_bords)
    tab_face_bord_int=1;
  if (noms_des_bords.size()!=0)
    {
      if (avec_les_bords)
        {
          Cerr<<" the option avec_les_bords is incompatible with the option avec_certains_bords"<<finl;
          exit();
        }
      ToDo_Kokkos("critical");
      for (int b=0; b<noms_des_bords.size(); b++)
        {
          const Frontiere& fr=domaine_vf.frontiere_dis(domaine_vf.rang_frontiere(noms_des_bords[b])).frontiere();
          int deb=fr.num_premiere_face();
          int fin=deb+fr.nb_faces();
          for (int f=deb; f<fin; f++)
            tab_face_bord_int[f]=1;
        }
    }

  // on marque toutes les faces que l'on veut mettre dans le domaine
  ParserView parser_condition_faces(condition_faces);
  parser_condition_faces.parseString();
  CDoubleTabView xv = domaine_vf.xv().view_ro();
  CIntTabView face_voisin = domaine_vf.face_voisins().view_ro();
  CIntArrView face_bord_int = tab_face_bord_int.view_ro();
  IntArrView marq = tab_marq.view_rw();
  Kokkos::parallel_reduce(start_gpu_timer(__KERNEL_NAME__), nb_faces, KOKKOS_LAMBDA(const int fac, int& local_nb_t)
  {
    int elem0=face_voisin(fac,0);
    int elem1=face_voisin(fac,1);
    int val0=-1;
    if (elem0!=-1) val0=marq_elem(elem0);
    int val1=-1;
    if (elem1!=-1) val1=marq_elem(elem1);

    if (val0!=val1)
      {
        if (((val0*val1==0)&&(val0+val1==1))||((val0==1)&&(face_bord_int[fac]==1))||((val1==1)&&(face_bord_int[fac]==1)))
          //if (domaine_test.chercher_elements(xv(fac,0),xv(fac,1),xv(fac,2))==0)
          {
            int threadId = parser_condition_faces.acquire();
            parser_condition_faces.setVar(0,xv(fac,0),threadId);
            parser_condition_faces.setVar(1,xv(fac,1),threadId);
            if (dim==3)
              parser_condition_faces.setVar(2,xv(fac,2),threadId);
            double res=parser_condition_faces.eval(threadId);
            if (std::fabs(res)>1e-5)
              if (marq[fac]!=-1)  // pas un joint, ou on est le proprietaire
                {
                  marq[fac]=1;
                  local_nb_t++;
                }
            parser_condition_faces.release(threadId);
          }
      }
  }, nb_t);
  end_gpu_timer(__KERNEL_NAME__);

  for (auto && nom : groupes_faces)
    {
      const ArrOfInt& arr = domaine_vf.domaine().groupe_faces(nom).get_indices_faces();
      for (int i = 0; i < arr.size_array(); i++)
        if (marq[arr[i]] != 1)
          marq[arr[i]] = 1, nb_t++;
    }

  ArrOfInt tab_indices(nb_faces);
  IntArrView indices = tab_indices.view_wo();
  Kokkos::parallel_scan(start_gpu_timer(__KERNEL_NAME__), nb_faces, KOKKOS_LAMBDA(const int fac, int& update, const bool final)
  {
    if (marq[fac] == 1)
      {
        if (final) indices(fac) = update;
        update++;
      }
  });
  end_gpu_timer(__KERNEL_NAME__);



  Cerr<<"Number of elements of the new domain : "<<nb_t<<finl;
  CIntTabView face_sommets = domaine_vf.face_sommets().view_ro();
  int nb_sommet_face=(int)face_sommets.extent(1);
  IntTab& tab_les_elems=domaine_surfacique.les_elems();
  tab_les_elems.resize(nb_t,nb_sommet_face);
  CDoubleTabView coord = domaine_surfacique.les_sommets().view_ro();
  IntTabView les_elems = tab_les_elems.view_wo();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_faces, KOKKOS_LAMBDA(const int fac)
  {
    if (marq[fac] == 1)
      {
        int nb = indices(fac);
        double normal[3] = {0., 0., 0.} , normal_b[3], point0b[3], point1b[3], point2b[3];
        int el1 = face_voisin(fac, 0);
        if (el1 == -1) el1 = face_voisin(fac, 1);
        if (marq_elem(el1) != 1)
          el1 = face_voisin(fac, 1);
        for (int i = 0; i < dim; i++)
          normal[i] = xv(fac, i) - xp(el1, i);
        for (int s = 0; s < nb_sommet_face; s++)
          les_elems(nb, s) = face_sommets(fac, s);
        // on calcule la normale
        if (nb_sommet_face > 1)
          {
            if (dim == 3)
              {
                for (int i = 0; i < 3; i++)
                  {
                    point0b[i] = coord(les_elems(nb, 0), i);
                    point1b[i] = coord(les_elems(nb, 1), i);
                    point2b[i] = coord(les_elems(nb, 2), i);
                  }
                calcul_normal(point0b, point1b, point2b, normal_b);
                double dot_product=0;
                for (int i = 0; i < dim; i++)
                  dot_product+=normal[i]*normal_b[i];
                if (dot_product < 0)
                  {
                    // si normal a l'envers on inverse les deux sommets
                    les_elems(nb, 1) = face_sommets(fac, 2);
                    les_elems(nb, 2) = face_sommets(fac, 1);
                  }
              }
            else
              {
                for (int i = 0; i < nb_sommet_face; i++)
                  point0b[i] = coord(les_elems(nb, 1), i) - coord(les_elems(nb, 0), i);
                double produit = normal[0] * point0b[1] - normal[1] * point0b[0];
                if (produit < 0)
                  {
                    // si normal a l'envers on inverse les deux sommets
                    les_elems(nb, 0) = face_sommets(fac, 1);
                    les_elems(nb, 1) = face_sommets(fac, 0);
                  }
              }
          }
      }
  });
  end_gpu_timer(__KERNEL_NAME__);
}
