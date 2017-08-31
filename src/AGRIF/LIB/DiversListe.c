/******************************************************************************/
/*                                                                            */
/*     CONV (converter) for Agrif (Adaptive Grid Refinement In Fortran)       */
/*                                                                            */
/* Copyright or   or Copr. Laurent Debreu (Laurent.Debreu@imag.fr)            */
/*                        Cyril Mazauric (Cyril_Mazauric@yahoo.fr)            */
/* This software is governed by the CeCILL-C license under French law and     */
/* abiding by the rules of distribution of free software.  You can  use,      */
/* modify and/ or redistribute the software under the terms of the CeCILL-C   */
/* license as circulated by CEA, CNRS and INRIA at the following URL          */
/* "http://www.cecill.info".                                                  */
/*                                                                            */
/* As a counterpart to the access to the source code and  rights to copy,     */
/* modify and redistribute granted by the license, users are provided only    */
/* with a limited warranty  and the software's author,  the holder of the     */
/* economic rights,  and the successive licensors  have only  limited         */
/* liability.                                                                 */
/*                                                                            */
/* In this respect, the user's attention is drawn to the risks associated     */
/* with loading,  using,  modifying and/or developing or reproducing the      */
/* software by the user in light of its specific status of free software,     */
/* that may mean  that it is complicated to manipulate,  and  that  also      */
/* therefore means  that it is reserved for developers  and  experienced      */
/* professionals having in-depth computer knowledge. Users are therefore      */
/* encouraged to load and test the software's suitability as regards their    */
/* requirements in conditions enabling the security of their systems and/or   */
/* data to be ensured and,  more generally, to use and operate it in the      */
/* same conditions as regards security.                                       */
/*                                                                            */
/* The fact that you are presently reading this means that you have had       */
/* knowledge of the CeCILL-C license and that you accept its terms.           */
/******************************************************************************/
/* version 1.7                                                                */
/******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "decl.h"

/******************************************************************************/
/*                           Add_Common_var_1                                 */
/******************************************************************************/
/*  This subroutines is used to add the variable defined in common in the     */
/*     List_Common_Var                                                        */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void Add_Common_var_1()
{
   listvar *newvar;
   listvar *newvar2;
   variable *newvariable;
   listdim *dims;
   char listdimension[LONG_C];
   char ligne[LONG_C];
   int out;

   if ( firstpass == 1 )
   {

   newvar = (listvar *)malloc(sizeof(listvar));
   newvariable = (variable *)malloc(sizeof(variable));
   /*                                                                         */
   Init_Variable(newvariable);
   /*                                                                         */
   strcpy(newvariable->v_nomvar,commonvar);
   Save_Length(commonvar,4);
   strcpy(newvariable->v_commonname,commonblockname);
   Save_Length(commonblockname,7);
   strcpy(newvariable->v_modulename,curmodulename);
   Save_Length(curmodulename,6);
   strcpy(newvariable->v_subroutinename,subroutinename);
   Save_Length(subroutinename,11);
   newvariable->v_positioninblock= positioninblock;
   newvariable->v_common=1;
   strcpy(newvariable->v_commoninfile,mainfile);
   Save_Length(mainfile,10);

   newvar->var = newvariable;

   if ( commondim )
   {
      newvariable->v_dimension=commondim;
      newvariable->v_dimensiongiven=1;
      newvariable->v_nbdim=num_dims(commondim);
      /* Creation of the string for the dimension of this variable            */
      dimsempty = 1;
      strcpy(listdimension,"");

      if ( commondim )
      {
         dims = commondim;
         while (dims)
         {
            if ( strcasecmp(dims->dim.first,"") ||
                 strcasecmp(dims->dim.last,""))  dimsempty = 0;
            sprintf(ligne,"%s:%s",dims->dim.first,dims->dim.last);
            strcat(listdimension,ligne);
            if ( dims->suiv ) strcat(listdimension,",");
            dims = dims->suiv;
         }
         if ( dimsempty == 1 ) newvariable->v_dimsempty=1;
      }
      strcpy(newvariable->v_readedlistdimension,listdimension);
      Save_Length(listdimension,15);
   }


   newvar->suiv = NULL;

   if ( !List_Common_Var )
   {
      List_Common_Var = newvar;
   }
   else
   {
      newvar2 = List_Common_Var;
      out = 0 ;
      while ( newvar2 && out == 0 )
      {
         if ( !strcasecmp(newvar2->var->v_nomvar,commonvar) &&
              !strcasecmp(newvar2->var->v_commonname,commonblockname) &&
              !strcasecmp(newvar2->var->v_subroutinename,subroutinename)
                          ) out = 1 ;
         else newvar2 = newvar2->suiv;
      }
      if ( out == 0 )
      {
         newvar->suiv = List_Common_Var;
         List_Common_Var = newvar;
      }
      else
      {
         free(newvar);
      }
   }
   }
}

/******************************************************************************/
/*                           Addtolistnom                                     */
/******************************************************************************/
/* This subroutine is used to add a variable to the list                      */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
listnom *Addtolistnom(char *nom, listnom *listin,int value)
{
   listnom *newnom;
   listnom *parcours;
   int out;

   newnom=(listnom *) malloc (sizeof (listnom));
   strcpy(newnom->o_nom,nom);
   Save_Length(nom,23);
   newnom->o_val = value;
   newnom->suiv = NULL;

   if ( !listin ) listin = newnom;
   else
   {
      parcours = listin;
      out = 0 ;
      while ( parcours && out == 0 )
      {
         if ( !strcasecmp(parcours->o_nom,nom) ) out = 1 ;
         else parcours=parcours->suiv;
      }
      if ( out == 0 )
      {
          newnom->suiv = listin;
          listin = newnom;
      }
      else
      {
         free(newnom);
      }
   }
   return listin;
}

/******************************************************************************/
/*                           Addtolistname                                    */
/******************************************************************************/
/* This subroutine is used to add a        variable to the list               */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ glob +--->+ glob +--->+ glob +--->+ glob +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
listname *Addtolistname(char *nom,listname *input)
{
   listname *newnom;
   listname *parcours;
   int out;

   if ( !input )
   {
      newnom=(listname *) malloc (sizeof (listname));
      strcpy(newnom->n_name,nom);
      Save_Length(nom,20);
      newnom->suiv = NULL;
      input = newnom;
   }
   else
   {
      parcours = input;
      out = 0 ;
      while ( parcours && out == 0 )
      {
         if ( !strcasecmp(parcours->n_name,nom) ) out = 1;
         else parcours=parcours->suiv;
      }
      if ( out == 0 )
      {
         newnom=(listname *) malloc (sizeof (listname));
         strcpy(newnom->n_name,nom);
         Save_Length(nom,20);
         newnom->suiv = input;
         input = newnom;
      }
   }
   return input;
}

/******************************************************************************/
/*                    ModuleIsDefineInInputFile                               */
/******************************************************************************/
/* This subroutine is used to know if the module is defined in the input file */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/******************************************************************************/
int ModuleIsDefineInInputFile(char *name)
{
   listnom *newnom;
   int out;

   out = 0;
   if ( listofmodules )
   {
      newnom = listofmodules;
      while( newnom && out == 0 )
      {
         if ( !strcasecmp(newnom->o_nom,name) ) out = 1 ;
         else newnom=newnom->suiv;
      }
   }
   return out;
}

/******************************************************************************/
/*                      Addmoduletothelisttmp                                 */
/******************************************************************************/
/* This subroutine is used to add a record to a list of struct                */
/* listusemodule                                                              */
/******************************************************************************/
/*                                                                            */
/*       subroutine sub ... USE mod1 ===> insert in list                      */
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ list +--->+ list +--->+ list +--->+ list +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*       list =  listofmoduletmp                                              */
/*                                                                            */
/******************************************************************************/
void Addmoduletothelisttmp(char *name)
{
  listusemodule *newmodule;
  listusemodule *parcours;
  int out;

  if ( !listofmoduletmp)
  {
    newmodule =(listusemodule *)malloc(sizeof(listusemodule));
    strcpy(newmodule->u_usemodule,name);
    Save_Length(name,16);
    strcpy(newmodule->u_cursubroutine,subroutinename);
    Save_Length(subroutinename,18);
    newmodule->suiv = NULL;
    listofmoduletmp = newmodule ;
  }
  else
  {
    parcours = listofmoduletmp;
    out = 0;
    while( parcours && out == 0 )
    {
       if ( !strcasecmp(parcours->u_usemodule,name) ) out = 1;
       else parcours = parcours->suiv;
    }
    if ( out == 0 )
    {
       newmodule =(listusemodule *)malloc(sizeof(listusemodule));
       strcpy(newmodule->u_usemodule,name);
       Save_Length(name,16);
       strcpy(newmodule->u_cursubroutine,subroutinename);
       Save_Length(subroutinename,18);
       newmodule->suiv = listofmoduletmp;
       listofmoduletmp = newmodule;
    }
  }
}

/******************************************************************************/
/*                          Add_NameOfModule_1                                */
/******************************************************************************/
/* This subroutine is used to add a        variable to the list               */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ glob +--->+ glob +--->+ glob +--->+ glob +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void Add_NameOfModule_1(char *nom)
{
   listnom *newnom;

   if ( firstpass == 1 )
   {
      newnom=(listnom *) malloc (sizeof (listnom));
      strcpy(newnom->o_nom,nom);
      Save_Length(nom,23);
      newnom->suiv = List_NameOfModule;
      List_NameOfModule = newnom;
   }
}

/******************************************************************************/
/*                          Add_NameOfCommon_1                                */
/******************************************************************************/
/* This subroutine is used to add a        variable to the list               */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ glob +--->+ glob +--->+ glob +--->+ glob +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void Add_NameOfCommon_1(char *nom,char *cursubroutinename)
{
   listnom *newnom;
   listnom *parcours;

   if ( firstpass == 1 )
   {
      parcours = List_NameOfCommon;
      while ( parcours && strcasecmp(parcours->o_nom,nom) )
                                                      parcours = parcours->suiv;
      if ( !parcours )
      {
         newnom=(listnom *) malloc (sizeof (listnom));
         strcpy(newnom->o_nom,nom);
         strcpy(newnom->o_subroutinename,cursubroutinename);
         Save_Length(nom,23);
         newnom->suiv = List_NameOfCommon;
         List_NameOfCommon = newnom;
      }
   }
}

/******************************************************************************/
/*                     Add_CouplePointed_Var_1                                */
/******************************************************************************/
/* Firstpass 1                                                                */
/* We should complete the listvarpointtovar                                   */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void Add_CouplePointed_Var_1(char *namemodule,listcouple *couple)
{
   listvarpointtovar *pointtmp;

   if ( firstpass == 1 )
   {
      /* we should complete the List_CouplePointed_Var                        */
      pointtmp=(listvarpointtovar *)malloc(sizeof(listvarpointtovar));
      strcpy(pointtmp->t_usemodule,namemodule);
      Save_Length(namemodule,28);
      strcpy(pointtmp->t_cursubroutine,subroutinename);
      Save_Length(subroutinename,29);
      pointtmp->t_couple = couple;
      if ( List_CouplePointed_Var )
      {
         pointtmp->suiv = List_CouplePointed_Var;
         List_CouplePointed_Var = pointtmp;
      }
      else
      {
         pointtmp->suiv = NULL;
         List_CouplePointed_Var = pointtmp;
      }
   }
}

/******************************************************************************/
/*                           Add_Include_1                                    */
/******************************************************************************/
/* This subroutine is used to add a record to a list of struct                */
/*  List_Include                                                              */
/******************************************************************************/
/*                                                                            */
/*       subroutine sub ... include mod1 ===> insert in list                  */
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ list +--->+ list +--->+ list +--->+ list +             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/*       list =  List_Include                                                 */
/*                                                                            */
/******************************************************************************/
void Add_Include_1(char *name)
{
  listusemodule *newinclude;

  if ( firstpass == 1 )
  {
  newinclude =(listusemodule *)malloc(sizeof(listusemodule));
  strcpy(newinclude->u_usemodule,name);
  Save_Length(name,16);
  strcpy(newinclude->u_cursubroutine,subroutinename);
  Save_Length(subroutinename,18);
  newinclude->suiv = NULL;

  if ( !List_Include)
  {
     List_Include  = newinclude ;
  }
  else
  {
    newinclude->suiv = List_Include;
    List_Include = newinclude;
  }
  }
}

/******************************************************************************/
/*                     Add_ImplicitNoneSubroutine_1                           */
/******************************************************************************/
/* This subroutine is used to add a record to a list of struct                */
/******************************************************************************/
/*                                                                            */
/*                                                                            */
/******************************************************************************/
void Add_ImplicitNoneSubroutine_1()
{

  if ( firstpass == 1 )
  {
     List_ImplicitNoneSubroutine = Addtolistname(subroutinename,
                                                   List_ImplicitNoneSubroutine);
  }
}


/******************************************************************************/
/*                        WriteIncludeDeclaration                             */
/******************************************************************************/
/* Firstpass 0                                                                */
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void WriteIncludeDeclaration()
{
  listusemodule *newinclude;

  newinclude = List_Include;
  fprintf(fortran_out,"\n");
  while ( newinclude )
  {
     if ( !strcasecmp(newinclude->u_cursubroutine,subroutinename) )
     {
        fprintf(fortran_out,"      INCLUDE %s \n",newinclude->u_usemodule);
     }
     newinclude = newinclude ->suiv;
  }
}

/******************************************************************************/
/*                          Add_Save_Var_1                                    */
/******************************************************************************/
/* This subroutine is used to add a record to List_Save_Var                   */
/******************************************************************************/
/*        _______     _______     _______     _______     _______             */
/*       +      +    +      +    +      +    +      +    +      +             */
/*       + NEW  +--->+ Save +--->+ Save +--->+ Save +--->+  Save+             */
/*       +______+    +______+    +______+    +______+    +______+             */
/*                                                                            */
/******************************************************************************/
void Add_Save_Var_1 (char *name,listdim *d)
{
  listvar *newvar;
  listdim *dims;
  char ligne[LONG_C];
  char listdimension[LONG_C];

  if ( firstpass == 1 )
  {
     newvar=(listvar *)malloc(sizeof(listvar));
     newvar->var=(variable *)malloc(sizeof(variable));
     /*                                                                       */
     Init_Variable(newvar->var);
     /*                                                                       */
     newvar->var->v_save=1;
     strcpy(newvar->var->v_nomvar,name);
     Save_Length(name,4);
     strcpy(newvar->var->v_modulename,curmodulename);
     Save_Length(curmodulename,6);
     strcpy(newvar->var->v_subroutinename,subroutinename);
     Save_Length(subroutinename,11);
     strcpy(newvar->var->v_commoninfile,mainfile);
     Save_Length(mainfile,10);

     newvar->var->v_dimension=d;
     /* Creation of the string for the dimension of this variable             */
     dimsempty = 1;

     if ( d )
     {
        newvar->var->v_dimensiongiven=1;
        dims = d;
        while (dims)
        {
           if ( strcasecmp(dims->dim.first,"") || strcasecmp(dims->dim.last,""))
                                                                  dimsempty = 0;
           sprintf(ligne,"%s:%s",dims->dim.first,dims->dim.last);
           strcat(listdimension,ligne);
           if ( dims->suiv )
           {
              strcat(listdimension,",");
           }
           dims = dims->suiv;
        }
        if ( dimsempty == 1 ) newvar->var->v_dimsempty=1;
     }

/*     strcpy(newvar->var->v_readedlistdimension,listdimension);
     Save_Length(listdimension,15);*/
     /*                                                                       */
     newvar->suiv = NULL;

     if ( !List_Save_Var )
     {
        List_Save_Var  = newvar ;
     }
     else
     {
        newvar->suiv = List_Save_Var;
        List_Save_Var = newvar;
     }
  }
}

void Add_Save_Var_dcl_1 (listvar *var)
{
  listvar *newvar;
  listvar *parcours;

  if ( firstpass == 1 )
  {
     parcours = var;
     while ( parcours )
     {
        newvar=(listvar *)malloc(sizeof(listvar));
        newvar->var=(variable *)malloc(sizeof(variable));
        /*                                                                    */
        Init_Variable(newvar->var);
        /*                                                                    */
        newvar->var->v_save=1;
        strcpy(newvar->var->v_nomvar,parcours->var->v_nomvar);
        strcpy(newvar->var->v_modulename,curmodulename);
        Save_Length(curmodulename,6);
        strcpy(newvar->var->v_subroutinename,subroutinename);
        Save_Length(subroutinename,11);
        strcpy(newvar->var->v_commoninfile,mainfile);
        Save_Length(mainfile,10);
        /*                                                                    */
        strcpy(newvar->var->v_readedlistdimension,
             parcours->var->v_readedlistdimension);
        newvar->var->v_nbdim = parcours->var->v_nbdim;
        newvar->var->v_dimension = parcours->var->v_dimension;
        /*                                                                    */
        newvar->var->v_dimensiongiven=parcours->var->v_dimensiongiven;
        /*                                                                    */
        newvar->suiv = NULL;

        if ( !List_Save_Var ) List_Save_Var  = newvar ;
        else
        {
           newvar->suiv = List_Save_Var;
           List_Save_Var = newvar;
        }
        parcours = parcours->suiv;
     }
  }
}
