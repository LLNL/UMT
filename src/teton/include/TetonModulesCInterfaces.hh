//--------------------------------------------------------------------------//
// TetonModulesCInterfaces.hh
//--------------------------------------------------------------------------//

#pragma once

// C versions of Teton Fortran modules interfaces
// These are for internal use only.
extern "C"
{
conduit::Node *teton_get_datastore_cptr();

#if !defined(TETON_ENABLE_MINIAPP_BUILD)
conduit::Node *teton_conduitcheckpoint_get_cptr();
void teton_conduitcheckpoint_prep_for_load();
void teton_conduitcheckpoint_data_loaded();
void teton_conduitcheckpoint_external_data_loaded();
void teton_conduitcheckpoint_prep_for_save();
void teton_conduitcheckpoint_teardown();
#endif

// Wrapped functions in mods/QuadratureList_mod.F90
void *teton_quadraturelist_getquadlist();
int teton_quadraturelist_getnumberofcommsets(void *);
int teton_quadraturelist_getnumberofsets(void *);
int teton_quadraturelist_getnumberofgtasets(void *);
int teton_quadraturelist_getnumberofanglesets(void *);
int teton_quadraturelist_getnumberofgroupsets(void *);
int teton_quadraturelist_getnumberofzonesets(void *);
int teton_quadraturelist_getnumberofhyperdomains(void *, int);

void *teton_quadraturelist_getquad(void *, int);
int teton_quadrature_getnumberofenergygroups(void *);
int teton_quadrature_getnumberofangles(void *);

void *teton_size_getmeshsize();
int teton_size_getnumberofzones(void *);
int teton_size_getnumberofcorners(void *);
int teton_size_getnumberofcommneighbors(void *);
}
