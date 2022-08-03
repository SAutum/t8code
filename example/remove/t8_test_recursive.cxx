/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <iostream>
#include <t8.h>
#include <t8_eclass.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_vec.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_forest/t8_forest_private.h>
#include "t8_cmesh/t8_cmesh_testcases.h"
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_types.h>

static int
t8_adapt_callback_remove (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements, 
                          t8_element_t * elements[])
{
  int                 num_children;
  int                 child_id;
  num_children = ts->t8_element_num_children (elements[0]);
  child_id     = ts->t8_element_child_id (elements[0]);
  if (num_children-1 != child_id) {
    return -2;
  }
  return 0;
}

static int
t8_adapt_callback_coarse (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements, 
                          t8_element_t * elements[])
{
  if (is_family) {
    return -1;
  }
  return 0;
}

static int
t8_adapt_callback_refine (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements, 
                          t8_element_t * elements[])
{
  int                 child_id, level;

  level    = ts->t8_element_level (elements[0]);
  child_id = ts->t8_element_child_id (elements[0]);

  if (child_id == 0 && level < 6) {
      return 1;
  }
  return 0;
}

static t8_forest_t
t8_adapt_forest (t8_forest_t forest_from, t8_forest_adapt_t adapt_fn,
                 int recursive)
{
  t8_forest_t         forest_new;

  t8_forest_init (&forest_new);
  t8_forest_set_adapt (forest_new, forest_from, adapt_fn, recursive);
  //t8_forest_set_partition (forest_new, NULL, 0);
  t8_forest_commit (forest_new);

  return forest_new;
}

void
t8_test_elements_remove (int cmesh_id)
{
  int                 min_level;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;
  t8_forest_t         forest_base;
  t8_scheme_cxx_t    *scheme;

  scheme = t8_scheme_new_default_cxx ();

  /* Construct a cmesh */
  //cmesh = t8_test_create_cmesh (cmesh_id);
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, sc_MPI_COMM_WORLD, 0, 0, 0);

  /* Compute the first level, such that no process is empty */
  min_level = t8_forest_min_nonempty_level (cmesh, scheme);

  forest = t8_forest_new_uniform (cmesh, scheme, min_level, 0, sc_MPI_COMM_WORLD);

  t8_debugf("[IL] start refine\n");
  forest = t8_adapt_forest (forest, t8_adapt_callback_refine, 1);
  t8_forest_write_vtk (forest, "/home/ioannis/VBshare/paraview_export/t8_test_rec_refine");

  t8_debugf("[IL] start remove\n"); 
  forest = t8_adapt_forest (forest, t8_adapt_callback_remove, 0);
  t8_forest_write_vtk (forest, "/home/ioannis/VBshare/paraview_export/t8_test_rec_remove");

  t8_debugf("[IL] start coarse\n");
  forest = t8_adapt_forest (forest, t8_adapt_callback_coarse, 1);
  t8_forest_write_vtk (forest, "/home/ioannis/VBshare/paraview_export/t8_test_rec_coarse");


  t8_debugf("[IL] Recursive done\n");

  SC_CHECK_ABORT (t8_forest_no_overlap (forest),
            "The forest has overlapping elements.");

  t8_cmesh_ref (cmesh);
  t8_scheme_cxx_ref (scheme);
  forest_base = t8_forest_new_uniform (cmesh, scheme, min_level, 0, sc_MPI_COMM_WORLD);
  SC_CHECK_ABORT (t8_forest_is_equal (forest, forest_base),
            "The forest are not equal.");
  
  t8_forest_unref (&forest);
  t8_forest_unref (&forest_base);
}

void
t8_test_cmesh_elements_remove_all ()
{
  int bigmesh_id;
  bigmesh_id = t8_get_number_of_comm_only_cmesh_testcases () +
               t8_get_number_of_new_hypercube_cmesh_testcases () +
               t8_get_number_of_new_empty_cmesh_testcases () +
               t8_get_number_of_new_from_class_cmesh_testcases () +
               t8_get_number_of_new_hypercube_hybrid_cmesh_testcases () +
               t8_get_number_of_new_periodic_cmesh_testcases ();
  /* Test all cmeshes over all different inputs we get through their id */
  for (int cmesh_id = 0; cmesh_id < t8_get_number_of_all_testcases (); cmesh_id++) {
    /* Skip all t8_test_create_new_bigmesh_cmesh since they are without geometry */
    if (cmesh_id < bigmesh_id || 
        cmesh_id >= bigmesh_id + t8_get_number_of_new_bigmesh_cmesh_testcases ()) {
        if (cmesh_id < 66) {
            t8_test_elements_remove (cmesh_id);
        }
    }
  }
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         mpic;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  mpic = sc_MPI_COMM_WORLD;
  sc_init (mpic, 1, 1, NULL, SC_LP_PRODUCTION);
  t8_init (SC_LP_DEFAULT);

  //t8_test_cmesh_elements_remove_all ();
  t8_test_elements_remove (0);
  
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}