/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

/*
 * This example demonstrates the use of the t8code library to generate a p4est
 * forest from a native cmesh. The cmesh is generated from a simple geometry and
 * then refined, partitioned and balanced. The forest is then output to vtu
 * files.
 *
 * Usage: t8_cmesh_p4est_examples <configuration> <level>
 *        possible configurations:
 *        o drop      Refinement on a 5-trees geometry with an inner hole.
 *        o bowtie    Refinement on a 2-tree bowtie domain.
 */

#include <t8.h>                                     /* General t8code header, always include this. */
#include <t8_cmesh.h>                               /* cmesh definition and basic interface. */
#include <t8_cmesh/t8_cmesh_examples.h>             /* A collection of exemplary cmeshes, including the p4est wrapper*/
#include <t8_forest/t8_forest_general.h>            /* forest definition and general interface. */
#include <t8_forest/t8_forest_io.h>                 /* forest io interface. */
#include <t8_forest/t8_forest.c>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h> /* linear geometry of the cmesh */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */
#include <t8_vec.h>                                 /* Basic operations on 3D vectors. */
#include <t8_forest/t8_forest_cxx.cxx>
#include <t8_cmesh_vtk_writer.h>
#include <sc_shmem.h>

/* Define the configuration of the geometry */
typedef enum
{
  P4EST_CONFIG_NULL,
  P4EST_CONFIG_DROP,
  P4EST_CONFIG_BOWTIE,

} simple_config_t;

/* Define the adapt data structure */
struct adapt_data
{
  int level_max;
};

/* Define the adapt callback function. If the tree is not at the maximum level
and it is the first tree, then it will be refined. */
int
adapt_callback (t8_forest_t forest, t8_forest_t forest_from,
              t8_locidx_t which_tree, t8_locidx_t lelement_id,
                         t8_eclass_scheme_c *ts, const int is_family,
                         const int num_elements, t8_element_t *elements[])
{
  const int *adapt_data = (const int *) t8_forest_get_user_data (forest);
  const int level_max = adapt_data[0];
  const int level = ts->t8_element_level (elements[0]);
  if (which_tree == 0 && level < level_max)
  {
    return 1;
  };

  return 0;
}

/* Define the adapt function */
t8_forest_t
adapt_forest (t8_forest_t forest, int max_level)
{
  t8_forest_t forest_adapt;
  struct adapt_data adapt_data = {
    max_level
  };
  T8_ASSERT (t8_forest_is_committed (forest));
  forest_adapt = t8_forest_new_adapt (forest, adapt_callback, 1, 0, &adapt_data);
  return forest_adapt;
}

/* Define the partition function */
t8_forest_t
partition_forest (t8_forest_t forest)
{
  t8_forest_t forest_partition;
  T8_ASSERT (t8_forest_is_committed (forest));
  t8_forest_init (&forest_partition);
  t8_forest_set_partition (forest_partition, forest, 0);
  t8_forest_commit (forest_partition);

  return forest_partition;
}

/* Define the balance function */
t8_forest_t
balance_forest (t8_forest_t forest)
{
  t8_forest_t forest_balance;
  T8_ASSERT (t8_forest_is_committed (forest));
  t8_forest_init (&forest_balance);
  t8_forest_set_balance (forest_balance, forest, 0);
  t8_forest_set_ghost (forest_balance, 1, T8_GHOST_FACES);
  t8_forest_commit (forest_balance);

  return forest_balance;
}

/* Define the cmesh generation function for the bowtie geometry */
t8_cmesh_t
t8_bowtie (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;

  double vertices[24] = {
    -1.4142,  0,      0,
    -0.7071, -0.7071, 0,
    -0.7071,  0.7071, 0,
     0,       0,      0,
     0.7071, -0.7071, 0,
     1.4142,  0,      0,
     0,       0,      0,
     0.7071,  0.7071, 0,
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_geometry_c *linear_geom = t8_geometry_linear_new (2);
  t8_cmesh_register_geometry (cmesh, linear_geom); /* Use linear geometry */

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_QUAD);

  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 4);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 12, 4);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/* Define the cmesh generation function for the drop geometry */
t8_cmesh_t
t8_drop (sc_MPI_Comm comm)
{
  t8_cmesh_t cmesh;

  double vertices[120] = {
    0, 0, 0,
    1, 0, 0,
    0, 1, 0,
    1, 1, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1,

    1, 1, 0,
    3, 1, 0,
    1, 2, 0,
    3, 2, 0,
    1, 1, 1,
    2, 1, 1,
    1, 2, 1,
    2, 2, 1,

    2, 1, 1,
    3, 1, 0,
    2, 2, 1,
    3, 2, 0,
    2, 1, 2,
    3, 1, 3,
    2, 2, 2,
    3, 2, 3,

    1, 1, 2,
    2, 1, 2,
    1, 2, 2,
    2, 2, 2,
    0, 1, 3,
    3, 1, 3,
    0, 2, 3,
    3, 2, 3,

    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 1,
    0, 0, 3,
    1, 0, 2,
    0, 1, 3,
    1, 1, 2,
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_geometry_c *linear_geom = t8_geometry_linear_new (5);
  t8_cmesh_register_geometry (cmesh, linear_geom); /* Use linear geometry */

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_HEX);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_HEX);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_HEX);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_HEX);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_HEX);

  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 8);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 24, 8);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 48, 8);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 72, 8);
  t8_cmesh_set_tree_vertices (cmesh, 4, vertices + 96, 8);

  t8_cmesh_set_join (cmesh, 0, 4, 5, 4, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 4, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 5, 1, 0);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
}

/* Main function */
int
main (int argc, char **argv)
{
  int                  mpiret;
  int                  wrongusage;
  const char           *usage;
  simple_config_t      config;
  t8_scheme_cxx_t      *scheme;
  t8_cmesh_t             cmesh;
  int                  init_level, refine_level;
  int                  global_num_trees;

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);
  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);
  /* We will use MPI_COMM_WORLD as a communicator. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  usage =
  "Arguments: <configuration> <level>\n"
  "   Configuration can be any of\n"
  "   Level controls the maximum depth of refinement\n";

  /* The prefix for our output files. */
  const char *prefix = argv[1];

  /* Check the command line arguments. */
  wrongusage = 0;
  config = P4EST_CONFIG_NULL;
  if (!wrongusage && argc < 3) {
    wrongusage = 1;
  }
  if (!wrongusage) {
    if (!strcmp (argv[1], "drop")) {
      config = P4EST_CONFIG_DROP;
    }
    else if (!strcmp (argv[1], "bowtie")) {
      config = P4EST_CONFIG_BOWTIE;
    }
    else {
      wrongusage = 1;
    }
  }

  if (wrongusage) {
    t8_global_errorf (usage);
    sc_abort_collective ("Usage error");
  }

  /* Set the initial and refinement level */
  init_level = 2;
  refine_level = atoi(argv[2]);

  /* create p4est connectivity and forest structures. These connectivies
  are legacies with no native built-in cmesh generation functions*/
  if (config == P4EST_CONFIG_BOWTIE) {
    cmesh = t8_bowtie (comm);
  }
  else if (config == P4EST_CONFIG_DROP) {
    cmesh = t8_drop (comm);
  }
  else{
    sc_abort_collective ("Usage error");
  }

  /* Output the meshes to vtu files. */
  scheme = t8_scheme_new_default_cxx ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, init_level, 0, comm);

  /* Compute local and global number of trees. */
  global_num_trees = t8_cmesh_get_num_trees (cmesh);
  t8_global_productionf (" Global number of trees:\t%i\n", global_num_trees);

  /* adapt, partition and balance the forest */
  forest = adapt_forest (forest, refine_level);
  forest = partition_forest(forest);
  forest = balance_forest(forest);

  /* Output the forest to vtu files. */
  t8_forest_write_vtk (forest, prefix);
  t8_global_productionf (" Generating vtk files.\n");
  /*
   * Clean-up
   */
  /* Deallocate the cmeshes */

  t8_cmesh_destroy (&cmesh);
  t8_global_productionf (" Deallocate cmeshes.\n");

  t8_forest_unref (&forest);
  t8_global_productionf (" Destroyed forests.\n");

  /* Finalize the sc library */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
