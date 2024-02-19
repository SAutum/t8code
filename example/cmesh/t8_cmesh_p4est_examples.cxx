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
 * This example demonstrates how to create a forest from a p4est connectivity
 * and refine it to a certain level. The forest is then adapted, partitioned,
 * and balanced. The forest is then written to vtk files.
 *
 * Usage: t8_cmesh_p4est_examples <configuration> <level>
 *        possible configurations:
 *        o unit      Refinement on the unit square.
 *        o brick     Refinement on a regular forest of octrees.
 *        o three     Refinement on a forest with three trees.
 *        o evil      Check second round of refinement with np=5 level=7
 *        o evil3     Check second round of refinement on three trees
 *        o pillow    Refinement on a 2-tree pillow-shaped domain.
 *        o moebius   Refinement on a 5-tree Moebius band.
 *        o star      Refinement on a 6-tree star shaped domain.
 *        o cubed     Refinement on a 6-tree cubed sphere surface.
 *        o disk      Refinement on a 5-tree spherical standard disk.
 *        o xdisk     Refinement on a 5-tree spherical disk periodic in x.
 *        o ydisk     Refinement on a 5-tree spherical disk periodic in y.
 *        o pdisk     Refinement on a 5-tree spherical disk, periodic b.c.
 *        o periodic  Refinement on the unit square with all-periodic b.c.
 *        o rotwrap   Refinement on the unit square with weird periodic b.c.
 *        o circle    Refinement on a 6-tree donut-like circle.
 *        o icosahedron   Refinement on the icosahedron sphere with geometry.
 *        o shell2d       Refinement on a 2d shell with geometry.
 *        o disk2d        Refinement on a 2d disk with geometry.
 *        o sphere2d      Refinement on a 6-tree sphere surface with geometry.
 */

#include <t8.h>                                     /* General t8code header, always include this. */
#include <t8_cmesh.h>                               /* cmesh definition and basic interface. */
#include <t8_forest/t8_forest_general.h>            /* forest definition and basic interface. */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */
#include <t8_cmesh_vtk_writer.h>                    /* write file in vtu file */
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h> /* linear geometry of the cmesh */
#include <t8_forest/t8_forest.c>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <p4est_connectivity.h>

typedef enum
{
  P4EST_CONFIG_NULL,
  P4EST_CONFIG_UNIT,
  P4EST_CONFIG_BRICK,
  P4EST_CONFIG_THREE,
  P4EST_CONFIG_EVIL,
  P4EST_CONFIG_EVIL3,
  P4EST_CONFIG_PILLOW,
  P4EST_CONFIG_MOEBIUS,
  P4EST_CONFIG_STAR,
  P4EST_CONFIG_CUBED,
  P4EST_CONFIG_DISK,
  P4EST_CONFIG_XDISK,
  P4EST_CONFIG_YDISK,
  P4EST_CONFIG_PDISK,
  P4EST_CONFIG_PERIODIC,
  P4EST_CONFIG_ROTWRAP,
  P4EST_CONFIG_ICOSAHEDRON,
  P4EST_CONFIG_SHELL2D,
  P4EST_CONFIG_DISK2D,
  P4EST_CONFIG_SPHERE2D,
  P4EST_CONFIG_LAST
}

/* Define adapt data */
simple_config_t;
struct adapt_data
{
  int level_max;
};

/* Define adapt callback */
int
adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                         t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
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

/* Define adapt forest. The callback function is used to adapt the forest */
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

/* Define partition forest */
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

/* Define balance forest */
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

/* Main function */
int
main (int argc, char **argv)
{
  int                  mpiret;
  int                  wrongusage;
  const char           *usage;
  p4est_connectivity_t *connectivity;
  simple_config_t      config;
  t8_scheme_cxx_t      *scheme;
  t8_cmesh_t             cmesh;
  int                  init_level, refine_level;
  int                 nbrick_x = 1, nbrick_y = 1;


  /*
   * Initialization.
   */

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
  /* Initialize p4est library. if not initialized, p4est will be initialized
  automatically by t8_forest_new_from_p4est. However, an error will be thrown
  during finalization. */
  p4est_init (NULL, SC_LP_ESSENTIAL);

  usage =
  "Arguments: <configuration> <level>\n"
  "   Configuration can be any of\n"
  "   Level controls the maximum depth of refinement\n";

  /* The prefix for our output files. */
  const char *prefix = argv[1];

  wrongusage = 0;
  config = P4EST_CONFIG_NULL;
  if (!wrongusage && argc < 3) {
    wrongusage = 1;
  }
  if (!wrongusage) {
    if (!strcmp (argv[1], "unit")) {
      config = P4EST_CONFIG_UNIT;
    }
    else if (!strcmp (argv[1], "brick")) {
      config = P4EST_CONFIG_BRICK;
    }
    else if (!strcmp (argv[1], "three")) {
      config = P4EST_CONFIG_THREE;
    }
    else if (!strcmp (argv[1], "evil")) {
      config = P4EST_CONFIG_EVIL;
    }
    else if (!strcmp (argv[1], "evil3")) {
      config = P4EST_CONFIG_EVIL3;
    }
    else if (!strcmp (argv[1], "pillow")) {
      config = P4EST_CONFIG_PILLOW;
    }
    else if (!strcmp (argv[1], "moebius")) {
      config = P4EST_CONFIG_MOEBIUS;
    }
    else if (!strcmp (argv[1], "star")) {
      config = P4EST_CONFIG_STAR;
    }
    else if (!strcmp (argv[1], "cubed")) {
      config = P4EST_CONFIG_CUBED;
    }
    else if (!strcmp (argv[1], "disk")) {
      config = P4EST_CONFIG_DISK;
    }
    else if (!strcmp (argv[1], "xdisk")) {
      config = P4EST_CONFIG_XDISK;
    }
    else if (!strcmp (argv[1], "ydisk")) {
      config = P4EST_CONFIG_YDISK;
    }
    else if (!strcmp (argv[1], "pdisk")) {
      config = P4EST_CONFIG_PDISK;
    }
    else if (!strcmp (argv[1], "periodic")) {
      config = P4EST_CONFIG_PERIODIC;
    }
    else if (!strcmp (argv[1], "rotwrap")) {
      config = P4EST_CONFIG_ROTWRAP;
    }
    else if (!strcmp (argv[1], "icosahedron")) {
      config = P4EST_CONFIG_ICOSAHEDRON;
    }
    else if (!strcmp (argv[1], "shell2d")) {
      config = P4EST_CONFIG_SHELL2D;
    }
    else if (!strcmp (argv[1], "disk2d")) {
      config = P4EST_CONFIG_DISK2D;
    }
    else if (!strcmp (argv[1], "sphere2d")) {
      config = P4EST_CONFIG_SPHERE2D;
    }
    else {
      wrongusage = 1;
    }
  }

  if (wrongusage) {
    t8_global_errorf (usage);
    sc_abort_collective ("Usage error");
  }

  init_level = 3;
  refine_level = atoi(argv[2]);

  /* create p4est connectivity and forest structures. These connectivies
  are legacies with no native built-in cmesh generation functions*/
  if (config == P4EST_CONFIG_BRICK) {
    nbrick_x = argc > 3 ? atoi (argv[3]) : 3;
    nbrick_y = argc > 4 ? atoi (argv[4]) : 2;
    connectivity = p4est_connectivity_new_brick (nbrick_x, nbrick_y, 0, 0);
  }
  else if (config == P4EST_CONFIG_THREE || config == P4EST_CONFIG_EVIL3) {
    connectivity = p4est_connectivity_new_corner ();
  }
  else if (config == P4EST_CONFIG_PILLOW) {
    connectivity = p4est_connectivity_new_pillow ();
  }
  else if (config == P4EST_CONFIG_MOEBIUS) {
    connectivity = p4est_connectivity_new_moebius ();
  }
  else if (config == P4EST_CONFIG_STAR) {
    connectivity = p4est_connectivity_new_star ();
  }
  else if (config == P4EST_CONFIG_CUBED) {
    connectivity = p4est_connectivity_new_cubed ();
  }
  else if (config == P4EST_CONFIG_DISK) {
    connectivity = p4est_connectivity_new_disk (0, 0);
  }
  else if (config == P4EST_CONFIG_XDISK) {
    connectivity = p4est_connectivity_new_disk (1, 0);
  }
  else if (config == P4EST_CONFIG_YDISK) {
    connectivity = p4est_connectivity_new_disk (0, 1);
  }
  else if (config == P4EST_CONFIG_PDISK) {
    connectivity = p4est_connectivity_new_disk (1, 1);
  }
  else if (config == P4EST_CONFIG_PERIODIC) {
    connectivity = p4est_connectivity_new_periodic ();
  }
  else if (config == P4EST_CONFIG_ROTWRAP) {
    connectivity = p4est_connectivity_new_rotwrap ();
  }
  else if (config == P4EST_CONFIG_ICOSAHEDRON) {
    connectivity = p4est_connectivity_new_icosahedron ();
  }
  else if (config == P4EST_CONFIG_SHELL2D) {
    connectivity = p4est_connectivity_new_shell2d ();
  }
  else if (config == P4EST_CONFIG_DISK2D) {
    connectivity = p4est_connectivity_new_disk2d ();
  }
  else if (config == P4EST_CONFIG_SPHERE2D) {
    connectivity = p4est_connectivity_new_cubed ();
  }
  else{
    connectivity = p4est_connectivity_new_unitsquare ();
  }

  /* Create the cmesh from the p4est connectivity */
  cmesh = t8_cmesh_new_from_p4est (connectivity, comm, 0);
  p4est_connectivity_destroy (connectivity);

  /* Create the default refinement scheme */
  scheme = t8_scheme_new_default_cxx ();

  /* Create the forest from the cmesh */
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, init_level, 0, comm);

  /* Adapt, partition, and balance the forest */
  forest = adapt_forest (forest, refine_level);
  forest = partition_forest(forest);
  forest = balance_forest(forest);

  /* Write the forest to vtk files */
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
