/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

  Copyright (C) 2023 the developers

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

/* Show-case several cmesh examples without curvilinear geometries. */

/*
 * Usage: t8_cmesh_nogeo <configuration> <level>
 *        possible configurations:
 *        o hypercube           Constructs a hypercube forest
 *        o hypercube_pad       Creates a hypercube mesh with custom boundaries
 *        o hypercube_hybrid    Generates a mixed-element mesh consisting of
 *                              tets, prisms, and hexes, forming a hypercube
 *        o periodic_tri        Generates a unit square mesh of two triangles,
 *                              with periodicity in both x and y directions.
 *        o periodic_hybrid     Creates a unit square mesh combining quads and
 *                              triangles, featuring periodicity in x and y
 *                              directions.
 *        o periodic_line_more_trees  Constructs a line mesh with three trees,
 *                              introducing periodicity to create a continuous
 *                              line effect.
 *        o bigmesh             Allows for the creation of a large mesh
 *                              consisting of a specified number of trees of the
 *                              same type, tailored to the given element class.
 */

#include <t8.h>                                     /* General t8code header, always include this. */
#include <t8_cmesh.h>                               /* Cmesh definition and basic interface. */
#include <t8_forest/t8_forest_general.h>            /* Forest definition and basic interface. */
#include <t8_forest/t8_forest_geometrical.h>        /* Forest-related geometry operations. */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* Default refinement scheme. */
#include <t8_cmesh_vtk_writer.h>                    /* Write file in vtu file */
#include <t8_forest/t8_forest_io.h>
#include <t8_cmesh/t8_cmesh_examples.h>

/* The following enum is used to select the configuration of the mesh. */
typedef enum
{
  T8_CONFIG_NULL,
  T8_CONFIG_HYPERCUBE,
  T8_CONFIG_HYPERCUBE_PAD,
  T8_CONFIG_HYPERCUBE_HYBRID,
  T8_CONFIG_PERIODIC_TRI,
  T8_CONFIG_PERIODIC_HYBRID,
  T8_CONFIG_PERIODIC_LINE_MORE_TREES,
  T8_CONFIG_LAST
} t8_cmesh_config_t;

static void
t8_write_forest_to_vtu (t8_forest_t forest, const char *prefix)
{
  const t8_locidx_t num_elements = t8_forest_get_local_num_elements (forest);

  /* We need to allocate a new array to store the data on their own.
   * These arrays have one entry per local element. */
  double *diameters = T8_ALLOC (double, num_elements);

  /* The number of user defined data fields to write. */
  const int num_data = 1;

  /* For each user defined data field we need one t8_vtk_data_field_t variable. */
  t8_vtk_data_field_t vtk_data[num_data];
  /* Set the type of this variable. Since we have one value per element, we pick T8_VTK_SCALAR. */
  vtk_data[0].type = T8_VTK_SCALAR;
  /* The name of the field as should be written to the file. */
  strcpy (vtk_data[0].description, "diameter");
  vtk_data[0].data = diameters;

  /* Get the number of trees that have elements of this process. */
  const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees (forest);

  /* Loop over all local trees in the forest. */
  for (t8_locidx_t itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
    const t8_locidx_t num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);

    /* Loop over all local elements in the tree and compute diameter estimate. */
    for (t8_locidx_t ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
      const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);
      diameters[current_index] = t8_forest_element_diam (forest, itree, element);
    }
  }

  {
    /* Write user defined data to vtu file. */
    const int write_treeid = 1;
    const int write_mpirank = 1;
    const int write_level = 1;
    const int write_element_id = 1;
    const int write_ghosts = 0;
    t8_forest_write_vtk_ext (forest, prefix, write_treeid, write_mpirank, write_level, write_element_id, write_ghosts,
                             0, 0, num_data, vtk_data);
  }

  T8_FREE (diameters);
}

int
main (int argc, char **argv)
{
  const char         *usage;
  int                 wrongusage;
  t8_cmesh_config_t     config;
  t8_cmesh_t cmesh;

  /*
   * Initialization.
   */

  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  int mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);

  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_PRODUCTION);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);

  /* We will use MPI_COMM_WORLD as a communicator. */
  sc_MPI_Comm comm = sc_MPI_COMM_WORLD;

  /* process command line arguments */
  usage =
    "Arguments: <configuration>\n"
    "   Configuration can be any of\n"
    "      hypercube|hypercube_pad|hypercube_hybrid|periodic_tri|\n"
    "         periodic_hybrid|periodic_line_more_trees\n"
    "   Level controls the maximum depth of refinement\n";
  wrongusage = 0;
  config = T8_CONFIG_NULL;
  if (!wrongusage && argc != 2) {
    wrongusage = 1;
  }
  if (!wrongusage) {
    if (!strcmp (argv[1], "hypercube")) {
      config = T8_CONFIG_HYPERCUBE;
    }
    else if (!strcmp (argv[1], "hypercube_pad")) {
      config = T8_CONFIG_HYPERCUBE_PAD;
    }
    else if (!strcmp (argv[1], "hypercube_hybrid")) {
      config = T8_CONFIG_HYPERCUBE_HYBRID;
    }
    else if (!strcmp (argv[1], "periodic_tri")) {
      config = T8_CONFIG_PERIODIC_TRI;
    }
    else if (!strcmp (argv[1], "periodic_hybrid")) {
      config = T8_CONFIG_PERIODIC_HYBRID;
    }
    else if (!strcmp (argv[1], "periodic_line_more_trees")) {
      config = T8_CONFIG_PERIODIC_LINE_MORE_TREES;
    }
    else {
      wrongusage = 1;
    }
  }
  if (wrongusage) {
    t8_global_errorf("Usage error\n");
  }

  if (config == T8_CONFIG_HYPERCUBE) {
    /*
     * Creation of a hypercube mesh and storing it to disk.
     */
    cmesh = t8_cmesh_new_hypercube(T8_ECLASS_HEX, sc_MPI_COMM_WORLD, 0, 0, 0);
  }
  else if (config == T8_CONFIG_HYPERCUBE_PAD) {
    /*
     * Creation of a hypercube mesh with custom boundaries and storing it to disk.
     */
    const double boundary_coords[24] = { 1, 0, 0, 4, 0, 0, 0, 6, 0, 5, 5, 0, -1, -2, 8, 9, 0, 10, 0, 8, 9, 10, 10, 10 };
    cmesh = t8_cmesh_new_hypercube_pad (T8_ECLASS_HEX, sc_MPI_COMM_WORLD, boundary_coords, 3, 3, 3, 0);
  }
  else if (config == T8_CONFIG_HYPERCUBE_HYBRID) {
    /*
     * Creation of a mixed-element mesh consisting of tets, prisms, and hexes, forming a hypercube and storing it to disk.
     */
    cmesh = t8_cmesh_new_hypercube_hybrid (sc_MPI_COMM_WORLD, 0, 0);
  }
  else if (config == T8_CONFIG_PERIODIC_TRI) {
    /*
     * Creation of a unit square mesh of two triangles, with periodicity in both x and y directions and storing it to disk.
     */
    cmesh = t8_cmesh_new_periodic_tri (sc_MPI_COMM_WORLD);
  }
  else if (config == T8_CONFIG_PERIODIC_HYBRID) {
    /*
     * Creation of a unit square mesh combining quads and triangles, featuring periodicity in x and y directions and storing it to disk.
     */
    cmesh = t8_cmesh_new_periodic_hybrid (sc_MPI_COMM_WORLD);
  }
  else if (config == T8_CONFIG_PERIODIC_LINE_MORE_TREES) {
    /*
     * Creation of a line mesh with three trees, introducing periodicity to create a continuous line effect and storing it to disk.
     */
    cmesh = t8_cmesh_new_periodic_line_more_trees (sc_MPI_COMM_WORLD);
  }

  const int uniform_level = 5;
  t8_forest_t forest = t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), uniform_level, 0, comm);

  /*
   * Creation of several meshes and storing them to disk.
   */

  const char *prefix_cmesh = "t8_non_geo_cmesh";
  const char *prefix_forest = "t8_non_geo_forest";

  t8_cmesh_vtk_write_file (cmesh, prefix_cmesh);
  t8_global_productionf ("Wrote t8_%s_cmesh.\n", argv[1]);

  t8_write_forest_to_vtu (forest, prefix_forest);
  t8_global_productionf ("Wrote t8_%s_forest.\n\n", argv[1]);

  t8_forest_unref (&forest);

  t8_global_productionf ("Done!\n");

  /* Finalize the sc library */
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
