#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_vec.h>

T8_EXTERN_C_BEGIN ();

struct t8_adapt_data
{
  double              midpoint[3];
  double              radius;
  double              with_ring;
};

int
t8_adapt_callback_refine (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c *ts,
                          const int is_family,
                          const int num_elements, t8_element_t *elements[])
{
  const struct t8_adapt_data *adapt_data =
    (const struct t8_adapt_data *) t8_forest_get_user_data (forest);

  double              centroid[3];
  double              dist;

  T8_ASSERT (adapt_data != NULL);

  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  dist = t8_vec_dist (adapt_data->midpoint, centroid);

  if (dist < adapt_data->radius + adapt_data->with_ring) {
    return 1;
  }
  return 0;
}

int
t8_adapt_callback_remove (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c *ts,
                          const int is_family,
                          const int num_elements, t8_element_t *elements[])
{
  const struct t8_adapt_data *adapt_data =
    (const struct t8_adapt_data *) t8_forest_get_user_data (forest);

  double              centroid[3];
  double              dist;

  T8_ASSERT (adapt_data != NULL);

  t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);
  dist = t8_vec_dist (adapt_data->midpoint, centroid);

  if (dist < adapt_data->radius) {
    return -2;
  }
  return 0;
}

int
t8_adapt_callback_coarse (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c *ts,
                          const int is_family,
                          const int num_elements, t8_element_t *elements[])
{
  /* coarse every possible element */
  if (is_family == 1) {
    return -1;
  }
  return 0;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_MPI_Comm         comm;
  t8_cmesh_t          cmesh;
  t8_forest_t         forest;

  const int           level = 3;

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  comm = sc_MPI_COMM_WORLD;

  // Build cmesh and uniform forest.
  cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
  forest =
    t8_forest_new_uniform (cmesh, t8_scheme_new_default_cxx (), level, 0,
                           comm);
  T8_ASSERT (t8_forest_is_committed (forest));

  struct t8_adapt_data adapt_data = { {0.5, 0.5, 0}, 0.3, 0.1 };

  t8_forest_write_vtk (forest, "t8_example_0_base");

  t8_debugf ("[IL] REFINE \n");
  forest =
    t8_forest_new_adapt (forest, t8_adapt_callback_refine, 0, 0, &adapt_data);
  t8_forest_write_vtk (forest, "t8_example_1_refine");

  t8_debugf ("[IL] REMOVE \n");
  forest =
    t8_forest_new_adapt (forest, t8_adapt_callback_remove, 0, 0, &adapt_data);
  t8_forest_write_vtk (forest, "t8_example_2_remove");

  t8_debugf ("[IL] COARSE \n");
  forest =
    t8_forest_new_adapt (forest, t8_adapt_callback_coarse, 0, 0, &adapt_data);
  t8_forest_write_vtk (forest, "t8_example_3_coarse");

  t8_forest_unref (&forest);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}

T8_EXTERN_C_END ();
