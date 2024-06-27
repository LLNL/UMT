#include <cmath>

#include "TetonInterface.hh"
#include "TetonSurfaceTallies.hh"
#include "dbc_macros.h"

#include "conduit/conduit_relay.hpp"

namespace TetonSurfaceTallies
{

void dumpTallyToJson(const conduit::Node &blueprint, const conduit::Node &options, int mpi_rank)
{
   // The tally definition is split into two places.
   // The SURFACE information is in blueprint.
   // The other details of the tally (shape, groups, frame, etc.) live in options.
   if (blueprint.has_path("teton/surface_edits"))
   {
      conduit::Node tally_node;
      // Update these numbers if the format changes:
      tally_node["fileFormat/version"] = "0.6.0";
      tally_node["fileFormat/major"] = 0;
      tally_node["fileFormat/minor"] = 6;
      tally_node["fileFormat/patch"] = 0;
      tally_node["runInfo/code/name"] = "teton"; // TODO read in host code name?
      std::string conduit_base_path = "xray/surfaceTallies/";

      const double *group_bounds = options.fetch_existing("quadrature/gnu").as_double_ptr();

      int nanglebin_teton = -1;
      teton_getnumanglebins(&nanglebin_teton);

      std::string tally_file_name = "radenergyflux.tal";
      if (options.has_path("surface_edits_output_filename"))
      {
         tally_file_name = options.fetch_existing("surface_edits_output_filename").as_string();
      }

      const conduit::Node &surface_edit_options_all = options["surface_edits"];
      const conduit::Node &surface_edit_blueprint_all = blueprint["teton/surface_edits"];
      conduit::NodeConstIterator surface_edit_blueprint_it = surface_edit_blueprint_all.children();
      // Loop over edit surfaces and dump
      while (surface_edit_blueprint_it.has_next())
      {
         const conduit::Node &surface_info = surface_edit_blueprint_it.next(); // surface info
         std::string tally_name = surface_info.name();
         const conduit::Node &surface_edit_option = surface_edit_options_all[tally_name]; // options for tallying

         std::string conduit_tally_path = conduit_base_path + tally_name;

         // Assemble the input for writing to file:
         const bool labFrame = surface_edit_option["transform_to_lab_frame"].as_int();
         const bool tShift = surface_edit_option["apply_time_shift"].as_int();
         const bool compInc = surface_edit_option.fetch_existing("calculate_incident").as_int();
         const bool dumpErr = surface_edit_option.fetch_existing("calculate_error_metrics").as_int();

         // if not integrating over angles, then this is equal to the number of polar levels
         const bool integrate_over_angles = surface_edit_option["integrate_over_angles"].as_int();
         int nanglebin = 1;
         if (!integrate_over_angles)
            teton_getnumanglebins(&nanglebin);

         // if not integrating over all groups, then this is equal to the number of groups
         const bool integrate_over_all_groups = surface_edit_option["integrate_over_all_groups"].as_int();
         int ngrp = 1;
         if (!integrate_over_all_groups)
            ngrp = options["quadrature/num_groups"].as_int();

         const double *tbinbnds = surface_edit_option["time_bin_boundaries"].as_double_ptr();
         const int ntimebin = surface_edit_option["time_bin_boundaries"].dtype().number_of_elements() - 1;
         const int size3d = nanglebin * ngrp * ntimebin;

         const double *tally_values_esc = blueprint["fields/" + tally_name + "_tallies/values"].as_double_ptr();
         const double *tally_values_inc = nullptr;
         if (compInc)
         {
            tally_values_inc = blueprint["fields/" + tally_name + "_tallies_incident/values"].as_double_ptr();
         }

         tally_node[conduit_tally_path + "/info/labFrame"] = labFrame;
         tally_node[conduit_tally_path + "/info/isPointSource"] = tShift;
         // All our tallies are "piecewise constant" and "extensive"
         tally_node[conduit_tally_path + "/info/interpolation"] = "Piecewise_Constant";
         tally_node[conduit_tally_path + "/info/intensive"] = false;

         const double scale_tally = surface_edit_option.fetch_existing("scale_tally").as_double();
         std::string energy_unit = "jerk";
         if (fabs(scale_tally - 1.e4) < 1.e-6)
         {
            energy_unit = "Terg";
         }
         else if (fabs(scale_tally - 1.) > 1.e-10)
         {
            energy_unit = std::to_string(scale_tally) + " jerks";
         }
         // tally_node[conduit_tally_path+"/instantaneousEnergy/units"] = energy_unit;
         // tally_node[conduit_tally_path+"/finalCumulativeEnergy/units"] = energy_unit;
         // tally_node.print();
         tally_node[conduit_tally_path + "/instantaneousEnergy/units"] = "Terg";
         tally_node[conduit_tally_path + "/finalCumulativeEnergy/units"] = "Terg";
         tally_node.print();

         int idim = 0;
         int idimfinal = 0;
         std::vector<int> shape;
         std::vector<int> shape_final;
         conduit::Node time_dim_info;
         if (ntimebin > 1)
         {
            time_dim_info["units"] = "sh";
            time_dim_info["type"] = "time step";
            time_dim_info["bounds"].set(tbinbnds, ntimebin + 1);
            time_dim_info["stride"] = ngrp * nanglebin;
            tally_node[conduit_tally_path + "/instantaneousEnergy/dim" + std::to_string(idim)] = time_dim_info;
            // The finalCumulativeEnergy tallies are integrated over time, so there is no time dimension.
            shape.push_back(ntimebin);
            idim++;
         }
         if (nanglebin > 1)
         {
            conduit::Node angle_dim_info;
            angle_dim_info["units"] = "cos(theta)";
            angle_dim_info["type"] = "angle bin";

            TETON_VERIFY_C(mpi_rank,
                           nanglebin == nanglebin_teton,
                           "nanglebin must match nanglebin_teton for surface tallies");

            std::vector<double> abinbnds(nanglebin + 1, 0.);
            teton_getanglebins(&nanglebin, abinbnds.data());
            angle_dim_info["bounds"].set(abinbnds.data(), nanglebin + 1);

            angle_dim_info["stride"] = ngrp;
            tally_node[conduit_tally_path + "/instantaneousEnergy/dim" + std::to_string(idim)] = angle_dim_info;
            tally_node[conduit_tally_path + "/finalCumulativeEnergy/dim" + std::to_string(idimfinal)] = angle_dim_info;
            shape.push_back(nanglebin);
            shape_final.push_back(nanglebin);
            idim++;
            idimfinal++;
         }
         if (ngrp > 1)
         {
            conduit::Node grp_dim_info;
            grp_dim_info["units"] = "keV";
            grp_dim_info["type"] = "group bound";
            grp_dim_info["bounds"].set(group_bounds, ngrp + 1);
            grp_dim_info["stride"] = 1;
            tally_node[conduit_tally_path + "/instantaneousEnergy/dim" + std::to_string(idim)] = grp_dim_info;
            tally_node[conduit_tally_path + "/finalCumulativeEnergy/dim" + std::to_string(idimfinal)] = grp_dim_info;
            shape.push_back(ngrp);
            shape_final.push_back(ngrp);
            idim++;
            idimfinal++;
         }
         tally_node[conduit_tally_path + "/instantaneousEnergy/numDimension"] = idim;
         tally_node[conduit_tally_path + "/finalCumulativeEnergy/numDimension"] = idimfinal;
         if (idim > 0)
         {
            tally_node[conduit_tally_path + "/instantaneousEnergy/shape"].set(shape.data(), idim);
            // idimfinal should always be less than or equal to idim:
            if (idimfinal > 0)
            {
               tally_node[conduit_tally_path + "/finalCumulativeEnergy/shape"].set(shape_final.data(), idimfinal);
            }
         }

         if (compInc)
         {
            // Everything is the same for the incident tallies except the data itself:
            tally_node[conduit_tally_path + "/instantaneousEnergyIncident"].set(
               tally_node[conduit_tally_path + "/instantaneousEnergy"]);
            tally_node[conduit_tally_path + "/finalCumulativeEnergyIncident"].set(
               tally_node[conduit_tally_path + "/finalCumulativeEnergy"]);

            tally_node[conduit_tally_path + "/instantaneousEnergyIncident/data"].set_external_node(
               blueprint["fields/" + tally_name + "_tallies_incident/values"]);

            // Recompute the data sums for the incident tallies:
            double datasum_incident = 0.;
            for (int i = 0; i < size3d; i++)
            {
               datasum_incident += tally_values_inc[i];
            }
            tally_node[conduit_tally_path + "/instantaneousEnergyIncident/sumOfData"] = datasum_incident;
            tally_node[conduit_tally_path + "/finalCumulativeEnergyIncident/sumOfData"] = datasum_incident;
         }

         // Set the exiting tally data down here so that the above lines of code
         //    don't make an unnecessary copy of it into the incident tally storage
         tally_node[conduit_tally_path + "/instantaneousEnergy/data"].set_external_node(
            blueprint["fields/" + tally_name + "_tallies/values"]);
         double datasum = 0.;
         for (int i = 0; i < size3d; i++)
         {
            datasum += tally_values_esc[i];
         }
         tally_node[conduit_tally_path + "/instantaneousEnergy/sumOfData"] = datasum;
         tally_node[conduit_tally_path + "/finalCumulativeEnergy/sumOfData"] = datasum;

         // Integrate over time to get the final cumulative exiting tally:
         const int ngrpang = ngrp * nanglebin;
         {
            std::vector<double> cumulativeEnergy(ngrpang, 0.);
            int igrpang = 0;
            for (int i = 0; i < size3d; i++)
            {
               cumulativeEnergy[igrpang] += tally_values_esc[i];
               igrpang++;
               if (igrpang == ngrpang)
               {
                  igrpang = 0;
               }
            }
            tally_node[conduit_tally_path + "/finalCumulativeEnergy/data"].set(cumulativeEnergy.data(), ngrpang);
         }
         // Integrate over time to get the final cumulative incident tally:
         if (compInc)
         {
            std::vector<double> cumulativeEnergy(ngrpang, 0.);
            int igrpang = 0;
            for (int i = 0; i < size3d; i++)
            {
               cumulativeEnergy[igrpang] += tally_values_inc[i];
               igrpang++;
               if (igrpang == ngrpang)
               {
                  igrpang = 0;
               }
            }
            tally_node[conduit_tally_path + "/finalCumulativeEnergyIncident/data"].set(cumulativeEnergy.data(),
                                                                                       ngrpang);
         }

         if (tShift && dumpErr)
         {
            // First convert raw error estimates to normalized error estimates:
            std::vector<double> integratedTally(ntimebin, 0.);

            // Get the normalization factors (tallies integrated over group and angle)
            const int stride = nanglebin * ngrp;
            int index1d;
            for (int itimebin = 0; itimebin < ntimebin; itimebin++)
            {
               integratedTally[itimebin] += 1.e-12; // Avoid dividing by zero
               for (int istride = 0; istride < stride; istride++)
               {
                  index1d = itimebin * stride + istride;
                  integratedTally[itimebin] += tally_values_esc[index1d];
               } // end loop over group and angle indices
            }    // end loop over time bin index

            std::vector<double> normalizedErrEstShift(ntimebin, 0.);
            std::vector<double> normalizedErrEstSrcSize(ntimebin, 0.);
            const double *error_est_shift = blueprint["fields/" + tally_name + "_error_est_shift/values"]
                                               .as_double_ptr();
            const double *error_est_src_size = blueprint["fields/" + tally_name + "_error_est_src_size/values"]
                                                  .as_double_ptr();
            for (int itimebin = 0; itimebin < ntimebin; itimebin++)
            {
               normalizedErrEstShift[itimebin] = error_est_shift[itimebin] / integratedTally[itimebin];
               normalizedErrEstSrcSize[itimebin] = error_est_src_size[itimebin] / integratedTally[itimebin];
            }

            // Note that we only store the exiting error estimates
            tally_node[conduit_tally_path + "/displacementParallel/units"] = "cm";
            tally_node[conduit_tally_path + "/displacementParallel/numDimension"] = (ntimebin > 1) ? 1 : 0;
            if (ntimebin > 1)
            {
               int shape_displacement[1] = {ntimebin};
               tally_node[conduit_tally_path + "/displacementParallel/shape"].set(shape_displacement, 1);
               tally_node[conduit_tally_path + "/displacementParallel/dim0"] = time_dim_info;
            }
            // Dimensions and units are the same for both error estimates:
            tally_node[conduit_tally_path + "/displacementPerp"].set(
               tally_node[conduit_tally_path + "/displacementParallel"]);

            tally_node[conduit_tally_path + "/displacementPerp/data"].set(normalizedErrEstShift.data(), ntimebin);
            tally_node[conduit_tally_path + "/displacementParallel/data"].set(normalizedErrEstSrcSize.data(), ntimebin);

            double sum_perp = 0.;
            double sum_parallel = 0.;
            for (int itimebin = 0; itimebin < ntimebin; itimebin++)
            {
               sum_perp += normalizedErrEstShift[itimebin];
               sum_parallel += normalizedErrEstSrcSize[itimebin];
            }
            tally_node[conduit_tally_path + "/displacementPerp/sumOfData"] = sum_perp;
            tally_node[conduit_tally_path + "/displacementParallel/sumOfData"] = sum_parallel;
         }
      }

      conduit::relay::io::save(tally_node, tally_file_name, "json");
   }

} // end dumpTallyToJson()

} // end namespace TetonSurfaceTallies
