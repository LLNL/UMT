#include "TetonSources.hh"
#include "dbc_macros.h"

#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_relay.hpp"
#include "conduit/conduit_relay_config.h"
#include "conduit/conduit_relay_mpi_io_blueprint.hpp"

void TetonSourceManager::AddPointSourceFromTally(int zone_index,
                                                 std::string tally_file_name,
                                                 std::string tally_name,
                                                 double multiplier)
{
   conduit::Node source_node, source_node_full;
   conduit::relay::io::load_merged(tally_file_name, "json", source_node_full);
   source_node.set_external_node(
      source_node_full.fetch_existing("xray/surfaceTallies/" + tally_name + "/instantaneousEnergy"));
   TETON_VERIFY_C(mRank, source_node["numDimension"].as_int64() == 3, "Replayed point source must have 3D data.");

   const long *shape = source_node["shape"].as_int64_ptr();
   int ntimebins = shape[0];
   int nanglebins = shape[1];
   int ngroups_tally = shape[2];
   int total_size = ntimebins * nanglebins * ngroups_tally;

   // int ngroups = options["quadrature/num_groups"].as_int();
   // int npolar = options["quadrature/npolar"].as_int();
   // TETON_VERIFY_C(mRank,
   //                options["quadrature/qtype"].as_int() == 2,
   //                "Angle-dependent sources require product quadrature (qtype = 2).");
   // TETON_VERIFY_C(mRank,
   //                ngroups_tally == ngroups,
   //                "Cannot interpolate between source and problem group structure yet.");
   // TETON_VERIFY_C(mRank,
   //                nanglebins == 2 * npolar,
   //                "Cannot interpolate between source and problem polar bin boundaries yet.");
   // TETON_VERIFY_C(mRank, nanglebins == 4, "TODO remove hardcoding for 4 angle bins ");

   TETON_VERIFY_C(mRank, source_node["dim0/units"].as_string() == "sh", "Time unit of tally must be in shakes.");

   const double *time_bin_bounds = source_node["dim0/bounds"].as_double_ptr();
   const double *angle_bin_bounds = source_node["dim1/bounds"].as_double_ptr();

   double factor = multiplier;
   double energy_unit_conversion = 1.0;
   if (source_node["units"].as_string() == "Terg")
   {
      energy_unit_conversion = 1.0e-4;
   }
   factor *= energy_unit_conversion;

   // Convert extensive tally (units of energy) to intensive tally (units of energy per time per angle unit)
   const double *source_profile_integrated = source_node["data"].as_double_ptr();
   std::vector<double> source_profile(total_size);
   int counter = 0;
   for (int timebin = 0; timebin < ntimebins; timebin++)
   {
      for (int anglebin = 0; anglebin < nanglebins; anglebin++)
      {
         double angle_bin_width = angle_bin_bounds[anglebin + 1] - angle_bin_bounds[anglebin];
         for (int group = 0; group < ngroups_tally; group++)
         {
            source_profile[counter] = source_profile_integrated[counter] * factor / angle_bin_width;
            counter++;
         }
      }
   }

   m_source_list.push_back(
      new PointSource(nanglebins, ngroups_tally, zone_index, ntimebins, time_bin_bounds, &source_profile[0]));
}

PointSource::PointSource(int nangles,
                         int ngroups,
                         int zone_index,
                         int ntimebins,
                         const double *timebinbounds,
                         const double *source_profile)
   : TetonSource(nangles, zone_index > 0 ? 1 : 0, ngroups),
     m_profile(),
     m_time_bin_bounds(),
     m_num_time_bins(ntimebins)
{
   if (zone_index <= 0)
      return;

   m_zone_list[0] = zone_index;

   int shape_size = nangles * ngroups;
   m_time_bin_bounds.resize(ntimebins + 1);
   m_profile.resize(shape_size * ntimebins);

   // Copy over data:
   int index = 0;
   for (int timebin = 0; timebin < ntimebins; timebin++)
   {
      m_time_bin_bounds[timebin] = timebinbounds[timebin];
      for (int i = 0; i < shape_size; i++)
      {
         m_profile[index] = source_profile[index];
         index++;
      }
   }
   m_time_bin_bounds[ntimebins] = timebinbounds[ntimebins];
}

void PointSource::UpdateSourceValues(double time, double dtrad)
{
   // If this rank has no source zones, return
   if (m_num_srczones < 1)
      return;

   m_time = time;

   int shape_size = m_num_groups * m_num_angles;

   // TODO, something to interpolate between angle bin differences?

   double time_prev = time - dtrad;

   // Zero out source:
   for (int i = 0; i < shape_size; i++)
   {
      m_src_values[i] = 0.;
   }

   // Before the first time value or after the last time value,
   //   source is just zero.
   if (time < m_time_bin_bounds[0] || time_prev > m_time_bin_bounds[m_num_time_bins])
   {
      return;
   }

   // Otherwise, loop through the time bins and accumulate source contributions:
   for (int timebin = 0; timebin < m_num_time_bins; timebin++)
   {
      // as timebin goes up, the time bin in the diagrams below moves to the right
      //      [ source profile's time bin ] -->

      // If time bin lower bound is beyond this time step, we're done
      //                  [ source profile's time bin ]
      //   [ time step ]
      if (m_time_bin_bounds[timebin] > time)
         break;

      // If time bin upper bound < time_prev, keep going
      //  [ source profile's time bin ]
      //                                 [ time step ]
      if (m_time_bin_bounds[timebin + 1] < time_prev)
         continue;

      // If we're here, then the current time bin overlaps with the time step.

      double factor = 0.;

      // else, we want a fraction of the energy from this time bin:
      const double source_time_bin_width = m_time_bin_bounds[timebin + 1] - m_time_bin_bounds[timebin];
      if (m_time_bin_bounds[timebin] < time_prev && m_time_bin_bounds[timebin + 1] > time)
      {
         // the time step falls entirely within the time step.
         //  [ source profile's time bin ]
         //     [ time step ]
         factor = dtrad / source_time_bin_width;
      }
      else if (m_time_bin_bounds[timebin] < time_prev)
      // equivalent to else if (m_time_bin_bounds[timebin] < time_prev && m_time_bin_bounds[timebin+1] <= time)
      {
         //  [ source profile's time bin ]
         //                          [ time step ]
         factor = (m_time_bin_bounds[timebin + 1] - time_prev) / source_time_bin_width;
      }
      else if (m_time_bin_bounds[timebin + 1] > time)
      // equivalent to else if (m_time_bin_bounds[timebin] >= time_prev && m_time_bin_bounds[timebin+1] > time)
      {
         //        [ source profile's time bin ]
         //  [ time step ]
         factor = (time - m_time_bin_bounds[timebin]) / source_time_bin_width;
      }
      else
      // equivalent to else if (m_time_bin_bounds[timebin] >= time_prev && m_time_bin_bounds[timebin+1] <= time)
      {
         // factor = 1 if the source time bin is entirely within the time step
         //       [ source profile's time bin ]
         //   [ time step                                ]
         //   (We want to add the contribution of the entire time bin to this time step's source)
         factor = 1.0;
      }

      if (factor > 1.0)
      {
         std::cerr << "Factor should not ever be larger than 1" << std::endl;
         exit(1);
      }

      const double *data_start = m_profile.data() + timebin * shape_size;
      for (int i = 0; i < shape_size; i++)
      {
         m_src_values[i] += factor * data_start[i];
      }
   }

   // Divide by dtrad to convert from energy to power:
   for (int i = 0; i < shape_size; i++)
   {
      m_src_values[i] /= dtrad;
   }
}
