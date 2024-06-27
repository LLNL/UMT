#include "TetonTesting.hh"
#include "TetonUtilities.hh"

#include "conduit/conduit_blueprint.hpp"
#include "conduit/conduit_blueprint_mesh.hpp"
#include "conduit/conduit_relay.hpp"
#include "conduit/conduit_relay_mpi.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

namespace Teton
{

namespace testing
{

namespace detail
{

//-----------------------------------------------------------------------------
std::vector<std::string> split(const std::string &str)
{
   std::vector<std::string> retval;
   std::istringstream f(str);
   std::string s;
   while (getline(f, s, '/'))
   {
      retval.push_back(s);
   }
   return retval;
}

//-----------------------------------------------------------------------------
std::vector<std::string> filter(const std::vector<std::string> &words, const std::vector<std::string> &exclude)
{
   std::vector<std::string> retval;
   for (const auto &w : words)
   {
      if (std::find(exclude.begin(), exclude.end(), w) == exclude.end())
         retval.push_back(w);
   }
   return retval;
}

//-----------------------------------------------------------------------------
std::string join(const std::vector<std::string> &words, const std::string &delim)
{
   std::string retval;
   for (size_t i = 0; i < words.size(); i++)
   {
      if (i > 0)
         retval += delim;
      retval += words[i];
   }
   return retval;
}

//-----------------------------------------------------------------------------
void diff_info_to_html_helper(std::ostream &os,
                              const conduit::Node &baseline,
                              const conduit::Node &current,
                              const conduit::Node &info)
{
   conduit::Node opts;
   opts["num_elements_threshold"] = 50;
   opts["num_children_threshold"] = 10000;

   const std::vector<std::string> exclude{"children", "diff"};
   for (conduit::index_t i = 0; i < info.number_of_children(); i++)
   {
      const conduit::Node &n = info[i];
      if (n.name() == "errors")
      {
         auto path = join(filter(split(n.parent()->path()), exclude), "/");

         os << "<tr>\n";
         os << "<td>" << path << "</td>\n";

         os << "<td>";
         if (baseline.has_path(path))
         {
            const conduit::Node &obj = baseline.fetch_existing(path);
            obj.to_summary_string_stream(os, opts);

            if (obj.dtype().number_of_elements() > 1)
               os << "<br><i>len=" << obj.dtype().number_of_elements() << "</i>";
         }
         os << "</td>\n";

         os << "<td>";
         if (current.has_path(path))
         {
            const conduit::Node &obj = current.fetch_existing(path);
            obj.to_summary_string_stream(os, opts);

            if (obj.dtype().number_of_elements() > 1)
               os << "<br><i>len=" << obj.dtype().number_of_elements() << "</i>";
         }
         os << "</td>\n";

         os << "</tr>\n";
      }
      else if (n.number_of_children() > 0)
      {
         diff_info_to_html_helper(os, baseline, current, n);
      }
#if 0
      // Enable this if we want to see nodes that are the same.
      else
      {
         auto path = join(filter(split(n.parent()->path()), exclude), "/");

         os << "<tr>\n";
         os << "<td>" << path << "</td>\n";

         os << "<td colspan=\"2\">match</td>\n";

         os << "</tr>\n";
      }
#endif
   }
}

//---------------------------------------------------------------------------
void diff_info_to_html(MPI_Comm comm,
                       std::ostream &os,
                       const conduit::Node &baseline,
                       const conduit::Node &current,
                       const conduit::Node &info,
                       int cycle)
{
   // Table style
   const char *style = R"(
<style>
#diffs {
  font-family: monospace;
  border-collapse: collapse;
  width: 100%;
}

#diffs td, #diffs th {
  border: 1px solid #ddd;
  padding: 8px;
}

#diffs tr:nth-child(even){background-color: #f2f2f2;}

#diffs tr:hover {background-color: #ddd;}

#diffs th {
  padding-top: 12px;
  padding-bottom: 12px;
  text-align: left;
  background-color: #04AA6D;
  color: white;
}

</style>
)";

   int rank = 0;
   MPI_Comm_rank(comm, &rank);

   // Write real diffs out as a table.
   os << style << "\n";
   if (cycle > 0)
   {
      std::stringstream prev;
      prev << "diff." << rank << "." << (cycle - 1) << ".html";
      os << "<a href=\"" << prev.str() << "\">Prev</a>&nbsp;&nbsp;&nbsp;&nbsp;";
   }
   std::stringstream next;
   next << "diff." << rank << "." << (cycle + 1) << ".html";
   os << "<a href=\"" << next.str() << "\">Next</a><br>\n";

   os << "<p>Cycle " << cycle << "</p><br>\n";
   os << "<table id=\"diffs\">\n";
   os << "<tr><td><b>Path</b></td><td><b>Baseline</b></td><td><b>Current</b></td></tr>\n";
   diff_info_to_html_helper(os, baseline, current, info);
   os << "</table>\n";
}

//-----------------------------------------------------------------------------
/**
 @return True if the node is the same as the baseline node.
 */
bool compare_baseline(MPI_Comm comm,
                      const conduit::Node &baseline,
                      const conduit::Node &current,
                      conduit::Node &info,
                      int cycle,
                      bool forceFile = true)
{
   const double tolerance = 1.e-6;

   int rank = 0;
   MPI_Comm_rank(comm, &rank);

   bool equal;
   if (forceFile)
   {
      // Sometimes Node::diff lies about data arrays being equal when they are not.
      // Save the results node n out to a file first and then read it back in.
      // That seems to work around the issue.
      std::stringstream ss;
      ss << "conduit_tmp_result" << rank << "." << cycle << ".yaml";
      std::string tmp_file(ss.str());
      conduit::relay::io::save(current, tmp_file, "yaml");
      conduit::Node n_fromfile;
      conduit::relay::io::load(tmp_file, "yaml", n_fromfile);
      conduit::utils::remove_file(tmp_file);
      // Node::diff returns true if the nodes are different. We want not different.
      equal = !baseline.diff(n_fromfile, info, tolerance, true);
   }
   else
   {
      // Node::diff returns true if the nodes are different. We want not different.
      equal = !baseline.diff(current, info, tolerance, true);
   }

   return equal;
}

//-----------------------------------------------------------------------------
void write_info(MPI_Comm comm,
                const conduit::Node &baseline,
                const conduit::Node &current,
                const conduit::Node &info,
                int cycle)
{
   int rank = 0;
   MPI_Comm_rank(comm, &rank);

   conduit::Node opts;
   opts["num_elements_threshold"] = 20;
   opts["num_children_threshold"] = 10000;

   // Write the info as YAML.
   std::stringstream ss;
   ss << "diff." << rank << ".yaml";
   std::ofstream f;
   f.open(ss.str().c_str());
   info.to_summary_string_stream(f, opts);
   f.close();

   // Write the important info as HTML.
   std::stringstream ss2;
   ss2 << "diff." << rank << "." << cycle << ".html";
   std::ofstream html;
   html.open(ss2.str().c_str());
   html << "<html>\n";
   html << "  <head>\n";
   html << "    <title>diff</title>\n";
   html << "  </head>\n";
   html << "  <body bgcolor=\"white\">\n";
   diff_info_to_html(comm, html, baseline, current, info, cycle);
   html << "  </body>\n";
   html << "</html>\n";
   html.close();
}

} // namespace detail

bool test(const conduit::Node &n, const std::string &fileBase, int cycle, bool make, MPI_Comm comm)
{
   int rank = 0;
   MPI_Comm_rank(comm, &rank);

   if (n.has_child("fields"))
   {
      // Test that the fields contain "good" values.
      const conduit::Node &fields = n["fields"];
      for (conduit::index_t i = 0; i < fields.number_of_children(); i++)
      {
         const conduit::Node &values = fields[i]["values"];
         if (values.number_of_children() > 0)
         {
            for (conduit::index_t comp = 0; comp < values.number_of_children(); comp++)
            {
               utilities::scan_field_values(rank, values[comp]);
            }
         }
         else
         {
            utilities::scan_field_values(rank, values);
         }
      }
   }

   // Make baseline filename.
   std::string baselineFilename;
   if (getenv("TETON_TESTING_BASELINE_DIR") != nullptr)
   {
      std::vector<std::string> s{getenv("TETON_TESTING_BASELINE_DIR"), "baseline" + fileBase};
      baselineFilename = detail::join(s, "/");
   }
   else
   {
      baselineFilename = "baseline" + fileBase;
   }

   bool rv = true;
   if (make)
      conduit::relay::io::save(n, baselineFilename, "yaml");
   else
   {
      // Make current filename.
      std::string filename;
      if (getenv("TETON_TESTING_CURRENT_DIR") != nullptr)
      {
         std::vector<std::string> s{getenv("TETON_TESTING_CURRENT_DIR"), "current" + fileBase};
         filename = detail::join(s, "/");
      }
      else
      {
         filename = "current" + fileBase;
      }

      // Save the current node.
      conduit::relay::io::save(n, filename, "yaml");

      // Load the baseline node.
      conduit::Node baseline;
      conduit::relay::io::load(baselineFilename, "yaml", baseline);

      // Compare current with the baseline and populate diff info.
      conduit::Node info;
      rv = detail::compare_baseline(comm, baseline, n, info, cycle);

      // If there were differences, write the info.
      if (!rv)
         detail::write_info(comm, baseline, n, info, cycle);
   }
   return rv;
}

} // namespace testing

} // namespace Teton
