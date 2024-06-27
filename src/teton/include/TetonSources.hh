#ifndef TETON_SOURCES_HH__
#define TETON_SOURCES_HH__

#include "TetonInterface.hh"
#include <vector>

class TetonSource
{
  public:
   TetonSource(int nangles, int nsrczones, int ngroups)
      : m_num_angles(nangles),
        m_num_srczones(nsrczones),
        m_num_groups(ngroups),
        m_src_values(),
        m_zone_list(),
        m_time(-1.)
   {
      if (m_num_srczones > 0)
      {
         m_src_values.resize(nangles * nsrczones * ngroups);
         m_zone_list.resize(nsrczones);
      }
   }

   virtual ~TetonSource() = default;

   virtual void UpdateSourceValues(double time, double dtrad) = 0;

   const double *GetSourceValues() const
   {
      return m_src_values.data();
   }
   const int *GetZoneList() const
   {
      return m_zone_list.data();
   }

   int GetNumAngles() const
   {
      return m_num_angles;
   }
   int GetNumZones() const
   {
      return m_num_srczones;
   }
   int GetNumGroups() const
   {
      return m_num_groups;
   }

  protected:
   const int m_num_angles;
   const int m_num_srczones;
   const int m_num_groups;

   // nangles x nsrczones x ngroups array with the source values with this time step
   //    groups varies fastest, angles varies slowest
   //    Units should be energy per time
   std::vector<double> m_src_values;

   // length nsrczones list of Fortran indices corresponding to the zones for m_src_values
   // (indices should be between 1 and the number of zones in the local domain, inclusive)
   std::vector<int> m_zone_list;

   // Current value of the time
   double m_time;
};

class TetonSourceManager
{
  public:
   TetonSourceManager(int rank = -1) : mRank(rank), m_source_list()
   {
   }

   void AddSource(TetonSource *source)
   {
      m_source_list.push_back(source);
   }

   void AddPointSourceFromTally(int teton_zone_index,
                                std::string tally_file_name,
                                std::string tally_name,
                                double multiplier = 1.0);

   void UpdateSources(double time, double dtrad)
   {
      for (auto source : m_source_list)
      {
         source->UpdateSourceValues(time, dtrad);
      }
   }

   void UpdatePsiWithSources() const
   {
      for (auto source : m_source_list)
      {
         const int num_zones = source->GetNumZones();
         if (num_zones < 1)
            continue; // this rank has no source zones
         const int num_angles = source->GetNumAngles();
         const int num_groups = source->GetNumGroups();
         // Some call to Teton:
         teton_appendsourcetopsi(&num_zones,
                                 &num_angles,
                                 &num_groups,
                                 source->GetZoneList(),
                                 source->GetSourceValues());
      }
   }

   void SetRank(int rank)
   {
      mRank = rank;
   }

   ~TetonSourceManager()
   {
      for (auto source : m_source_list)
      {
         delete source;
      }
      m_source_list.clear();
   }

  protected:
   int mRank;
   std::vector<TetonSource *> m_source_list;
};

class PointSource : public TetonSource
{
   // Class for an extensive point source
  public:
   // timevals is of size ntimes
   // source_profile is of size ntimes x ngroups x nangles
   //
   // @note The zone_index value should be -1 if the zone is not owned by the current MPI rank.
   PointSource(int nangles,
               int ngroups,
               int zone_index, // Which zone is the point source in? TODO, convert from coordinate to zone index
               int ntimebins,
               const double *time_bin_bounds,
               const double *source_profile); // units of energy

   virtual void UpdateSourceValues(double time, double dtrad);

  protected:
   // TODO, do we want to make copies of these, or just let these be const pointers?
   std::vector<double> m_profile;         // size m_num_times x m_num_groups x m_num_angles
   std::vector<double> m_time_bin_bounds; // length m_num_time_bins+1
   int m_num_time_bins;                   //  Number of time bins
};

#endif // TETON_SOURCES_HH__
