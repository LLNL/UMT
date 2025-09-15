//--------------------------------------------------------------------------//
// TetonNDAccessor.hh
//
// This class provides an accessor for multidimensional data that can be
// expressed as mcarrays in Blueprint.
//
//--------------------------------------------------------------------------//

#ifndef __TETON_NDACCESSOR_HH__
#define __TETON_NDACCESSOR_HH__

#include "conduit/conduit_node.hpp"
#include <sstream>
#include <stdexcept>

namespace Teton
{

namespace utilities
{

// Store a named dimension
struct NDDimension
{
   std::string name;
   conduit::index_t size{0};
};

/**
 \brief Iterate over dimensions such that each mcarray component is visited once.
 \param dims The dimensions.
 \param func A function to call on each relevant index.
 */
template <typename Functor> void iterate_dimensions(const std::vector<NDDimension> &dims, Functor &&func)
{
   std::vector<conduit::index_t> idx(dims.size(), 0);

   if (dims.size() == 1)
   {
      func(0, idx);
   }
   else if (dims.size() == 2)
   {
      for (conduit::index_t i = 0; i < dims[1].size; i++)
      {
         idx[1] = i;
         func(i, idx);
      }
   }
   else if (dims.size() == 3)
   {
      conduit::index_t c = 0;
      for (conduit::index_t j = 0; j < dims[2].size; j++)
      {
         for (conduit::index_t i = 0; i < dims[1].size; i++, c++)
         {
            idx[1] = i;
            idx[2] = j;
            func(c, idx);
         }
      }
   }
}

/**
 \brief This class creates mcarrays for multidimensional data.
 */
class NDAccessor
{
  public:
   /**
    \brief Constructor
    \param v The "values" node for the Blueprint field.
    \param d Dimension information.
    \param inter Whether the data are interleaved (x0,y0,z0,x1,y1,z1,...)
    */
   NDAccessor(conduit::Node &v, const NDDimension &d, bool inter = false) : values(v), dims(), interleave(inter)
   {
      dims.push_back(d);
   }

   /**
    \brief Constructor
    \param v The "values" node for the Blueprint field.
    \param d A vector of dimension information.
    \param inter Whether the data are interleaved (x0,y0,z0,x1,y1,z1,...)
    */
   NDAccessor(conduit::Node &v, const std::vector<NDDimension> &d, bool inter = false)
      : values(v),
        dims(d),
        interleave(inter)
   {
   }

   /**
    \brief Constructor
    \param v The "values" node for the Blueprint field.
    \param d An initialization vector of dimension information.
    \param inter Whether the data are interleaved (x0,y0,z0,x1,y1,z1,...)
    */
   NDAccessor(conduit::Node &v, const std::initializer_list<NDDimension> &d, bool inter = false)
      : values(v),
        dims(d),
        interleave(inter)
   {
   }

   /**
    \brief Return the number of dimensions.
    \return The number of dimensions.
    */
   size_t number_of_dimensions() const
   {
      return dims.size();
   }

   /**
    \brief Return the total number of elements for all components.
    \return The number of elements for all components.
    */
   conduit::index_t total_number_of_elements() const
   {
      conduit::index_t sz = 1;
      for (const auto &d : dims)
         sz *= static_cast<conduit::index_t>(d.size);
      return sz;
   }

   /**
    \brief Return the number of array elements in the data's fastest dimension.
    \return The number of array elements.
    */
   conduit::index_t number_of_elements() const
   {
      return (values.number_of_children() > 0) ? values[0].dtype().number_of_elements()
                                               : values.dtype().number_of_elements();
   }

   /**
    \brief Associate a data pointer with the Conduit node(s). If there are multiple
           conventions or dimensions then mcarray components will be created.

    \param data A pointer to the buffer that contains the data.
    */
   template <typename T> void set_external(const T *data)
   {
      // TODO: This could be more general to handle more dimensions. We only need up to 3D for Teton.
      auto data_ptr = const_cast<T *>(data);
      constexpr auto elemSizeBytes = sizeof(T);
      conduit::index_t size = dims[0].size;
      if (dims.size() == 1)
      {
         values.set_external(data_ptr, size);
      }
      else
      {
         // Compute stride for the last dimensions for the interleaved case.
         // We intentionally omit the first dimension.
         conduit::index_t stride = 1;
         for (size_t i = 1; i < dims.size(); i++)
            stride *= dims[i].size;

         const auto _this = this;
         iterate_dimensions(dims,
                            [&](conduit::index_t component, const std::vector<conduit::index_t> &idx)
         {
            conduit::Node &comp = values[_this->index_to_name(idx)];
            if (interleave)
            {
               // Data are interleaved [comp0][comp1]...[comp0][comp1]...
               comp.set_external(data_ptr, size, component * elemSizeBytes, stride * elemSizeBytes);
            }
            else
            {
               // Data are contiguous [all comp 0][all comp 1]...
               comp.set_external(data_ptr, size, component * (size * elemSizeBytes));
            }
         });
      }
   }

   /**
    \brief Return the data for the specified array index.

    \param idx The index of the element whose data will be returned.

    \return The data at the index \a idx.
    */
   double operator()(const std::vector<conduit::index_t> &idx) const
   {
      double retval{};
#ifndef _NDEBUG
      if (idx.size() != dims.size())
      {
         throw std::invalid_argument("idx,dim size mismatch");
      }
#endif
      if (values.number_of_children() > 0)
      {
         // Fetch the right mcarray component and get the right value out.
         // TODO: this node lookup is a potential performance issue.
         conduit::Node &comp = values.fetch_existing(index_to_name(idx));
         auto acc = comp.as_double_accessor();
         retval = acc[idx[0]];
      }
      else
      {
         // All data are in a single array. Make an index.
         auto size_factor = [&](int dim)
         {
            conduit::index_t s = 1;
            for (int i = 0; i < dim; i++)
               s *= dims[i].size;
            return s;
         };
         conduit::index_t index = 0;
         if (interleave)
         {
            // 2D
            index = idx[0] * dims[1].size + idx[1];
         }
         else
         {
            for (size_t i = 0; i < idx.size(); i++)
               index += idx[i] * size_factor(i);
         }
         auto acc = values.as_double_accessor();
#ifndef _NDEBUG
         if (index >= acc.number_of_elements())
         {
            std::stringstream ss;
            ss << "Out of bounds index " << index << " in " << values.path()
               << ". number_of_elements=" << acc.number_of_elements() << " idx={";
            for (const auto ival : idx)
               ss << ival << ", ";
            ss << "}\n";
            throw std::range_error(ss.str());
         }
#endif
         retval = acc[index];
      }
      return retval;
   }

   /**
    \brief Turn the data, however it is organized, into a linear array in \a dest.
           Data are copied so all values from component 0 are contiguous, followed
           by component 1, and so on (double[ncomp][nelem]).

    \param[out] The destination array that contains the contiguous data.
    */
   void to_contiguous(double *dest) const
   {
      if (values.number_of_children() > 0)
      {
         double *dptr = dest;
         const auto _this = this;
         iterate_dimensions(dims,
                            [&](conduit::index_t /*c*/, const std::vector<conduit::index_t> &idx)
         {
            const conduit::Node &comp = values.fetch_existing(_this->index_to_name(idx));
            auto acc = comp.as_double_accessor();
            auto n = acc.number_of_elements();
            for (conduit::index_t i = 0; i < n; i++)
               *dptr++ = acc[i];
         });
      }
      else
      {
         auto acc = values.as_double_accessor();
         auto n = acc.number_of_elements();
         for (conduit::index_t i = 0; i < n; i++)
            dest[i] = acc[i];
      }
   }

   /**
    \brief Print the dimensions to a stream.
    \param os The stream to use for printing.
    */
   void print(std::ostream &os) const
   {
      os << "{";
      for (const auto &d : dims)
         os << "{" << d.name << ", " << d.size << "},";
      os << "}";
   }

  private:
   /**
    \brief Produces an mcarray name given a tuple of dimension indices.
    \param idx A tuple of indices within the set of dimensions provided to
               the object. The first dimension is assumed to vary the fastest
               so it is ignored.
    \return The name of the mcarray component.
    */
   std::string index_to_name(const std::vector<conduit::index_t> &idx) const
   {
      std::stringstream ss;
      // Dimension 0 skipped on purpose.
      for (size_t i = 1; i < dims.size(); i++)
      {
         if (i > 1)
            ss << "_";
         ss << dims[i].name << idx[i];
      }
      return ss.str();
   }

  private:
   conduit::Node &values;
   std::vector<NDDimension> dims;
   bool interleave{false};
};

} // namespace utilities

} // namespace Teton
#endif
