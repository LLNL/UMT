# This file generates repetative access functions for the iteration controls.
# Maybe this indicates a bad design.


template = """!***********************************************************************
!   Created:  08/2021, TAB
!
!   GetDefault{upname}
!
!   Get default tolerance controls, for users to query before Teton setup
!   in case they need it for their parsers, etc.
!
!***********************************************************************

   subroutine GetDefault{upname}(value) &
                        BIND(C,NAME="teton_get_default_{lowname}_internal")

   USE ISO_C_BINDING
   use default_iter_controls_mod
   use Size_mod

   implicit none

   ! Input
   {vartype}, intent(out) :: value

   value = {lowname}
   return
   end subroutine GetDefault{upname}
"""

names = [ "outer_temp_reltol",
       "outer_intensity_reltol",
       "grey_reltol",
       "incident_flux_reltol",
       "inner_nl_reltol",
       "outer_max_it",
       "grey_max_it",
       "incident_flux_max_it",
       "inner_nl_max_it",
       ]

def to_camel_case(snake_str):
   components = snake_str.split('_')
   # We capitalize the first letter of each component except the first one
   # with the 'title' method and join them together.
   return components[0] + ''.join(x.title() for x in components[1:])

def to_all_camel_case(snake_str):
   components = snake_str.split('_')
   # We capitalize the first letter of each component except the first one
   # with the 'title' method and join them together.
   return ''.join(x.title() for x in components)

for fname in names:
   if "_reltol" in fname:
      vartype = "real(C_DOUBLE)"
      cvar = "double"
   else:
      vartype = "integer(C_INT)"
      cvar = "int"

   cc = to_all_camel_case(fname).replace("Reltol","RelTol")
   fd = {'upname': cc,
         'lowname': fname,
         'vartype': vartype,
         'cvar': cvar,
         }

   f = open("GetDefault{upname}.F90".format(**fd),'w')
   f.write(template.format(**fd))
   f.close()
   print "  inline {cvar} teton_get_default_{lowname}(){{\n    {cvar} var{{}};\n    teton_get_default_{lowname}_internal(&var);\n    return var;\n  }}".format(**fd)
   print "  void teton_get_default_{lowname}_internal(const {cvar}* var);\n".format(**fd)
