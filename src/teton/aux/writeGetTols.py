# This file generates repetative access functions for the iteration controls.
# Maybe this indicates a bad design.


template = """!***********************************************************************
!   Created:  08/2021, TAB
!
!   Get{upname}{uptype}
!
!   Gets one of the iteration control values.
!
!***********************************************************************

   subroutine Get{upname}{uptype}(value) &
                        BIND(C,NAME="teton_get_{lowname}_{lowtype}_internal")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   {vartype}, intent(out) :: value

   ! Local variables
   type(IterControl), pointer  :: {control}Control => NULL()

   {control}Control => getIterationControl(IterControls,"{control}")

   value = {controlVar}({control}Control)

   return
   end subroutine Get{upname}{uptype}
"""

for n, c in [('FluxExchange','incidentFlux'),
      ('Grey','grey'),
      ('Nonlinear','nonLinear'),
      ('RadEnergyDensity','intensity'),
      ('Temperature','temperature'),
   ]:
   for t in ("RelTol","MaxIts"):
      if t == "RelTol":
         vartype = "real(C_DOUBLE)"
         controlVar = "getEpsilonPoint"
         cvar = "double"
      else:
         if c == "intensity":
            continue
         vartype = "integer(C_INT)"
         controlVar = "getMaxNumberOfIterations"
         cvar = "int"


      fd = {'upname': n,
            'uptype': t,
            'lowname': n.lower(),
            'lowtype': t.lower(),
            'control': c,
            'controlVar': controlVar,
            'vartype': vartype,
            'cvar': cvar,
            }

      f = open("Get{upname}{uptype}.F90".format(**fd),'w')
      f.write(template.format(**fd))
      f.close()
      print "void teton_get_{lowname}_{lowtype}_internal(const {cvar}* {uptype});".format(**fd)
      print "inline {cvar} teton_get_{lowname}_{lowtype}(){{ {cvar} v{{0}}; teton_get_{lowname}_{lowtype}_internal(&v); return v;}} ".format(**fd)

