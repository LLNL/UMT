# This file generates repetative access functions for the iteration controls.
# Maybe this indicates a bad design.


template = """!***********************************************************************
!   Created:  08/2021, TAB
!
!   Adjust{upname}{uptype}
!
!   Sets one of the iteration control values.
!
!***********************************************************************

   subroutine Adjust{upname}{uptype}(value) &
                        BIND(C,NAME="teton_adjust_{lowname}_{lowtype}")

   USE ISO_C_BINDING
   use iter_control_mod
   use iter_control_list_mod
   use Size_mod

   implicit none

   ! Input
   {vartype}, intent(in) :: value

   ! Local variables
   type(IterControl), pointer  :: {control}Control => NULL()

   {control}Control => getIterationControl(IterControls,"{control}")

   call setControls({control}Control, {controlVar}=value)

   return
   end subroutine Adjust{upname}{uptype}
"""

for n, c in [('FluxExchange','incidentFlux'),
      #('Grey','grey'), # DIABLED FOR NOW.  grey was customized for sweep/iteration differences.
      ('Nonlinear','nonLinear'),
      ('RadEnergyDensity','intensity'),
      ('Temperature','temperature'),
   ]:
   for t in ("RelTol","MaxIts"):
      if t == "RelTol":
         vartype = "real(C_DOUBLE)"
         controlVar = "epsilonPoint"
         cvar = "double"
      else:
         if c == "intensity":
            continue
         vartype = "integer(C_INT)"
         controlVar = "maxNumberOfIterations"
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

      f = open("Adjust{upname}{uptype}.F90".format(**fd),'w')
      f.write(template.format(**fd))
      f.close()
      print "void teton_adjust_{lowname}_{lowtype}(const {cvar}* {uptype});".format(**fd)

