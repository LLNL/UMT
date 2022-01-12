Teton's use of macros to facilitate creating overloaded functions
==================================================================

Creating a family of overloaded Function functions or subroutines to support different data types or array ranks can
be tedious in Fortran, as the language does not support templating and this can result in copying essentially the same
function implementation many times.

Teton uses the code preprocessor and a few utility macros dubbed 'Fortran Templating Macros (FTM)' to simulate this
functionality and provide overloaded functions in some of its modules such as:
* `src/teton/mods/MemoryAllocator_mod.F90.templates` has a family of overloaded allocate/deallocate calls
   supporting real and integer arrays with one to three ranks.
* `src/teton/gpu/OMPWrappers.F90.templates` has a family of overloaded subroutines supporting arrays of
   real, int, and derived type with one to three ranks.

Example
----------------------------------

We now show an example of using macros to create two overloaded versions of a function that will accept
a double or integer 1d array.

### Step 1

Place the function implementation in a separate source file, for example `foo.F90.implementation`.  Replace the data type
of the parameter you want to overload on with a `DEFINE` that will later be replaced by the preprocessor.

You will need the preprocessor to create a unique function name for each instance, so one way to accomplish this is to
generate the name using a macro that will append the parameter type onto your function name for each.

A family of `FTM_CAT` concatenation macros is included in `src/teton/include/ftm_utils.h`.  This is an
exmaple of template code you would write.
```fortran
subroutine FTM_CAT( _, foo, FTM_TYPE ) ( var )
   FTM_TYPE, pointer, contiguous, dimension(:), intent(inout) :: var

   ! < body of function >

end subroutine FTM_CAT(_, foo, FTM_TYPE ) ( var )
```

### Step 2

Insert the function body into your source code.  Define the data type for each instance of the function you want to create.
In C++, this would be like explicitly instantiating the templates.

```fortran
#define FTM_TYPE real
#include "MemoryAllocator_mod.F90.templates"
#undef FTM_TYPE

#define FTM_TYPE integer
#include "foo.f90.implementation
#undef FTM_TYPE
```

The preprocessor will substitute the defines and produce code like the following:

```fortran
subroutine foo_real ( var )
   real, pointer, contiguous, dimension(:), intent(inout) :: var

   ! < body of function >

end subroutine foo_real ( var )

subroutine foo_integer ( var )
   integer, pointer, contiguous, dimension(:), intent(inout) :: var

   ! < body of function >

end subroutine foo_integer ( var )
```
