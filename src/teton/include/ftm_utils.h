!-----------------------------------------------------------------------------
! Defines utility macros used by the Fortran Templating Macros
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Utility macros for concatenating up to four macro parameters with specified
! delimiter
!-----------------------------------------------------------------------------

#if defined(TETON_STRICT_FPP_MODE)

#define FTM_PASTE(a) a
#define FTM_ADD_DELIM(a) a ## _
#define FTM_APPEND(a,b) a ## b

#define FTM_CAT2_INTER(a,b) a ## b
#define FTM_CAT2(a,b) FTM_CAT2_INTER(FTM_ADD_DELIM(a),b)

#define FTM_CAT3_INTER(a,b,c) a ## b ## c
#define FTM_CAT3(a,b,c) FTM_CAT3_INTER(FTM_APPEND(a,_),FTM_APPEND(b,_),c)

#define FTM_CAT4_INTER(a,b,c,d) a ## b ## c ## d
#define FTM_CAT4(a,b,c,d) FTM_CAT4_INTER(FTM_APPEND(a,_), FTM_APPEND(b,_),FTM_APPEND(c,_),d)

#elif defined(__PGI)

#define FTM_PASTE(a)a
#define FTM_ADD_DELIM(a)FTM_PASTE(a)##_
#define FTM_APPEND(a,b)a##b
#define FTM_CAT2_INTER(a,b)a##b
#define FTM_CAT2(a,b)FTM_CAT2_INTER(FTM_ADD_DELIM(a),b)
#define FTM_CAT3_INTER(a,b,c)a##b##c
#define FTM_CAT3(a,b,c)FTM_CAT3_INTER(FTM_APPEND(a,_),FTM_APPEND(b,_),c)
#define FTM_CAT4_INTER(a,b,c,d)a##b##c##d
#define FTM_CAT4(a,b,c,d)FTM_CAT4_INTER(FTM_APPEND(a,_),FTM_APPEND(b,_),FTM_APPEND(c,_),d)

#else

#define FTM_PASTE(a)a
#define FTM_ADD_DELIM(a)FTM_PASTE(a)_
#define FTM_CAT2(a,b)FTM_ADD_DELIM(a)b
#define FTM_CAT3(a,b,c)FTM_ADD_DELIM(a)FTM_ADD_DELIM(b)c
#define FTM_CAT4(a,b,c,d)FTM_ADD_DELIM(a)FTM_ADD_DELIM(b)FTM_ADD_DELIM(c)d

#endif
