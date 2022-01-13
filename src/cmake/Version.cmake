
include(cmake/GetGitRevisionDescription.cmake)

get_git_head_revision(GIT_REFSPEC GIT_SHA1 ALLOW_LOOKING_ABOVE_CMAKE_SOURCE_DIR)

set( TETON_GIT_SHA ${GIT_SHA1} )
