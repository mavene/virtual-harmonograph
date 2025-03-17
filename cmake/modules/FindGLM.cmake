# FindGLM - attempts to locate the glm matrix/vector library.
#
# This module defines the following variables (on success):
# GLM_INCLUDE_DIRS - where to find glm/glm.hpp
# GLM_FOUND - if the library was successfully located
#
# It is trying a few standard installation locations, but can be customized
# with the following variables:
# GLM_ROOT - root directory of a glm installation
# Headers are expected to be found in either:
# <GLM_ROOT>/glm/glm.hpp OR
# <GLM_ROOT>/include/glm/glm.hpp
# This variable can either be a cmake or environment
# variable. Note however that changing the value
# of the environment varible will NOT result in
# re-running the header search and therefore NOT
# adjust the variables set by this module.
#=============================================================================
# Copyright 2012 Carsten Neumann
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
# License text for the above reference.)
# default search dirs

set(_glm_HEADER_SEARCH_DIRS
"/usr/include"
"/usr/local/include"
"${CMAKE_SOURCE_DIR}/include")

# Check environment for root search directory
set( _glm_ENV_ROOT "$ENV{GLM_ROOT}" )
if( NOT GLM_ROOT AND _glm_ENV_ROOT )
	set(GLM_ROOT "${_glm_ENV_ROOT}")
endif() #NOT GLM_ROOT_DIR AND _glm_ENV_ROOT

# Put user specified location at beginning of search
if( GLM_ROOT )
	set(_glm_HEADER_SEARCH_DIRS "${GLM_ROOT}"
	"${GLM_ROOT}/include"
	${_glm_HEADER_SEARCH_DIRS})
endif() # GLM_ROOT

# Search for the header
FIND_PATH(GLM_INCLUDE_DIR "glm/glm.hpp"
PATHS ${_glm_HEADER_SEARCH_DIRS} )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GLM DEFAULT_MSG
GLM_INCLUDE_DIR)