%module(package="alps") transform

%include "numpy.i"
%include "transform.i"
%include "attribute.i"

%{
#define SWIG_FILE_WITH_INIT

#include <alps/transform/fourier.hpp>
%}

/*
 * %init
 * %{
 *     // Ensures that numpy is set up properly
 *     import_array();
 * %}
 */

%include <alps/transform/common.hpp>
%include <alps/transform/fourier.hpp>

