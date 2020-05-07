#pragma once

#include <QtCore/qglobal.h>

#ifndef BUILD_STATIC
# if defined(VDRCQTLIB_LIB)
#  define VDRCQTLIB_EXPORT Q_DECL_EXPORT
# else
#  define VDRCQTLIB_EXPORT Q_DECL_IMPORT
# endif
#else
# define VDRCQTLIB_EXPORT
#endif
