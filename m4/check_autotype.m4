# Check for type specifier 'auto'

AC_DEFUN([AX_CXX_HAVE_AUTOTYPE],
  [AC_CACHE_CHECK(
    [for automatically typed variables],
    ax_cv_cxx_have_autotype,
    [dnl
      AC_LANG_PUSH([C++])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
        [
        ],
        [
          [int a = 1;]
          [auto b = a;]
        ]
        )],
        [ax_cv_cxx_have_autotype=1],
        [ax_cv_cxx_have_autotype=0]
      )
    AC_LANG_POP([C++])])
    if test "$ax_cv_cxx_autotype" == 1
    then
      AC_DEFINE(HAVE_CXX_CREF,
        1,
        [Define if functional defines the std::cref class.])
    fi
  ])

