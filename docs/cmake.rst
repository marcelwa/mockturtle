CMake integration and migration
===============================

mockturtle requires CMake 3.25 or newer and C++17. It can be embedded using
``add_subdirectory`` or FetchContent, or installed and consumed using
``find_package(mockturtle CONFIG REQUIRED)``. Both modes provide the same targets:

.. list-table::
   :header-rows: 1

   * - Target
     - Requirements
   * - ``mockturtle::mockturtle``
     - Headers and compile requirements, without compiled backends
   * - ``mockturtle::sat``
     - Base plus ABC SAT, threading, and Nauty when enabled
   * - ``mockturtle::esop``
     - Base plus ABC ESOP
   * - ``mockturtle::all``
     - SAT, ESOP, and explicitly enabled optional integrations

For example, an application constructing and reading networks can link only the
base, whereas an application performing SAT-based equivalence checking needs SAT:

.. code-block:: cmake

   find_package(mockturtle CONFIG REQUIRED)
   target_link_libraries(network_reader PRIVATE mockturtle::mockturtle)
   target_link_libraries(equivalence_checker PRIVATE mockturtle::sat)

Link both ``mockturtle::sat`` and ``mockturtle::esop`` when both backends are used.
Including an algorithm header does not normally require its archive; calling an
algorithm that uses a backend does. Public templates still need bundled dependency
headers and compile definitions even when no archive is linked.

Breaking changes
----------------

The bare ``mockturtle`` target now means the lightweight base. Existing consumers
that need the previous aggregate behavior must use ``mockturtle::all``, or choose
the narrower components above. Generic targets such as ``kitty``, ``fmt``,
``percy``, ``bill``, ``libabcsat`` and ``libabcesop`` are no longer created or
adopted from a parent project. Namespaced targets beginning with an underscore in
the installed export are internal implementation details, not consumer APIs.

Set your own target's ``CXX_STANDARD`` or compile features to request a standard
above C++17; ``MOCKTURTLE_CXX_STANDARD`` has been removed. mockturtle does not change
parent warning flags, ``BUILD_SHARED_LIBS``, or compile-command export preferences.
Warnings apply only to mockturtle's own example, test, and experiment executables.

``MOCKTURTLE_BUILD_EXAMPLES`` defaults to on for standalone builds and off when
embedded. Tests and experiments remain opt-in. The two bundled backend libraries
are always static and position independent, even with ``BUILD_SHARED_LIBS=ON``.
They are excluded from ordinary builds until a consumer links them, unless
installation is enabled.

Installing
----------

``MOCKTURTLE_INSTALL`` defaults to on for standalone builds and off when embedded.
Enabling it intentionally builds both archives so the installed package has all
components. For example:

.. code-block:: console

   cmake -S . -B build -DMOCKTURTLE_BUILD_EXAMPLES=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/mockturtle
   cmake --build build --config Release --parallel
   cmake --install build --config Release

The installed package is relocatable and contains bundled dependency headers under
``include/mockturtle/dependencies`` with their original relative layouts and
licenses. Consumers should use target requirements instead of adding dependency
include directories themselves. The package is unversioned: do not request a
version in ``find_package``. Separately installed optional libraries must also be
available to consumers of a package built with those integrations enabled.

Optional integrations and optimization
--------------------------------------

``MOCKTURTLE_ENABLE_NAUTY`` enables the bundled Nauty archive. It requires a native
Unix environment with a shell and C compiler; Windows and cross compilation fail
explicitly. Generated Nauty headers are included in the installed package. Nauty follows the
SAT component because Percy exact synthesis uses it. Host-specific popcount
detection is disabled to keep installed archives portable.

``MOCKTURTLE_ENABLE_ABC`` requires an ABC archive built according to
https://github.com/lsils/abc-staticlib. Set ``MockturtleABC_ROOT`` to its prefix or
``MockturtleABC_LIBRARY`` to the archive. The old hard-coded ``lib/abc_static``
location is no longer used. ``BILL_Z3`` enables Z3; set ``MockturtleZ3_ROOT`` or
``MockturtleZ3_INCLUDE_DIR`` and ``MockturtleZ3_LIBRARY`` to locate it. The legacy
``BILL_Z3_INCLUDE_PATH`` and ``BILL_Z3_LIBRARY_PATH`` are accepted as discovery
hints. ``MOCKTURTLE_ENABLE_MATPLOTLIB`` requires Python development files and NumPy.
Link ``mockturtle::all`` to use these optional integrations. Discovery happens
only when the corresponding option is enabled and repeats at package consumption
time, without embedding build-machine paths in the export.

Backend archives use hidden visibility, position-independent code, and function
and data sections in optimized configurations. Consumers control final-link
section collection and exported symbols. ``MOCKTURTLE_ENABLE_IPO`` is opt-in and
fails at configuration if the compiler/toolchain does not support it. Installed
IPO archives require a compatible consumer compiler and linker.

``MOCKTURTLE_ENABLE_COVERAGE`` and ``MOCKTURTLE_ENABLE_ASAN`` instrument mockturtle's
own archives and executables with GCC or Clang. Required runtime link flags follow
instrumented archives to their consumers; consumer source compilation flags are
not modified. Configure consumer instrumentation separately when needed.
