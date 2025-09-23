find_package(ROOT REQUIRED COMPONENTS Core Hist Physics)
add_executable(analysis analysis.cc)
target_link_libraries(analysis PUBLIC HepMC3::HepMC3)

set_target_properties(analysis PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Achilles HepMC3 Root Analysis executables
add_executable(analyze_hepmc ../achilles_hepmc3_root_analysis/analyze_hepmc.cpp)
target_link_libraries(analyze_hepmc PUBLIC HepMC3::HepMC3 ROOT::Core ROOT::Hist ROOT::Physics)
set_target_properties(analyze_hepmc PROPERTIES  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)


add_executable(analyze_hepmcCC ../achilles_hepmc3_root_analysis/analyze_hepmcCC.cpp)
target_link_libraries(analyze_hepmcCC PUBLIC HepMC3::HepMC3 ROOT::Core ROOT::Hist ROOT::Physics)
set_target_properties(analyze_hepmcCC PROPERTIES  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)


add_executable(analyze_hepmcNC ../achilles_hepmc3_root_analysis/analyze_hepmcNC.cpp)
target_link_libraries(analyze_hepmcNC PUBLIC HepMC3::HepMC3 ROOT::Core ROOT::Hist ROOT::Physics )
set_target_properties(analyze_hepmcNC PROPERTIES  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)


add_executable(outgoing_particle_multiplicity ../achilles_hepmc3_root_analysis/outgoing_particle_multiplicity.cpp)
target_link_libraries(outgoing_particle_multiplicity PUBLIC HepMC3::HepMC3 ROOT::Core ROOT::Hist ROOT::Physics)
set_target_properties(outgoing_particle_multiplicity PROPERTIES  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

