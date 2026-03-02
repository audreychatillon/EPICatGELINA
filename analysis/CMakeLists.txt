install(DIRECTORY ./ DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY} FILES_MATCHING PATTERN "*.h" )
add_library(NPUserAnalysis SHARED Analysis.h Analysis.cxx)
#non root analysis
target_link_libraries(NPUserAnalysis PUBLIC nptool )
# root analysis
#target_link_libraries(NPUserAnalysis PUBLIC ROOT::RIO ROOT::Tree ROOT::Physics ROOT::Imt ROOT::Hist ROOT::Gpad nptool )


