add_executable(gravitysim gravitysim.cpp)

target_compile_features(gravitysim PRIVATE cxx_std_20)

target_link_libraries(gravitysim PRIVATE gravitysim_library SimpleCanvas_library SDL2::SDL2)