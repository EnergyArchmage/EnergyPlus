add_executable (ga_demo ga_demo.cc)

target_link_libraries (ga_demo LINK_PUBLIC global_search)

set_property(TARGET ga_demo PROPERTY CXX_STANDARD 11)
set_property(TARGET ga_demo PROPERTY CXX_STANDARD_REQUIRED ON)

