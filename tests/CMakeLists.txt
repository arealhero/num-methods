include (GoogleTest)

enable_testing ()

add_executable (tests
  algebra.specs.cc
)

target_link_libraries (tests
  GTest::gtest_main
)

gtest_discover_tests (tests)

