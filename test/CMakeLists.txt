include(bench_add_test)

bench_add_test("test_oneshot" "oneshot.cpp" "${PROJECT_SOURCE_DIR}/qp_collection/qpoases_test/json/oneshot/")
bench_add_test("test_iterative" "iterative.cpp" "${PROJECT_SOURCE_DIR}/qp_collection/qpoases_test/json/iterative/")

bench_add_test("test_oneshot_profiling" "oneshot_profiling.cpp" "${PROJECT_SOURCE_DIR}/qp_collection/qpoases_test/json/oneshot/")
