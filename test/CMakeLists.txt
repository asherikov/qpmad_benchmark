include(bench_add_test)

bench_add_test("test_oneshot" "oneshot.cpp" "${PROJECT_SOURCE_DIR}/qp_collection/qpoases_test/json/oneshot/")
bench_add_test("test_oneshot_reserve" "oneshot_reserve.cpp" "${PROJECT_SOURCE_DIR}/qp_collection/qpoases_test/json/oneshot/")
bench_add_test("test_oneshot_profiling" "oneshot_profiling.cpp" "${PROJECT_SOURCE_DIR}/qp_collection/qpoases_test/json/oneshot/")

bench_add_test("test_iterative" "iterative.cpp" "${PROJECT_SOURCE_DIR}/qp_collection/qpoases_test/json/iterative/")
bench_add_test("test_iterative_sparse" "iterative_sparse.cpp" "${PROJECT_SOURCE_DIR}/qp_collection/qpoases_test/json/iterative/")
bench_add_test("test_iterative_reserve" "iterative_reserve.cpp" "${PROJECT_SOURCE_DIR}/qp_collection/qpoases_test/json/iterative/")

