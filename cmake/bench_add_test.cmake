function(bench_add_test     TEST_NAME   TEST_SOURCES    DIR_PROBLEMS)
    set(TGT_NAME   "TGT_${TEST_NAME}")
    add_executable(${TGT_NAME} "${TEST_SOURCES}")
    set_target_properties(${TGT_NAME} PROPERTIES OUTPUT_NAME "${TEST_NAME}")

    target_link_libraries(${TGT_NAME} "eiquadprog")

    target_link_libraries(${TGT_NAME} "ariles2_visitor_rapidjson")
    target_link_libraries(${TGT_NAME} "ariles2_visitor_octave")
    target_link_libraries(${TGT_NAME} ${Boost_FILESYSTEM_LIBRARIES})

    target_link_libraries(${TGT_NAME} "qpmad")

    target_link_libraries(${TGT_NAME} "qpOASES")
    target_include_directories(${TGT_NAME} PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/qpOASES/include/>
    )

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} ${DIR_PROBLEMS})
endfunction()
