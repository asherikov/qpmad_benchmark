MAKE_FLAGS?=-j7


CMAKE_DIR=./cmake/
BUILD_DIR=./build/
ROOT_DIR=../


TARGETS?=all


#----------------------------------------------
# Cleaning
#----------------------------------------------

clean:
	rm -Rf build;


#----------------------------------------------
# Generic targets
#----------------------------------------------
BUILD_SUBDIR=${BUILD_DIR}/

build:
	mkdir -p ${BUILD_SUBDIR};
	cd ${BUILD_SUBDIR}; cmake 	${EXTRA_CMAKE_PARAM} \
								${ROOT_DIR};
	cd ${BUILD_SUBDIR}; ${MAKE} ${MAKE_FLAGS} ${TARGETS}

build-tests: build
	cd ${BUILD_SUBDIR}; env ${TEST_ENV} ctest --verbose
#----------------------------------------------


#----------------------------------------------
# test
#----------------------------------------------
test: clean
	${MAKE} build-tests
#----------------------------------------------


#----------------------------------------------
# other
#----------------------------------------------

update:
	git submodule update --init --recursive

format:
	find -l ./src -iname "*.cpp" "*.h" | xargs clang-format80 -verbose -i

.PHONY: clean cmake build
