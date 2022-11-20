MAKE_FLAGS?=-j1


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
	cd ${BUILD_SUBDIR}; env ${TEST_ENV} ctest --verbose -j 1
#----------------------------------------------


#----------------------------------------------
# test
#----------------------------------------------
test: clean
	${MAKE} build-tests

plots:
	octave --no-history --silent build/test/oneshot.m
	octave --no-history --silent build/test/iterative.m
#	octave --no-history --silent build/test/iterative_sparse.m
#----------------------------------------------


#----------------------------------------------
# other
#----------------------------------------------

update:
	git submodule update --init --recursive

format:
	find -l ./src -iname "*.cpp" "*.h" | xargs clang-format80 -verbose -i

.PHONY: clean cmake build
