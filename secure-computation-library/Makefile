BUILD_DIR=build
DEBUG_BUILD_DIR=${BUILD_DIR}/Debug
RELEASE_BUILD_DIR=${BUILD_DIR}/Release

default: debug

test: debug
	cd ${DEBUG_BUILD_DIR} && make && ./tests.x --use-colour no

coverage: debug
	cd ${DEBUG_BUILD_DIR} && make coverage

debug: _cmake_debug
	cd ${DEBUG_BUILD_DIR} && make

documentation: _cmake_release
	cd ${RELEASE_BUILD_DIR} && make documentation

release: _cmake_release
	cd ${RELEASE_BUILD_DIR} && make

check_formatting:
	./scripts/check_formatting.sh

_cmake_release: _create_build_dir
	cmake -DCMAKE_BUILD_TYPE=Release -B ${RELEASE_BUILD_DIR}

_cmake_debug: _create_build_dir
	cmake -DCMAKE_BUILD_TYPE=Debug -B ${DEBUG_BUILD_DIR}

_create_build_dir:
	mkdir -p ${BUILD_DIR}

.PHONY: default
	test
	coverage
	debug
	release
	_create_build_dir
	_cmake_release
	_cmake_debug
