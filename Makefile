CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic

TARGETS := build/h_tree build/zst_dme build/bst_dme build/pt_dme build/hp_dme build/gh_tree build/salt build/sllt build/case_gen
COMMON_OBJ := build/cts_common.o
DME_OBJ := build/dme_common.o
SALT_INCLUDES := -Ithird_party -Ithird_party/salt -Ithird_party/flute3
SALT_SRCS := \
	third_party/salt/salt.cpp \
	third_party/salt/base/eval.cpp \
	third_party/salt/base/flute.cpp \
	third_party/salt/base/mst.cpp \
	third_party/salt/base/net.cpp \
	third_party/salt/base/rsa.cpp \
	third_party/salt/base/tree.cpp \
	third_party/salt/refine/cancel_intersect.cpp \
	third_party/salt/refine/flip.cpp \
	third_party/salt/refine/remove_redundant_coincident.cpp \
	third_party/salt/refine/substitute.cpp \
	third_party/salt/refine/u_shift.cpp \
	third_party/salt/utils/log.cpp \
	third_party/flute3/flute.cpp \
	third_party/flute3/flute_lut_stub.cpp
SALT_OBJ := $(patsubst third_party/%.cpp,build/third_party/%.o,$(SALT_SRCS))

.PHONY: all clean cases check test experiments

all: $(TARGETS)

build:
	@mkdir -p build

build/%.o: %.cpp | build
	$(CXX) $(CXXFLAGS) -I. -c -o $@ $<

build/cts_common.o: cts_common.cpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) -I. -c -o $@ $<

build/dme_common.o: dme_common.cpp dme_common.hpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) -I. -c -o $@ $<

build/h_tree.o: h_tree.cpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) -I. -c -o $@ $<

build/zst_dme.o: zst_dme.cpp dme_common.hpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) -I. -c -o $@ $<

build/bst_dme.o: bst_dme.cpp dme_common.hpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) -I. -c -o $@ $<

build/case_gen.o: case_gen.cpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) -I. -c -o $@ $<

build/pt_dme.o: pt_dme.cpp pt_dme.hpp dme_common.hpp cts_common.hpp third_party/delaunator.hpp | build
	$(CXX) $(CXXFLAGS) -pthread -I. -c -o $@ $<

build/pt_dme_test_core.o: pt_dme.cpp pt_dme.hpp dme_common.hpp cts_common.hpp third_party/delaunator.hpp | build
	$(CXX) $(CXXFLAGS) -pthread -DPT_DME_NO_MAIN -I. -c -o $@ $<

build/pt_dme_tests.o: pt_dme_tests.cpp pt_dme.hpp | build
	$(CXX) $(CXXFLAGS) -pthread -I. -c -o $@ $<

build/hp_dme.o: hp_dme.cpp hp_dme.hpp dme_common.hpp cts_common.hpp third_party/delaunator.hpp | build
	$(CXX) $(CXXFLAGS) -pthread -I. -c -o $@ $<

build/hp_dme_test_core.o: hp_dme.cpp hp_dme.hpp dme_common.hpp cts_common.hpp third_party/delaunator.hpp | build
	$(CXX) $(CXXFLAGS) -pthread -DHP_DME_NO_MAIN -I. -c -o $@ $<

build/hp_dme_tests.o: hp_dme_tests.cpp hp_dme.hpp | build
	$(CXX) $(CXXFLAGS) -pthread -I. -c -o $@ $<

build/gh_tree.o: gh_tree.cpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) -I. -c -o $@ $<

build/gh_tree_test_core.o: gh_tree.cpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) -DGH_TREE_NO_MAIN -I. -c -o $@ $<

build/gh_tree_tests.o: gh_tree_tests.cpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) -I. -c -o $@ $<

build/third_party/%.o: third_party/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(SALT_INCLUDES) -c -o $@ $<

build/bst_dme_core.o: bst_dme.cpp bst_dme.hpp cts_common.hpp dme_common.hpp | build
	$(CXX) $(CXXFLAGS) -DBST_DME_NO_MAIN -I. -c -o $@ $<

build/salt_common.o: salt_common.cpp salt_common.hpp bst_dme.hpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) $(SALT_INCLUDES) -I. -c -o $@ $<

build/salt.o: salt.cpp salt_common.hpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) $(SALT_INCLUDES) -I. -c -o $@ $<

build/sllt.o: sllt.cpp salt_common.hpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) $(SALT_INCLUDES) -I. -c -o $@ $<

build/salt_tests.o: salt_tests.cpp salt_common.hpp cts_common.hpp | build
	$(CXX) $(CXXFLAGS) $(SALT_INCLUDES) -I. -c -o $@ $<

build/h_tree: build/h_tree.o $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -o $@ $^

build/zst_dme: build/zst_dme.o $(DME_OBJ) $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -o $@ $^

build/bst_dme: build/bst_dme.o $(DME_OBJ) $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -o $@ $^

build/case_gen: build/case_gen.o $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -o $@ $^

build/pt_dme: build/pt_dme.o $(DME_OBJ) $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -pthread -o $@ $^

build/pt_dme_tests: build/pt_dme_tests.o build/pt_dme_test_core.o $(DME_OBJ) $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -pthread -o $@ $^

build/hp_dme: build/hp_dme.o $(DME_OBJ) $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -pthread -o $@ $^

build/hp_dme_tests: build/hp_dme_tests.o build/hp_dme_test_core.o $(DME_OBJ) $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -pthread -o $@ $^

build/gh_tree: build/gh_tree.o $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -o $@ $^

build/gh_tree_tests: build/gh_tree_tests.o build/gh_tree_test_core.o $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -o $@ $^

build/salt: build/salt.o build/salt_common.o build/bst_dme_core.o $(SALT_OBJ) $(DME_OBJ) $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -o $@ $^

build/sllt: build/sllt.o build/salt_common.o build/bst_dme_core.o $(SALT_OBJ) $(DME_OBJ) $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -o $@ $^

build/salt_tests: build/salt_tests.o build/salt_common.o build/bst_dme_core.o $(SALT_OBJ) $(DME_OBJ) $(COMMON_OBJ) | build
	$(CXX) $(CXXFLAGS) -o $@ $^

cases: $(TARGETS)
	./build/case_gen generate-all --force

check: $(TARGETS)
	@printf 'built binaries:\n'
	@ls -l $(TARGETS)

test: build/pt_dme_tests build/hp_dme_tests build/gh_tree_tests build/salt_tests
	./build/pt_dme_tests
	./build/hp_dme_tests
	./build/gh_tree_tests
	./build/salt_tests

experiments: $(TARGETS)
	python3 scripts/run_experiments.py

clean:
	rm -rf $(TARGETS) build/pt_dme_tests build/hp_dme_tests build/gh_tree_tests build/salt_tests build/*.o build/third_party
