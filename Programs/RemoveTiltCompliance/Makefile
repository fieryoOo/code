BIN1 := RemoveTilt
BIN2 := RemoveCompliance
BIN3 := RemoveTiltCompliance
BIN4 := RemoveTiltComplianceList
BIN5 := NoiseDirectionality
BIN6 := NoiseDirectionalityList
BINall := $(BIN1) $(BIN2) $(BIN3) $(BIN4) ${BIN5} ${BIN6}

all : $(BINall)

CC := g++
cflags := -std=c++11 -O3 -fopenmp -lSacRec -lfftw3

LIBS = -lstdc++ -O3 -fopenmp -lSacRec -lfftw3
OBJS= 

# --- main rules --- #
define make-bin
$(1) : $(OBJS) $(1).o
	$(CC) $$^ -o $$@ $(LIBS)
endef
$(foreach bin,$(BINall),$(eval $(call make-bin,$(bin))))

# --- .cpp rules with dependencies assembled by gcc --- #
define make-cpp
$(shell $(CC) $(cflags) -MM -MT $(patsubst %.cpp,%.o,$(1)) $(1) | tr -d '\\\n' | awk '{print}' )
	$(CC) $(cflags) -c $$< -o $$@
endef
MOD_DIRS = .
FSRC = $(foreach moddir,$(MOD_DIRS) .,$(wildcard $(moddir)/*.cpp))
$(foreach fcpp,$(FSRC),$(eval $(call make-cpp,$(fcpp))))
