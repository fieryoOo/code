BIN1 = SAC_chhdr
BIN2 = SAC_filter
BIN3 = SAC_addall
BIN4 = SAC_dump
BIN5 = SAC_envelope
BIN6 = SAC_FTAN
BIN7 = SAC_mul
BIN8 = SAC_norm
BIN9 = SAC_Correlation
BIN10 = SAC_CrossCorr

LIBS = -L${HOME}/usr/lib -lfftw3 -lSacRec

cflags = -lstdc++ -std=c++11 -O3

#OBJS = SacRec.o
OBJS = 

BINALL = $(BIN1) $(BIN2) ${BIN3} $(BIN4) $(BIN5) $(BIN6) $(BIN7) $(BIN8) $(BIN9) $(BIN10)
all : $(BINALL)

define makebin
$(1) : $(OBJS) $(1).cpp
	g++ $(LIBS) $(cflags) -o $$@ $$^
endef

$(foreach bin,$(BINALL),$(eval $(call makebin,$(bin))))

%.o : %.cpp
	g++ $(cflags) -c $^

.PHONY : clean
clean:
	rm -f $(OBJS) $(BINALL)

.PHONY : install
install:
	cp $(BIN1) /home/yeti4009/bin
