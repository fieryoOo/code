EXECUTABLE = RGB2GREY

NVCC = nvcc
cflags = -O3 

LDLIBS = `pkg-config opencv --cflags --libs`

FOBJS = driver.o image_IO.o rgb2grey.o

$(EXECUTABLE) : $(FOBJS)
	$(NVCC) $(LDLIBS) $(cflags) -o $@ $^

%.o : %.cu
	$(NVCC) $(cflags) -c $<

clean :
	rm -f $(EXECUTABLE) core.* $(FOBJS)
