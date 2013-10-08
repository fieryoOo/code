EXECUTABLE = BLUR

NVCC = nvcc
cflags =

LDLIBS = `pkg-config opencv --cflags --libs`

FOBJS = driver.o image_IO.o blur.o

$(EXECUTABLE) : $(FOBJS)
	$(NVCC) $(LDLIBS) $(cflags) -o $@ $^

%.o : %.cu
	$(NVCC) $(cflags) -c $<

clean :
	rm -f $(EXECUTABLE) core.* $(FOBJS)
