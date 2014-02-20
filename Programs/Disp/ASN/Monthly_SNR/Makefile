BIN = Monthly_SNR_from_Disp

FJOBS = Filter.o get_snr.o driver.o

cflags = -O3

FFTW3FLAG = -lfftw3 -lm

$(BIN) : $(FJOBS)
	g++ $(cflags) $(FFTW3FLAG) -o $@ $^

%.o : %.c
	g++ $(cflags) -c $^

clean : 
	rm -f $(BIN) $(FJOBS)
