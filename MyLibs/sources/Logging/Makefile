BIN=Listing_two

OBJs=Listing_eight.o  Listing_five.o  Listing_four.o  Listing_nine.o  Listing_seven.o  Listing_six.o  Listing_ten.o  Listing_three.o  Listing_two.o

all : $(BIN)

$(BIN) : $(OBJs)
	g++ -std=c++0x -o $@ $^

%.o : %.cpp
	g++ -std=c++0x -c $<
