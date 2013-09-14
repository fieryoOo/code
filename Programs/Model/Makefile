#
#BIN=mhr_grvel_predict_earth_v3_cv_for_sm
BIN = mhr_grvel_predict_earth_v2
 
FFLAGS =  -e -static
FC=f77
CFLAGS = 
 
DBG = 
 
$(BIN) : mhr_grvel_predict_earth_v2.o azidl.o spline.o  project.o new_path.o integr.o
	$(FC) $(FFLAGS) $(CFLAGS)  mhr_grvel_predict_earth_v2.o  azidl.o spline.o project.o new_path.o integr.o -o $(BIN) $(LDFLAGS) $(LDLIBS)
