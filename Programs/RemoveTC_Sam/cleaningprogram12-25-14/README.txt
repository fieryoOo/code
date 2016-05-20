This program, developed at Brown University by Sam Bell, Don Forsyth, and Youyi Ruan, removes water wave noise from OBS data.

For this program to run, Fortran G77, SAC, and Python must be installed.  (For a different Fortran compiler, a few modifications may be necessary.)  The program contains pre-compiled Fortran scripts that may or may not run without compiling.

To start the program, load the data into the main directory and stored in folders named after their stations.  Data files should be in 24-hour blocks, with the file name beginning with the decimal day (with the form 20??.???).  They must also be decimated to five samples per second and given the endings Z.dec.SAC, 1.dec.SAC, 2.dec.SAC, or H.dec3.SAC.  (Note the dec3 for the DPG files.)  Instrument response corrections should also be applied.

To load the files, either copy all of them into the main directory and run fileData.py, or modify fileDataFrom.py by specifying the directory to copy the data from, and then run it.  (These scripts may need to be modified for different file naming schemes.)  

Example files for BB330 (Blanco), G21B (WHOI-ARRA), and J23B (WHOI-Keck) are included.

Seismic data should be in displacement, not velocity or acceleration.

Pressure data should be differential, not absolute pressure.  (Absolute pressure data should work in theory but probably would require some modifications to the technique.)

Each of the Fortran scripts must be compiled, and the tcsh and Python scripts must have had chmod +x applied to them.  Example code for compiling Fortran scripts:
g77 fortranscript.f -o fortranscript -lm /usr/local/sac/lib/libsacio.a

You must also load the station depths into the depths.csv file, being careful not to modify its formatting.

Either create your own event-free window files for each day manually or load the days you want into events/daylist.txt, and run the tcsh script events/findwindows.  You will need to modify events/findwindows to specify the file path.  The events/findwindows script is imperfect, and hand-selected windows are advised for high-accuracy work.  To run the automated program, you must install the IRIS FetchEvent program.  If you aren't working with Cascadia data, you must modify events/findwindows to change the coordinates of the focus region.  Window files are of the form events/estring.20??.???.txt.  Examples are included.

To create hand windows, copy the relevant estring files into a special event directory in events/handwindows, and modify the files there.  If there are any estring files in a special event directory in the handwindows directory, the program will automatically read those files instead of the ones in the main events directory.

To start the program, load the four input files for each station into inputfiles.txt, using code of the form: 
ls BB330/2013.015*dec*SAC > inputfiles.txt

Then run one of the four Python cleaning scripts.  The fastClean scripts do not make the plots of the transfer function, but they do run significantly faster.  The tilt first scripts should be used on stations with significantly larger tilt than compliance noise (typically the case for the older station designs), and the compliance first scripts should be used wherever there tilt does not totally swamp the compliance (typically the case for the newer compact Trillium stations).

The new files should have the following endings:
Z.not.SAC--tilt removed
Z.dep.SAC--pressure removed
Z.nop.SAC--tilt and pressure removed

Some data files may not contain a whole day.  This program may not handle those files properly.

Contact Sam Bell at Samuel_Bell@Brown.EDU if you have any questions.
