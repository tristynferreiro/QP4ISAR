# Quick-look Processor
This folder contains the initial QLP MATLAB script which was then used to create the QLP MATLAB app.

## QLP scripts
> All scripts require additional functions which are available in the algorithms/ and helper functions/ folders. These need to be copied into the same path as the QLP scripts for them to work.
Additionally, the measured data files are not freely availble, you will need to input your own MATLAB structure data to use the scripts.
The scripts are setup to use **high resolution range profiles** obtained from **stepped-frequency waveform** data.

- [QuicklookProcessor.m](https://github.com/tristynferreiro/QP4ISAR/blob/main/src/Quick-look%20Processor/QuicklookProcessor.m) is the final QLP MATLAB script. It makes use of only the final selected algorithms: [HaywoodRA.m](https://github.com/tristynferreiro/QP4ISAR/blob/main/src/algorithms/Haywood%20RA/HaywoodRA.m) and [YuanAF.m](https://github.com/tristynferreiro/QP4ISAR/blob/main/src/algorithms/Multiple%20DS%20AF/YuanAF.m). This scrip outputs a .avi video format that may not be viewable on all computers (the file type can be changed within the script however this can lead to degraded quality).

- [videoViewer.m](https://github.com/tristynferreiro/QP4ISAR/blob/main/src/Quick-look%20Processor/videoViewer.m) is a simple MATLAB script that can be used to view the videos produced by the QLP.

## /gui 
This contains all the source code used to make the MATLAB app. The app is still in beta and some error-handling and possibly other functionality can still be added.

