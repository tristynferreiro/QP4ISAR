# Simple ISAR Simulator Files
This folder contains all the files and versions used to create the final simulator, which was used to generate simulated ISAR data and images. It uses a stepped frequency pulsed waveform and allows scatterer co-ordinates, amplitude and object motion to be configured.

> All scripts require additional functions which are available in the [algorithms/](https://github.com/tristynferreiro/QP4ISAR/tree/main/src/algorithms) and [helper functions/](https://github.com/tristynferreiro/QP4ISAR/tree/main/src/helper%20functions) folders. These need to be copied into the same path as the QLP scripts for them to work.

#### ISAR_Simulator.m
This script is the **only** the ISAR simulator :

1. Customizable Scatterer Amplitudes: In addition to x-y coordinates, you can now set scatterer amplitudes either manually or by utilizing a Gaussian-like distribution.
2. Scatterer Configuration Display: The simulator plots the scatterer configuration.
3. User-Friendly Command Line Interface (CLI): An easy-to-use CLI is provided, allowing the user to define the rotation and translation motion parameters of the object.
4. Saves the simulated profiles to a .csv file.

There are many other parameters that can be set by changing the code. The CLI was added with the idea of making repetitive testing simpler and is very simplistic in nature.

#### ISAR_Simulator_with_RA_and_AF.m
This script is similar to ISAR_Simulator.m, however it also uses the range-alignment and autofocus algorithms to streamline testing. The algorithms to use on the data can be selected in the CLI.

## Simulator Development/
This folder contains both the original simulator files and their revised versions. The revisions were made to introduce additional features and address previous issues.

## Matlab2Tikz
The [matlab2tikz](http://www.mathworks.com/matlabcentral/fileexchange/22022-matlab2tikz-matlab2tikz?download=true) package was used to save the MATLAB plots in LaTeX-compatible file formats
> When dealing with imagesc() plots, follow these steps:
>1. Run matlab2tikz() and save the resulting .tex file.
2. Configure the MATLAB plot to remove all titles, axes, etc.
3. Save only the image as a PNG.
4. Update the .tex file to use the newly created PNG instead of the original one generated by MATLAB.

>This process is necessary because the original PNG saved by matlab2tikz may be compressed, which reduces the quality of the image. While there is a built-in function for producing the TikZ code for the entire image within the .tex file, it may be very slow to compile in LaTeX."

