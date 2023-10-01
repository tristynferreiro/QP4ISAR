
This is how the plots were saved to LaTeX friendly code for easy ploting and formatting. 
Note: for imagesc() plots, you need to run the matlab2tikx() and save the .tex file. Then you need to configure the matlab plot and remove all titles,axes, etc and then save only the image as a png. Then update the .tex to use that png instead of the original png that it creates. This is because it compresses the png which causes defocusing of the image. There is a built in function that produces the tikz code for the entire image in the tex file (rather than producing a png) but this takes very long to compile and is not ideal.

**Package origin: [matlab2tikz](http://www.mathworks.com/matlabcentral/fileexchange/22022-matlab2tikz-matlab2tikz?download=true)**

