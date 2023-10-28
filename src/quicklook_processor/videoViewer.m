% %% View videos produced by the QLP
vidview = VideoReader('ISARmovie_dataset3.avi');

% Loop through the frames and display them
while hasFrame(vidview)
    frame = readFrame(vidview);
    imshow(frame); % Display the frame

    pause(0.5);
end
