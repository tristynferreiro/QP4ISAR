# Quick-look Processor for ISAR imaging
A QLP is a radar data processing tool that produces a focused ISAR image from measured HRR profiles. It was developed for use in validating experiment setups in the field to ensure the collection of high-quality data. The project involved research into ISAR imaging and the implementation of low computation ISAR image processing algorithms. These algorithms were used in the design and implementation of a QLP, which was validated using multiple measured datasets.

>The report detailing this process is available: [report](./docs/Final_Report.pdf/).

The general idea is that HRR ISAR radar data is processed into frames, each frame is motion compensated using a Range Alignment and Autofocus algorithm and then the frames are collated into a video for easy viewing.
![Final QLP Design](./docs/QLPDesign.pdf/).

## /app
This folder contains the MATLAB package that can be installed into MATLAB for easy use. Screenshots of the app's GUI are used to describe the app's design and functionality below.
![QLP GUI system Interactions](./docs/UXAppDesign.pdf/)

## /src
This folder contains all MATLAB source code that was used in the development of the QLP. Various versions of the source code are available along with the final implementations. 

## /tests
This folder contains MATLAB scripts that were used to test various aspects of the QLP.

## /docs
This folder contains the final project report and project poster. Additionally, it contains images used in the readme files in this repo.
