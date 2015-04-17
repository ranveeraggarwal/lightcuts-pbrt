Lightcuts
====

Abhinav Gupta
Ranveer Aggarwal

http://www.ranveeraggarwal.com/projects/lightcuts/

##Introduction
We are implementing the paper, Lightcuts by Bruce Walter, Sebastian Fernandez, Adam Arbree, Kavita Bala, Mike Donikian and Don Greenberg. The link for the paper is: http://www.cs.cornell.edu/~kb/projects/lightcuts/  

##Implementation
We are using pbrt-v2, which we have forked from https://github.com/mmp/pbrt-v2 and its README can be found in pbrt-README.txt.

##Installation
After cloning this repository or downloading the source code, run
	
	make

in the src/ folder and generate the binaries.

##Running
In the root directory, run

	./src/bin/pbrt scenes/<pbrtfilename> --outputfile result.tga

##Behind the Scenes
We implemented our own light tree data structure, and then we used error metrics to prune the tree from top to bottom and find light cuts. The tree construction was bottom up.
We wrote a Python script to generate random light sources and used them to compare the time taken while rendering using Vanilla PBRT and our implementation.