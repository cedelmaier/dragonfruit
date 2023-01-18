# Masterclass 1

## Before you start...
Before you start, it turns out there are a couple of things that make this process much, much easier.

First, acquire the entire dragonfruit github repository from me from online. This makes it so that you have access to some of the files we will need for the entire project. This is slightly different depending on if you are on Windows or Mac, and so figure out how to get a git client for your system, and then download this repository!

Second, you probably want to put all of this into a directory somewhere on your computer that you can get to easily. There are going to be some small instances here where you need to issue commands to PyMOL's CLI (command line interface) and so having everything in one place just makes this much, much easier.

## Introduction to building peptides for fun and glory
Finding the correct information online to build a peptide sequence from scratch is hard. However, I've gone through the trouble and have a couple of versions that work it seems.

### PyMOL fab
PyMOL has a built-in utility named fab that you can use to build a peptide from a provided sequence. This is useful as you can just hand it the amino acid single letter sequence, and it seems to work quite a bit better than the builder that is GUI accessible.

	fab LEEIQGKVKKLEEQVKSL, Cdc12ah, resi=1, chain=A, ss=1

This builds an AH sequence corresponding with a chain identifier (useful for CHARMM-GUI and later MD simulations). You also want to verify that you have the correct sequence with the sequence viewer.

	Display -> Sequence

This just does part 1 of what we need, as this create an **uncapped** peptide sequence. So, we need to go in an identify where the caps are, and if we can modify them.

### Capping a protein in PyMOL

To change the capping sequence, I also used PyMOL. There are slightly different things you need to do for both the N-term and the C-term. For this exact example, we will be adding NH2 (amide?) groups to the caps of the protein.

The N-term residue can be selected by clicking on the corresponding letter in the sequence displayed by PyMOL. This selects all the atoms in the residue. We don't want that, we just want to pick one residue at a time, so change your Mouse to something that can edit.

	Mouse -> 3 Button Editing
	
We can then click on the nitrogen atom (blue in my case) near the N-term of the amino acid. This can be found by looking around for the C_alpha atom, as the nitrogen should be attached to it. It should get a white wrap looking selection sphere around the point. Then, in PyMOL, also click on the 'Builder' option (upper right for me) in the group of (Reset, Zoom, Orient...).

	Builder
	
This will bring up a dialog box. Then, you can run the 'Fix H' command on the Nitrogen atom to make it a neutral group cap (NH2 in this case). Another hydrogen should appear attached to the nitrogen.

The C-term residue is slightly more complicated. We once again want to find the C_alpha atom of the C-term, and then look for the adjacent carboxyl group (another carbon atom, probably has an oxygen double bonded to it). We then need to add the capping to this portion of the sequence. In my case, it is missing a hydrogen atom, and so if you click 'Fix H' it will add a hydrogen atom. If we want to add an NH2 to the end, what we then have to do is change this hydrogen to a nitrogen atom, which can be done by selecting the newly added hydrogen and clicking the 'N' in the builder menu. Depending how how PyMOL feels, it might add the 2 hydrogen atoms automatically for this group.

Congratulations! You have now built a capped peptide that we can use for simulations! To save, export as a PDB file for further use.

	File -> Export Molecule
	PDB Options
	Save as a PDB file (not the default type!)
	
### PyMOL for charged caps

Charged caps require some slight modifications to the above. For example, for an NH3+ group on the N-term, and an O- group on the C-term, you need to do the following.

N-term you need to change the charge (via the builder menu) to +1 on the nitrogen atom, which then should automatically place the 3 hydrogen atoms it needs nearby.

C-term you need to star the same way by adding an H (Fix H) to the carbon that has the double bound oxygen at the end. Then, using the builder, you can change the hydrogen to an oxygen. This will automatically place a hydrogen atom, which then we need to get rid of. Click on the oxygen and assign a charge of -1, which should remove the extra hydrogen.

Congratulations! You have now placed an NH3+ cap on the N-term, and an O- cap on the C-term. This same workflow can be used to change the caps to really anything you want, but really the only other cap that might be of interest is a neutral C-term cap of OH, instead of O-.

### Using CHARMM-GUI to cap the protein

**This is honestly the method that I would use, as it automatically determines the charged caps and places them in PDB format that is understandable by CHARMM-GUI.**

This is a workaround I figured out, since sometimes it can turn out to be a huge pain to cap the sequence ourselves using PyMOL. You can simply build the sequence in PyMOL, but without caps. Then, use this version in the CHARMM-GUI when prompted. If you start going through the CHARMM-GUI sequence for building a membrane, you can change the terminal caps on your protein, and just make sure that they are alright. On the page that selects terminal cappings, you can populate this with the different termini that are possible (at least in CHARMM-GUI), and look into what they are. For instance, the NTER option is NH3+, while CTER is O-. For a neutral capping, you can use NNEU and CT2, which cap with NH2 on both.

### ChimeraX

You can also use ChimeraX to build a peptide sequence. You can start this process with the Builder tool.

	Tools -> Structure Editing

This will bring up a dialog box where you can entire the peptide sequence. Do this now, using the following sequence (it will bring up another box for setting the phi and psi angles). In the following dialog, put in that you want an idealized alpha-helix, and it should do the rest. As a note, you probably want to put something in the Chain ID, as that will be important later. You can then save the structure of this as a PDB file.


### Orienting your structure with PyMOL

Chances are that you want this structure in a particular orientation for use in CHARMM-GUI. This means that you should *probably* orient it in PyMOL first, as then you can translate it up and down later, but getting the right rotation and transormation and orientation in CHARMM-GUI is a matter of trial and error. Which nobody wants to do.

So, you're going to want two helper python files, which are both in this directory -- axes.py and transform_by_camera_rotation.py. Then, in PyMOL, one has to run these first before using them.

	run axes.py
	run transform_by_camera_rotation.py
	
Then, build the axes on the PyMOL window (honestly, I don't know why PyMOL doesn't come with axes as part of the program, yet here we are).

	axes
	
You should see some axes show up in the bottom left of PyMOL. We can now move on to 'correctly' orienting the structure for use in CHARMM-GUI later.

First, you will notice the 'Orient' button in the PyMOL window. You can click on this, which will orient the molecule so that it lays 'flat' with respect to the axes. So, we can maniuplate the molecule here so that it lays along a specified set of axes, which can be combined with some of the functionality that CHARMM-GUI gives us later to properly position the molecule. Move the molecule so your two favorite residues are pointed **towards the camera**.

	Orient structure so that residues are towards your face.
	transform_by_camera_rotation
	
This will reset the axes in the corner so that the BLUE arrow is pointed towards you, along with the residues of interest. Hooray, you have now oriented the file!


## Preparing the peptide + membrane system with CHARMM-GUI.

Our next task is to prepare the combined peptide+membrane system in CHARMM-GUI. To start this task, you need to figure out which of your structures you are going to use in the PDB file.

### Importing the peptide to CHARMM-GUI

Go to charmm-gui.org and go into the membrane builder.

	Input Generator -> Membrane Builder -> Bilayer builder
	
The first page gives you an option of uploading a PDB file. Depending on your choice here, the next couple of steps might be different.

If you choose the uncapped protein, you **have** to choose a terminal capping in the next section(s). First, proceed through the page where you determine which residues you are going to import (do them all). Next, you get some PDB manipulation options. If, for instance, you wanted the neutral caps for your uncapped protein, you would select

	Terminal group patching
	NNEU, CT2 (nh2, nh2)

For the standard terminal patchings (charged), you would choose

	Terminal group patching
	NTER, CTER (NH3+, O-)

You can choose which options you want here to your heart's content. In the end, you should **always** go look at your PDB file in either the terminal or PyMOL, as you want to double check the terminal patching options to make sure they make sense. Note that PyMOL does try to automatically assign bonds to the atoms, and will sometimes get it wrong, but I would always check to see if you have the right hydrogen atoms, the right total charge, etc. You can check your protein by downloading the PDB file from the top options; there are several things you can either view online (view structure) or download (the PDB files).

### Orienting the peptide with CHARMM-GUI

**THIS IS NOT AN IDEAL WAY TO DO THIS AS YOU WILL SEE IN A MOMENT, AS IT REQUIRES TRIAL AND ERROR**

Now you have to choose how you are going to place your peptide in the membrane system. For this example, we are going to be placing the peptide so that the hydrophobic residues lay flat along the surface (or near the surface) of the bilayer. This means we have to give it some geometry options.

First we set the alignment of two residues to be **perpendicular** to the membrane because we have no other option. Don't worry, we will be able to change this in a moment (**this is where things might go horribly wrong**).

	Align a Vector (Two Atoms) Along Z
	Choose by clicking: ILE4 and LEU11

Then, we can shift this to actually lay flat in the positioning options. For our particular setup, we want to set it up thusly.

	Rotate Molecule respect to the Y axis 270 Degree
	Translate Molecule along Z axis 15 Angstrom
	
You can adapt these as you see fit. You can also use the first option to rotate the helix in the plane of the membrane, which can be useful to create several replicate configurations. When you click on the 'Next Step' button, you will then be taken to the membrane builder step.

**!!! BEFORE YOU DO ANYTHING ELSE YOU NEED TO CHECK YOUR STRUCTURE !!!**

This is the huge caveat that I've alluded to several times by now. Depending on how the peptide was oriented when you wrote it out, you might come up with some 'creative' orientations in this step, so we need to check that. Check this by clicking on the 'view structure' link, and then examining your system to make sure that things are pointed in the correct direction. Note that, at least on my system, +z and -z are reversed to start in their online viewer. Not sure why.

### Selecting your membrane
Now you can actually build your membrane! First, we want to select how much water to have above/below. I like to give it quite a bit, but for our case today the default number will work (22.5 angstroms).

Now you have to guess how big your system is going to be (in angstroms). For my peptide, around 110 x 110 angstroms worked, although larger peptides should be taking up more space. This is just a guess, as the algorithm later will fill the area with the proper number of lipid molecules.

Now, for the main portion of this page, you need to select the lipids/whatever you are building your system out of. My basic membrane is a 3:1 ratio of DOPC:PLPI lipids (shoutout to Brandy), but you can really do whatever you want. The hard part is if you want PIP2, then you have to figure out what arcane combination of peptides model the correct tail-state of the lipids, along with having the extra charge on the head groups. Also, you set the upper and lower leaflet ratios separately, so you're going to want to do that. Once you've built your system, you go back up towards where your initial system guess it, and click 'Show the system info'. This will bring up information for the next step.

Click 'Next step'. It might complain about the number of lipids mismatch, or things like that. Just ignore it for now, as the peptide takes up some of this area.

### Ions (salt matters)
On the next page you can select your salt type and concentration for the system. They don't tell you what the concentration units are, but turns out it is just pure molar. So, 0.05 is 50 mM, which is what I use. So, set this to 50 mM KCl, and click the 'Calculate Solvent Composition' to reset this, and make sure you have the correct number of atoms.

### Finishing the process
At this point, you're going to click through the next several pages while CHARMM-GUI figures out how to pack the lipids/water/ions/peptide for you. You can keep clicking on the 'Component PDB (view structure)' link to see what the system is doing which is pretty fun. Some of these steps can take a LONG time.

The only additional thing we need to do is tell the system what kind of thermodynamic ensemble to use, as well as the temperature to use. On one of the final pages, click the GROMACS checkbox (you can also see the other types of simulation frameworks to use, or just click them all to see what happens!).

You also need to set the temperature. I use 300K, but use whatever you feel like.

Once you click 'Next' it will take a LONG time to generate the results, so sit back, relax, go do something else, whatever, and then we will pick up next time, which is actually equilibrating and running your simulation!

### Checking on jobs
If you want to close the window, you can do so, and come back later to harvest the results. Go into CHARMM-GUI and click on your 'User Profile' link. This will take you to a page where you can go look at the Job IDs. Click on that.

This will take you to a page where you can see the list of jobs that were run by CHARMM-GUI, including what should be the one we are using. Once this changes status to Done, you can click on the Job ID link to go to the results.

### The results
Once you have reached the result page, you need to download the *.tgz file with the files for your simulation in it. There is usually a link in the upper right of the page for downloading the file, click that, and save it somewhere. Congratulations, you have finished the first masterclass!





	