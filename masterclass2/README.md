# Masterclass 2

## May the odds be ever in your favor
Hooray! You have successfully made it through the first masterclass (hunger games?). Before we start the next part, make sure that you have downloaded your .tgz file from CHARMM-GUI, as we will be using things from it for the next section. There are a couple of logistical things that will help, as well.

### Github access
Make sure that you have a git client running on your local machine, and are able to clone the main dragonfruit repository (where this file is located). On a mac (or linux box) command line, this would be the equivalent of the command

	git clone https://github.com/cedelmaier/dragonfruit.git
	
NOTE: This is for the HTTP version of the repository, which means that you can't push changes to it. If you want to have write access to the repository, you should talk to me. For the purposes of this masterclass, you only need to download things, you should never be uploading.

We are also going to be using git commands on the longleaf cluster, and downloading the repository there, so that we have access to some of my helper scripts.

### Command line interface
Our windows users need to have a command line interface client so that they can login to the longleaf cluster. I've used PuTTY in the past for this, which is what I would recommend. There is also talk about a native Windows 10/11 SSH/SCP client, but we need something that allows us to login to the longleaf cluster.

### Longleaf access
Make sure that you have access to the longleaf cluster at UNC! Moreover, you need to have access to the GPU partition(s), so that we can actually try to run some of these examples. If you can login (via terminal) to longleaf, you can then test to see if you have GPU access by issuing the following command.

	srun --ntasks=1 --cpus-per-task=8 --mem=60G --time=1:00:00 --partition=volta-gpu --gres=gpu:1 --qos=gpu_access --pty /bin/bash

This will log you into the volta-gpu partition in an *interactive* session. Type exit to get out of this session. If this didn't work, then you must not have the permission to run on the GPU partition at the moment, which is kind of a bummer.

## Running simulations with GROMACS
From masterclass 1 we have a peptide+membrane system (of your choice even!) that is now ready to be simulated. First up on our list is getting the simulation onto longleaf.

### Move simulation data to the cluster (longleaf)
First, we need to untar/unzip the file that we got from CHARMM-GUI, which has all of the information we need to run the simulations. Once this has been unzipped, there should be a `gromacs` directory that you can do into. We need everything in this directory for our simulation. For my setup on longleaf, I will use scp to move this to the `scratch` space.

Clusters usually have lots of shared resources for users. However, in many cases, this space is not unlimited, and so a scratch space is usually set up that is periodically wiped clean of data. This is what we will be using for the output of our simulations. You have your own personally scratch space on longleaf under `pine`, which is where you will want to change to as your current directory.

	cd /pine/scr/<first initial of username>/<next initial of username>/<username>
	cd /pine/scr/e/d/edelmaie
	
I have put mine as an example above. Then, make a directory, named whatever you want, to run this simulation in. Again, for example, mind might be something like the following.

	mkdir /pine/scr/e/d/edelmaie/20230120
	mkdir /pine/scr/e/d/edelmaie/20230120/gromacs_masterclass2

You can then use scp to move the data from your machine to the cluster (or the windows equivalent of scp). You should have a directory then with the following files somewhere on the `longleaf`.

	index.ndx
	README
	step5_input.gro
	step5_input.pdb
	step5_input.psf
	step6.0_minimization.mdp
	step6.1_equilibration.mdp
	step6.2_equilibration.mdp
	step6.3_equilibration.mdp
	step6.4_equilibration.mdp
	step6.5_equilibration.mdp
	step6.6_equilibration.mdp
	step7_production.mdp
	topol.top
	toppar/
	
Where `toppar/` is a directory that contains our forcefield files.

### Equilibrating the system
So, now that you have the files on longleaf, the only thing left to do is actually run them. There is a lot involved with running a GROMACS simulation, but I've put a lot of that into a couple of scripts that are going to be useful here. The first script you want to find is one called `run_equilibration.sh`, either in this directory, or in the main dragonfruit directory under `dragonfruit/gromacs/run_equilibration.sh`. If you have access to the volta-gpu partition at UNC, then you can equilibrate your system using the following command (copy run_equilibration.sh into the directory you setup earlier first).

	sbatch run_equilibration.sh
	
This will submit a job to the slurm scheduler to run your equilibration script on the volta-gpu partition. It will then take some amount of time to run.

If you want a more detailed explanation of what is going on, I would suggest looking at the gromacs tutorials on [mdtutorials](http://www.mdtutorials.com/), or the new introduction they put on the gromacs website [located here](https://tutorials.gromacs.org/md-intro-tutorial.html). These should take you through the whole process we are doing here, which involves equilibrating your system first, and then running some amount of 'production' time.

### Simulating our system
At last we will reveal ourselves to the Jedi. Congratulations, you have reached the simulation point of this exercise. You should have an equilibrated structure that is ready for simulating. MD simulations are usually measured in how many ns (or ps, depending on the output) of time that they have simulated. For the 'production' runs I put together a nice little script in the main dragonfruit code-base to set up a sequence of simulations. This is so that if one fails, then the process stops, and it is easy to recover from where it last correctly did a simulation via checkpoint files. If that doesn't make sense, at some point during this, it will.

First, you want to generate some run scripts in the directory where you have your equilibrated files. I wrote a script that **should** do this for you (lots of asterisks here, as there are some caveats, such as when things change on the cluster with new versions of gromacs). This file can be found in the dragonfruit repository as `md_drivers/create_run_scripts.py`. This is a python file that generates the run scripts for your with the proper cluster topology and everything (within reason). First, you are going to want to load python into your cluster environment.

	module purge
	module load python/3.9.6
	
You can then see the options for `create_run_scripts.py` by typing the following command.

	python3 <repositorylocation>/md_drivers/create_run_scripts.py -h
	
For today, we are going to use it on the longleaf cluster at UNC without PLUMED (covered in an upcoming class). This means that we are using the gromacs engine, gpu runtype, volta-gpu partition, and longleaf cluster. Don't worry about the other options for now, as since I wrote it, I also have it so that it will configure runs at the Flatiron Institute. We want to simulate 100 ns worth of time, which means that we have to set the start, end, and stride values appropriately, in this case each 'number/block' here runs for 1 ns of lab time, and so we need 100 total, with a stride of 10 ns at a time in a single block (I know this is confusing, but just look at the command in a moment). Usually we get around 60 ns/day for these simulations, which is a bummer, as GPU computing should be faster, but there are some issues in the way of that. In the end, you really just need to run this command (run in the directory where your equilibrated simulation exists).

	python3 <repo>/md_drivers/create_run_scripts.py --engine gromacs --runtype gpu --partition volta-gpu --cluster longleaf --start 1 --end 100 --stride 10 --nnodes 1 --ntmpi 3 --ntomp 8 --ngpu 1 --gfile step7_production.mdp --nsday 60.0

You should see some files named stuff like `run_prod_stride10_v1_10.sh` and so on. Each of these files runs 10 ns worth of lab time (in 10 different checkpoints), and then submits the next file in the series to the slurm scheduling system. All you have to do now is run the following.

	sbatch run_prod_stride10_v1_10.sh
	
And that's it! You should come back in around a day to see how your simulations are doing.

### Monitoring your simulations while running
This is cool and all, but you also probably want to monitor the system while it is running. We can ask the slurm system this information with a handy command that tells us what all is running (or waiting to run) at the moment.

	squeue -a -u <username> -S 'i' 
	
I put this in my ~/.bashrc file (along with some other helper aliases) so that you can check on the status of simulations. In the masterclass2 directory, I do have a copy of my bashrc file so that you can edit yours to look the same, add helper functions, etc. The one in question here is the `showsq` function, which does this automatically. You can query what you have running/waiting to run with

	showsq
	
and then check when things will start running if they are waiting with

	showsq --start
	
isn't that neat?




	