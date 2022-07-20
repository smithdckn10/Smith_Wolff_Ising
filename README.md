# Smith_Wolff_Ising
Very similar to First_Ising, except instead of Metropolis-Hastings Monte Carlo sweeps, this one uses the Wolff cluster algorithm

I included two plots made using the wolff alg, these are the png files. It seems to be in agreement with the MH MC method for the most part, at least to me.


tester.cpp is just to test and check the cluster algorithm, so it outputs the lattice to the terminal. 
Starting at lines 338 and 363 in the functions output and artOut, you can remove the comments to output the lattice positions into the terminal as well.
output() is for showing the spins of the sites (not very useful), but artOut shows how the cluster grows.
If a site is in the lattice, it is represented by '*'. Neighboring sites need to be checked, so if the site is unchecked its represented with 'o'. 
If a site has been checked and is not included in the cluster, it is represented by 'x'. Sites that are not adjacent to the cluster are represented with '|'.
It should be a depth first wolff alg, but the sites are not actually changed until all the checks have been done. The equation to check if a site is included
is simplified for the Ising model, because the first thing the alg does is check whether or not the neighbor has the same spin as the initial site. If it 
doesn't it' automatically set to 'x' or "not included". Also, it should just do one cluster update, but this can be changed at the top of the program (I 
wouldn't recommend this, since pretty much every step of the algorithm outputs the lattice, which makes it easy to follow for small clusters and a single
update). The temperature (J) is pretty low just to keep the cluster size small, but this can also be changed. 

