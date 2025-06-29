# Photometric Analysis Algorithms

Least-squares prewhitening and isochrone analysis algorithms used in a paper on beta Cepheids in the cluster NGC 4755. More information will come as it is published.

## High-Level Overview
Together the code in these files forms a comprehensive analysis methodology for both TESS light curves and isochrones on CMD diagrams. A high-level overview of the files is as follows.

* `iso-err.js`: Calculates the uncertainty in isochrone metallicity, extinction, and distance modulus. Needs a sources file exported from [Astromancer](https://astromancer.skynet.unc.edu/cluster) (`ISOCHRONE_FILE` in the code), and individual cluster members' metallicity standard deviation (`METALS_FILE` in the code). This rendition of the code can work with a [VizieR](https://vizier.cds.unistra.fr/) output file with semicolon separators and individual cluster members' metallicity standard deviation as the fifth element in each row.
* `isochrone-show.js`: Outputs a comparison of two different isochrones based on the input files. Needs the graph data export files from Astromancer.
* `prewhiten.js`: Frequency analysis using the Lomb-Scargle periodogram and prewhitening. Needs an input CSV file with `time,flux` rows representing a light curve of the star.
* `compare-freqs.js`: Compare the frequency tables exported from `prewhiten.js`.

## Sample workflow

Before running the code, don't forget to install `nodeplotplib`: `npm i nodeplotplib`. Also, create the directory `/freq-results` for storing the frequency results. Other recommended directory names for input files are set in the code.

A general workflow involving all of the analysis methods here might look like so:

1. Obtain TESS time-series photometry of some target stars in CSV format.
2. Run the prewhitening code (`prewhiten.js`) on all of the light curves from step 1.
3. Compare frequency tables with `compare-freqs.js`
4. Fit an isochrone to a cluster using Astromancer.
5. Export sources, graph data, and literature graph data from step 3.
6. Get metallicity standard deviation of cluster members.
7. Compute isochrone uncertainties with `iso-err.js`.
8. Look at a visual comparison of isochrones with `isochrone-show.js`.

The routine above is just shown for reference, and it is optimized for simplicity and convienience. In reality, many verifications, other checks, and further research would need to be done to scientifically interpret the results.
