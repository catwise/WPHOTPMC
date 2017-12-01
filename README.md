# WPHOTPMC
Pipeline module(program) that carries out the profile(PSF) and aperture photometry.It has a direct legacy(design and underlying code) from the ALLSKY processing, known as “WPHOT”, and the ALLWISE/NEOWISE processing, known as “WPHOTpm”. For CatWISE, the modifications feature the use of UnWISE epoch-mosaics as “single frames”, carrying out profile fitting photometry and source characterization(WPRO) using these epoch frames.

The modifications and ([REDACTED]) porting were carried out by T. Jarrett.

The original authors of WPHOT/PM were Ken Marsh and Tom Jarrett.See the powerpoint presentation “Wphotpmc.Jarrett.pdf” for additional details of the module and performance testing; the presentation is appended to the end of this document.

The code lives here on[REDACTED]
`[REDACTED]`


The binary and ancillary files (e.g., the header file; namelist) are located in
`[REDACTED]`

## To compile:
`./make_gfortran WPHotpmc` where you are using my csh script to run the compile and library linkages.

## How to run WPHOTpmc
The following directory is where you can play with wphot:
[REDACTED]

I suggest you try running there first, and then if that works, just copy the entire directory.

"Testcase" to your own sandbox and start your battery of tests. Read on...

Ok, to run you will be using the csh script: `run_wphotpm.csh`

Have a look at the script, you can see how it points at the various files that are needed, including
* MDET
* frames
* coadds
* psfs, as well as setting the modes and file root names, etc.

You can change things, of course, but I strongly recommend you know what you are doing before modifying the script.Understand the data directory structure that is in play here.

So, for example:
```
cd [REDACTED]/testcase
./run_wphotpm.csh
```

run will take about 930 seconds to run.So be patient.You can look at the horrid verbiage spewing from the program in the file `out`

Meanwhile, the output tables (aka mdex table) is in the directory `output`.

See if this works.If it does, you are ready for the next (tougher) step, which is to compare the mdex table with ALLWISE and assess the run/test performed.

I've written a cross-match fortran program that parses the ALLWISE and unwise tables, and does a coordinate match, and computes some simple differences in photometry and astrometry.

Makes some simple ascii tables that have the results.

Then you will need your fav plotting package to view; e.g., using python.

In my case, I use PGPLOT and WIP, but they do not exist on[REDACTED], so I copy the relevant tables and run my own with my mac.

# WPHotpmc
T.Jarret (July 20, 2017)

## WPHotpmc
WPHotpmc is the direct successor to WPHotpm(ALLWISE) and WPHOT(ALLSKY).

It uses UnWISE epoch-mosaics as “single frames”, carrying out profile fiQng photometry and source characterizaSon(WPRO).

The modifications and([REDACTED]) porSng were carried out by T.Jarre.

This presentation introduces the key aspects of Wphotpmc, and provides a “tutorial” on how to run and interpret the output.

## WPHotpm Background
Photometry and source characterizaSon of the ALLWISE pipeline (IPAC)

### Primary components
* Data manager (inputs & allocations; field division)
* Pre-processing(iniSal measurements, use by WPRO)
* ** WPRO** (profile fiQng photometry and characterization)
* Merging(field division)
* WAPP(coadd or fulldepth mosaic photometry )
* Output(meta; mdex)

#### Inputs
* PSFs
* mdet
* frames
* coadds
* Namelist
* headers WPHotpmc

#### Notes
* Instead of a stack of single exposure frames, WPRO works on a stack of unWISE epoch-coadds.
* (and some other small fiddles, such as the frame “masks” ).

### Inputs
* PSF & PSFunc(one for each band)
* Mdet(input list of target sources)
* “Frames” (unWISE epoch mosaics)
* “coadd” (full depth mosaics)
* Namelist(one time set parameters)
* Headers, misc

### Outputs
* Meta & house-keeping
* QA diagnostics
* MEP(binary table) … To be, or not to be
* Mdex(the big kahuna)
* 300 columns …
  * Name
  * coords
  * WPRO
  * WAPP
  * statistics & vari
  * ext.src(gmag)
  * pm biz
  
## Take WPHot for a walkies …
`./run_wphotpm.csh`

[ADD IMAGE1]

## WPHotpmc run script

[ADD IMAGE2]
[ADD IMAGE3]

## WPHotpmc run …

[ADD IMAGE4]
[ADD IMAGE5]

## WPHotpmc run …

[ADD IMAGE 6]

## Photometry Validation
* W`N`rchi2 distribution & error model
* Auto-Compare WPRO to WAPP
* Cross-compare with ALLWISE
* Cross-compare with IRAC

## Validation
* W`N`rchi2 distribution & error model

[ADD IMAGE 7]

## Validation
* W`N`rchi2 distribution & error model

[ADD IMAGE 8]

## Validation
* W`N`rchi2 distribution & error model

[ADD IMAGE 9]

## Validation
* W`N`rchi2 distribution & error model

[ADD IMAGE 10]

## Validation
* Auto-Compare WPRO to WAPP

[ADD IMAGE 11]

## Validation
* Cross-compare with ALLWISE::WPRO

[ADD IMAGE 12]

## Validation
* Cross-compare with ALLWISE::WAPP

[ADD IMAGE 13]

## Validation
* Cross-compare with IRAC

(TBD)

## Validation
* Astrometry : vs.ALLWISE

[ADD IMAGE 14]

## Validation
* Astrometry : PM

(TBD)

## Issues
* **Bad PSFs**: w`N`rchi2 distribution(badly) deviates from unity(but you have some wiggle room by scaling the UNCs)
* **Bad PSFs**: poor photometry performance RELATIVE to WAPP and to ALLSKY ; also, flux is lost(!!!)
* **Up vs.Down scans(epochs)**: astrometric implications
* **Error model**: external validation needed(although internal validation looks plausible)

## Issues: PSFs
WPRO is sensitive to the **cores** of PSFs.If you get the core wrong, WPRO will perform badly.Even if you think you have made the world’s greatest PSFs(because of your awesome brilliance), you will be smacked down by WPRO.Hence, the lesson here is this:
  ALWAYS TEST YOUR PSFs with WPHOTPMC


What we know so far (w/ somewhat limited testing w/ cosmos)
  [REDACTED]

* PSFs constructed from unWISE mosaics (epoch “frames” and full depth) perform badly
* PSFs model using gaussians perform badly
* ALLWISE PSFs work pretty well(surprising since they are not rotated to the cosmos field orientation)
* NEOWISE PSFs also work pretty well
* Weighted combo of cryo + postcryo works well(but get the weighting wrong, say too much cryo, and the performance goes to hell); rotated OK also
