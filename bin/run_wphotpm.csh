#! /bin/tcsh -f

#  shell script that runs WPHotpmc
#  T. Jarrett ; June 25, 2017


###  these will not change during processing (one-time change only)

set wphot = /Users/jarrett/wphotpm/bin/WPHotpmc

set namelist = /Users/jarrett/wphotpm/bin/nl.WPHot_unwise


#### MDET (this will change every time)

set mdetfile =  ../1497p015_ac51-mdet-raw-1234-2.400.tbl  

###### Frame parameters  (file names here will likely change)

set flist =  ../frames.list_unwise.tbl

set imageid = cosmos

set ftype = u

set fint = img

set func = std

set fcov = n

####  PSF parameters  (likely will change, depending on scan angles, etc)

set psfdir = ../psflib/jwf/cosmos-rotated 

set psfname = unWISE


####  Coadd parameters  (coadd name will change with every set)

set cname = ../fulldepth/unwise-1497p015

set ctype = u


#### QA

set outQAdir = QA
mkdir -p $outQAdir

#### output 

set outname = mdex.tbl
set metaname = meta.tbl

########################################  run the command line
echo ' '

echo $wphot -namlis  $namelist \
-mdettab $mdetfile \
-mdex mdet \
-imageid $imageid  \
-ifile $flist \
-level $ftype -int $fint -unc $func -cov $fcov  \
-calpsfdir $psfdir \
-calbname $psfname \
-calgridszX 1  -calgridszY 1  \
-coadd  $cname \
-clevel $ctype \
-qadir $outQAdir \
-ofile $outname \
-meta $metaname


echo ' '
$wphot -namlis  $namelist \
-mdettab $mdetfile \
-mdex mdet \
-imageid $imageid  \
-ifile $flist \
-level $ftype -int $fint -unc $func -cov $fcov  \
-calpsfdir $psfdir \
-calbname $psfname \
-calgridszX 1  -calgridszY 1  \
-coadd  $cname \
-clevel $ctype \
-qadir $outQAdir \
-ofile $outname \
-meta $metaname -v  >& out



