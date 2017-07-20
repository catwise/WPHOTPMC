IFORT = 1
#G95 = 1

include $(WISE_MAKEINC)/stdheader.mk
include $(WISE_MAKEINC)/stdmacros.mk

EXECLIST = WPHotpm wphotpm_galmag
INCLIST = WPHotpm.header.tbl \
          WPHotpm.header.multiframe.4bandcryo.tbl \
          WPHotpm.header.multiframe.postcryo.tbl \
          WPHotpm.header.singleframe.4bandcryo.tbl \
          WPHotpm.header.singleframe.postcryo.tbl

FSUBS = most_wise_pm.f unspew_pm.f func_pm.f dfunc_pm.f chkvar.f polygrid.f wpro_v6.f medsort.f estmode.f keyhead.f line_fit.f varworks.f getFITSkey.f divider.f CDmatrix.f GALcphot.f CoaddPhotGAL.f check.f loadXSC.f MergeWPRO.f framefunc.f Noise_Formal.f PosCoMoment.f pure.f  FrameStarMask.f PosMoment.f BIGcphot.f CoaddAnnulus.f CoaddNearby.f  CoaddPhot.f  CoaddStarMask.f photom_error.f cphot.f wfits.f QA_stats.f act_wise.f AnnulusStats.f bilint.f CircPhot.f dfunc.f find_neighbors.f func.f matinv.f most_wise.f unspew.f headpar.f PSFnames.f WhatLoc.f Phot_Wrapper.f GetApCorr.f LoadApCorr.f WhatPos.f TripNearby.f readFhead.f SigClipWmean.f HistStats.f MetWr.f VecStats.f annstat.f centroid.f coordstuff.f countfile.f dcon.f elphot.f find_peak.f imoment.f initmeta.f interp.f isovec.f isovec2.f loadtable.f mannuli.f MergeFluxes.f merge_mags.f parg.f pinput.f readfits.f rot_moment.f smoment.f sortstuff.f statstuff.f stringstuff.f twod_ellip2.f twod_ellip_new.f wimage.f writeTABLE.f

ALL = bin

FOPT = -O3

include $(WISE_MAKEINC)/stdtargets.mk

