        integer*4 nBadBlend, Nok2a, Nok2b, nAllZero(3)      !  JWF B30130
        integer*4 maxsteps, nPMNaN, MeanObsEpochType, nbmax !  JWF B30117/B30529
        integer*4 nWgtFail, nBlendSwap, nPMsolns, Nwpros    !  JWF B30228
        integer*4 max_mep_lines,mwrote,msources             !  CJG B30228
        integer*4 nSingMatFrmFunc, nTossedItAll, BlendType  !  JWF B30319
        integer*4 MaxWarnPMNaN, MaxWarnBlendSwap            !  JWF B30502/B30520
        integer*4 BGtype(4), nFallbackMMM                   !  JWF B30604
        logical*4 PMsubtract, PMinitADB, TossZeroFlux       !  JWF B30126
        logical*4 GenJD0posns, QuitEarly, WarnZeroFlux      !  JWF B30126
        logical*4 poscon0, InitPMflux0, Init        integer*4 nBadBlend, Nok2a, Nok2b, nAllZero(3)      !  JWF B30130
        integer*4 maxsteps, nPMNaN, MeanObsEpochType, nbmax !  JWF B30117/B30529
        integer*4 nWgtFail, nBlendSwap, nPMsolns, Nwpros    !  JWF B30228
        integer*4 max_mep_lines,mwrote,msources             !  CJG B30228
        integer*4 nSingMatFrmFunc, nTossedItAll, BlendType  !  JWF B30319
        integer*4 MaxWarnPMNaN, MaxWarnBlendSwap            !  JWF B30502/B30520
        integer*4 BGtype(4), nFallbackMMM                   !  JWF B30604
        logical*4 PMsubtract, PMinitADB, TossZeroFlux       !  JWF B30126
        logical*4 GenJD0posns, QuitEarly, WarnZeroFlux      !  JWF B30126
        logical*4 poscon0, InitPMflux0, InitPMposn0         !  JWF B30312
        logical*4 UseNonPMposns, WarnRunaway, CorrPMerr     !  JWF B30329/B30512
        logical*4 next_file                                 !  CJG B30228
        logical*4 WarnAlloc8err, tdbg, DumPvec, WrtFSimg    !  JWF B30419/B30524/B30823
        logical*4 SingleFrame, postcryo                     !  JWF B31121
        real*4    PMfac, minPMsig, ftol_npm, ftol_pm        !  JWF B30312
        real*4    RblendFac, mLogQfac, pminFac, PMstepFac   !  JWF B30402/B30507/B30524
        real*4    BGtrimFrac, SkThresh                      !  JWF B30604
        real*8    MJD0, TargetRA(10), TargetDec(10)         !  JWF B30227
        real*4    SrcSubSNR(4), WarnBigFrameSig             !  TPC B30605/JWF B30711
        real*4    PMminSub, PMrchi2Ratio                    !  JWF B91109
        integer*4 nPMsub                                    !  JWF B91109
        common /JWFcom/ MJD0, nBadBlend, PMfac,             !  JWF B21207
     +                  minPMsig, Nok2a, Nok2b,             !  JWF B21211
     +                  maxsteps, nAllZero, nPMNaN,         !  JWF B30130
     +                  PMsubtract, PMinitADB,              !  JWF B30207
     +                  MeanObsEpochType, nWgtFail,         !  JWF B30208
     +                  nBlendSwap, TossZeroFlux,           !  JWF B30221
     +                  GenJD0posns, nPMsolns,              !  JWF B30227
     +                  TargetRA, TargetDec, QuitEarly,     !  JWF B30227
     +                  Nwpros, WarnZeroFlux, poscon0,      !  JWF B30308
     +                  InitPMflux0, InitPMposn0,           !  JWF B30312
     +                  ftol_npm, ftol_pm, UseNonPMposns,   !  JWF B30312
     +                  nSingMatFrmFunc, nTossedItAll,      !  JWF B30312
     +                  max_mep_lines,mwrote,msources,      !  CJG B30228
     +                  next_file, WarnRunaway, mLogQfac,   !  CJG B30228/JWF B30402
     +                  RblendFac, WarnAlloc8err,           !  JWF B30402
     +                  BlendType, tdbg, MaxWarnPMNaN,      !  JWF B30501
     +                  pminFac, MaxWarnBlendSwap,          !  JWF B30507/B30520
     +                  CorrPMerr, PMstepFac, DumPvec,      !  JWF B30521
     +                  nbmax, BGtype, BGtrimFrac,          !  JWF B30529/B30604
     +                  SkThresh, SrcSubSNR, nFallbackMMM,  !  JWF/TPC B30606
     +                  WarnBigFrameSig, WrtFSimg,          !  JWF B30711/B30823
     +                  SingleFrame, postcryo,              !  JWF B31121
     +                  PMminSub, PMrchi2Ratio, nPMsub      !  JWF B91109
     PMposn0         !  JWF B30312
        logical*4 UseNonPMposns, WarnRunaway, CorrPMerr     !  JWF B30329/B30512
        logical*4 next_file                                 !  CJG B30228
        logical*4 WarnAlloc8err, tdbg, DumPvec, WrtFSimg    !  JWF B30419/B30524/B30823
        logical*4 SingleFrame, postcryo                     !  JWF B31121
        real*4    PMfac, minPMsig, ftol_npm, ftol_pm        !  JWF B30312
        real*4    RblendFac, mLogQfac, pminFac, PMstepFac   !  JWF B30402/B30507/B30524
        real*4    BGtrimFrac, SkThresh                      !  JWF B30604
        real*8    MJD0, TargetRA(10), TargetDec(10)         !  JWF B30227
        real*4    SrcSubSNR(4), WarnBigFrameSig             !  TPC B30605/JWF B30711
        common /JWFcom/ MJD0, nBadBlend, PMfac,             !  JWF B21207
     +                  minPMsig, Nok2a, Nok2b,             !  JWF B21211
     +                  maxsteps, nAllZero, nPMNaN,         !  JWF B30130
     +                  PMsubtract, PMinitADB,              !  JWF B30207
     +                  MeanObsEpochType, nWgtFail,         !  JWF B30208
     +                  nBlendSwap, TossZeroFlux,           !  JWF B30221
     +                  GenJD0posns, nPMsolns,              !  JWF B30227
     +                  TargetRA, TargetDec, QuitEarly,     !  JWF B30227
     +                  Nwpros, WarnZeroFlux, poscon0,      !  JWF B30308
     +                  InitPMflux0, InitPMposn0,           !  JWF B30312
     +                  ftol_npm, ftol_pm, UseNonPMposns,   !  JWF B30312
     +                  nSingMatFrmFunc, nTossedItAll,      !  JWF B30312
     +                  max_mep_lines,mwrote,msources,      !  CJG B30228
     +                  next_file, WarnRunaway, mLogQfac,   !  CJG B30228/JWF B30402
     +                  RblendFac, WarnAlloc8err,           !  JWF B30402
     +                  BlendType, tdbg, MaxWarnPMNaN,      !  JWF B30501
     +                  pminFac, MaxWarnBlendSwap,          !  JWF B30507/B30520
     +                  CorrPMerr, PMstepFac, DumPvec,      !  JWF B30521
     +                  nbmax, BGtype, BGtrimFrac,          !  JWF B30529/B30604
     +                  SkThresh, SrcSubSNR, nFallbackMMM,  !  JWF/TPC B30606
     +                  WarnBigFrameSig, WrtFSimg,          !  JWF B30711/B30823
     +                  SingleFrame, postcryo               !  JWF B31121
