! useVfm: read/write files at the original RAMS vfm format (.true.)
!         or read/write binary files (.false.). VFM format converts
!         floating point data into characters, allowing file portability
!         among different machines, which was required at the time
!         RAMS was written. Currently, floating point data is portable
!         across machines.
!         This flag applyes only to files produced by BRAMS modes MkSfc or MkVFile
!         and consumed by BRAMS mode Initial. Currently, files produced by
!         BRAMS mode Initial are unnafected by this flag (still VFM format).
     
logical, parameter :: useVfm=.true.

