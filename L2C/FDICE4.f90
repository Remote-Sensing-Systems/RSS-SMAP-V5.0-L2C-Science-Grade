!     ice4 MASK HARD WIRED
!     contains iceberg alley

      SUBROUTINE FDICE4(IMON,XLAT,XLON, ICE)
         use external_files_l2c_module
         use dir_paths_module
         implicit none
         CHARACTER(1) AMASK(45,180,12)
         CHARACTER(LEN=200) FILENAME
         INTEGER*4 IFAC(8)
         INTEGER*4 ISTART
         INTEGER*4 ILON,ILAT,JLON,IMON
         INTEGER*4 IBIT,III,ICE
         REAL*4 XLAT,XLON

         DATA IFAC/128,64,32,16,8,4,2,1/

         ISTART = 1

         call get_dir_paths(table_dir, data_dir)
         filename = trim(table_dir) // '/' // icefile

         IF(ISTART.EQ.1) THEN
            ISTART=0
            OPEN(UNIT=3,FILE=FILENAME,STATUS='OLD',ACTION='READ',FORM='unformatted',access='stream')
            READ(3) AMASK
            CLOSE(3)
         ENDIF

         ILAT=1+NINT(XLAT+89.5)
         ILON=1+NINT(XLON- 0.5)
         IF(ILAT.LT.1) ILAT=1
         IF(ILON.LT.1) ILON=1
         IF(ILAT.GT.180) ILAT=180
         IF(ILON.GT.360) ILON=360

         JLON=1 + INT((ILON-1)/8)
         IBIT=ILON-8*(JLON-1)
         III=INT(ICHAR(AMASK(JLON,ILAT,IMON))/IFAC(IBIT))
         ICE=1
         IF(2*INT(III/2).EQ.III) ICE=0

         RETURN
      END SUBROUTINE FDICE4
