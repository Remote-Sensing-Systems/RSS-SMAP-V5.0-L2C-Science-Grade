subroutine check_orbit(jorbit,  ibad)
use dir_paths_module
! check if jorbit is in badorbit list (ibad=1) or in bad ancillary orbit list (ibad=2)

implicit none
integer(4), intent(in)               :: jorbit
integer(4), intent(out)              :: ibad

character(len=*), parameter        :: filename_bad     = 'bad_orbits.txt'
character(len=*), parameter        :: filename_bad_anc = 'bad_ancillaries.txt'


integer(4), save                     :: istart = 1

integer(4), save                     :: nbad=0
integer(4), dimension(10000)         :: iorbit_bad=-1

integer(4), save                     :: nbad_anc=0
integer(4), dimension(10000)         :: iorbit_bad_anc=-1


integer(4)                           :: kbad

call get_dir_paths(table_dir, data_dir)
    if (istart==1) then
    
        istart=0
        
        ! read in  bad orbit data 
        write(*,*) ' reading bad orbit list'       
        open(3,file=trim(table_dir) // '/' // filename_bad,status='old')  
        read(3,'(i5)') nbad
        
        if(nbad.gt.10000) stop 'pgm stopped, nbad.gt.10000'
        
        do kbad=1,nbad
            read(3,'(i5)') iorbit_bad(kbad)
        enddo
        
        close(3)
        write(*,*) ' nbad = ',nbad
         
        ! read in  bad ancillary orbit data 
        write(*,*) ' reading bad ancillary data list'
        open(3,file=trim(table_dir) // '/' // filename_bad_anc,status='old')  
        read(3,'(i5)') nbad_anc
        
        if(nbad_anc.gt.10000) stop 'pgm stopped, nbad_anc.gt.10000'
        
        do kbad=1,nbad_anc
            read(3,'(i5)') iorbit_bad_anc(kbad)
        enddo
        
        close(3) 
        write(*,*) ' nbad_anc = ',nbad_anc

        
    endif


   ibad=0  ! orbit OK

    do kbad=1,nbad
        if(jorbit == iorbit_bad(kbad)) then
        ibad=1
        return
        endif
    enddo

    do kbad=1,nbad_anc
        if(jorbit == iorbit_bad_anc(kbad)) then
        ibad=2
        return
        endif
    enddo


return
end subroutine check_orbit