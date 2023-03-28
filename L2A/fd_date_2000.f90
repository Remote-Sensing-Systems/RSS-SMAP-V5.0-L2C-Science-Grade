      subroutine fd_date_2000(time_2000, secyr,lyear,idayjl,imon,idaymo,secdy)
        implicit none
!     time_2000= sec from begin of 2000
!     isecyr= sec from begin of lyear
!     lyear= year, full for digits, i.e., 1987
!     idayjl= julian dat of year
!     imon= month (1 to 12)                   
!     idaymo= day of month (1 to 31)
!     isecdy= seconds of day (0 to 86399)
      
        real(8) time_2000,secyr,secdy,xtime 
        integer(4) lyear,idayjl,imon,idaymo

        integer(4) idaytot,idaybg,ileap,jmon,idayfx(12,2) 
                                                 
      data idayfx/1,32,60,91,121,152,182,213,244,274,305,335, &          
                  1,32,61,92,122,153,183,214,245,275,306,336/  
     
      xtime=time_2000 + 410227200.d0 !convert to time87 
      
        if(xtime.lt.0 .or. xtime.gt.3534451200.d0)  stop 'time oob in fd_date_2000'  !3534451200 is jan 1 2099       
  
      idaytot=1 + int(xtime/86400.d0)
      lyear=1987 + int((idaytot-1)/365)                        !kyear may be 1 year too big
      idaybg= 1 + 365*(lyear-1987) + int((lyear-1985)/4)       !begin day of kyear
      if(idaytot.lt.idaybg) lyear=lyear-1
  
      secyr=xtime-31536000.d0*(lyear-1987)-86400.d0*int((lyear-1985)/4)
                                                  
      ileap=1                                                           
      if(lyear.eq.4*int(lyear/4)) ileap=2                               
      idayjl=1+int(secyr/86400.d0)                                           
      do jmon=2,12                                                  
      imon=jmon-1                                                       
      if(idayfx(jmon,ileap).gt.idayjl) goto 210  
        enddo                        
      imon=12                                                         
      210 continue                                                          
      idaymo=idayjl-idayfx(imon,ileap)+1                                  
      secdy=secyr-(idayjl-1)*86400.d0                                         
      
      return                                                            
      end subroutine fd_date_2000                                                               
