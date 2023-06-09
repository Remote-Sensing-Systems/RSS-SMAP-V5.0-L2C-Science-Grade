subroutine stokes_to_vh( ta)
   implicit none
   real(4) ta(2),tav,tah
   tav=0.5*(ta(1)+ta(2))
   tah=0.5*(ta(1)-ta(2))
   ta(1)=tav
   ta(2)=tah
   return
end


subroutine vh_to_stokes( ta)
   implicit none
   real(4) ta(2),ta1,ta2
   ta1=ta(1)+ta(2)
   ta2=ta(1)-ta(2)
   ta(1)=ta1
   ta(2)=ta2
   return
end


subroutine vh_to_stokes8( ta)
   implicit none
   real(8) ta(2),ta1,ta2
   ta1=ta(1)+ta(2)
   ta2=ta(1)-ta(2)
   ta(1)=ta1
   ta(2)=ta2
   return
end
