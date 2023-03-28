! module for salinity retrieval
! adapted to SMAP
! TM 05/04/2015
! changed 03/23/2018:
! use new real(8) GNU BRENT routine as in Aquarius V5.0


module sss_module
   implicit none
   public
   save

   real(4), parameter, private         ::      freq_aq = 1.413 ! Aquarius/SMAP center frequency


! for chi2 minimization to determine sss
   real(8), parameter          :: sss_min=0.d0, sss_max=45.d0, eps=1.0D-10, tol=1.0D-8, zeta=1.0d0


   real(4), dimension(2)       :: ytb0 ! global variable accessed by chi2
   real(4)                     :: y_sst, y_tht ! global variable accessed by chi2
   real(4),dimension(2)        :: ydel_tb0


   real(4), parameter, dimension(2)    :: del_tb0 =  (/1.0, 1.0/)   ! V-pol and H-pol equal weight
! error estimate for tbsur0


contains


   subroutine sss_MLE(tht,sst,surtb_mea,   sss, chisq_sss, iflag_sss)
      use l2c_module_smap
      implicit none

      real(4),    intent(in )  :: tht,sst,surtb_mea(2)

      real(4),    intent(out)  :: sss, chisq_sss
      integer(4), intent(out)  :: iflag_sss

      real(8)                  :: x_min, y_min, b0, b1, eps1, tol1

      real(8)                  ::  bsss_0, bsss_min, bsss_max, ychi2, xv, yv, xvsv
      real(8), parameter       ::  very_large_number = 1.0D15, bstep=0.1D0
      integer(4)               ::  istep, mstep, jm


      ! save global variables to be accesses through chi2
      ydel_tb0 = del_tb0
      y_tht =tht
      y_sst =sst
      ytb0  =surtb_mea

      ! get first estimate for minimum by brute force method running through [0,45] in 0.1 steps
      mstep = nint((sss_max-sss_min)/bstep)
      jm=-1
      ychi2 = very_large_number
      do istep=0,mstep
         xv = sss_min+istep*bstep
         yv = chi2(xv)
         if (yv<ychi2) then
            jm=istep
            xvsv=xv
            ychi2=yv
         endif
      enddo

      if (jm==-1) then
         write(*,*) ' warning. could not find 1st guess for minimum'
         iflag_sss=1
         chisq_sss=missing_val4
         sss=missing_val4
         return
      endif

      bsss_0=xvsv
      bsss_min=xvsv-zeta
      bsss_max=xvsv+zeta
      if (bsss_min<sss_min) bsss_min=sss_min
      if (bsss_max>sss_max) bsss_max=sss_max


      ! BRENT search (from GNU software)

      b0=bsss_min
      b1=bsss_max
      eps1=eps
      tol1=tol

      y_min = local_min(b0, b1, eps1, tol1, chi2, x_min)

      iflag_sss=0 ! o.k.
      sss=x_min
      chisq_sss=y_min

      return
   end subroutine sss_MLE


   real(8) function chi2(sss)
! use V and H pol

!     sos evaluation
!     fixed sst
!     variable sss

      implicit none
      real(8)                    :: sss
      real(4), dimension(2)      :: ze0, ztb0
      real(4)                    :: zsst, ztht, surtep, zsss

      zsst=y_sst
      ztht=y_tht
      surtep=(zsst+273.15)
      zsss=sss

      call fdem0_meissner_wentz(freq_aq,ztht,zsst,zsss, ze0)
      ztb0 = ze0*surtep
      chi2 = ((ytb0(1)-ztb0(1))/ydel_tb0(1))**2 + ((ytb0(2)-ztb0(2))/ydel_tb0(2))**2

      return
   end function chi2


   function local_min ( a, b, eps, t, f, x )
!*****************************************************************************80
!
!! LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much slower
!    than that for a Fibonacci search.  If F has a continuous second
!    derivative which is positive at the minimum (which is not at A or
!    B), then convergence is superlinear, and usually of the order of
!    about 1.324....
!
!    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
!    F is never evaluated at two points closer than TOL.
!
!    If F is a unimodal function and the computed values of F are always
!    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
!    LOCAL_MIN approximates the abscissa of the global minimum of F on the
!    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.
!
!    If F is not unimodal, then LOCAL_MIN may approximate a local, but
!    perhaps non-global, minimum to the same accuracy.
!
!    Thanks to Jonathan Eggleston for pointing out a correction to the
!    golden section step, 01 July 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 July 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, real ( kind = 8 ) EPS, a positive relative error tolerance.
!    EPS should be no smaller than twice the relative machine precision,
!    and preferably not much less than the square root of the relative
!    machine precision.
!
!    Input, real ( kind = 8 ) T, a positive absolute error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose local minimum is being sought.
!
!    Output, real ( kind = 8 ) X, the estimated value of an abscissa
!    for which F attains a local minimum value in [A,B].
!
!    Output, real ( kind = 8 ) LOCAL_MIN, the value F(X).
!
      implicit none

      real ( kind = 8 ) a
      real ( kind = 8 ) b
      real ( kind = 8 ) c
      real ( kind = 8 ) d
      real ( kind = 8 ) e
      real ( kind = 8 ) eps
      real ( kind = 8 ) f
      real ( kind = 8 ) fu
      real ( kind = 8 ) fv
      real ( kind = 8 ) fw
      real ( kind = 8 ) fx
      real ( kind = 8 ) local_min
      real ( kind = 8 ) m
      real ( kind = 8 ) p
      real ( kind = 8 ) q
      real ( kind = 8 ) r
      real ( kind = 8 ) sa
      real ( kind = 8 ) sb
      real ( kind = 8 ) t
      real ( kind = 8 ) t2
      real ( kind = 8 ) tol
      real ( kind = 8 ) u
      real ( kind = 8 ) v
      real ( kind = 8 ) w
      real ( kind = 8 ) x
!
!  C is the square of the inverse of the golden ratio.
!
      c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

      sa = a
      sb = b
      x = sa + c * ( b - a )
      w = x
      v = w
      e = 0.0D+00
      fx = f ( x )
      fw = fx
      fv = fw

      do

         m = 0.5D+00 * ( sa + sb )
         tol = eps * abs ( x ) + t
         t2 = 2.0D+00 * tol
!
!  Check the stopping criterion.
!
         if ( abs ( x - m ) <= t2 - 0.5D+00 * ( sb - sa ) ) then
            exit
         end if
!
!  Fit a parabola.
!
         r = 0.0D+00
         q = r
         p = q

         if ( tol < abs ( e ) ) then

            r = ( x - w ) * ( fx - fv )
            q = ( x - v ) * ( fx - fw )
            p = ( x - v ) * q - ( x - w ) * r
            q = 2.0D+00 * ( q - r )

            if ( 0.0D+00 < q ) then
               p = - p
            end if

            q = abs ( q )

            r = e
            e = d

         end if

         if ( abs ( p ) < abs ( 0.5D+00 * q * r ) .and. &
            q * ( sa - x ) < p .and. &
            p < q * ( sb - x ) ) then
!
!  Take the parabolic interpolation step.
!
            d = p / q
            u = x + d
!
!  F must not be evaluated too close to A or B.
!
            if ( ( u - sa ) < t2 .or. ( sb - u ) < t2 ) then

               if ( x < m ) then
                  d = tol
               else
                  d = - tol
               end if

            end if
!
!  A golden-section step.
!
         else

            if ( x < m ) then
               e = sb - x
            else
               e = sa - x
            end if

            d = c * e

         end if
!
!  F must not be evaluated too close to X.
!
         if ( tol <= abs ( d ) ) then
            u = x + d
         else if ( 0.0D+00 < d ) then
            u = x + tol
         else
            u = x - tol
         end if

         fu = f ( u )
!
!  Update A, B, V, W, and X.
!
         if ( fu <= fx ) then

            if ( u < x ) then
               sb = x
            else
               sa = x
            end if

            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu

         else

            if ( u < x ) then
               sa = u
            else
               sb = u
            end if

            if ( fu <= fw .or. w == x ) then
               v = w
               fv = fw
               w = u
               fw = fu
            else if ( fu <= fv .or. v == x .or. v == w ) then
               v = u
               fv = fu
            end if

         end if

      end do

      local_min = fx

      return
   end function local_min




end module sss_module
