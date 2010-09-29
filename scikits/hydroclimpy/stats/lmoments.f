C***********************************************************************
C*                                                                     *
C*  FORTRAN CODE WRITTEN FOR INCLUSION IN IBM RESEARCH REPORT RC20525, *
C*  'FORTRAN ROUTINES FOR USE WITH THE METHOD OF L-MOMENTS, VERSION 3' *
C*                                                                     *
C*  J. R. M. HOSKING                                                   *
C*  IBM RESEARCH DIVISION                                              *
C*  T. J. WATSON RESEARCH CENTER                                       *
C*  YORKTOWN HEIGHTS                                                   *
C*  NEW YORK 10598, U.S.A.                                             *
C*                                                                     *
C*  VERSION 3     AUGUST 1996                                          *
C*                                                                     *
C***********************************************************************
C

      subroutine samlmr(x,n,xmom,nmom,a,b,ierr)

C  SAMPLE L-MOMENTS OF A DATA ARRAY
C
C  PARAMETERS OF ROUTINE:
C  x      * input* array of length n. contains the data, in ascending
C                  order.
C  n      * input* number of data values
C  xmom   *output* array of length nmom. on exit, contains the sample
C                  L-moments L-1, L-2, T-3, T-4, ... .
C  nmom   * input* number of l-moments to be found. at most max(n,20).
C  a      * input* ) parameters of plotting
C  b      * input* ) position (see below)
C  ierr   * error code : ierr =  0  : no error
C                        ierr = -1  : invalid parameter nmom
C                        ierr = -2  : invalid paramters a, b
C                        ierr = -3  : All data values are equal	
C
C  For unbiased estimates (of the lambda's) set a=b=zero. otherwise,
C  plotting-position estimators are used, based on the plotting position
C  (j+a)/(n+b)  for the j'th smallest of n observations. for example,
C  a=-0.35d0 and b=0.0d0 yields the estimators recommended by
C  Hosking et al. (1985, technometrics) for the gev distribution.
C
C 
C
      implicit double precision (a-h,o-z)
      double precision x(n),xmom(nmom),sum(20)

cf2py double precision dimension(n),intent(in) :: x
cf2py integer depend(x) :: n = len(x)
cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
cf2py integer intent(in) :: nmom
cf2py double precision intent(in) :: a
cf2py double precision intent(in) :: b
cf2py integer intent(out) :: ierr

      data zero/0d0/,one/1d0/

      ierr=0
      if(nmom.gt.20.or.nmom.gt.n)goto 1000
      do 10 j=1,nmom
   10   sum(j)=zero
      if(a.eq.zero.and.b.eq.zero)goto 50
      if(a.le.-one.or.a.ge.b)goto 1010
C
C         Plotting-position estimates of PWM's
C
      do 30 i=1,n
        ppos=(i+a)/(n+b)
        term=x(i)
        sum(1)=sum(1)+term
        do 20 j=2,nmom
          term=term*ppos
   20     sum(j)=sum(j)+term
   30 continue
      do 40 j=1,nmom
   40   sum(j)=sum(j)/n
      goto 100
c
c     Unbiased estimates of PWM's
c
   50 do 70 i=1,n
        z=i
        term=x(i)
        sum(1)=sum(1)+term
        do 60 j=2,nmom
          z=z-one
          term=term*z
   60     sum(j)=sum(j)+term
   70 continue
      y=n
      z=n
      sum(1)=sum(1)/z
      do 80 j=2,nmom
        y=y-one
        z=z*y
   80   sum(j)=sum(j)/z
c
c         l-moments
c
  100 k=nmom
      p0=one
      if(nmom-nmom/2*2.eq.1)p0=-one
      do 120 kk=2,nmom
        ak=k
        p0=-p0
        p=p0
        temp=p*sum(1)
        do 110 i=1,k-1
          ai=i
          p=-p*(ak+ai-one)*(ak-ai)/(ai*ai)
  110     temp=temp+p*sum(i+1)
        sum(k)=temp
  120   k=k-1
      xmom(1)=sum(1)
      if(nmom.eq.1)return
      xmom(2)=sum(2)
      if(sum(2).eq.zero)goto 1020
      if(nmom.eq.2)return
      do 130 k=3,nmom
  130   xmom(k)=sum(k)/sum(2)
      return

C     Parameter nmom invalid
 1000 ierr=-1
      return
C     Plotting parameters invalid
 1010 ierr=-2
      return
C     All data values are equal
 1020 ierr=-3
      return
c
      end

C=======================================================================

      subroutine samlmu(x,n,xmom,nmom,ierr)
C
C  SAMPLE L-MOMENTS OF A DATA ARRAY
C
C  PARAMETERS OF ROUTINE:
C  x      * input* array of length n. contains the data, in ascending
C                  order.
C  n      * input* number of data values
C  xmom   *output* array of length nmom. contains the sample l-moments,
C                  stored as described below.
C  nmom   * input* number of l-moments to be found. at most 100.
C  ierr   *output* error code: ierr =  0 : no error
C                              ierr = -1 : invalid paramter nmom
C                              ierr = -3 : All data values are equal
C
      implicit double precision (a-h,o-z)
      parameter (maxmom=100)
      double precision x(n),xmom(nmom),coef(2,maxmom)
      data zero/0d0/,one/1d0/,two/2d0/
cf2py double precision dimension(n),intent(in) :: x
cf2py integer depend(x) :: n = len(x)
cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
cf2py integer intent(in) :: nmom
cf2py integer intent(out) :: ierr

C
      if(nmom.gt.maxmom)goto 1000
      dn=n
      do 10 j=1,nmom
   10   xmom(j)=zero
      if(nmom.le.2)goto 100
C
C     Unbiased estimates of L-moments -- the 'do 30' loop
C     recursively calculates discrete Legendre polynomials, via
C     eq.(9) of Neuman and Schonbach (1974, int.j.num.meth.eng.)
C
      do 20 j=3,nmom
        temp=one/dfloat((j-1)*(n-j+1))
        coef(1,j)=dfloat(j+j-3)*temp
        coef(2,j)=dfloat((j-2)*(n+j-2))*temp
   20 continue
      temp=-dn-one
      const=one/(dn-one)
      nhalf=n/2
      do 40 i=1,nhalf
        temp=temp+two
        xi=x(i)
        xii=x(n+1-i)
        termp=xi+xii
        termn=xi-xii
        xmom(1)=xmom(1)+termp
        s1=one
        s=temp*const
        xmom(2)=xmom(2)+s*termn
        do 30 j=3,nmom,2
          s2=s1
          s1=s
          s=coef(1,j)*temp*s1-coef(2,j)*s2
          xmom(j)=xmom(j)+s*termp
          if(j.eq.nmom)goto 30
          jj=j+1
          s2=s1
          s1=s
          s=coef(1,jj)*temp*s1-coef(2,jj)*s2
          xmom(jj)=xmom(jj)+s*termn
   30   continue
   40 continue
      if(n.eq.nhalf+nhalf)goto 60
      term=x(nhalf+1)
      s=one
      xmom(1)=xmom(1)+term
      do 50 j=3,nmom,2
        s=-coef(2,j)*s
        xmom(j)=xmom(j)+s*term
   50 continue
C
C         L-moment ratios
C
   60 continue
      xmom(1)=xmom(1)/dn
      if(xmom(2).eq.zero)goto 1010
      do 70 j=3,nmom
   70   xmom(j)=xmom(j)/xmom(2)
      xmom(2)=xmom(2)/dn
      return
c
c         at most two l-moments
c
  100 continue
      sum1=zero
      sum2=zero
      temp=-dn+one
      do 110 i=1,n
        sum1=sum1+x(i)
        sum2=sum2+x(i)*temp
        temp=temp+two
  110 continue
      xmom(1)=sum1/dn
      if(nmom.eq.1)return
      xmom(2)=sum2/(dn*(dn-one))
      return
C
C     Parameter nmom invalid
 1000 ierr=-1
      return
C     All data values equa;
 1010 ierr=-3
      do 1020 j=2,nmom
 1020   xmom(j)=zero
      return
C
      end

C=======================================================================

      subroutine sampwm(x,n,xmom,nmom,a,b,kind,ierr)
C
C  PROBABILITY WEIGHTED MOMENTS OF A DATA ARRAY
C
C  PARAMETERS OF ROUTINE:
C  x      * input* array of length n. contains the data, in ascending
C                  order.
C  n      * input* number of data values
C  xmom   *output* array of length nmom. on exit, contains the sample
C                  probability weighted moments. xmom(i) contains
C                  alpha-sub-(i-1) or beta-sub-(i-1).
C  nmom   * input* number of probability weighted moments to be found.
C                  at most max(n,20).
C  a      * input* ) parameters of plotting
C  b      * input* ) position (see below)
C  kind   * input* specifies which kind of pwm's are to be found.
C                  1  alpha-sub-r = e ( x (1-f(x))**r )
C                  2  beta -sub-r = e ( x f(x)**r )
C  ierr   *output* error code: ierr =  0 : no error
C                              ierr = -1 : invalid paramter nmom
C                              ierr = -2 : invalid plotting positions
C                              ierr = -4 : invalid parameter kind
C 
C  For unbiased estimates set a and b equal to zero. otherwise,
C  plotting-position estimators are used, based on the plotting position
C  (j+a)/(n+b)  for the j'th smallest of n observations. for example,
C  a=-0.35d0 and b=0.0d0 yields the estimators recommended by
C  hosking et al. (1985, technometrics) for the gev distribution.
C
      implicit double precision (a-h,o-z)
      double precision x(n),xmom(nmom)
Cf2py double precision dimension(n),intent(in) :: x
Cf2py integer depend(x) :: n = len(x)
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py double precision intent(in) :: a
Cf2py double precision intent(in) :: b
Cf2py integer intent(in) :: kind
Cf2py integer intent(out) :: ierr

      data zero/0d0/,one/1d0/
      if(nmom.gt.20.or.nmom.gt.n)goto 1000
      if(kind.ne.1.and.kind.ne.2)goto 1010
      do 10 j=1,nmom
   10   xmom(j)=zero
      dn=n
      if(a.eq.zero.and.b.eq.zero)goto 50
      if(a.le.-one.or.a.ge.b)goto 1020
C
C     Plotting-position estimates of PWM's
C
      do 30 i=1,n
        ppos=(i+a)/(n+b)
        if(kind.eq.1)ppos=one-ppos
        term=x(i)
        xmom(1)=xmom(1)+term
        do 20 j=2,nmom
          term=term*ppos
   20     xmom(j)=xmom(j)+term
   30 continue
      do 40 j=1,nmom
   40   xmom(j)=xmom(j)/dn
      return
C
C     Unbiased estimates of pwm's
C
   50 do 70 i=1,n
        di=i
        weight=one/dn
        xmom(1)=xmom(1)+weight*x(i)
        do 60 j=2,nmom
          dj=j-one
          if(kind.eq.1)weight=weight*(dn-di-dj+one)/(dn-dj)
          if(kind.eq.2)weight=weight*(di-dj)/(dn-dj)
   60     xmom(j)=xmom(j)+weight*x(i)
   70 continue
      return
C
 1000 ierr=-1
      return
 1010 ierr=-4
      return
 1020 ierr=-2
      return
C
      end


C=============================================================================
C 
      double precision function cdfwak(x,para)
C
C  CUMULATIVE DISTRIBUTION FUNCTION OF THE WAKEBY DISTRIBUTION
C
C  OTHER ROUTINES USED: QUAWAK
C
C  METHOD: the equation x=g(z), where g(z) is the Wakeby quantile
C  expressed as a function of z=-log(1-f), is solved using Halley's
C  method (the 2nd-order analogue of Newton-Raphson iteration).
C
C  NOTE : This routine has been modified for an easier integration with
C         Numpy/Scipy. It assumes a location parameter of 0 and a scale
C         parameter of 1.
C         Checking the validity of the input parameters is delegated to
C         Numpy/Scipy.
C
C  PARAMETERS OF ROUTINE:
C  x      * input* input value
C  para   * input* array of length 3. contains the parameters of the
C                  distribution, in the order beta, gamma, delta.

      implicit double precision (a-h,o-z)
      double precision para(3)
Cf2py double precision dimension(3),intent(in) :: para
Cf2py double precision intent(in) :: x
Cf2py integer intent(out) :: ierr
      data zero/0d0/,half/0.5d0/,one/1d0/
      data P1/0.1D0/,P7/0.7D0/,P99/0.99D0/
C
C     eps,maxit control the test for convergence of the iteration
C     zincmx is the largest permitted iterative step
C     zmult controls what happens when the iteration steps below zero
C     ufl should be chosen so that dexp(ufl) just does not cause
C     underflow
      data eps/1d-8/,maxit/20/,zincmx/3d0/,zmult/0.2d0/
      data ufl/-170d0/
C
      b=para(1)
      c=para(2)
      d=para(3)
C
      ierr=0
      cdfwak=zero
      if(x.le.zero)return
C
C     Test for special cases
      if(b.eq.zero.and.c.eq.zero.and.d.eq.zero)goto 100
      if(c.eq.zero)goto 110
C
C     General case
      cdfwak=one
      if(d.lt.zero.and.x.ge.one/b-c/d)return
C
C     Initial values for iteration:
C     * If x is in the lowest decile of the distribution, start at z=0
C       (F=0);
C     * If x is in the highest percentile of the distribution,
C       starting value is obtained from asymptotic form of the
C       distribution for large z (F near 1);
C     * Otherwise start at z=0.7 (close to F=0.5).
C
      z=p7
      if(x.lt.quawak(p1,para))z=zero
      if(x.lt.quawak(p99,para))goto 10
      if(d.lt.zero)z=dlog((x-one/b)*d/c+one)/d
      if(d.eq.zero)z=(x-one/b)/c
      if(d.gt.zero)z=dlog(x*d/c+one)/d
   10 continue
C
C     Halley's method, with modifications:
C     * If Halley iteration would move in wrong direction
C       (temp.le.zero), use ordinary Newton-Raphson instead;
C     * If step goes too far (zinc.gt.zincmx or znew.le.zero),
C       limit its length.
      do 30 it=1,maxit
        eb=zero
        bz=-b*z
        if(bz.ge.ufl)eb=dexp(bz)
        gb=z
        if(dabs(b).gt.eps)gb=(one-eb)/b
        ed=dexp(d*z)
        gd=-z
        if(dabs(d).gt.eps)gd=(one-ed)/d
        xest=gb-c*gd
        func=x-xest
        deriv1=eb+c*ed
        deriv2=-b*eb+c*d*ed
        temp=deriv1+half*func*deriv2/deriv1
        if(temp.le.zero)temp=deriv1
        zinc=func/temp
        if(zinc.gt.zincmx)zinc=zincmx
        znew=z+zinc
        if(znew.le.zero)goto 20
        z=znew
        if(dabs(zinc).le.eps)goto 200
        goto 30
   20   z=z*zmult
   30 continue
C
C     Not converged
      goto 300
C
      ierr=0
C     Special case b=c=d=0: Wakeby is exponential
  100 continue
      z=x
      goto 200
C
C     Special case c=0: wakeby is generalized Pareto, bounded above
  110 continue
      cdfwak=one
      if(x.ge.one/b)return
      z=-dlog(one/b-x)
      goto 200
C
C     Convert z value to probability
  200 cdfwak=one
      if(-z.lt.ufl)return
      cdfwak=one-dexp(-z)
      return
C
  300 cdfwak=zero/zero
      return
      end



      double precision function quawak(f,para)
C
C  QUANTILE FUNCTION OF THE WAKEBY DISTRIBUTION
C
C  NOTE: This routine has been modified for an easier integration 
C        with Numpy/Scipy.  It assumes a location parameter of 0 
C        and a scale parameter of 1.
C        This routine is to be used internally by cdfwak.
C        Therefore, no validity tests are perforemd.
C        The number of parameters has been reduced to 3, assuming
C        0 for the location parameter and 1 for the scale parameter.
      implicit double precision (a-h,o-z)
      double precision para(3)
      data zero/0d0/,one/1d0/
C
C     ufl should be chosen so that exp(ufl) just does not cause
C     underflow
      data ufl/-170d0/
      b=para(1)
      c=para(2)
      d=para(3)
C
      if(f.le.zero.or.f.ge.one)goto 10
      z=-dlog(one-f)
      y1=z
      if(b.eq.zero)goto 5
      temp=-b*z
      if(temp.lt.ufl)y1=one/b
      if(temp.ge.ufl)y1=(one-dexp(temp))/b
    5 continue
      y2=z
      if(d.ne.zero)y2=(one-dexp(d*y2))/(-d)
      quawak=y1+c*y2
      return
C
   10 if(f.eq.zero)goto 20
      if(f.eq.one)goto 30
      goto 1010
   20 quawak=zero
      return
   30 if(d.gt.zero)goto 1010
      if(d.lt.zero)quawak=one/b-c/d
      if(d.eq.zero.and.c.gt.zero)goto 1010
      if(d.eq.zero.and.c.eq.zero)quawak=one/b
      return
C
 1010 quawak=zero/zero
      return
      end


C=============================================================================
C L-moments estimations
C
      subroutine lmrexp(para,xmom,nmom,ierr)
C
C  L-MOMENT RATIOS FOR THE EXPONENTIAL DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  para   * input* array of length 2. contains the parameters of the
C                  distribution, in the order xi, alpha (location,
C                  scale).
C  xmom   *output* array of length nmom. on exit, contains the l-moments
C                  lambda-1, lambda-2, tau-3, tau-4, ... .
C  nmom   * input* number of l-moments to be found. at most 20
C  ierr   *output* error code : ierr =  0 : no error
C                               ierr = -1 : invalid parameter nmom
C                               ierr = -2 : invalid input parameters
C
      implicit double precision (a-h,o-z)
      double precision para(2),xmom(nmom)
      data zero/0d0/,half/0.5d0/,two/2d0/
Cf2py double precision dimension(2),intent(in) :: para
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py integer intent(out) :: ierr
      a=para(2)
      if(a.le.zero)goto 1000
      if(nmom.gt.20)goto 1010
      xmom(1)=para(1)+a
      if(nmom.eq.1)return
      xmom(2)=half*a
      if(nmom.eq.2)return
      do 10 j=3,nmom
   10   xmom(j)=two/dfloat(j*(j-1))
      return
C
C     Invalid input parameters
 1000 ierr=-2
      RETURN
C     Parameter nmom too large
 1010 ierr=-1
      return
C
      end


      subroutine lmrgam(para,xmom,nmom,ierr)
C
C  L-MOMENT RATIOS FOR THE GAMMA DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  para   * input* array of length 3. contains the parameters of the
C                  distribution, in the order loc, beta, alpha (scale,shape).
C  xmom   *output* array of length nmom. on exit, contains up to 4 of
C                  the l-moments lambda-1, lambda-2, tau-3, tau-4.
C  nmom   * input* number of l-moments to be found. at most 4.
C  ierr   *output* error code : ierr =  0 : no error
C                               ierr = -1 : invalid parameter nmom
C                               ierr = -2 : invalid input parameters
C
C  other routines used: dlgama
C
      implicit double precision (a-h,o-z)
      double precision para(3),xmom(nmom)
Cf2py double precision dimension(3),intent(in) :: para
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py integer intent(out) :: ierr

      data zero/0d0/,half/0.5d0/,one/1d0/
C
C         const is 1/sqrt(pi)
C
      data const/0.56418 95835 47756 287d0/
C
C         coefficients of rational-function approximations
C         a0 is 1/sqrt(3*pi)
C         c0 is tau-4 for the normal distribution
C
      data a0      / 0.32573501d+00/
      data a1,a2,a3/ 0.16869150d+00, 0.78327243d-01,-0.29120539d-02/
      data b1,b2   / 0.46697102d+00, 0.24255406d+00/
      data c0      / 0.12260172d+00/
      data c1,c2,c3/ 0.53730130d-01, 0.43384378d-01, 0.11101277d-01/
      data d1,d2   / 0.18324466d+00, 0.20166036d+00/
      data e1,e2,e3/ 0.23807576d+01, 0.15931792d+01, 0.11618371d+00/
      data f1,f2,f3/ 0.51533299d+01, 0.71425260d+01, 0.19745056d+01/
      data g1,g2,g3/ 0.21235833d+01, 0.41670213d+01, 0.31925299d+01/
      data h1,h2,h3/ 0.90551443d+01, 0.26649995d+02, 0.26193668d+02/
C
      alpha=para(3)
      beta=para(2)
      if(alpha.le.zero.or.beta.le.zero)goto 1000
      if(nmom.gt.4)goto 1010
C
C     lambda-1
C
      xmom(1)=alpha*beta
      if(nmom.eq.1)return
C
C     lambda-2
C
      xmom(2)=beta*const*dexp(dlgama(alpha+half)-dlgama(alpha))
      if(nmom.eq.2)return
C
C     higher moments
C
      if(alpha.lt.one)goto 10
      z=one/alpha
      xmom(3)=dsqrt(z)*(((a3*z+a2)*z+a1)*z+a0)/((b2*z+b1)*z+one)
      if(nmom.eq.3)return
      xmom(4)=(((c3*z+c2)*z+c1)*z+c0)/((d2*z+d1)*z+one)
      if(nmom.gt.4)ierr=-1
      return
C
   10 z=alpha
      xmom(3)=(((e3*z+e2)*z+e1)*z+one)/(((f3*z+f2)*z+f1)*z+one)
      if(nmom.eq.3)return
      xmom(4)=(((g3*z+g2)*z+g1)*z+one)/(((h3*z+h2)*z+h1)*z+one)
      return
C
C     Invalid input parameters
 1000 ierr=-2
      RETURN
C     Parameter nmom too large
 1010 ierr=-1
      return
C
      end


      subroutine lmrgev(para,xmom,nmom,ierr)
C
C  L-MOMENT RATIOS FOR THE GENERALIZED EXTREME-VALUE DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  para   * input* array of length 3. contains the parameters of the
C                  distribution, in the order xi, alpha, k (location,
C                  scale, shape).
C  xmom   *output* array of length nmom. on exit, contains the l-moments
C                  lambda-1, lambda-2, tau-3, tau-4, ... .
C  nmom   * input* number of l-moments to be found. at most 20.
C  ierr   *output* error code : ierr =  0 : no error
C                               ierr = -1 : invalid parameter nmom
C                               ierr = -2 : invalid input parameters
C
C
C  OTHER ROUTINES USED: DLGAMA
C
      implicit double precision (a-h,o-z)
      double precision para(3),xmom(nmom),zmom(20)
Cf2py double precision dimension(3),intent(in) :: para
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py integer intent(out) :: ierr
      data zero/0d0/,one/1d0/,two/2d0/,three/3d0/,four/4d0/,six/6d0/
C
C         array zmom contains the l-moment ratios of the standard
C         Gumbel distribution (xi=0, alpha=1).
C         zmom(1) is Euler's constant, zmom(2) is log(2).
C
      data zmom/
     *  0.57721 56649 01532 861D 0,  0.69314 71805 59945 309D 0,
     *  0.16992 50014 42312 363D 0,  0.15037 49927 88438 185D 0,
     *  0.55868 35005 77583 138D-1,  0.58110 02399 99710 876D-1,
     *  0.27624 25842 97309 125D-1,  0.30556 37665 79053 126D-1,
     *  0.16465 02822 58328 802D-1,  0.18784 66242 98170 912D-1,
     *  0.10932 82150 63027 148D-1,  0.12697 31266 76329 530D-1,
     *  0.77898 28180 57231 804D-2,  0.91483 61796 21999 726D-2,
     *  0.58333 23893 28363 588D-2,  0.69010 42875 90348 154D-2,
     *  0.45326 79701 80679 549D-2,  0.53891 68113 26595 459D-2,
     *  0.36240 77677 72368 790D-2,  0.43238 76086 05538 096D-2/
C
C         small is used to test whether k is effectively zero
C
      data small/1d-6/
C
      u=para(1)
      a=para(2)
      g=para(3)
      if(a.le.zero.or.g.le.-one)goto 1000
      if(nmom.gt.20)goto 1010
C
C     test for k=0
      if(dabs(g).gt.small)goto 5
      xmom(1)=u
      if(nmom.eq.1)return
      xmom(2)=a*zmom(2)
      if(nmom.eq.2)return
      do 2 i=3,nmom
    2   xmom(i)=zmom(i)
      return
    5 continue
C
C     First 2 moments
      gam=dexp(dlgama(one+g))
      xmom(1)=u+a*(one-gam)/g
      if(nmom.eq.1)return
      xx2=one-two**(-g)
      xmom(2)=a*xx2*gam/g
      if(nmom.eq.2)return
C
C     Higher moments
      z0=one
      do 50 j=3,nmom
        dj=j
        beta=(one-dj**(-g))/xx2
        z0=z0*(four*dj-six)/dj
        z=z0*three*(dj-one)/(dj+one)
        sum=z0*beta-z
        if(j.eq.3)goto 40
        do 30 i=2,j-2
          di=i
          z=z*(di+di+one)*(dj-di)/((di+di-one)*(dj+di))
          sum=sum-z*xmom(i+1)
   30   continue
   40   xmom(j)=sum
   50 continue
      return
C
C     Invalid input parameters
 1000 ierr=-2
      return
C     Parameter nmom too large
 1010 ierr=-1
      return
C
      end



      subroutine lmrglo(para,xmom,nmom,ierr)
C  L-MOMENT RATIOS FOR THE GENERALIZED LOGISTIC DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  para   * input* array of length 3. contains the parameters of the
C                  distribution, in the order xi, alpha, k (location,
C                  scale, shape).
C  xmom   *output* array of length nmom. on exit, contains the l-moments
C                  lambda-1, lambda-2, tau-3, tau-4, ... .
C  nmom   * input* number of l-moments to be found. at most 20.
C  ierr   *output* error code : ierr =  0 : no error
C                               ierr = -1 : invalid parameter nmom
C                               ierr = -2 : invalid input parameters
C
      implicit double precision (a-h,o-z)
      double precision para(3),xmom(nmom),z(10,20)
Cf2py double precision dimension(3),intent(in) :: para
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py integer intent(out) :: ierr
      data zero/0d0/,one/1d0/
      data pi/3.141592653589793238d0/
C
C     small is used to decide whether to approximate the first 2
C     L-moments by a power-series expansion when g is near zero.
C     c1,c2 are coefficients of this power-series expansion.
C     c1 is pi**2/6, c2 is 7*pi**4/360.
C
      data small/1d-4/
      data c1,c2/
     *  0.16449 34066 84822 644d 1,  0.18940 65658 99449 184d 1/
C
C     z-array contains coefficients of the representations of
C     L-moment ratios as polynomials in the shape parameter k
C
      data z(1,3)/1d0/
      data (z(i, 4),i=1, 2)/
     *  0.16666 66666 66666 667d 0,  0.83333 33333 33333 333d 0/
      data (z(i, 5),i=1, 2)/
     *  0.41666 66666 66666 667d 0,  0.58333 33333 33333 333d 0/
      data (z(i, 6),i=1, 3)/
     *  0.66666 66666 66666 667d-1,  0.58333 33333 33333 333d 0,
     *  0.35000 00000 00000 000d 0/
      data (z(i, 7),i=1, 3)/
     *  0.23333 33333 33333 333d 0,  0.58333 33333 33333 333d 0,
     *  0.18333 33333 33333 333d 0/
      data (z(i, 8),i=1, 4)/
     *  0.35714 28571 42857 143d-1,  0.42083 33333 33333 333d 0,
     *  0.45833 33333 33333 333d 0,  0.85119 04761 90476 190d-1/
      data (z(i, 9),i=1, 4)/
     *  0.15099 20634 92063 492d 0,  0.51562 50000 00000 000d 0,
     *  0.29791 66666 66666 667d 0,  0.35466 26984 12698 413d-1/
      data (z(i,10),i=1, 5)/
     *  0.22222 22222 22222 222d-1,  0.31889 32980 59964 727d 0,
     *  0.47997 68518 51851 852d 0,  0.16550 92592 59259 259d 0,
     *  0.13398 36860 67019 400d-1/
      data (z(i,11),i=1, 5)/
     *  0.10650 79365 07936 508d 0,  0.44766 31393 29805 996d 0,
     *  0.36081 01851 85185 185d 0,  0.80390 21164 02116 402d-1,
     *  0.46285 27336 86067 019d-2/
      data (z(i,12),i=1, 6)/
     *  0.15151 51515 15151 515d-1,  0.25131 61375 66137 566d 0,
     *  0.46969 52160 49382 716d 0,  0.22765 04629 62962 963d 0,
     *  0.34713 95502 64550 265d-1,  0.14727 13243 54657 688d-2/
      data (z(i,13),i=1, 6)/
     *  0.79569 50456 95045 695d-1,  0.38976 59465 02057 613d 0,
     *  0.39291 73096 70781 893d 0,  0.12381 31062 61022 928d 0,
     *  0.13499 87139 91769 547d-1,  0.43426 15974 56041 900d-3/
      data (z(i,14),i=1, 7)/
     *  0.10989 01098 90109 890d-1,  0.20413 29966 32996 633d 0,
     *  0.44773 66255 14403 292d 0,  0.27305 34428 27748 383d 0,
     *  0.59191 74382 71604 938d-1,  0.47768 77572 01646 091d-2,
     *  0.11930 26366 63747 775d-3/
      data (z(i,15),i=1, 7)/
     *  0.61934 52050 59490 774d-1,  0.34203 17593 92870 504d 0,
     *  0.40701 37051 73427 396d 0,  0.16218 91928 06752 331d 0,
     *  0.25249 21002 35155 791d-1,  0.15509 34276 62872 107d-2,
     *  0.30677 82085 63922 850d-4/
      data (z(i,16),i=1, 8)/
     *  0.83333 33333 33333 333d-2,  0.16976 83649 02293 474d 0,
     *  0.42219 12828 68366 202d 0,  0.30542 71728 94620 811d 0,
     *  0.84082 79399 72285 210d-1,  0.97243 57914 46208 113d-2,
     *  0.46528 02829 88616 322d-3,  0.74138 06706 96146 887d-5/
      data (z(i,17),i=1, 8)/
     *  0.49716 60284 16028 416d-1,  0.30276 58385 89871 328d 0,
     *  0.41047 33000 89185 506d 0,  0.19483 90265 03251 764d 0,
     *  0.38659 80637 04648 526d-1,  0.34139 94076 42897 226d-2,
     *  0.12974 16173 71825 705d-3,  0.16899 11822 91033 482d-5/
      data (z(i,18),i=1, 9)/
     *  0.65359 47712 41830 065d-2,  0.14387 48475 95085 690d 0,
     *  0.39643 28537 10259 464d 0,  0.32808 41807 20899 471d 0,
     *  0.10797 13931 65194 318d 0,  0.15965 33699 32077 769d-1,
     *  0.11012 77375 69143 819d-2,  0.33798 23645 82066 963d-4,
     *  0.36449 07853 33601 627d-6/
      data (z(i,19),i=1, 9)/
     *  0.40878 45705 49276 431d-1,  0.27024 42907 25441 519d 0,
     *  0.40759 95245 14551 521d 0,  0.22211 14264 89320 008d 0,
     *  0.52846 38846 29533 398d-1,  0.59829 82392 72872 761d-2,
     *  0.32859 39655 65898 436d-3,  0.82617 91134 22830 354d-5,
     *  0.74603 37711 50646 605d-7/
      data (z(i,20),i=1,10)/
     *  0.52631 57894 73684 211d-2,  0.12381 76557 53054 913d 0,
     *  0.37185 92914 44794 917d 0,  0.34356 87476 70189 607d 0,
     *  0.13019 86628 12524 058d 0,  0.23147 43648 99477 023d-1,
     *  0.20519 25194 79869 981d-2,  0.91205 82581 07571 930d-4,
     *  0.19023 86116 43414 884d-5,  0.14528 02606 97757 497d-7/
C
      u=para(1)
      a=para(2)
      g=para(3)
      if(a.le.zero.or.dabs(g).ge.one)goto 1000
      if(nmom.gt.20)goto 1010
C
C     first 2 moments
      gg=g*g
      alam1=-g*(c1+gg*c2)
      alam2=one+gg*(c1+gg*c2)
      if(dabs(g).gt.small)alam2=g*pi/dsin(g*pi)
      if(dabs(g).gt.small)alam1=(one-alam2)/g
      xmom(1)=u+a*alam1
      if(nmom.eq.1)return
      xmom(2)=a*alam2
      if(nmom.eq.2)return
C
C     higher moments
      do 20 m=3,nmom
        kmax=m/2
        sum=z(kmax,m)
        do 10 k=kmax-1,1,-1
   10     sum=sum*gg+z(k,m)
        if(m.ne.m/2*2)sum=-g*sum
        xmom(m)=sum
   20 continue
      return
C
C     Invalid input parameters
 1000 ierr=-2
      RETURN
C     Parameter nmom too large
 1010 ierr=-1
      return
C
      end



      subroutine lmrgno(para,xmom,nmom,ierr)
C
C  L-MOMENT RATIOS FOR THE GENERALIZED NORMAL DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  para   * input* array of length 3. contains the parameters of the
C                  distribution, in the order xi, alpha, k (location,
C                  scale, shape).
C  xmom   *output* array of length nmom. on exit, contains the l-moments
C                  lambda-1, lambda-2, tau-3, tau-4, ... .
C  nmom   * input* number of l-moments to be found. at most 20.
C  ierr   *output* error code : ierr =  0 : no error
C                               ierr = -1 : invalid parameter nmom
C                               ierr = -2 : invalid input parameters
C                               ierr =-100: convergence failed 
C 
C  other routines used: derf
C
      implicit double precision (a-h,o-z)
      double precision para(3),xmom(nmom),est(20),estx(20),sum(20),
     *  zmom(20)
Cf2py double precision dimension(3),intent(in) :: para
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py integer intent(out) :: ierr
      data zero/0d0/,half/0.5d0/,one/1d0/
C
C     array zmom contains l-moments of the standard normal dist.
      data zmom/
     *  0D0,   0.56418 95835 47756 287D 0,
     *  0D0,   0.12260 17195 40890 947D 0,
     *  0D0,   0.43661 15389 50024 944D-1,
     *  0D0,   0.21843 13603 32508 776D-1,
     *  0D0,   0.12963 50158 01507 746D-1,
     *  0D0,   0.85296 21241 91705 402D-2,
     *  0D0,   0.60138 90151 79323 333D-2,
     *  0D0,   0.44555 82586 47650 150D-2,
     *  0D0,   0.34264 32435 78076 985D-2,
     *  0D0,   0.27126 79630 48139 365D-2/
C
C     rrt2 is 1/sqrt(2), rrtpi is 1/sqrt(pi)
      data rrt2 /0.70710 67811 86547 524d0/
      data rrtpi/0.56418 95835 47756 287d0/
C
C     range,eps,maxit control the iterative procedure for numerical
C     integration
      data range/5d0/,eps/1d-8/,maxit/10/
C
      u=para(1)
      a=para(2)
      g=para(3)
      if(a.le.zero)goto 1000
      if(nmom.gt.20)goto 1010
C
C     test for k=0
      if(dabs(g).gt.eps)goto 5
      xmom(1)=u
      if(nmom.eq.1)return
      xmom(2)=a*zmom(2)
      if(nmom.eq.2)return
      do 2 i=3,nmom
    2   xmom(i)=zmom(i)
      return
    5 continue
C
C     Lambda-1
      egg=dexp(half*g*g)
      alam1=(one-egg)/g
      xmom(1)=u+a*alam1
      if(nmom.eq.1)return
C
C     Lambda-2
      alam2=egg*derf(half*g)/g
      xmom(2)=a*alam2
      if(nmom.eq.2)return
C
C     Higher moments. The integral defining Lambda-r is evaluated
C     by iterative application of the trapezium rule.
C 
C     - Initial estimate, using 16 ordinates  (the 'do 20' loop
C       calculates legendre polynomials recursively)
      cc=-g*rrt2
      xmin=cc-range
      xmax=cc+range
      do 10 m=3,nmom
   10   sum(m)=zero
      n=16
      xinc=(xmax-xmin)/n
      do 30 i=1,n-1
        x=xmin+i*xinc
        e=dexp(-((x-cc)**2))
        d=derf(x)
        p1=one
        p=d
        do 20 m=3,nmom
          c1=m+m-3
          c2=m-2
          c3=m-1
          p2=p1
          p1=p
          p=(c1*d*p1-c2*p2)/c3
   20     sum(m)=sum(m)+e*p
   30   continue
      do 40 m=3,nmom
   40   est(m)=sum(m)*xinc
C
C     - double the number of ordinates until converged
      do 90 it=1,maxit
        do 50 m=3,nmom
   50     estx(m)=est(m)
        n=n*2
        xinc=(xmax-xmin)/n
        do 70 i=1,n-1,2
          x=xmin+i*xinc
          e=dexp(-((x-cc)**2))
          d=derf(x)
          p1=one
          p=d
          do 60 m=3,nmom
            c1=m+m-3
            c2=m-2
            c3=m-1
            p2=p1
            p1=p
            p=(c1*d*p1-c2*p2)/c3
   60       sum(m)=sum(m)+e*p
   70   continue
C
C     --- test for convergence
        notcgd=0
        do 80 m=nmom,3,-1
          est(m)=sum(m)*xinc
          if(dabs(est(m)-estx(m)).gt.eps*dabs(est(m)))notcgd=m
   80   continue
        if(notcgd.eq.0)goto 100
   90 continue
C
      ierr=-100-(notcgd-1)
  100 continue
      const=-dexp(cc*cc)*rrtpi/(alam2*g)
      do 110 m=3,nmom
  110   xmom(m)=const*est(m)
      return
C
C     Invalid input parameters
 1000 ierr=-2
      return
C     Parameter nmom too large
 1010 ierr=-1
      return
C
      end



      subroutine lmrgpa(para,xmom,nmom,ierr)
C  L-MOMENT RATIOS FOR THE GENERALIZED PARETO DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  para   * input* array of length 3. contains the parameters of the
C                  distribution, in the order xi, alpha, k (location,
C                  scale, shape).
C                  NOTE : the parameter k is the OPPOSITE of the original
C                  version (for compatibility w/ scipy)
C  xmom   *output* array of length nmom. on exit, contains the l-moments
C                  lambda-1, lambda-2, tau-3, tau-4, ... .
C  nmom   * input* number of l-moments to be found. at most 20.
C  ierr   *output* error code : ierr =  0 : no error
C                               ierr = -1 : invalid parameter nmom
C                               ierr = -2 : invalid input parameters
C
      implicit double precision (a-h,o-z)
      double precision para(3),xmom(nmom)
Cf2py double precision dimension(3),intent(in) :: para
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py integer intent(out) :: ierr
      data zero/0d0/,one/1d0/,two/2d0/
C
      u=para(1)
      a=para(2)
      g=-para(3)
      if(a.le.zero.or.g.le.-one)goto 1000
      if(nmom.gt.20)goto 1010
C
C     Lambda-1
      y=one/(one+g)
      xmom(1)=u+a*y
      if(nmom.eq.1)return
C
C     Lambda-2
      y=y/(two+g)
      xmom(2)=a*y
      if(nmom.eq.2)return
C
C     Higher moments
      y=one
      do 10 m=3,nmom
        am=m-two
        y=y*(am-g)/(m+g)
        xmom(m)=y
   10 continue
      return
C
C     Invalid input parameters
 1000 ierr=-2
      RETURN
C     Parameter nmom too large
 1010 ierr=-1
      return
C
      end



      subroutine lmrgum(para,xmom,nmom,ierr)
C  L-MOMENT RATIOS FOR THE GUMBEL DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  para   * input* array of length 2. contains the parameters of the
C                  distribution, in the order xi, alpha (location,
C                  scale).
C  xmom   *output* array of length nmom. on exit, contains the l-moments
C                  lambda-1, lambda-2, tau-3, tau-4, ... .
C  nmom   * input* number of l-moments to be found. at most 20.
C  ierr   *output* error code : ierr =  0 : no error
C                               ierr = -1 : invalid parameter nmom
C                               ierr = -2 : invalid input parameters
C
C
      implicit double precision (a-h,o-z)
      double precision para(2),xmom(nmom),zmom(20)
Cf2py double precision dimension(2),intent(in) :: para
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py integer intent(out) :: ierr
      data zero/0d0/
C
C     array zmom contains the L-moment ratios of the standard
C     Gumbel distribution (xi=0, alpha=1).
C     zmom(1) is Euler's constant, zmom(2) is log(2).
C
      data zmom/
     *  0.57721 56649 01532 861D 0,  0.69314 71805 59945 309D 0,
     *  0.16992 50014 42312 363D 0,  0.15037 49927 88438 185D 0,
     *  0.55868 35005 77583 138D-1,  0.58110 02399 99710 876D-1,
     *  0.27624 25842 97309 125D-1,  0.30556 37665 79053 126D-1,
     *  0.16465 02822 58328 802D-1,  0.18784 66242 98170 912D-1,
     *  0.10932 82150 63027 148D-1,  0.12697 31266 76329 530D-1,
     *  0.77898 28180 57231 804D-2,  0.91483 61796 21999 726D-2,
     *  0.58333 23893 28363 588D-2,  0.69010 42875 90348 154D-2,
     *  0.45326 79701 80679 549D-2,  0.53891 68113 26595 459D-2,
     *  0.36240 77677 72368 790D-2,  0.43238 76086 05538 096D-2/
C
      a=para(2)
      if(a.le.zero)goto 1000
      if(nmom.gt.20)goto 1010
      xmom(1)=para(1)+a*zmom(1)
      if(nmom.eq.1)return
      xmom(2)=a*zmom(2)
      if(nmom.eq.2)return
      do 10 j=3,nmom
   10   xmom(j)=zmom(j)
      return
C
C     Invalid input parameters
 1000 ierr=-2
      RETURN
C     Parameter nmom too large
 1010 ierr=-1
      return
C
      end



      subroutine lmrkap(para,xmom,nmom,ierr)
C  L-MOMENT RATIOS FOR THE KAPPA DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  para   * input* array of length 4. contains the parameters of the
C                  distribution, in the order xi, alpha, k, h.
C  xmom   *output* array of length nmom. on exit, contains the l-moments
C                  lambda-1, lambda-2, tau-3, tau-4, ... .
C  nmom   * input* number of l-moments to be found. at most 20.
C  ierr   *output* error code : ierr =  0 : no error
C                               ierr = -1 : invalid parameter nmom
C                               ierr = -2 : invalid input parameters
C                               ierr =-99 : calculation broke down
C
C  Other routines used: dlgama,digamd
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION PARA(4),XMOM(NMOM),BETA(20)
Cf2py double precision dimension(4),intent(in) :: para
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py integer intent(out) :: ierr
      data zero/0d0/,half/0.5d0/,one/1d0/,three/3d0/,four/4d0/,six/6d0/
C
C     eu  is Euler's constant
C
      data eu/0.577215664901532861d0/
C
C     small is used to test whether h is effectively zero
C     ofl should be chosen so that exp(ofl) just does not cause
C     overflow
C
      data small/1d-8/,ofl/170d0/
C
      u=para(1)
      a=para(2)
      g=para(3)
      h=para(4)
C
C     test for feasible parameters
C
      if(a.le.zero)goto 1000
      if(g.le.-one)goto 1000
      if(h.lt.zero.and.g*h.le.-one)goto 1000
      if(nmom.gt.20)goto 1010
C
C     Calculate functions occurring in the pwm's beta-sub-r
      dlgam=dlgama(one+g)
      icase=1
      if(h.gt.zero)icase=3
      if(dabs(h).lt.small)icase=2
      if(g.eq.zero)icase=icase+3
      goto(10,30,50,70,90,110),icase
C
C   - case h<0, g nonzero
   10 do 20 ir=1,nmom
        r=ir
        arg=dlgam+dlgama(-r/h-g)-dlgama(-r/h)-g*dlog(-h)
        if(dabs(arg).gt.ofl)goto 1020
   20   beta(ir)=dexp(arg)
      goto 130
C
C   - case h small, g nonzero
   30 do 40 ir=1,nmom
        r=ir
   40   beta(ir)=dexp(dlgam-g*dlog(r))*(one-half*h*g*(one+g)/r)
      goto 130
C
C   - case h>0, g nonzero
   50 do 60 ir=1,nmom
        r=ir
        arg=dlgam+dlgama(one+r/h)-dlgama(one+g+r/h)-g*dlog(h)
        if(dabs(arg).gt.ofl)goto 1020
   60   beta(ir)=dexp(arg)
      goto 130
C
C   - case h<0, g=0
   70 do 80 ir=1,nmom
        r=ir
   80   beta(ir)=eu+dlog(-h)+digamd(-r/h)
      goto 130
C
C   - case h small, g=0
   90 do 100 ir=1,nmom
        r=ir
  100   beta(ir)=eu+dlog(r)
      goto 130
C
C   - case h>0, g=0
  110 do 120 ir=1,nmom
        r=ir
  120   beta(ir)=eu+dlog(h)+digamd(one+r/h)
      goto 130
C
C     Lambda-1
  130 continue
      if(g.eq.zero)xmom(1)=u+a*beta(1)
      if(g.ne.zero)xmom(1)=u+a*(one-beta(1))/g
      if(nmom.eq.1)return
C
C     Lambda-2
      alam2=beta(2)-beta(1)
      if(g.eq.zero)xmom(2)=a*alam2
      if(g.ne.zero)xmom(2)=a*alam2/(-g)
      if(nmom.eq.2)return
C
C     Higher moments
      z0=one
      do 170 j=3,nmom
        dj=j
        z0=z0*(four*dj-six)/dj
        z=z0*three*(dj-one)/(dj+one)
        sum=z0*(beta(j)-beta(1))/alam2-z
        if(j.eq.3)goto 160
        do 150 i=2,j-2
          di=i
          z=z*(di+di+one)*(dj-di)/((di+di-one)*(dj+di))
          sum=sum-z*xmom(i+1)
  150   continue
  160   xmom(j)=sum
  170 continue
      return
C
C     Invalid input parameters
 1000 ierr=-2
      return
C     Parameter nmom too large
 1010 ierr=-1
      return
C     Calculations of L-moments have broken down
 1020 ierr=-99
      return
C
      end



      subroutine lmrnor(para,xmom,nmom,ierr)
C  L-MOMENT RATIOS FOR THE NORMAL DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  para   * input* array of length 2. contains the parameters of the
C                  distribution, in the order mu,sigma (location,scale).
C  xmom   *output* array of length nmom. on exit, contains the l-moments
C                  lambda-1, lambda-2, tau-3, tau-4, ... .
C  nmom   * input* number of l-moments to be found. at most 20.
C  ierr   *output* error code : ierr =  0 : no error
C                               ierr = -1 : invalid parameter nmom
C                               ierr = -2 : invalid input parameters
C
      implicit double precision (a-h,o-z)
      double precision para(2),xmom(nmom),zmom(20)
Cf2py double precision dimension(2),intent(in) :: para
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py integer intent(out) :: ierr
      data zero/0d0/
C
C     array zmom contains l-moments of the standard normal dist.
C
      data zmom/
     *  0D0,   0.56418 95835 47756 287D 0,
     *  0D0,   0.12260 17195 40890 947D 0,
     *  0D0,   0.43661 15389 50024 944D-1,
     *  0D0,   0.21843 13603 32508 776D-1,
     *  0D0,   0.12963 50158 01507 746D-1,
     *  0D0,   0.85296 21241 91705 402D-2,
     *  0D0,   0.60138 90151 79323 333D-2,
     *  0D0,   0.44555 82586 47650 150D-2,
     *  0D0,   0.34264 32435 78076 985D-2,
     *  0D0,   0.27126 79630 48139 365D-2/
C
      if(para(2).le.zero)goto 1000
      if(nmom.gt.20)goto 1010
      xmom(1)=para(1)
      if(nmom.eq.1)return
      xmom(2)=para(2)*zmom(2)
      if(nmom.eq.2)return
      do 10 m=3,nmom
   10   xmom(m)=zmom(m)
      return
C
C     Invalid input parameters
 1000 ierr=-2
      return
C     Parameter nmom too large
 1010 ierr=-1
      return
C
      end



      subroutine lmrpe3(para,xmom,nmom,ierr)
C  L-MOMENT RATIOS FOR THE PEARSON TYPE 3 DISTRIBUTION
C
C  PARAMETERS OF ROUTINE:
C  para   * input* array of length 3. contains the parameters of the
C                  distribution, in the order mu, sigma, gamma (mean,
C                  s.d., skewness).
C  xmom   *output* array of length nmom. on exit, contains up to 4 of
C                  the l-moments lambda-1, lambda-2, tau-3, tau-4.
C  nmom   * input* number of l-moments to be found. at most 4.
C  ierr   *output* error code : ierr =  0 : no error
C                               ierr = -1 : invalid parameter nmom
C                               ierr = -2 : invalid input parameters
C
C  Other routines used: dlgama
C
      implicit double precision (a-h,o-z)
      double precision para(3),xmom(nmom)
Cf2py double precision dimension(3),intent(in) :: para
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py integer intent(out) :: ierr
      data zero/0d0/,half/0.5d0/,one/1d0/,four/4d0/
C
C     small is used to test whether skewness is effectively zero
      data small/1d-6/
C
C     const is 1/sqrt(pi)
      data const/0.56418 95835 47756 287d0/
C
C     Coefficients of rational-function approximations
C     A0 is 1/sqrt(3*pi)
C     C0 is tau-4 for the normal distribution
      data a0      / 0.32573501d+00/
      data a1,a2,a3/ 0.16869150d+00, 0.78327243d-01,-0.29120539d-02/
      data b1,b2   / 0.46697102d+00, 0.24255406d+00/
      data c0      / 0.12260172d+00/
      data c1,c2,c3/ 0.53730130d-01, 0.43384378d-01, 0.11101277d-01/
      data d1,d2   / 0.18324466d+00, 0.20166036d+00/
      data e1,e2,e3/ 0.23807576d+01, 0.15931792d+01, 0.11618371d+00/
      data f1,f2,f3/ 0.51533299d+01, 0.71425260d+01, 0.19745056d+01/
      data g1,g2,g3/ 0.21235833d+01, 0.41670213d+01, 0.31925299d+01/
      data h1,h2,h3/ 0.90551443d+01, 0.26649995d+02, 0.26193668d+02/
C
      sd=para(2)
      if(sd.le.zero)goto 1000
      if(nmom.gt.4)goto 1010
C
C     Lambda-1
      xmom(1)=para(1)
      if(nmom.eq.1)return
C
C     Lambda-2
      gamma=para(3)
      if(dabs(gamma).lt.small)goto 20
      alpha=four/(gamma*gamma)
      beta=dabs(half*sd*gamma)
      alam2=const*dexp(dlgama(alpha+half)-dlgama(alpha))
      xmom(2)=alam2*beta
      if(nmom.eq.2)return
C
C     Higher moments
      if(alpha.lt.one)goto 10
      z=one/alpha
      xmom(3)=dsqrt(z)*(((a3*z+a2)*z+a1)*z+a0)/((b2*z+b1)*z+one)
      if(gamma.lt.zero)xmom(3)=-xmom(3)
      if(nmom.eq.3)return
      xmom(4)=(((c3*z+c2)*z+c1)*z+c0)/((d2*z+d1)*z+one)
      return
C
   10 z=alpha
      xmom(3)=(((e3*z+e2)*z+e1)*z+one)/(((f3*z+f2)*z+f1)*z+one)
      if(gamma.lt.zero)xmom(3)=-xmom(3)
      if(nmom.eq.3)return
      xmom(4)=(((g3*z+g2)*z+g1)*z+one)/(((h3*z+h2)*z+h1)*z+one)
      if(nmom.gt.4)ierr=-1
      return
C
C     Case of zero skewness
   20 xmom(1)=para(1)
      if(nmom.eq.1)return
      xmom(2)=const*para(2)
      if(nmom.eq.2)return
      xmom(3)=0
      if(nmom.eq.3)return
      xmom(4)=c0
      if(nmom.gt.4)ierr=-1
      return
C
C     Invalid input parameters
 1000 ierr=-2
      return
C     Parameter nmom too large
 1010 ierr=-1
      return
C
      end



      subroutine lmrwak(para,xmom,nmom,ierr)
C  L-MOMENT RATIOS FOR THE WAKEBY DISTRIBUTION
C
C  parameters of routine:
C  para   * input* array of length 5. contains the parameters of the
C                  distribution, in the order xi, alpha, beta, gamma,
C                  delta.
C  xmom   *output* array of length nmom. on exit, contains the l-moments
C                  lambda-1, lambda-2, tau-3, tau-4, ... .
C  nmom   * input* number of l-moments to be found. at most 20.
C  ierr   *output* error code : ierr =  0 : no error
C                               ierr = -1 : invalid parameter nmom
C                               ierr = -2 : invalid input parameters
C
      implicit double precision (a-h,o-z)
      double precision para(5),xmom(nmom)
Cf2py double precision dimension(5),intent(in) :: para
Cf2py double precision dimension(nmom),intent(out),depend(nmom) :: xmom
Cf2py integer intent(in) :: nmom
Cf2py integer intent(out) :: ierr
      data zero/0d0/,one/1d0/,two/2d0/
C
      xi=para(1)
      a=para(2)
      b=para(3)
      c=para(4)
      d=para(5)
C
C     Test for valid parameters
C
      if(d.ge.one)goto 1000
      if(b+d.le.zero.and.(b.ne.zero.or.c.ne.zero.or.d.ne.zero))goto 1000
      if(a.eq.zero.and.b.ne.zero)goto 1000
      if(c.eq.zero.and.d.ne.zero)goto 1000
      if(c.lt.zero)goto 1000
      if(a+c.lt.zero)goto 1000
      if(a.eq.zero.and.c.eq.zero)goto 1000
      if(nmom.gt.20)goto 1010
C
C     Lambda-1
      y=a/(one+b)
      z=c/(one-d)
      xmom(1)=xi+y+z
      if(nmom.eq.1)return
C
C     Lambda-2
      y=y/(two+b)
      z=z/(two-d)
      alam2=y+z
      xmom(2)=alam2
      if(nmom.eq.2)return
C
C     Higher moments
      do 10 m=3,nmom
        am=m
        y=y*(am-two-b)/(am+b)
        z=z*(am-two+d)/(am-d)
        xmom(m)=(y+z)/alam2
   10 continue
      return
C
C     Invalid input parameters
 1000 ierr=-2
      return
C     Parameter nmom too large
 1010 ierr=-1
      return
C
      end



C=============================================================================
C Parameter estimation from L-moments
C
      subroutine pelkap(xmom,para,ifail)
C  Parameter estimation via L-moments for the Kappa distribution
C
C  Parameters of routine:
C  xmom   * input* array of length 4. contains the l-moments lambda-1,
C                  lambda-2, tau-3, tau-4.
C  para   *output* array of length 4. on exit, contains the parameters
C                  in the order xi, alpha, k, h.
C  ifail  *output* fail flag. on exit, it is set as follows.
C                  0  successful exit
C                  1  l-moments invalid
C                  2  (tau-3, tau-4) lies above the generalized-logistic
C                     line (suggests that l-moments are not consistent
C                     with any kappa distribution with h.gt.-1)
C                  3  iteration failed to converge
C                  4  unable to make progress from current point in
C                     iteration
C                  5  iteration encountered numerical difficulties -
C                     overflow would have been likely to occur
C                  6  iteration for h and k converged, but overflow
C                     would have occurred when calculating xi and alpha
C
C  NB: Parameters are sometimes not uniquely defined by the first 4
C  L-moments. In such cases the routine returns the solution for which
C  the h parameter is largest.
C
C  Other routines used: dlgama,digamd
C
C  The shape parameters k and h are estimated using Newton-Raphson
C  iteration on the relationship between (tau-3,tau-4) and (k,h).
C  the convergence criterion is that tau-3 and tau-4 calculated from
C  the estimated values of k and h should differ by less than 'eps'
C  from the values supplied in array xmom.
C
      implicit double precision (a-h,o-z)
      double precision xmom(4),para(4)
Cf2py double precision dimension(4),intent(in) :: xmom
Cf2py double precision dimension(4),intent(out) :: para
Cf2py integer, intent(out) :: ifail

      data zero/0d0/,half/0.5d0/,one/1d0/,two/2d0/,three/3d0/,four/4d0/
      data five/5d0/,six/6d0/,twelve/12d0/,twenty/20d0/,thirty/30d0/
      data p725/0.725d0/,p8/0.8d0/
C 
C     eps,maxit control the test for convergence of NR iteration.
C     maxsr is the max. nb. of steplength reductions per iteration
C     hstart is the starting value for h
C     big is used to initialize the criterion function
C     oflexp is such that dexp(oflexp) just does not cause overflow
C     oflgam is such that dexp(dlgama(oflgam)) just does not cause
C     overflow
      data eps/1d-6/,maxit/20/,maxsr/10/,hstart/1.001d0/,big/10d0/
      data oflexp/170d0/,oflgam/53d0/
C
      t3=xmom(3)
      t4=xmom(4)
      do 10 i=1,4
   10   para(i)=zero
C
C     Test for feasibility
      if(xmom(2).le.zero)goto 1000
      if(dabs(t3).ge.one.or.dabs(t4).ge.one)goto 1000
      if(t4.le.(five*t3*t3-one)/four)goto 1000
      if(t4.ge.(five*t3*t3+one)/six )goto 1010

C     Set starting values for n-r iteration:
C     g is chosen to give the correct value of tau-3 on the
C     assumption that h=1 (i.e. a generalized pareto fit) -
C     but h is actually set to 1.001 to avoid numerical
C     difficulties which can sometimes arise when h=1 exactly
      g=(one-three*t3)/(one+t3)
      h=hstart
      z=g+h*p725
      xdist=big

C     Start of Newton-Raphson iteration
      do 100 it=1,maxit
C        reduce steplength until we are nearer to the required
C        values of tau-3 and tau-4 than we were at the previous step
        do 40 i=1,maxsr
C         - calculate current tau-3 and tau-4
C           notation:
C           u.    - ratios of gamma functions which occur in the pwm's
C                   beta-sub-r
C           alam. - l-moments (apart from a location and scale shift)
C           tau.  - l-moment ratios
C
          if(g.gt.oflgam)goto 1020
          if(h.gt.zero)goto 20
C         h <= 0.
          u1=dexp(dlgama(  -one/h-g)-dlgama(  -one/h+one))
          u2=dexp(dlgama(  -two/h-g)-dlgama(  -two/h+one))
          u3=dexp(dlgama(-three/h-g)-dlgama(-three/h+one))
          u4=dexp(dlgama( -four/h-g)-dlgama( -four/h+one))
          goto 30
C         h > 0.
   20     u1=dexp(dlgama(  one/h)-dlgama(  one/h+one+g))
          u2=dexp(dlgama(  two/h)-dlgama(  two/h+one+g))
          u3=dexp(dlgama(three/h)-dlgama(three/h+one+g))
          u4=dexp(dlgama( four/h)-dlgama( four/h+one+g))
   30     continue
          alam2=u1-two*u2
          alam3=-u1+six*u2-six*u3
          alam4=u1-twelve*u2+thirty*u3-twenty*u4
          if(alam2.eq.zero)goto 1020
          tau3=alam3/alam2
          tau4=alam4/alam2
          e1=tau3-t3
          e2=tau4-t4
C         - if nearer than before, exit this loop
          dist=dmax1(dabs(e1),dabs(e2))
          if(dist.lt.xdist)goto 50
C         - Otherwise, halve the steplength and try again
          del1=half*del1
          del2=half*del2
          g=xg-del1
          h=xh-del2
   40   continue
C
C       Too many steplength reductions
        ifail=4
        return

C       Test for convergence
   50   continue
        if(dist.lt.eps)goto 110
C
C       not converged: calculate next step
C       notation:
C       u1g  - derivative of u1 w.r.t. g
C       dl2g - derivative of alam2 w.r.t. g
C       d..  - matrix of derivatives of tau-3 and tau-4 w.r.t. g and h
C       h..  - inverse of derivative matrix
C       del. - steplength
        xg=g
        xh=h
        xz=z
        xdist=dist
        rhh=one/(h*h)
        if(h.gt.zero)goto 60
C       h <= 0.
        u1g=-u1*digamd(  -one/h-g)
        u2g=-u2*digamd(  -two/h-g)
        u3g=-u3*digamd(-three/h-g)
        u4g=-u4*digamd( -four/h-g)
        u1h=      rhh*(-u1g-u1*digamd(  -one/h+one))
        u2h=  two*rhh*(-u2g-u2*digamd(  -two/h+one))
        u3h=three*rhh*(-u3g-u3*digamd(-three/h+one))
        u4h= four*rhh*(-u4g-u4*digamd( -four/h+one))
        goto 70
C       h > 0.
   60   u1g=-u1*digamd(  one/h+one+g)
        u2g=-u2*digamd(  two/h+one+g)
        u3g=-u3*digamd(three/h+one+g)
        u4g=-u4*digamd( four/h+one+g)
        u1h=      rhh*(-u1g-u1*digamd(  one/h))
        u2h=  two*rhh*(-u2g-u2*digamd(  two/h))
        u3h=three*rhh*(-u3g-u3*digamd(three/h))
        u4h= four*rhh*(-u4g-u4*digamd( four/h))
   70   continue
        dl2g=u1g-two*u2g
        dl2h=u1h-two*u2h
        dl3g=-u1g+six*u2g-six*u3g
        dl3h=-u1h+six*u2h-six*u3h
        dl4g=u1g-twelve*u2g+thirty*u3g-twenty*u4g
        dl4h=u1h-twelve*u2h+thirty*u3h-twenty*u4h
        d11=(dl3g-tau3*dl2g)/alam2
        d12=(dl3h-tau3*dl2h)/alam2
        d21=(dl4g-tau4*dl2g)/alam2
        d22=(dl4h-tau4*dl2h)/alam2
        det=d11*d22-d12*d21
        h11= d22/det
        h12=-d12/det
        h21=-d21/det
        h22= d11/det
        del1=e1*h11+e2*h12
        del2=e1*h21+e2*h22
C
C       Take next n-r step
        g=xg-del1
        h=xh-del2
        z=g+h*p725
C       Reduce step if g and h are outside the parameter space
        factor=one
        if(g.le.-one)factor=p8*(xg+one)/del1
        if(h.le.-one)factor=dmin1(factor,p8*(xh+one)/del2)
        if(z.le.-one)factor=dmin1(factor,p8*(xz+one)/(xz-z))
        if(h.le.zero.and.g*h.le.-one)
     *  factor=dmin1(factor,p8*(xg*xh+one)/(xg*xh-g*h))
        if(factor.eq.one)goto 80
        del1=del1*factor
        del2=del2*factor
        g=xg-del1
        h=xh-del2
        z=g+h*p725
   80   continue
C
C       end of newton-raphson iteration
  100 continue
C
C     Not converged
      ifail=3
      return
C
C     Converged
  110 ifail=0
      para(4)=h
      para(3)=g
      temp=dlgama(one+g)
      if(temp.gt.oflexp)goto 1030
      gam=dexp(temp)
      temp=(one+g)*dlog(dabs(h))
      if(temp.gt.oflexp)goto 1030
      hh=dexp(temp)
      para(2)=xmom(2)*g*hh/(alam2*gam)
      para(1)=xmom(1)-para(2)/g*(one-gam*u1/hh)
      return
C
 1000 ifail=1
      return
 1010 ifail=2
      return
 1020 ifail=5
      return
 1030 ifail=6
      return
c
      end


C=============================================================================
C
C Helpers
C

      double precision function dlgama(x)
C
C  LOGARITHM OF GAMMA FUNCTION
C
C  BASED ON ALGORITHM ACM291, COMMUN. ASSOC. COMPUT. MACH. (1966)
C
      implicit double precision (a-h,o-z)
      data small,crit,big,toobig/1d-7,13d0,1d9,2d36/
C
C     c0 is 0.5*log(2*pi)
C     c1...c7 are the coeffts of the asymptotic expansion of dlgama
C
      data c0,c1,c2,c3,c4,c5,c6,c7/
     *   0.91893 85332 04672 742d 0,  0.83333 33333 33333 333d-1,
     *  -0.27777 77777 77777 778d-2,  0.79365 07936 50793 651d-3,
     *  -0.59523 80952 38095 238d-3,  0.84175 08417 50841 751d-3,
     *  -0.19175 26917 52691 753d-2,  0.64102 56410 25641 026d-2/
C
C    S1 is -(Euler's constant), S2 is pi**2/12
C
      data s1/-0.57721 56649 01532 861d 0/
      data s2/ 0.82246 70334 24113 218d 0/
C
      data zero/0d0/,half/0.5d0/,one/1d0/,two/2d0/
      dlgama=zero
      if(x.le.zero)goto 1000
      if(x.gt.toobig)goto 1000
C
C     Use small-x approximation if x is near 0, 1 or 2
C
      if(dabs(x-two).gt.small)goto 10
      dlgama=dlog(x-one)
      xx=x-two
      goto 20
   10 if(dabs(x-one).gt.small)goto 30
      xx=x-one
   20 dlgama=dlgama+xx*(s1+xx*s2)
      return
   30 if(x.gt.small)goto 40
      dlgama=-dlog(x)+s1*x
      return
C
C     reduce to dlgama(x+n) where x+n.ge.crit
C
   40 sum1=zero
      y=x
      if(y.ge.crit)goto 60
      z=one
   50 z=z*y
      y=y+one
      if(y.lt.crit)goto 50
      sum1=sum1-dlog(z)
C
C     Use asymptotic expansion if y.ge.crit
C
   60 sum1=sum1+(y-half)*dlog(y)-y+c0
      sum2=zero
      if(y.ge.big)goto 70
      z=one/(y*y)
      sum2=((((((c7*z+c6)*z+c5)*z+c4)*z+c3)*z+c2)*z+c1)/y
   70 dlgama=sum1+sum2
      return
C
 1000 write(6,7000)x
      return
C
 7000 format(' *** error *** routine dlgama :',
     *  ' argument out of range :',d24.16)
      end



      double precision function derf(x)
C  ERROR FUNCTION
C
C  Based on Algorithm 5666, J.F.Hart et al. (1968) 'Computer
C  Approximations'
C
C  Accurate to 15 decimal places
C
      implicit double precision (a-h, o-z)
      data zero/0d0/,one/1d0/,two/2d0/,three/3d0/,four/4d0/,p65/0.65d0/
C
C     Coefficients of rational-function approximation
      data p0,p1,p2,p3,p4,p5,p6/
     *  0.22020 68679 12376 1D3,    0.22121 35961 69931 1D3,
     *  0.11207 92914 97870 9D3,    0.33912 86607 83830 0D2,
     *  0.63739 62203 53165 0D1,    0.70038 30644 43688 1D0,
     *  0.35262 49659 98910 9D-1/
      data q0,q1,q2,q3,q4,q5,q6,q7/
     *  0.44041 37358 24752 2D3,   0.79382 65125 19948 4D3,
     *  0.63733 36333 78831 1D3,   0.29656 42487 79673 7D3,
     *  0.86780 73220 29460 8D2,   0.16064 17757 92069 5D2,
     *  0.17556 67163 18264 2D1,   0.88388 34764 83184 4D-1/
C
C     c1 is sqrt(2), c2 is sqrt(2/pi)
C     big is the point at which derf=1 to machine precision
C
      data c1/1.4142 13562 37309 5d0/
      data c2/7.9788 45608 02865 4d-1/
      data big/6.25d0/,crit/5d0/
C
      derf=zero
      if(x.eq.zero)return
      xx=dabs(x)
      if(xx.gt.big)goto 20
      expntl=dexp(-x*x)
      zz=dabs(x*c1)
      if(xx.gt.crit)goto 10
      derf=expntl*((((((p6*zz+p5)*zz+p4)*zz+p3)*zz+p2)*zz+p1)*zz+p0)/
     *  (((((((q7*zz+q6)*zz+q5)*zz+q4)*zz+q3)*zz+q2)*zz+q1)*zz+q0)
      if(x.gt.zero)derf=one-two*derf
      if(x.lt.zero)derf=two*derf-one
      return
C
   10 derf=expntl*c2/(zz+one/(zz+two/(zz+three/(zz+four/(zz+p65)))))
      if(x.gt.zero)derf=one-derf
      if(x.lt.zero)derf=derf-one
      return
C
   20 derf=one
      if(x.lt.zero)derf=-one
      return
      end


      double precision function digamd(x)
C  Digamma function (Euler's Psi function) - the first derivative of
C  log(gamma(x))
C
C  based on Algorithm AS103, Appl. Statist. (1976) vol.25 no.3
C
      implicit double precision (a-h,o-z)
      data zero/0d0/,half/0.5d0/,one/1d0/
      data small/1d-9/,crit/13d0/
C
C     C1...C7 are the coeffts of the asymptotic expansion of digamd
C     D1 is  -(Euler's constant)
C
      data c1,c2,c3,c4,c5,c6,c7,d1/
     *  0.83333 33333 33333 333D-1,  -0.83333 33333 33333 333D-2,
     *  0.39682 53968 25396 825D-2,  -0.41666 66666 66666 666D-2,
     *  0.75757 57575 75757 575D-2,  -0.21092 79609 27960 928D-1,
     *  0.83333 33333 33333 333D-1,  -0.57721 56649 01532 861D 0/
      digamd=zero
      if(x.le.zero)goto 1000
C
C     Use small-x approximation if x.le.small
      if(x.gt.small)goto 10
      digamd=d1-one/x
      return
C
C     Reduce to digamd(x+n) where x+n.ge.crit
   10 y=x
   20 if(y.ge.crit)goto 30
      digamd=digamd-one/y
      y=y+one
      goto 20
C
C     Use asymptotic expansion if y.ge.crit
C
   30 digamd=digamd+dlog(y)-half/y
      y=one/(y*y)
      sum=((((((c7*y+c6)*y+c5)*y+c4)*y+c3)*y+c2)*y+c1)*y
      digamd=digamd-sum
      return
C
 1000 write(6,7000)x
      return
c
 7000 format(' *** error *** routine digamd :',
     *  ' argument out of range :',d24.16)
      end
