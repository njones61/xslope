c   Version date 4 May 2005
c
      program Seep2D

C***********************************************************************
C                                                                      *
C                                                                      *
C                             S E E P 2 D                              *
C                            (Version 3.0)                             *
C                                                                      *
C              Seepage Analysis Developed by Fred Tracy                *
C            Engineering Research and Development Center               *
C                         Vicksburg, MS 39180                          *
C                     E-Mail: tracyf@wes.army.mil                      *
C              World-Wide-Web: http://hlnet.wes.army.mil               *
C                                                                      *
C            GMS Interface Developed by R. Jeffrey Davis               *
C                 Engineering Computer Graphics Lab                    *
C                       E-Mail: jeff@byu.edu                           *
C              World-Wide-Web: http://www.ecgl.byu.edu                 *
C                                                                      *
C***********************************************************************
C
C This program is furnished by the Government and is accepted and used
C by the recipient with the express understanding that the United States
C Government makes no warranties, expressed or implied, concerning the
C accuracy, completeness, reliability, usability, or suitability for any
C particular purpose of the information and data contained in this
C program or furnished in connection therewith, and the United States
C shall be under no liability whatsoever to any person by reason of any
C use made thereof.  The program belongs to the Government.  Therefore,
C the recipient further agrees not to assert any proprietary rights
C therein or to represent this program to anyone as other than a
C Government program.
C
C
      use DFPORT  ! emrl jig
      include 'seep.inc'
      implicit real * 8 (a - h, o - z)
      character*80 argv,flname,path
      character*10 sptype
      character*4 fltype
      integer*2 argc,status
      integer isup
c-------ECGL JIG
	logical usegeo
	character*10 stayopen
	common/emrl/stayopen
	integer iu
c-------ECGL JIG

      common  / comm / ibuf(20), ist, ieod
      character*4 ibuf

c-------ECGL JIG
	usegeo = .false.
c-------ECGL JIG

      print 100
  100 format( / )
      print *,'Entering seepage analysis'

      isup=19
      argc=1
c      call getarg(argc,argv,status)
c      if(status .eq. -1) argv = ' '
      numargs = iargc()
	if (numargs.eq.2) then
	  call getarg(2,stayopen)
	else
	  stayopen = '1'
	endif
	if (numargs.gt.0) call getarg(1,argv)
      flname = argv
	if(flname(1:14).eq.'-getArraySizes') then
	  iu = 77
	  open(unit=iu,file='arraySizes.txt',status='replace',
     +       action='write',iostat=ierror)
	  if(ierror.ne.0) stop
	  write(iu,*) 'MaxNode     ', MXNODS
        write(iu,*) 'MaxElement  ', MXELES
        write(iu,*) 'MaxBandwidth', MXBNDW
        write(iu,*) 'MaxMaterials', MXMATS
c690     format(1x,'(a)',i10)
	  stop
      else if(flname(1:1) .eq. ' ') then
666     print *, 'Enter the name of the Seep2D super file'
        read (*,'(a)') flname
      endif
      call getpath(flname,path)
      open(unit=isup,file=flname,status='old',err=666)
      read(isup, 667)sptype
667   format(a)
668   format(a,a)
      if (sptype(1:7) .ne. 'SEEPSUP') then
        print *, 'This file is not a GMS Seep2D superfile'
c-------ECGL JIG
        call stopfile
c-------ECGL JIG
        stop
      endif

669   read(isup,668,end=672) fltype,flname
c-------ECGL JIG
      if((fltype(1:4).ne.'SEEP').and.(fltype(1:4).ne.'ODAT')
     1.and.(fltype(1:4).ne.'OGEO').and.(fltype(1:4).ne.'DSET'))
     1go to 669
c-------ECGL JIG
      call setpath(path,flname)
      if(fltype(1:4).eq.'SEEP') then
        open(15,file=flname,err=670,form='formatted',status='old')
      elseif(fltype(1:4).eq.'ODAT') then
        open(16,file=flname,err=670,form='formatted',
     &       status='unknown')
      elseif(fltype(1:4).eq.'OGEO') then
        open(21,file=flname,err=670,form='formatted',
     &       status='unknown')
	  usegeo = .true.
      elseif(fltype(1:4).eq.'DSET') then
        open(22,file=flname,err=670,form='unformatted',status='unknown')
      else
      endif
      go to 669
c
670   write(*,671) flname
671   format(1x,'Error: file cannot be found or opened =>',a80)
c-------ECGL JIG
        call stopfile
c-------ECGL JIG
      stop

672	if (usegeo .ne. .true.) then
	  open(21,status='scratch',form='formatted')
	endif
      call sepage()
c           
        call stopfile
      end
c ======================================================================
      subroutine sepage()
c
c
c     2D finite element steady-state seepage program
c
c     Call Fred tracy, (601) 634-4112, for questions.
c
c
c-------ECGL JIG
      use DFLIB
c-------ECGL JIG
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common fx(mxnods), flow(mxnods), x(mxnods), y(mxnods),
     &  xk1(mxmats), xk2(mxmats), ang(mxeles), hlast(mxnods),
     &  count(mxnods), nbc(mxnods), ndmp, numnp, numel, np(5,mxeles),
     &  lplx
      common  / param / datum, nummat, nflcd, hed, flonet
      common  / banarg / mband, numblk, r(mxnods), c(2*mxbndw,mxbndw),
     &                   nd, nd2
      common  / plt / iplt, line
      common / fnet / q, hmin, hmax, smin, smax,
     & iout(mxnods), n1(mxnods), n2(mxnods), iuncf
      common/byuout/velx(mxnods),vely(mxnods),
     &              gradx(mxnods),grady(mxnods),
     &              hedout(mxnods),ibndry(2,mxnods),
     &              floout(mxnods),flonod(mxnods),
     &              out(mxnods),outv(mxnods,3)
      common / unsat / uspar(2, mxmats), iuntyp
      common / scale / sck
      real * 4 out, outv
      character*80 hed
      character*4 lplx
      character*1 flonet
      integer idnew(mxnods)
      integer   iversion,iobty,imodule,isflt,jsflt,isflg,jsflg,iscl,
     &          iseep,ivec,inode,ielem,iname,its,iactts,istate,iend
      real * 4 time
      character fname*40
c
      nd = mxbndw
      nd2 = nd * 2
      maxban = nd
c
      line = 20000
      iplt = 0
c
      open(1,status='scratch',form='unformatted')
      open(2,status='scratch',form='unformatted')
      open(4,status='scratch',form='unformatted')
      open(17,status='scratch',form='unformatted')
c
c         Input parameters.
c
  100 read(15,'(a80)') hed

      read(15,131) numnp, numel, nummat, nflcd, lplx, datum,
     &             flonet, unitwt, iuntyp
      datums = datum
c
c         Error condition.
c
      if ((iuntyp.lt.0) .or. (iuntyp.gt.2)) then
        print*, 'Error condition for iuntyp.  The valid values',
     &          'are 0, 1, and 2.'
c-------ECGL JIG
        call stopfile
c-------ECGL JIG
        stop
      endif
c
      do 2 n = 1,numnp
        ibndry(1,n) = 0
        ibndry(2,n) = 0
 2    continue
c
      call steady (iquit)
c
      mxndid = 0
      do 3 n=1,numnp
        ibndry(1,n) = 0
 3    continue
      do 5 n=1,numel
        do 4 nn=1,4
          ibndry(1,np(nn,n)) = ibndry(1,np(nn,n)) + 1
          if (np(nn,n) .gt. mxndid) then
            mxndid = np(nn,n)
          endif
 4      continue
 5    continue
      nmnode = 0
      do 6 n=1,numnp
        if (ibndry(1,n) .ne. 0) then
          nmnode = nmnode + 1
        endif
 6    continue
c Write out a 2D mesh file for the solution mesh
      do 61 i = 1, numnp
        idnew(i) = -1
  61  continue
      do 63 i = 1, numel
        do 62 j = 1, 4
          idnew(np(j,i)) = 0
  62  continue
  63  continue
      icount = 0
      do 64 i = 1, numnp
        if (idnew(i) .eq. 0) then
          icount = icount + 1
          idnew(i) = icount
        endif
  64  continue
      write(21, '(a)') 'MESH2D'
      do 10 n = 1, numel
        if (np(3,n) .ne. np(4,n)) then
          write(21,1001)'e4q',n,idnew(np(1,n)),idnew(np(2,n)),
     &                  idnew(np(3,n)),idnew(np(4,n)),np(5,n)
        else
          write(21,1002)'e3t',n,idnew(np(1,n)),idnew(np(2,n)),
     &                  idnew(np(3,n)),np(5,n)
        endif
  10  continue
      do 20 n = 1, numnp
        if (idnew(n) .ne. -1) then
          write(21,1003)'nd',idnew(n),x(n),y(n),0.d0
      endif
  20  continue
c Write out the data set file 
      iversion = 3000
      iobty = 100
      imodule = 3
      isflt = 110
      jsflt = 4
      isflg = 120
      jsflg = 4
      iscl = 130
      ivec = 140
      iseep = 145
      inode = 170
      ielem = 180
      iname = 190
      its = 200
      iactts = 220
      istate = 0
      time = 0.d0
      iend = 210
      datum = datums
c
      write(22) iversion,iobty,imodule,isflt,jsflt,isflg,jsflg
      if ((flonet .eq. 'F') .and. (iquit .eq. 0)) then
        write(22) iseep
        write(22) inode,icount
        write(22) ielem,numel
        fname = 'flowlines'
        write(22) iname,fname
        write(22) its,istate,time
        j = 0
        do i=1,numnp
          if (idnew(i) .ne. -1) then
            j = j + 1
            out(j) = floout(i) 
          endif
        enddo
        write(22) (out(i),i=1,j)
        write(22) iend
      end if
c
      write(22) iseep
      write(22) inode,icount
      write(22) ielem,numel
      fname = 'flowrate'
      write(22) iname,fname
      write(22) its,istate,time
      j = 0
      do i=1,numnp
        if (idnew(i) .ne. -1) then
          j = j + 1
          out(j) = flonod(i) 
        endif
      enddo
c
c         Unscale the flowrate data.
c
      do i = 1, j
        out(i) = out(i) * sck
      end do
c
      write(22) (out(i),i=1,j)
      write(22) iend
      write(22) iscl
      write(22) inode,icount
      write(22) ielem,numel
      fname = 'pressure head'
      write(22) iname,fname
      write(22) its,istate,time
      j = 0
      do i=1,numnp
        if (idnew(i) .ne. -1) then
          j = j + 1
          out(j) = hedout(i)+datum-y(i) 
        endif
      enddo
      write(22) (out(i),i=1,j)
      write(22) iend
      write(22) iscl
      write(22) inode,icount
      write(22) ielem,numel
      fname = 'pore pressure'
      write(22) iname,fname
      write(22) its,istate,time
      j = 0
      do i=1,numnp
        if (idnew(i) .ne. -1) then
          j = j + 1
          out(j) = unitwt*(hedout(i)+datum-y(i)) 
        endif
      enddo
      write(22) (out(i),i=1,j)
      write(22) iend
      write(22) ivec
      write(22) inode,icount
      write(22) ielem,numel
      fname = 'velocity'
      write(22) iname,fname
      write(22) its,istate,time
      j = 0
      do i=1,numnp
        if (idnew(i) .ne. -1) then
          j = j + 1
          outv(j,1) = velx(i)
          outv(j,2) = vely(i) 
          outv(j,3) = 0.d0 
        endif
      enddo
c
c         Unscale the velocity data.
c
      do i = 1, j
        do k = 1, 3
          outv(i, k) = outv(i, k) * sck
        end do
      end do
c
      write(22) ((outv(i,k),k=1,3),i=1,j)
      write(22) iend
      write(22) iscl
      write(22) iactts,time
      write(22) inode,icount
      write(22) ielem,numel
      fname = 'total head'
      write(22) iname,fname
      write(22) its,istate,time
      j = 0
      do i=1,numnp
        if (idnew(i) .ne. -1) then
          j = j + 1
          out(j) = hedout(i)
        endif
      enddo
      write(22) (out(i),i=1,j)
      write(22) iend
c
c         Output the gradient.
c
      write (22) ivec
      write (22) inode, icount
      write (22) ielem, numel
      fname = 'gradient'
      write (22) iname, fname
      write (22) its, istate, time
      j = 0
      do i=1, numnp
        if (idnew(i) .ne. -1) then
          j = j + 1
          outv(j, 1) = gradx(i)
          outv(j, 2) = grady(i) 
          outv(j, 3) = 0.0d0 
        end if
      end do
      write (22) ((outv(i, k), k=1, 3), i = 1, j)
      write (22) iend
c
      write (*,*)'Seep2D terminated successfully'
c-------ECGL JIG
        call stopfile
c-------ECGL JIG
      stop
c
 1001 format(a4,1x,i5,1x,i5,1x,i5,1x,i5,1x,i5,1x,i5)
 1002 format(a4,1x,i5,1x,i5,1x,i5,1x,i5,1x,i5)
 1003 format(a2,1x,i5,1x,e16.4,1x,e16.4,1x,e16.4)
 1004 format(e12.5,1x,e12.5)
  131 format(4i5, 1x, a4, f10.0, 4x, a1, f10.0, i5, 2f10.0)
c
      end
c ======================================================================
      subroutine bansol
c
c
c     This subroutine solves the system of equations by
c     Gauss Elimination.
c
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common  / banarg / mm, numblk, b(mxnods), a(2*mxbndw,mxbndw),
     &                   nn, nh
c
      nl = nn + 1
      rewind 2
      rewind 4
      nb = 0
      go to 120
c
c     Reduce equations by blocks
c
c     1. Shift block of equations
c
  100 nb = nb + 1
      do 110 n = 1, nn
      nm = nn + n
      b(n) = b(nm)
      b(nm) = 0.d0
      do 110 m = 1, mm
      a(n,m) = a(nm,m)
  110 a(nm,m) = 0.d0
c
c     2. Read next block of equations into core
c
      if (numblk - nb) 120, 130, 120
  120 read(4) (b(n), (a(n,m), m = 1, mm), n = nl, nh)
      if (nb) 130, 100, 130
c
c     3. Reduce block of equations
c
  130 do 180 n = 1, nn
      if (a(n,1)) 140, 180, 140
  140 b(n) = b(n) / a(n,1)
      do 170 l = 2, mm
      if (a(n,l)) 150, 170, 150
  150 c = a(n,l) / a(n,1)
      i = n + l - 1
      j = 0
      do 160 k = l, mm
      j = j + 1
  160 a(i,j) = a(i,j) - c * a(n,k)
      b(i) = b(i) - a(n,l) * b(n)
      a(n,l) = c
  170 continue
  180 continue
c
c     4. Write block of reduced equations on tape 2
c
      if (numblk - nb) 190, 200, 190
  190 write(2) (b(n), (a(n,m), m = 2, mm), n = 1, nn)
      go to 100
c     back-substitution
  200 do 220 m = 1, nn
      n = nn + 1 - m
      do 210 k = 2, mm
      l = n + k - 1
  210 b(n) = b(n) - a(n,k) * b(l)
      nm = n + nn
      b(nm) = b(n)
  220 a(nm,nb) = b(n)
      nb = nb - 1
      if (nb) 230, 240, 230
  230 if (nb .eq. 1) then
        rewind 2
      else
        backspace 2
      endif
      read(2) (b(n), (a(n,m), m = 2, mm), n = 1, nn)
      if (nb .eq. 1) then
        rewind 2
      else
        backspace 2
      endif
      go to 200
c
c     Order unknowns in b array
c
  240 k = 0
      do 250 nb = 1, numblk
      do 250 n = 1, nn
      nm = n + nn
      k = k + 1
  250 b(k) = a(nm,nb)
c
      return
c
      end
c ======================================================================
      subroutine bcrvrs(iquit)
c
c
c         This subroutine reverses the boundary conditions for the
c         second FEM solution for stream function.
c
c
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common fx(mxnods), flow(mxnods), x(mxnods), y(mxnods),
     &  xk1(mxmats), xk2(mxmats), ang(mxeles), hlast(mxnods),
     &  count(mxnods), nbc(mxnods), ndmp, numnp, numel, np(5,mxeles),
     &  lplx
      common  / param / datum, nummat, nflcd, hed, flonet
      common / fnet / q, hmin, hmax, smin, smax,
     & iout(mxnods), n1(mxnods), n2(mxnods), iuncf
      common / unsat / uspar(2, mxmats), iuntyp
      character*80 hed
      character*4 lplx
      character*1 flonet
c
      iquit = 0
c
c         For the default unconfined flow option, flag the nodes no longer
c         used and set nbc = 1 at these nodes.  Also fix free surface boundary
c         conditions.
c
      if (iuntyp .eq. 0) then
c
        do n = 1, numnp
          if (nbc(n) .eq. 200) then
            nbc(n) = 0
            fx(n) = 0.d0
            flow(n) = 0.d0
          else if (nbc(n) .eq. 600) then
            nbc(n) = 2
          else
            nbc(n) = iout(n)
          end if
          count(n) = 0.d0
        end do
c
        do n = 1, numel
          do i = 1, 4
            node = np(i,n)
            count(node) = count(node) + 1.d0
          end do 
        end do
c
          do n = 1, numnp
            if (count(n).eq.0.d0) then
              nbc(n) = 1
              fx(n) = y(n)
            end if
          end do
c
      end if
c
c         Make sure the triangles are defined with
c         nodes 3 and 4 the same.
c
      do n = 1, numel
        na = np(3, n)
        nb = np(4, n)
        i = 1
        do while ((na .ne. nb) .and. (i .le. 3))
          i = i + 1
          nb = np(1, n)
          do j = 1, 3
            np(j, n) = np(j + 1, n)
          end do
          np(4, n) = nb
          na = np(3, n)
        end do
      end do
c
c         Compute smin.
c
      ymax =  - 1.d30
      do 200 n = 1, numnp
      ymax = dmax1(y(n), ymax)
  200 continue
      smin = ymax + q + 100.d0
c
c         Determine which line segments are on the border.
c
      do n = 1, numnp
        iout(n) = 0
      end do
c
      do 250 n = 1, numel
c
      ic = np(3,n) - np(4,n)
      if (ic) 220, 210, 220
  210 nod = 3
      go to 230
  220 nod = 4
c
  230 do 240 i = 1, nod
      ip1 = mod(i,nod) + 1
      im1 = mod(i+nod-2,nod) + 1
      nn = np(i,n)
      iout(nn) = np(ip1,n) - np(im1,n) + iout(nn)
  240 continue
c
  250 continue
c
      lseg = 0
c
      do 320 n = 1, numel
c
      ic = np(3,n) - np(4,n)
      if (ic) 270, 260, 270
  260 nod = 3
      go to 280
  270 nod = 4
  280 nn = np(1,n)
      io2 = iout(nn)
c
      do 310 i = 1, nod
      ip1 = mod(i,nod) + 1
      io1 = io2
      mm = np(i,n)
      nn = np(ip1,n)
      io2 = iout(nn)
      if (io1) 290, 310, 290
  290 if (io2) 300, 310, 300
  300 lseg = lseg + 1
      n1(lseg) = mm
      n2(lseg) = nn
  310 continue
c
  320 continue
c
c         Delete duplicate line segments.
c
      lsm1 = lseg - 1
c
      do 350 i = 1, lsm1
      ip1 = i + 1
      in1 = n1(i)
      in2 = n2(i)
      do 340 j = ip1, lseg
      jj = iabs(n1(j) - in2) + iabs(n2(j) - in1)
      if (jj) 330, 330, 340
  330 n1(i) = 0
      n2(i) = 0
      n1(j) = 0
      n2(j) = 0
  340 continue
  350 continue
c
      ls = 0
      do 370 l = 1, lseg
      if (n1(l)) 370, 370, 360
  360 ls = ls + 1
      n1(ls) = n1(l)
      n2(ls) = n2(l)
  370 continue
c
c         Sort the line segments in a counterclockwise direction.
c
      lsm1 = ls - 1
c
      do 390 i = 1, lsm1
      in = n2(i)
      ip1 = i + 1
      do 380 j = ip1, ls
      if (n1(j).ne.in) go to 380
      if (j.eq.ip1) go to 390
      in1 = n1(ip1)
      in2 = n2(ip1)
      n1(ip1) = n1(j)
      n2(ip1) = n2(j)
      n1(j) = in1
      n2(j) = in2
      go to 390
  380 continue
      print *,'Flow net option not available.'
      iquit = 1
      return
  390 continue
c
c         If nbc = -1 and flow = 0., set nbc = 0.
c
      do i = 1, numnp
        if ((nbc(i) .eq. -1) .and. (flow(i) .eq. 0.d0)) then
          nbc(i) = 0
        end if
      end do
c
c         Provide new boundary conditions.
c
      n = n1(1)
      fx(n) = smin
c
      do i = 2, ls
c
        n = n1(i)
        ip1 = i + 1
        if (ip1 .gt. ls) ip1 = 1
        np1 = n1(ip1)
        im1 = i - 1
        if (im1 .eq. 0) im1 = ls
        nm1 = n1(im1)
        im2 = i - 2
        if (im2 .lt. 1) im2 = im2 + ls
        nm2 = n1(im2)
        ds1 = dsqrt ((x(nm1) - x(nm2)) ** 2 + (y(nm1) - y(nm2)) ** 2)
        ds2 = dsqrt ((x(n) - x(nm1)) ** 2 + (y(n) - y(nm1)) ** 2)
        ds3 = dsqrt ((x(np1) - x(n)) ** 2 + (y(np1) - y(n)) ** 2)
        q1 = flow(nm1) * 2.d0 / (ds1 + ds2)
        q2 = flow(n) * 2.d0 / (ds2 + ds3)
c
c   If this is a single specified head / flow node between two
c   impervious nodes, make an adjustment.
c
        if ((nbc(nm1) .eq. 0) .and. (nbc(n) .ne. 0) .and.
     &     (nbc(np1) .eq. 0)) then
          bq = flow(n) * 0.5d0
c
c   If this is the second impervious node where a single specified
c   head / flow node is between two impervious node, make an
c   adjustment.
c
        else if ((nbc(nm2) .eq. 0) .and. (nbc(nm1) .ne. 0) .and.
     &     (nbc(n) .eq. 0)) then
          bq = flow(nm1) * 0.5d0
c
c   If this is the first specified head / flow node with an
c   impervious node previous to it, make an adjustment.
c
        else if ((nbc(nm1) .eq. 0) .and. (nbc(n) .ne. 0)) then
          bq = 0.d0
c
c   If this is the second specified head / flow node with an
c   impervious node two nodes back, make an adjustment.
c
        else if ((nbc(n) .ne. 0) .and. (nbc(nm1) .ne. 0) .and.
     &     (nbc(nm2) .eq. 0)) then
          bq = flow(nm1) + q2 * ds2 * 0.5d0
c
c   If this is the first impervious node with a specified 
c   head / flow node previous to it, make an adjustment.
c
        else if ((nbc(n) .eq. 0) .and. (nbc(nm1) .ne. 0)) then
          bq = 0.0d0
c
c   If this is the last specified head / flow node with an
c   impervious node ahead of it, make an adjustment.
c
        else if ((nbc(n) .ne. 0) .and. (nbc(np1) .eq. 0)) then
          bq = q1 * ds2 * 0.5d0 + flow(n)
c
        else
c
          bq = (q1 + q2) * ds2 * .5
c
        end if
c
        fx(n) = fx(nm1) + bq
c
      end do
c
      do i = 1, ls
        n = n1(i)
        nbc(n) = 1
      end do
c
      return
      end
c ======================================================================
      subroutine elflow
c
c     This subroutine calculates the flowrates in the
c     principal directions for each element.
c
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common fx(mxnods), flow(mxnods), x(mxnods), y(mxnods),
     &  xk1(mxmats), xk2(mxmats), ang(mxeles), hlast(mxnods),
     &  count(mxnods), nbc(mxnods), ndmp, numnp, numel, np(5,mxeles),
     &  lplx
      common  / banarg / mband, numblk, r(mxnods), c(2*mxbndw,mxbndw),
     &                   nd, nd2
      common / unsat / uspar(2, mxmats), iuntyp
      common  / plt / iplt, line
      common/byuout/velx(mxnods),vely(mxnods),
     &              gradx(mxnods),grady(mxnods),
     &              hedout(mxnods),ibndry(2,mxnods),
     &              floout(mxnods),flonod(mxnods),
     &              out(mxnods),outv(mxnods,3)
      common / scale / sck
      real * 4 out, outv
      integer icnt(mxnods)
      character*4 lplx
      dimension xx(4), yy(4), hh(4), grad(2), v(2)
c
      dr = datan2(1.d0,0.d0) / 90.d0
      nume = 0
c
      write(16,190)
      do n = 1, numnp
        velx(n) = 0.d0
        vely(n) = 0.d0
        icnt(n) = 0
      end do
c
      do 170 n = 1, numel
c
      xc = 0.d0
      yc = 0.d0
c
c         Determine whether the element is a triangle or quad.
c
      if (np(3, n) .eq. np(4, n)) then
        nod = 3
      else
        nod = 4
      end if
c
      if (iuntyp .gt. 0) go to 110
c
c         Compute free surface factor for steady-state problem.
c         Also compute twice the area of the element.
c
      ithrow = 1
      area2 = 0.0d0
      do i = 1, nod
        ino = np(i, n)
        if (r(ino) .gt. y(ino)) ithrow = 0
        ip1 = i + 1
        if (ip1 .gt. nod) ip1 = 1
        jno = np(ip1, n)
        area2 = (x(ino) - x(jno)) * (y(ino) + y(jno))
     &    + area2
      end do
      if ((ithrow .eq. 1) .or. (dabs (area2) .lt. 1.0d-20))
     &  go to 170
c
  110 nume = nume + 1
c
c         Save the new grid for the new flow net option
c         for unconfined flow problems.
c
      do 112 i = 1, 5
      np(i, nume) = np(i, n)
  112 continue
c
c     Calculate velocity.
c
      do 120 j = 1, nod
        nd1 = np(j,n)
        xx(j) = x(nd1)
        yy(j) = y(nd1)
        hh(j) = r(nd1)
        xc = xc + xx(j)
        yc = yc + yy(j)
  120 continue
      mat = np(5,n)
      thetd = ang(n)
c
      if (nod .eq. 4) then
        call grvel4 (n, xx, yy, hh, thetd, mat, grad, v)
      else
        call grvel3 (n, xx, yy, hh, thetd, mat, grad, v)
      end if
c
      g3 = dsqrt (grad(1) * grad(1) + grad(2) * grad(2))
      v3 = dsqrt (v(1) * v(1) + v(2) * v(2))
      if (dabs (v(2)) + dabs (v(1))) 130, 140, 130
  130 dir = datan2 (v(2), v(1)) / dr + thetd
      go to 150
  140 dir = 0.d0
c
c         Unscale the velocities.
c
  150 v1sc = v(1) * sck
      v2sc = v(2) * sck
      v3sc = v3 * sck
      write(16,180) n, v1sc, v2sc, thetd, v3sc, dir
c
      xc = xc / float (nod)
      yc = yc / float (nod)
      cc = dcos (dir * dr)
      ss = dsin (dir * dr)
      gx = g3 * cc
      gy = g3 * ss
      vx = v3 * cc
      vy = v3 * ss
      do nn = 1, nod
        inode = np(nn, n)
        gradx(inode) = gradx(inode) + gx
        grady(inode) = grady(inode) + gy
        velx(inode) = velx(inode) + vx
        vely(inode) = vely(inode) + vy
        icnt(inode) = icnt(inode) + 1
      end do
c
  160 format(2i5, 2(1pe12.4), 4i5, 2(1pe12.4))
      line = line + 10
  170 continue
c
      do n = 1, numnp
        rdenom = 1.0d0 / dble (icnt(n))
        if (icnt(n) .ne. 0) then
          gradx(n) = gradx(n) * rdenom
          grady(n) = grady(n) * rdenom
          velx(n) = velx(n) * rdenom
          vely(n) = vely(n) * rdenom
        end if
      end do
c
      numel = nume
  180 format(3x, i5, 4e12.3, e14.3)
  190 format(////// 30x, 'Element Flowrates'///
     & 5x, 'Elmt', 6x, 'V1', 10x, 'V2', 5x,
     & 'P-axis ang', 5x, 'Res V', 8x, 'Dir of V' // )
c
      return
      end
c ======================================================================
      subroutine exitpt
c
c
c         This subroutine computes the exit points.
c
c
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common fx(mxnods), flow(mxnods), x(mxnods), y(mxnods),
     &  xk1(mxmats), xk2(mxmats), ang(mxeles), frx(mxnods),
     &  count(mxnods), nbc(mxnods), ndmp, numnp, numel, np(5,mxeles),
     &  lplx
      common  / banarg / mband, numblk, r(mxnods), c(2*mxbndw,mxbndw),
     &                   nd, nd2
      common / fnet / q, hmin, hmax, smin, smax,
     & iout(mxnods), n11(mxnods), n22(mxnods), iuncf
      character*4 lplx
      dimension fry(mxnods)
      equivalence (c, fry)
      dimension xf(3), yf(3), xb(2), yb(2)
c
c         Find the exit point to be computed.
c
  100 do 130 n = 1, numel
c
      do 120 i = 1, 4
c
      ip1 = i + 1
      if (ip1.gt.4) ip1 = 1
      im1 = i - 1
      if (im1.lt.1) im1 = 4
      n1 = np(i,n)
      if (nbc(n1).ne.600) go to 120
      n2 = np(ip1,n)
      n3 = np(im1,n)
c
      if (nbc(n2).ne.200) go to 110
      nb = n1
      nt = n2
      ip =  - 1
      go to 140
c
  110 if (nbc(n3).ne.200) go to 120
      nb = n1
      nt = n3
      ip = 1
      go to 140
c
  120 continue
c
  130 continue
      go to 230
c
c         Define boundary line segment.
c
  140 xb(1) = x(nb)
      yb(1) = y(nb)
      xb(2) = x(nt)
      yb(2) = y(nt)
c
c         Determine next three points on free surface.
c
      ndd = nt
c
      do 200 k = 1, 3
c
      do 180 n = 1, numel
c
      do 170 i = 1, 4
c
      if (np(i,n).ne.ndd) go to 170
      iim1 = i
      iip1 = i
c
      if (ip.lt.1) go to 160
  150 iim1 = iim1 - 1
      if (iim1.lt.1) iim1 = 4
      if (iim1.eq.i) go to 180
      nd1 = np(iim1,n)
      if (nd1.eq.ndd) go to 150
      if ((count(nd1).gt.1.d29).or.(x(nd1).le.x(ndd))) go to 150
      go to 190
c
  160 iip1 = iip1 + 1
      if (iip1.gt.4) iip1 = 1
      if (iip1.eq.i) go to 180
      nd1 = np(iip1,n)
      if (nd1.eq.ndd) go to 160
      if ((count(nd1).gt.1.d29).or.(x(nd1).ge.x(ndd))) go to 160
      go to 190
c
  170 continue
c
  180 continue
c
      xexit = xb(1)
      yexit = yb(1)
      t = 0.d0
      go to 210
c
  190 xf(4-k) = frx(nd1)
      yf(4-k) = fry(nd1)
      ndd = nd1
c
  200 continue
c
c         Compute the extrapolated value of the exit point.
c
      call extrap(xf, yf, xb, yb, xexit, yexit, t)
 210  count(nt) = 0.d0
      frx(nt) = xexit
      fry(nt) = yexit
      nbc(nt) = 600
      count(nb) = 1.d30
      nbc(nb) = iout(nb)
      go to 100
c
  230 return
      end
c ======================================================================
      subroutine extrap (xf, yf, xb, yb, xexit, yexit, t)
c
c
c         This subroutine extrapolates the free surface
c         to find the exit point.  Three poins on the
c         free surface are used for modeling.
c
c
      implicit real * 8 (a - h, o - z)
      dimension xf(3), yf(3), xb(2), yb(2)
c
c         Compute the slope of the line segment intersecting
c         the surface of seepage.
c
      fm1 = (yf(2) - yf(1)) / (xf(2) - xf(1))
      fm2 = (yf(3) - yf(2)) / (xf(3) - xf(2))
      sd = (fm2 - fm1) / (xf(2) - xf(1))
      fm = (xf(3) - xf(2)) * sd + fm2
c
c         Compute intersection of the line segment
c         with the surface of seepage.
c
      dx = xb(2) - xb(1)
      dy = yb(2) - yb(1)
      a = dy - fm * dx
      if (a.eq.0.d0) go to 100
      b = (xb(1) - xf(3)) * fm + yf(3) - yb(1)
      t = b / a
      if ((t.lt.0.d0).or.(t.gt.1.d0)) go to 100
      go to 110
  100 t = 0.05d0
  110 xexit = xb(1) + dx * t
      yexit = yb(1) + dy * t
c
      return
      end
c ======================================================================
      subroutine flows
c
c
c     This subroutine computes the flows.
c
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common fx(mxnods), flow(mxnods), x(mxnods), y(mxnods),
     &  xk1(mxmats), xk2(mxmats), ang(mxeles), hlast(mxnods),
     &  count(mxnods), nbc(mxnods), ndmp, numnp, numel, np(5,mxeles),
     &  lplx
      common  / banarg / mm, numblk, r(mxnods), c(2*mxbndw,mxbndw),
     &                   nd, nd2
      character*4 lplx
c
      nl =  - nd + 1
      n2 = nd + 1
      do 100 n = 1, numnp
      flow(n) = 0.d0
  100 continue
c
      do 190 k = 1, numblk
c
      nl = nl + nd
      nb = nl + nd - 1
      kshift = nb - nd2
c
      read(1) ((c(n,m), n = n2, nd2), m = 1, mm)
c
      do 160 n = nl, nb
      if (n - numnp) 110, 110, 160
  110 if (nbc(n)) 120, 160, 130
  120 flow(n) = fx(n)
      go to 160
  130 ist = max0(n-mm+1,1)
c
      do 140 i = ist, n
      jc = n - i + 1
      ic = i - kshift
      flow(n) = c(ic,jc) * r(i) + flow(n)
  140 continue
c
      n3 = min0(numnp-n+1,mm)
      ic = n - kshift
      do 150 jc = 2, n3
      jd = jc + n - 1
      flow(n) = c(ic,jc) * r(jd) + flow(n)
  150 continue
c
  160 continue
      if (k - numblk) 170, 190, 190
c
  170 do 180 n = 1, nd
      nn = n + nd
      do 180 m = 1, mm
      c(n,m) = c(nn,m)
  180 continue
c
  190 continue
      return
      end
c ======================================================================
      subroutine frees
c
c
c         This subroutine computes the position of the free
c         surface.
c
c
c         Description of nbc values of 100 or greater.
c
c         100 - above the free surface (fs)
c         200 - on fs with flow not given
c         300 - on fs with flow given
c         400 - below fs with flow not given
c         500 - below fs with flow given
c         600 - exit point node
c
c
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common fx(mxnods), flow(mxnods), x(mxnods), y(mxnods),
     &  xk1(mxmats), xk2(mxmats), ang(mxeles), frx(mxnods),
     &  count(mxnods), nbc(mxnods), ndmp, numnp, numel, np(5,mxeles),
     &  lplx
      common  / banarg / mband, numblk, r(mxnods), c(2*mxbndw,mxbndw),
     &                   nd, nd2
      dimension fry(5*mxnods)
      equivalence (c, fry)
      dimension node(4)
c
c         Initialize data.
c
      do 100 i = 1, numnp
      count(i) = 1.d30
      frx(i) = 0.d0
      fry(i) = 0.d0
  100 continue
c
c         Consider each element.
c
      do 190 n = 1, numel
c
      pmax = 0.d0
c
      do 180 i = 1, 4
c
      ip1 = i + 1
      if (ip1.gt.4) ip1 = 1
      im1 = i - 1
      if (im1.lt.1) im1 = 4
      n1 = np(i,n)
      n2 = np(ip1,n)
      n3 = np(im1,n)
      ph1 = r(n1) - y(n1)
      ph2 = r(n2) - y(n2)
      ph3 = r(n3) - y(n3)
      pmax = dmax1(pmax,ph1)
c
      if (ph1) 120, 110, 180
c
c         Fix the headwater free surface node.
c
  110 if (nbc(n1).ne.1) go to 180
      frx(n1) = x(n1)
      fry(n1) = y(n1)
      count(n1) = 0.d0
      go to 180
c
c         Consider line segment n1-n2.
c
  120 if (ph2) 150, 130, 140
  130 if ((nbc(n1).ne.2).or.((nbc(n2).ne.1).and.
     &    (nbc(n2).ne.2))) go to 150
      nbc(n1) = 200
      nbc(n2) = 600
      frx(n2) = x(n2)
      fry(n2) = y(n2)
      count(n2) = 0.d0
      go to 180
  140 s = ph1 / (ph1 - ph2)
      dx = (x(n2) - x(n1)) * s
      dy = (y(n2) - y(n1)) * s
      rr = dx * dx + dy * dy
      if (rr.ge.count(n1)) go to 150
      frx(n1) = x(n1) + dx
      fry(n1) = y(n1) + dy
      count(n1) = rr
c
c         Consider line segment n1-n3.
c
  150 if (ph3) 180, 160, 170
  160 if ((nbc(n1).ne.2).or.((nbc(n3).ne.1).and.
     &   (nbc(n3).ne.2))) go to 180
      nbc(n1) = 200
      nbc(n3) = 600
      frx(n3) = x(n3)
      fry(n3) = y(n3)
      count(n3) = 0.d0
      go to 180
  170 s = ph1 / (ph1 - ph3)
      dx = (x(n3) - x(n1)) * s
      dy = (y(n3) - y(n1)) * s
      rr = dx * dx + dy * dy
      if (rr.ge.count(n1)) go to 180
      frx(n1) = x(n1) + dx
      fry(n1) = y(n1) + dy
      count(n1) = rr
c
  180 continue
c
  190 continue
c
c         Do the exit point computations.
c
      call exitpt
c
c         Set flags.
c
      do 300 i = 1, numnp
      nbi = nbc(i)
      if (count(i).gt.1.d29) go to 280
      eps = dabs(frx(i) - x(i)) + dabs(fry(i) - y(i))
      if ((eps.le.1.d-4).and.(nbi.eq.200)) nbi = 300
      if (nbi.eq.0) nbi = 200
      if (nbi.eq.1) nbi = 300
      if (nbi.eq.2) nbi = 200
      count(i) = 1.d0
      x(i) = frx(i)
      y(i) = fry(i)
      r(i) = y(i)
      go to 290
  280 count(i) = 0.d0
      nbi = 400
      if (nbc(i).ne.0) nbi = 500
      if (r(i).lt.y(i)) nbi = 100
  290 nbc(i) = nbi
  300 continue
c
c         Fix elements close to the free surface.
c
      do n = 1, numel
c
        iflag = 1
        i = 1
        do ii = 1, 4
           node(ii) = np(ii, n)
        end do
c
        do while ((iflag .eq. 1) .and. (i .le. 4))
c
          i = i + 1
          n1 = node(1)
          n2 = node(2)
          n3 = node(3)
          n4 = node(4)
c
c  Look for an element where nodes 1 and 3 are on the free surface,
c  node 2 is below the free surface, and node 4 is above the free
c  surface.
c
          if ((count(n1) .eq. 1.d0) .and. (count(n3) .eq. 1.d0) .and.
     &    (r(n2) .gt. y(n2)) .and. (r(n4) .lt. y(n4))) then
c
            iflag = 0
c
c         Replace node n4 with n3 everywhere.
c
            do nn = 1, numel
              do  ii = 1, 4
                if (np(ii, nn) .eq. n4) then
                  np(ii, nn) = n3
                end if
              end do
            end do
c
          end if
c
c         Rotate the nodes.
c
          nsav = node(1)
          do ii = 1, 3
            node(ii) = node(ii + 1)
          end do
          node(4) = nsav
c
        end do
c
      end do
c
      return
      end
c ======================================================================
      subroutine modify(ia)
c
c
c     This subroutine modifies the stiffness matrix for the
c     boundary condition of specified head and specified flow.
c
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common fx(mxnods), flow(mxnods), x(mxnods), y(mxnods),
     &  xk1(mxmats), xk2(mxmats), ang(mxeles), hlast(mxnods),
     &  count(mxnods), nbc(mxnods), ndmp, numnp, numel, np(5,mxeles),
     &  lplx
      common  / banarg / mm, numblk, r(mxnods), c(2*mxbndw,mxbndw),
     &                   nd, nd2
c
      do 100 i = 1, nd2
      r(i) = 0.d0
  100 continue
c
      kshift =  - nd
      n3 = nd + 1
c
      do 270 kk = 1, numblk
c
      read(1) ((c(n,m), n = n3, nd2), m = 1, mm)
c
      do 150 n = n3, nd2
      i = n + kshift
      if (i - numnp) 110, 110, 160
  110 if (nbc(i).le.0) go to 150
      if (nbc(i).ne.2) go to 120
      if (count(i).eq.0.d0) go to 150
      h = y(i)
      go to 130
  120 h = fx(i)
  130 if (ia.gt.1) h = 0.d0
      n1 = n + 1
      do 140 k = 2, mm
      n2 = n1 - k
      if (n2.le.0) go to 140
      r(n2) =  - c(n2,k) * h + r(n2)
      c(n2,k) = 0.d0
  140 continue
      c(n,1) = 1.d0
      r(n) = h
  150 continue
c
  160 if (kk - 1) 180, 180, 170
  170 write(4) (r(n), (c(n,m), m = 1, mm), n = 1, nd)
c
c     Shift block of equations.
c
  180 do 190 n = 1, nd
      nn = nd + n
      r(n) = r(nn)
      r(nn) = 0.d0
      do 190 m = 1, mm
      c(n,m) = c(nn,m)
      c(nn,m) = 0.d0
  190 continue
c
      kshift = kshift + nd
c
      do 260 n = 1, nd
      i = n + kshift
      if (i - numnp) 200, 200, 270
  200 if (ia.gt.1) r(n) =  - flow(i)
      if (nbc(i)) 210, 260, 220
  210 r(n) = fx(i) + r(n)
      go to 260
  220 if (nbc(i).ne.2) go to 230
      if (count(i).eq.0.d0) go to 260
      h = y(i)
      go to 240
  230 h = fx(i)
  240 if (ia.gt.1) h = 0.d0
      n1 = n + 1
      do 250 k = 2, mm
      m = n + k - 1
      if (m.gt.numnp) go to 250
      r(m) =  - c(n,k) * h + r(m)
      c(n,k) = 0.d0
  250 continue
      c(n,1) = 1.d0
      r(n) = h
  260 continue
c
  270 continue
c
      write(4) (r(n), (c(n,m), m = 1, mm), n = 1, nd)
c
      return
      end
c ======================================================================
      subroutine qdflow(ia, itime)
c
c
c         This subroutine calculates the individual element
c         stiffnesses and the assembled matrix for each iteration.
c
c
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common fx(mxnods), flow(mxnods), x(mxnods), y(mxnods),
     &  xk1(mxmats), xk2(mxmats), ang(mxeles), hlast(mxnods),
     &  count(mxnods), nbc(mxnods), ndmp, numnp, numel, np(5,mxeles),
     &  lplx
      common  / banarg / mband, numblk, r(mxnods), c(2*mxbndw,mxbndw),
     &                   nd, nd2
      common / unsat / uspar(2, mxmats), iuntyp
      dimension xk(2,2), gg(5), s(5,5)
      dimension fji(2,2), pnst(2,5), temp(2,2), pres(4,4)
      dimension wf(4), pt(4)
      character*4 lplx
      data wf(4), wf(3)/.34785485, .65214516/
      data pt(4), pt(3)/.86113631, .33998104/
c
      wf(1) = wf(4)
      wf(2) = wf(3)
      pt(1) =  - pt(4)
      pt(2) =  - pt(3)
      ndmp = 0
      pi = datan2(1.d0,0.d0) * 2.d0
      a1 = 1.d0
c
c     Set the pressure head at the integration points to some
c     positive number.
c
      do 110 i = 1, 4
      do 100 j = 1, 4
      pres(i, j) = 1.d0
  100 continue
  110 continue
c
c     Initialize the system of equations.
c
      if (ia .eq. 0) go to 170
      if ((iuntyp .eq. 0) .and. (ia .gt. 1)) go to 140
      numblk = 0
      nl =  - nd + 1
c
      do 130 ii = 1, nd2
      do 130 jj = 1, nd
  130 c(ii,jj) = 0.d0
      if (iuntyp .eq. 0) go to 160
c
  140 do 150 ii = 1, numnp
      flow(ii) = 0.d0
  150 continue
      if (iuntyp .eq. 0) go to 170
c
c     Compute and assemble the element stiffnesses.
c
c
c     Form stiffness matrix by blocks.
c
  160 numblk = numblk + 1
      nl = nl + nd
      nm = nl + nd - 1
      nh = nm + nd
      kshift = nl - 1
  170 rewind 17
c
      do 500 nn = 1, numel
c
      mat = iabs(np(5,nn))
c
      if (ia) 290, 290, 180
  180 read(17) s
c
      if ((iuntyp .eq. 0) .and. (ia. gt. 1)) go to 230
      if (np(5,nn)) 500, 500, 190
  190 do 210 i = 1, 4
      if (np(i,nn) - nl) 210, 200, 200
  200 if (np(i,nn) - nm) 220, 220, 210
  210 continue
      go to 500
c
  220 np(5,nn) =  - np(5,nn)
c
c         Compute presure heads at the integration points.
c
  230 pmin = 1.d30
      pmax =  - pmin
      do 260 nt = 1, 4
      ttt = pt(nt)
      do 250 ns = 1, 4
      sss = pt(ns)
      pp = 0.d0
      do 240 i = 1, 4
      si =  - i * i + i * 5 - 5
      ti = i / 3 * 2 - 1
      fni = (sss * si + 1.d0) * (ttt * ti + 1.d0) * 0.25d0
      nd1 = np(i,nn)
      pp = (hlast(nd1) - y(nd1)) * fni + pp
  240 continue
      pres(nt,ns) = pp
      pmin = dmin1(pmin,pp)
      pmax = dmax1(pmax,pp)
  250 continue
  260 continue
c
c         Check for all positive pressure heads.
c
      if (pmin.ge.0.d0) go to 450
c
c         Check for some negative pressure heads.
c
      if ((iuntyp .gt. 0) .or. (pmax .gt. 0.d0)) go to 290
c
      do 280 i = 1, 4
      do 270 j = 1, 4
      s(i,j) = s(i,j) * uspar(1, mat)
  270 continue
  280 continue
      go to 450
c
c     Calculate element stiffness.
c
  290 xk3 = ang(nn) * pi / 180.d0
      xk(1,1) = xk1(mat) * dcos(xk3) ** 2 + xk2(mat) * dsin(xk3) ** 2
      xk(2,2) = xk2(mat) * dcos(xk3) ** 2 + xk1(mat) * dsin(xk3) ** 2
      xk(2,1) = dsin(xk3) * dcos(xk3) * (xk1(mat) - xk2(mat))
      xk(1,2) = xk(2,1)
c
c         Nodify for flow net option.
c
      if (itime .eq. 2) then
        det = xk(1, 1) * xk (2, 2) - xk( 1, 2) * xk(1, 2)
        xk(1,1) = xk(1, 1) / det
        xk(2, 2) = xk(2, 2) / det
        xk(1, 2) = xk(1, 2) / det
        xk(2, 1) = xk(2, 1) / det
      end if
c
      do 300 ii = 1, 5
      gg(ii) = 0.d0
      do 300 jj = 1, 5
  300 s(ii,jj) = 0.d0
c
      do 410 nt = 1, 4
      ttt = pt(nt)
      do 410 ns = 1, 4
      sss = pt(ns)
c
      fsfact = fkrel (pres(nt, ns), mat)
c
c         Modify for flow net option.
c
      if (itime .eq. 2) then
        if (fsfact .gt. 0.0d0) then
          fsfact = 1.0d0 / fsfact
        else
          fsfact = 1.0d10
        end if
      end if
c
      do 310 i = 1, 4
      si =  - i * i + i * 5 - 5
      ti = i / 3 * 2 - 1
      pnst(1,i) = (ttt * ti + 1.d0) * si * 0.25d0
      pnst(2,i) = (sss * si + 1.d0) * ti * 0.25d0
  310 continue
      pnst(1,5) = (ttt * ttt - 1.d0) * sss * 2.d0
      pnst(2,5) = (sss * sss - 1.d0) * ttt * 2.d0
      do 320 i = 1, 2
      do 320 j = 1, 2
      fji(i,j) = 0.d0
  320 temp(i,j) = 0.d0
      if(lplx.eq.'AXSY') a1=0.d0
  340 do 360 i = 1, 4
      nd1 = np(i,nn)
      xn1 = x(nd1)
      yn1 = y(nd1)
      fji(1,1) = pnst(2,i) * yn1 + fji(1,1)
      fji(1,2) =  - pnst(1,i) * yn1 + fji(1,2)
      fji(2,1) =  - pnst(2,i) * xn1 + fji(2,1)
      fji(2,2) = pnst(1,i) * xn1 + fji(2,2)
      if (lplx.eq.'AXSY') then
        si =  - i * i + i * 5 - 5
        ti = i / 3 * 2 - 1
        fni = (sss * si + 1.d0) * (ttt * ti + 1.d0) * 0.25d0
        a1 = fni * xn1 + a1
      endif
  360 continue
      det = fji(1,1) * fji(2,2) - fji(1,2) * fji(2,1)
      if (det) 370, 550, 370
  370 do 380 i = 1, 2
      do 380 j = 1, 2
      tem1 = fji(j,i)
      do 380 k = 1, 2
      tem2 = xk(j,k) * tem1
      do 380 l = 1, 2
      temp(i,l) = fji(k,l) * tem2 + temp(i,l)
  380 continue
c
      qid = wf(ns) * wf(nt) * a1 * fsfact / det
      do 390 i = 1, 2
      do 390 j = 1, 2
      temp(i,j) = temp(i,j) * qid
  390 continue
c
      do 400 i = 1, 5
      do 400 j = 1, 2
      tem1 = pnst(j,i)
      do 400 k = 1, 2
      tem2 = temp(j,k) * tem1
      do 400 l = 1, i
      s(i,l) = pnst(k,l) * tem2 + s(i,l)
  400 continue
  410 continue
c
      do 420 i = 1, 4
      ip1 = i + 1
      do 420 j = ip1, 5
      s(i,j) = s(j,i)
  420 continue
c
c     Decompose element stiffness.
c
      do 430 ii = 1, 4
      gg(ii) = gg(ii) - s(ii,5) * gg(5) / s(5,5)
      do 430 jj = 1, 4
  430 s(ii,jj) = s(ii,jj) - s(ii,5) * s(5,jj) / s(5,5)
c
      if (ia) 440, 440, 450
  440 write(17) s
      go to 500
c
c     Assemble element stiffness matrix into the system.
c
  450 if ((iuntyp .eq. 0) .and. (ia .gt. 1)) go to 470
      do 460 ii = 1, 4
      l1 = np(ii,nn) - kshift
      do 460 jj = 1, 4
      l2 = np(jj,nn) - kshift
      k2 = l2 - l1 + 1
      if (k2.le.0) go to 460
      c(l1,k2) = c(l1,k2) + s(ii,jj)
  460 continue
      if ((iuntyp .eq. 0) .or. (ia .eq. 1)) go to 500
c
c         Compute flows.
c
  470 do 490 ii = 1, 4
      sum = 0.d0
      do 480 jj = 1, 4
      l2 = np(jj,nn)
      sum = s(ii,jj) * r(l2) + sum
  480 continue
      l1 = np(ii,nn)
      flow(l1) = flow(l1) + sum
  490 continue
c
  500 continue
c
      if (ia .eq. 0) go to 570
      if (ia .eq. 1) go to 510
      if ((iuntyp .eq. 0) .and. (ia .gt. 1)) go to 570
c
c     Store unmodified stiffness matrix.
c
  510 write(1) ((c(n,m), n = 1, nd), m = 1, mband)
c
      do 520 n = 1, nd
      k = n + nd
      do 520 m = 1, mband
      c(n,m) = c(k,m)
      c(k,m) = 0.d0
  520 continue
c
c     Check for last block.
c
      if (nm - numnp) 160, 530, 530
c
c     Correct material numbers.
c
  530 do 540 n = 1, numel
      np(5,n) = iabs(np(5,n))
  540 continue
      go to 570
c
c     Error condition.
c
  550 write (16, 560) nn
  560 format(/ ' Zero or negative area element', i8 /)
      ndmp = 1
      do i = 1, 4
        nd1 = np(i, nn)
        xn1 = x(nd1)
        yn1 = y(nd1)
        write (16, *) 'node, x, y', nd1, xn1, yn1
      end do
c
  570 return
      end
c ======================================================================
      subroutine steady (iquit)
c
c
c         This subroutine performs a steady-state seepage
c         analysis.
c
c
c-------ECGL JIG
      use DFLIB
c-------ECGL JIG
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common fx(mxnods), flow(mxnods), x(mxnods), y(mxnods),
     &  xk1(mxmats), xk2(mxmats), ang(mxeles), hlast(mxnods),
     &  count(mxnods), nbc(mxnods), ndmp, numnp, numel, np(5,mxeles),
     &  lplx
      common  / param / datum, nummat, nflcd, hed, flonet
      common  / banarg / mband, numblk, r(mxnods), c(2*mxbndw,mxbndw),
     &                   nd, nd2
      common  / plt / iplt, line
      common / fnet / q, hmin, hmax, smin, smax,
     & iout(mxnods), n1(mxnods), n2(mxnods), iuncf
      common / unsat / uspar(2, mxmats), iuntyp
      common/byuout/velx(mxnods),vely(mxnods),
     &              gradx(mxnods),grady(mxnods),
     &              hedout(mxnods),ibndry(2,mxnods),
     &              floout(mxnods),flonod(mxnods),
     &              out(mxnods),outv(mxnods,3)
      common / scale / sck
      real * 4 out, outv
      dimension nsmall(2),ysmall(2)
      dimension angls(mxmats)
      character*80 hed
      character*4 ibl,iast,lplx
      character*1 flonet
      data ibl/'    '/
      data iast/'*   '/
c
      maxban = nd
      noitr = 300
c
      do 100 n = 1, numnp
      hlast(n) = 1.d30
      count(n) = 0.d0
      hedout(n) = 0.d0
      floout(n) = 0.d0
      flonod(n) = 0.d0
  100 continue
      itime = 0
      if (flonet.eq.'F') itime = 1
      iuncf = 0
c
      if (lplx.eq.'PLNE') then
        write(16,870)
       else
        write(16,860)
      endif
  120 write(16,800) hed, numnp, numel, nummat, datum, iuntyp
c
c     Input material properties.
c
      do 122 n = 1, nummat
c
      read (15, 760) i, xk1(n), xk2(n), angls(n), uspar(1, n),
     &     uspar(2, n)
c
c         If zero's are given for the unsaturated flow data, keep
c         the old default.
c
      if ((iuntyp .eq. 0) .and. (uspar(1, n) .eq. 0.d0))
     &   uspar(1, n) = .001
c
c         Check van genuchten parameters for valid values.
c
      if (iuntyp .eq. 2) then
        if (uspar(1, n) .le. 0.d0) then
          write (*, '(a, f10.4, a)') 'alpha for van genuchten = ',
     &      uspar(1, n), ' must be > 0.'
c-------ECGL JIG
        call stopfile
c-------ECGL JIG
          stop
        end if
        if (uspar(2, n) .lt. 1.d0) then
          write (*, '(a, f10.4, a)') 'n fror van genuchten = ',
     &      uspar(2, n), ' must be > 1.d0'
c-------ECGL JIG
        call stopfile
c-------ECGL JIG
          stop
        end if
      end if
c
  122 continue
c
      write (16, 830) (n, xk1(n), xk2(n), angls(n), uspar(1, n),
     &      uspar(2, n), n = 1, nummat)
c
c         Compute the maximum xk value and use it to scale the k's.
c
      sck = - 1.0d0
      do i = 1, nummat
        sck = dmax1 (sck, xk1(i))
        sck = dmax1 (sck, xk2(i))
      end do
c
      sc = 1.0d0 / sck
      do i = 1, nummat
        xk1(i) = xk1(i) * sc
        xk2(i) = xk2(i) * sc
      end do
c
c     Read nodal information.
c
      print *,'Reading nodes'
      read(15,770) m, im, nbc(m), x(m), y(m), fx(m)
      ymin = y(m)
      ymax = ymin
      do 170 i = 1, numnp
        read(15, 770) n, in, nbc(n), x(n), y(n), fx(n)
        ymin = dmin1(y(n), ymin)
        ymax = dmax1(y(n), ymax)
        if (n - m) 130, 140, 140
  130   print 890, m, n
c-------ECGL JIG
        call stopfile
c-------ECGL JIG
        stop
  140   continue
        nmm = n - m
        xnmm = nmm
        dx = (x(n) - x(m)) / xnmm
        dy = (y(n) - y(m)) / xnmm
        df = (fx(n) - fx(m)) / xnmm
        mp1 = m + 1
        if (mp1.ge.n) go to 160
        nm1 = n - 1
        do 150 nn = mp1, nm1
          x(nn) = x(nn-1) + dx
          y(nn) = y(nn-1) + dy
          nbc(nn) = 0
          fx(nn) = 0.d0
          if (im.eq.0) go to 150
          nbc(nn) = nbc(m)
          fx(nn) = fx(nn-1) + df
  150     continue
  160   if (n.ge.numnp) go to 180
        m = n
        im = in
  170 continue
  180 write(16,840) (n, nbc(n), x(n), y(n), fx(n), n = 1, numnp)
      do 200 n = 1, numnp
        if (nbc(n) .eq. 1) fx(n) = fx(n) + datum
        if (nbc(n) .eq. 2) fx(n) = y(n)
  200 continue
c
c         Save nbc data for later use.
c
      do 202 i = 1, numnp
        iout(i) = nbc(i)
  202 continue
c
c     Read element information.
c
      print *,'Reading elements'
      read(15,780)m,(np(i,m), i = 1, 5)
      ang(m) = angls(np(5,m))
      k = 0
      do 220 i1 = 1, 4
        do 210 l1 = 1, 4
          kk = iabs(np(i1,m) - np(l1,m))
          if (k.ge.kk) go to 210
          k = kk
  210   continue
  220 continue
      if (numel.le.1) go to 310
      do 300 j = 1, numel
        read(15,780) n, (np(i,n), i = 1, 5)
        ang(n) = angls(np(5,n))
        if (n - m) 230, 240, 240
  230   print 900, m, n
c-------ECGL JIG
        call stopfile
c-------ECGL JIG
        stop
  240   continue
        do 260 i1 = 1, 4
          do 250 l1 = 1, 4
            kk = iabs(np(i1,n) - np(l1,n))
            if (k.ge.kk) go to 250
            k = kk
  250     continue
  260   continue
        mp1 = m + 1
        if (mp1.ge.n) go to 290
        nm1 = n - 1
        do 280 nn = mp1, nm1
          do 270 i = 1, 4
  270       np(i,nn) = np(i,nn-1) + 1
          np(5,nn) = np(5,m)
          ang(nn) = ang(m)
  280   continue
  290   if (n.ge.numel) go to 310
        m = n
  300 continue
  310 nbndry = 1
      do 6 n = 1,numel
        do 4 nn = 1,4 
          ibond = 1
          in1 = np(nn,n)
          if (nn .eq. 4) then
            in2 = np(1,n)
          else
            in2 = np(nn+1,n)
          endif
          do 8 m = 1, numel
            if (m .ne. n) then
              do 7 mm = 1, 4
                im1 = np(mm,m)
                if (mm .eq. 4) then
                  im2 = np(1,m)
                else
                  im2 = np(mm+1,m)
                endif
                if (in1.eq.im2 .and. in2.eq.im1) then
                  ibond = 0
                endif
 7            continue
            endif
 8        continue
          if (ibond .eq. 1 .and. in1 .ne. in2) then
            ibndry(1,nbndry) = in1
            ibndry(2,nbndry) = in2
            nbndry = nbndry + 1
          endif
 4      continue
 6    continue
      nbndry = nbndry-1
      write(16,851) (i, i = 1, 4)
      do 330 n = 1, numel
        write(16,852) n, (np(i,n), i = 1, 5), ang(n)
  330 continue
      mband = k + 1
c
c     Check to determine if bandwidth is too large.
c
      if (mband.gt.maxban) go to 750
c
c     Read flowrate cards.
c
      if (nflcd.le.0) go to 350
      write(16,810) (i, i = 1, 2)
      do 340 n = 1, nflcd
        read(15,790) i, j, flrt
        write(16,820) i, j, flrt
c
c         Scale the flowrate card.
c
        flrt = flrt / sck
c
        sij = dsqrt((x(j) - x(i)) ** 2 + (y(j) - y(i)) ** 2)
        if (lplx.eq.'PLNE') then
          fx(i) = fx(i) + 0.5d0 * sij * flrt
          fx(j) = fx(j) + 0.5d0 * sij * flrt
        else
          fx(i) = (x(i) * 2.d0 + x(j)) * sij * flrt / 6.d0 + fx(i)
          fx(j) = (x(j) * 2.d0 + x(i)) * sij * flrt / 6.d0 + fx(j)
        endif
        nbc(i) = -1
        nbc(j) = -1
  340 continue
c
  350 call qdflow(0, itime)
c
      if (itime .le. 1) then
        print *,'Solving for heads'
      else if (itime .eq. 2) then
        print *,'Solving for flow lines'
      endif
c
c         Check for the tailwater level at the
c         bottom of the exit face.  Start nbc = 2
c         nodes off with head = elevation if needed.
c
      do 360 l = 1, 2
        nsmall(l) = 0
        ysmall(l) = 1.d30
  360 continue
      hsmall = 1.d30
      do 380 n = 1, numnp
        if (nbc(n).ne.2) go to 380
        hsmall = dmin1(fx(n), hsmall)
        do 370 l = 1, 2
          if (ysmall(l).le.y(n)) go to 370
          ysmall(l) = y(n)
          nsmall(l) = n
          go to 380
  370   continue
  380 continue
      if (nsmall(1).eq.0) go to 400
c
      isum = 0
      do 390 n = 1, numnp
        if (nbc(n).ne.1) go to 390
        if (fx(n).eq.hsmall) isum = isum + 1
        if (isum.eq.2) go to 400
  390 continue
c
      nsm = nsmall(1)
      count(nsm) = 1.d0
      if (isum.eq.1) go to 400
      nsm = nsmall(2)
      if (nsm.ne.0) count(nsm) = 1.d0
c
c         Establish convergence tolerance.
c
  400 eps = (ymax - ymin) * .0001
      relax = 1.d0
c
c         Calculate heads, flows, and flowrates.
c
      do 490 k = 1, noitr
c
c         Calculate assembled unmodified stiffness matrix for first
c         iteration.  Otherwise calculate flows.
c
      rewind 1
      kk = k
      call qdflow(kk, itime)
      if (ndmp) 410, 410, 740
c
c         Modify for boundary conditions.
c
  410 rewind 1
      rewind 4
      call modify(kk)
c
c         Solve system of equations.
c
      call bansol
c
c         Compute flows.
c
      if (k.gt.1) go to 420
      rewind 1
      call flows
      if (itime.eq.2) go to 700
c
c         Check for convergence.
c
  420 ifirst = 1
      dhmax = 0.d0
      if (k .gt. 20) relax = .5
      if (k .gt. 40) relax = .2
      if (k .gt. 60) relax = .1
      if (k .gt. 80) relax = .05
      if (k .gt. 100) relax = .02
      if (k .gt. 120) relax = .01
      do 440 i = 1, numnp
        if (k.gt.1) r(i) = r(i) * relax + hlast(i)
        if (r(i).ge.y(i)) go to 430
        ifirst = 0
        go to 440
  430   dh = dabs(r(i) - hlast(i))
        dhmax = dmax1(dhmax,dh)
  440 continue
      if (ifirst.eq.1) go to 540
c
c         Reset hlast.
c
      do 450 i = 1, numnp
      hlast(i) = r(i)
  450 continue
c
c         Check surface of seepage nodes.
c
      do 480 i = 1, numnp
        if (nbc(i).ne.2) go to 480
        if (count(i)) 460, 460, 470
  460   if (r(i).lt.y(i)) go to 480
        count(i) = 1.d0
        r(i) = y(i)
        hlast(i) = y(i)
        go to 480
  470   if (flow(i).lt.0.d0) go to 480
        count(i) = 0.d0
  480 continue
c
      print *,'Iteration no. =', k, '    Delta h max. =', dhmax
      if (dhmax.gt.eps) iconv = 0
      if (dhmax.le.eps) iconv = iconv + 1
      if (iconv.eq.2) go to 500
c
  490 continue
c
c         Compute free surface.
c
  500 if (iuntyp .gt. 0) go to 540
      call frees
      iuncf = 1
c
c         Print results.
c
  540 write(16,990)
c
      hmin = 1.d20
      hmax =  - 1.d20
      do 550 n = 1, numnp
        rr = dmax1(r(n), y(n))
        hmin = dmin1(rr,hmin)
        hmax = dmax1(rr,hmax)
  550 continue
      qp = 0.d0
      delh = hmax - hmin
      iflag = 1
      do 560 n = 1, numnp
        if (r(n) .lt. y(n)) go to 570
  560 continue
      iflag = 0
  570 if ((iflag .eq. 0) .or. (iuntyp .gt. 0)) go to 650
  580 write(16,940)
      do 645 n = 1, numnp
c
c         Unscale flow.
c
        flowsc = flow(n) * sck
c
        izero = 1
        head = r(n) - datum
        perc = (r(n) - hmin) * 100.d0 / delh
        num = nbc(n) / 100
        go to (590, 600, 610, 620, 630, 600), num
  590   write(16,880) n, iast
        go to 640
  600   write(16,950) n, head, perc, ibl, iast, ibl, 
     &                x(n), y(n)
        go to 640
  610   write(16,960) n, head, perc, flowsc, ibl, iast, ibl,
     &                x(n), y(n)
        if (flow(n).gt.0.d0) qp = qp + flow(n)
        izero = 0
        go to 640
  620   write(16,930) n, head, perc, ibl, ibl, iast
        go to 640
  630   write(16,920) n, head, perc, flowsc, ibl, ibl, iast
        if (flow(n).gt.0.d0) qp = qp + flow(n)
        izero = 0
  640   if (izero .eq. 1) flow(n) = 0.d0
  645 continue
      go to 690
  650 write(16,970)
      do 685 n = 1, numnp
c
c         Unscale flow.
c
        flowsc = flow(n) * sck
c
        izero = 1
        head = r(n) - datum
        perc = (r(n) - hmin) * 100.d0 / delh
        if (nbc(n)) 660, 670, 660
  660   write(16,980) n, head, perc, flowsc
        if (flow(n).gt.0.d0) qp = qp + flow(n)
        izero = 0
        go to 680
  670   write(16,980) n, head, perc
  680   if (izero .eq. 1) flow(n) = 0.d0
  685 continue
  690 continue
c
      q = qp
c
c         Unscale flow sum.
c
        qpsc = qp * sck
c
      write(16,692) qpsc
  692 format( //// 28x, 'Flow = ', 1pe12.4)
c
c         Write post-processor file.
c
  700 continue
      line = line + 10
      do 730 n = 1, numnp
        xx = x(n)
        yy = y(n)
        hh = r(n) - datum
        pp = r(n) - y(n)
        perc = (r(n) - hmin) * 100.d0 / delh
        if (pp.ge.0.d0) go to 710
        if (iuntyp .gt. 0) go to 710
        hh = y(n) - datum
        pp = 0.d0
  710   continue
        if (flonet .eq. 'F') then
          if (itime .eq. 1) then
            hedout(n) = hh
            flonod(n) = flow(n)
          else
            floout(n) = hh
          endif
        else
          hedout(n) = hh
          flonod(n) = flow(n)
        endif
  720   format(2i5, 4(1pe12.4))
        line = line + 10
  730 continue
c
c         Calculate element flowrates.
c
      if (itime .eq. 2) return
      call elflow
c
c         Process flow net option.
c
      if (itime .eq. 0) return
      hmx = hmax - datum
      hmn = hmin - datum

  733 format(i5, 4(1pe12.4))
      line = line + 10
      call bcrvrs(iquit)
      if (iquit .eq. 1) return
      if (iuntyp .eq. 0) then
        do 734 n = 1, numnp
          hlast(n) = 1.d30
  734   continue
      end if
      datum = 0.d0
      itime = 2
      go to 350
c
      line = line + 10
  740 return
c
  750 write(16,910)
      return
c
  760 format(i5, 5f15.0)
  770 format(i5, i2, i3, 3f15.0)
  780 format(6i5)
  790 format(2i5, f10.0)
  800 format(////,1x,a80,/// ,
     & ' Number of nodal points------',i5 /
     & ' Number of elements----------',i5 /
     & ' Number of diff. materials---',i5 /
     & ' Elevation of datum----------',f10.3 /
     & ' Unsaturated flow option-----',i5 //)
  810 format( ////// 26x, 'Specified flowrates'
     &/// 19x, 2(2x, '#', i1), 16x, 'Flowrate' // )
  820 format(19x, 2i4, 14x, f10.3)
  830 format ( ////// 28x, 'Material Properties'
     & // 2x, 'Mat', 8x, 'K1', 11x, 'K2', 9x, 'Angle',
     & 8x, 'Uspar1', 7x, 'Uspar2' // (i5, 5(2x, e11.4)))
  840 format( ////// 25x, 'Node Point Information',
     &  /// 14x, 'Node', 2x, 'BC', 9x, 'X', 11x,
     & 'Y', 6x, 'Flow-head' // (13x, i4, 3x, i2, 3(2x,f10.2)))
  851 format ( ////// 26x, 'Element Information'
     &/// 9x, 'Elmt', 2x, 4(5x, '#', i1), 7x, 'Mat', 6x, 'Angle' /)
  852 format ((8x, i4, 3x, 4(3x, i4), 7x, i2, 6x, f6.1))
  860 format(////,25x, 'Axisymmetric flow problem',////)
  870 format(////,25x, 'Plane flow problem',////)
  880 format(8x, i4, 62x, a1)
  890 format('Error- data card out of order'
     &, ', missplaced node card ',  / , 'node card ',
     & i5, ' found before ', i5)
  900 format('Error- data card out of order'
     &, ', missplaced element card ',  / , 'element card ',
     & i5, ' found before ', i5)
  910 format(' Maximum bandwidth exceeds allowable' // )
  920 format(8x, i4, 5x, e13.4, 9x, f5.1, 1x, '%',
     & 7x, e12.4, 8x, a1, 5x, a1, 4x, a1)
  930 format(8x, i4, 5x, e13.4, 9x, f5.1, 1x, '%',
     & 27x, a1, 5x, a1, 4x, a1)
  940 format( ////// 25x, 'Nodal Flows and Heads',
     & 28x, 'Position of phreatic surface',  /
     &// 34x, 'Percentage of',  / 8x, 'Node', 9x,
     & 'Head', 9x, 'available head', 8x, 'Flow',
     & 11x, 'Above', 2x, 'On', 3x, 'Below', 8x,
     & 'x', 10x, 'y' // )
  950 format(8x, i4, 5x, e13.4, 9x, f5.1, 1x, '%',
     & 27x, a1, 5x, a1, 4x, a1, 3x, 2(4x,f7.2))
  960 format(8x, i4, 5x, e13.4, 9x, f5.1, 1x, '%',
     & 7x, e12.4, 8x, a1, 5x, a1, 4x, a1, 3x, 2(4x,f7.2))
  970 format( ////// 25x, 'Nodal Flows and Heads'
     &/// 39x, 'Percentage of' / 9x, 'Node', 12x,
     & 'Head', 9x, 'available head', 8x, 'Flow' // )
  980 format(8x, i4, 7x, e13.4, 9x, f5.1, 1x, '%', 7x, e12.4)
  990 format( ////// )
      end
c ======================================================================
      subroutine grvel3 (n, xx, yy, hh, angle, mat, grad, v)
c
c
c         This subroutine computes the gradient and discharge velocity
c         for a triangular element containing coordinates (xx, yy) in
c         the principal axes directions.
c
c
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common fx(mxnods), flow(mxnods), x(mxnods), y(mxnods),
     &  xk1(mxmats), xk2(mxmats), ang(mxeles), hlast(mxnods),
     &  count(mxnods), nbc(mxnods), ndmp, numnp, numel, np(5,mxeles),
     &  lplx
      common / unsat / uspar(2, mxmats), iuntyp
      dimension xk(2,2), tmp(2), grad(2), v(2)
      dimension xx(3), yy(3), hh(3)
c
      pi = datan2 (1.d0, 0.d0) * 2.d0
c
c         Calculate the gradient.
c
      do i = 1, 2
        grad(i) = 0.d0
        v(i) = 0.d0
      end do
c
      yc = 0.d0
      hc = 0.d0
c
      do j = 1, 3
        yn1 = yy(j)
        hn1 = hh(j)
        yc = yc + yn1
        hc = hc + hn1
      end do
c
c         Compute constants.
c
      a12 = xx(1) * yy(2) - xx(2) * yy(1)
      a23 = xx(2) * yy(3) - xx(3) * yy(2)
      a31 = xx(3) * yy(1) - xx(1) * yy(3)
      tmp(1) = (yy(2) - yy(3)) * hh(1) + (yy(3) - yy(1)) * hh(2)
     &       + (yy(1) - yy(2)) * hh(3)
      tmp(2) = (xx(3) - xx(2)) * hh(1) + (xx(1) - xx(3)) * hh(2)
     &       + (xx(2) - xx(1)) * hh(3)
c
      det = a12 + a23 + a31
      if (det. eq. 0.d0) then
        print*, 'Element', n, ' with zero area has been encountered.'
        print*, 'Velocities were set to zero.'
        print*, 'xx =', xx
        print*, 'yy =', yy
        return
      end if
c
      theta = angle * pi / 180.d0
      xk(1, 1) =  dcos (theta)
      xk(1, 2) =  dsin (theta)
      xk(2, 1) = - xk(1, 2)
      xk(2, 2) =  xk(1, 1)
c
      do k = 1, 2
        t1 = tmp(k) / det
        do i = 1, 2
          grad(i) = xk(i, k) * t1 + grad(i)
        end do
      end do
c
c         Compute the relative hydraulic conduvtivity at the center
c         of the element.
c
      phc = (hc - yc) / 3.0d0
      fkr = fkrel (phc, mat)
c
c         Compute the discharge velocity.
c
c
      xk11 = xk1(mat)
      xk22 = xk2(mat)
      v(1) = - grad(1) * xk11 * fkr
      v(2) = - grad(2) * xk22 * fkr
c
      return
      end
c
c ======================================================================
      subroutine grvel4 (n, xx, yy, hh, angle, mat, grad, v)
c
c
c         This subroutine computes the gradient and discharge velocity
c         for a quad element containing coordinates (xx, yy) in the
c         principal axes directions.
c
c
      implicit real * 8 (a - h, o - z)
      include 'seep.inc'
      common fx(mxnods), flow(mxnods), x(mxnods), y(mxnods),
     &  xk1(mxmats), xk2(mxmats), ang(mxeles), hlast(mxnods),
     &  count(mxnods), nbc(mxnods), ndmp, numnp, numel, np(5,mxeles),
     &  lplx
      common / unsat / uspar(2, mxmats), iuntyp
      dimension xk(2, 2), pnst(2, 4), tmp(2), grad(2), v(2), fji(2, 2)
      dimension xx(4), yy(4), hh(4)
c
c         Evaluate p at the center of the element.
c
      data pnst /-0.25d0, -0.25d0,  0.25d0, -0.25d0,
     &            0.25d0,  0.25d0, -0.25d0,  0.25d0/
c
      pi = datan2 (1.d0, 0.d0) * 2.d0
c
c         Calculate the gradient.
c
      do i = 1, 2  
        grad(i) = 0.d0
        v(i) = 0.d0
        tmp(i) = 0.d0
          do j = 1, 2
            fji(i, j) = 0.d0
          end do
        end do
c
      if (iuntyp .eq. 0) then
c
c         Check for very small elements.
c
        itim = 0
        amn = 1.d30
        amx = 0.d0
        do i = 1, 4
          ip1 = mod(i,4) + 1
          x1 = xx(i) - xx(ip1)
          y1 = yy(i) - yy(ip1)
          r1 = x1 * x1 + y1 * y1
          amx = dmax1 (r1, amx)
          if (r1 .ne. 0.d0) then
            amn = dmin1 (amn, r1)
          else
            itim = itim + 1
            if (itim .gt. 1) return
          end if
        end do
c
        if (itim .ne. 0) then
          if (amn / amx.lt.1.d - 3) return
        end if
c
      end if
c
      yc = 0.d0
      hc = 0.d0
c
      do j = 1, 4
        xn1 = xx(j)
        yn1 = yy(j)
        hn1 = hh(j)
        yc = yc + yn1
        hc = hc + hn1
        fji(1, 1) = pnst(2,j) * yn1 + fji(1, 1)
        fji(1, 2) =  - pnst(1, j) * yn1 + fji(1, 2)
        fji(2,1) =  - pnst(2, j) * xn1 + fji(2, 1)
        fji(2, 2) = pnst(1, j) * xn1 + fji(2, 2)
        do i = 1, 2
          tmp(i) = pnst(i, j) * hn1 + tmp(i)
        end do
      end do
c
      det = fji(1,1) * fji(2,2) - fji(1,2) * fji(2,1)
      if (det. eq. 0.d0) then
        print*, 'Element', n, ' with zero area has been encountered.'
        print*, 'Velocities were set to zero.'
        print*, 'xx =', xx
        print*, 'yy =', yy
        return
      end if
c
      theta = angle * pi / 180.d0
      xk(1, 1) =  dcos (theta)
      xk(1, 2) =  dsin (theta)
      xk(2, 1) = - xk(1, 2)
      xk(2, 2) =  xk(1, 1)
c
      do k = 1, 2
        t1 = tmp(k) / det
        do j = 1, 2
          t2 = fji(j, k) * t1
          do i = 1, 2
            grad(i) = xk(i, j) * t2 + grad(i)
          end do
        end do
      end do
c
c         Compute the relative hydraulic conduvtivity at the center
c         of the element.
c
      phc = (hc - yc) * 0.25d0
      fkr = fkrel (phc, mat)
c
c         Compute the discharge velocity.
c
      xk11 = xk1(mat)
      xk22 = xk2(mat)
      v(1) = - grad(1) * xk11 * fkr
      v(2) = - grad(2) * xk22 * fkr
c
      return
      end
c----------------------------------------------------------------
      function fkrel (ph, mat)
c
c
c   purpose:  this function computes the relative hydraulic
c             conductivity where the relative hydraulic
c             conductivity vs. pressure head curve is one of three
c             options:  step function, front, and van genuchten
c             model.
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
      implicit real * 8 (a - h, o - z)
c
c   Input variables:
c
      real * 8 ph     ! pressure head.
c
      integer mat     ! material type number.
c
c
c   Output variables:
c
      real * 8 fkrel  ! relative hydraulic conductivity.
c
c
c   Common block variables:
c
c   / unsat /
c
      real * 8 uspar  ! parameters used for modeling the relative
c                       hydraulic conductivity vs. pressure head
c                       curve.
c
      integer iuntyp  ! flag for determining which type of relative
c                       hydraulic conductivity vs. pressure head
c                       curve is to be used.
c                       0 - step function.
c                       1 - front.
c                       2 - van genuchten model.
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
      include 'seep.inc'
      common / unsat / uspar(2, mxmats), iuntyp
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
c         Step function.
c
      if (iuntyp .eq. 0) fkrel = fkrels (ph, mat)
c
c         Front.
c
      if (iuntyp .eq. 1) fkrel = fkrelf (ph, mat)
c
c         Van Genuchten model.
c
      if (iuntyp .eq. 2) fkrel = fkrelv (ph, mat)
c
      return
      end
c----------------------------------------------------------------
      function fkrelf (ph, mat)
c
c
c   Purpose:  This function computes the relative hydraulic
c             conductivity where the relative hydraulic
c             conductivity vs. pressure head curve is a small
c             value for all negative pressure heads until the
c             pressure head reaches a transition value at which
c             time the curve goes linearly to 1 as pressure head
c             goes to 0.  The curve has the shape of a front.
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
      implicit real * 8 (a - h, o - z)
c
c   Input variables:
c
      real * 8 ph     ! pressure head.
c
      integer mat     ! material type number.
c
c
c   Output variables:
c
      real * 8 fkrelf ! relative hydraulic conductivity from the front
c                       model.
c
c
c   Common block variables:
c
c   / unsat /
c
      real * 8 uspar  ! parameters used for modeling the relative
c                       hydraulic conductivity vs. pressure head
c                       curve.  for the front model,
c                       uspar(1) - small relative hydraulic
c                                  conductivity.
c                       uspar(2) - transition negative pressure head
c                                  value.
c
      integer iuntyp  ! flag for determining which type of relative
c                       hydraulic conductivity vs. pressure head
c                       curve is to be used.
c                       0 - step function.
c                       1 - front.
c                       2 - van genuchten model.
c
c
c   Intermediate variables:
c
      real * 8 fksmll ! small relative hydraulic conductivity.
      real * 8 phtran ! transition negative pressure head value.
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
      include 'seep.inc'
      common / unsat / uspar(2, mxmats), iuntyp
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
      fksmll = uspar(1, mat)
      phtran = uspar(2, mat)
c
c         Determine relative hydraulic conductivity.
c
      if (ph .ge. 0.d0) then
        fkrelf = 1.d0
      else 
        if (ph .gt. phtran) then
          fkrelf = (fksmll - 1.d0) * ph / phtran + 1.d0
        else
          fkrelf = fksmll
        endif
      endif
c
      return
      end
c----------------------------------------------------------------
      function fkrels (ph, mat)
c
c
c   Purpose:  This function computes the relative hydraulic
c             conductivity where the relative hydraulic
c             conductivity vs. pressure head curve is a step
c             function.  that is, the relative hydraulic
c             conductivity is one if the pressure head is positive
c             or zero and a small value otherwise.
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
      implicit real * 8 (a - h, o - z)
c
c   Input variables:
c
      real * 8 ph     ! pressure head.
c
      integer mat     ! material type number.
c
c
c   Output variables:
c
      real * 8 fkrels ! relative hydraulic conductivity from the step
c                       function model.
c
c
c   Common block variables:
c
c   / unsat /
c
      real * 8 uspar  ! parameters used for modeling the relative
c                       hydraulic conductivity vs. pressure head
c                       curve.  for the step function model,
c                       uspar(1) - small relative hydraulic
c                                  conductivity.
c
      integer iuntyp  ! flag for determining which type of relative
c                       hydraulic conductivity vs. pressure head
c                       curve is to be used.
c                       0 - step function.
c                       1 - front.
c                       2 - van genuchten model.
c
c
c   Intermediate variables:
c
      real * 8 fksmll ! small relative hydraulic conductivity.
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
      include 'seep.inc'
      common / unsat / uspar(2, mxmats), iuntyp
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
      fksmll = uspar(1, mat)
c
c         Determine relative hydraulic conductivity.
c
      if (ph .ge. 0.d0) then
        fkrels = 1.d0
      else 
        fkrels = fksmll
      endif
c
      return
      end
c----------------------------------------------------------------
      function fkrelv (ph, mat)
c
c
c   Purpose:  This function computes the relative hydraulic
c             conductivity according to the van genuchten model. 
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
      implicit real * 8 (a - h, o - z)
c
c   Input variables:
c
      real * 8 ph     ! pressure head.
c
      integer mat     ! material type number.
c
c
c   Output variables:
c
      real * 8 fkrelv ! relative hydraulic conductivity from the van
c                       genuchten model. 
c
c
c   Common block variables:
c
c   / unsat /
c
      real * 8 uspar  ! parameters used for modeling the relative
c                       hydraulic conductivity vs. pressure head
c                       curve.  for the van genuchten model,
c                       uspar(1) - parameter alpha > 0.
c                       uspar(2) - parameter n, where n > 1.
c
      integer iuntyp  ! flag for determining which type of relative
c                       hydraulic conductivity vs. pressure head
c                       curve is to be used.
c                       0 - step function.
c                       1 - front.
c                       2 - van genuchten model.
c
c
c   Intermediate variables:
c
      real * 8 alpha  ! parameter alpha.
      real * 8 fm     ! exponent m.
      real * 8 fn     ! parameter n.
      real * 8 seff   ! effective saturation.
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
      include 'seep.inc'
      common / unsat / uspar(2, mxmats), iuntyp
c
c - - - - - - - - - - - - - - - - - - - - - - - - -
c
      alpha = uspar(1, mat)
      fn = uspar(2, mat)
c
c         Determine relative hydraulic conductivity.
c
      if (ph .ge. 0.d0) then
        fkrelv = 1.d0
      else
        fm = 1.d0 - 1.d0 / fn
        seff = (1.d0 + (- alpha * ph) ** fn) ** (- fm)
        fkrelv = (1.d0 - (1.d0 - seff ** (1.d0 / fm)) ** fm) ** 2 *
     &           dsqrt (seff)
      endif
c
      return
      end
c ======================================================================
      subroutine getpath(name,path)
      implicit real * 8 (a - h, o - z)
      character*80 name,path
      integer*4 count

      count = 80
      path = name
200   if ((count .gt. 0) .and. (path(count:count) .ne. '\')) then
        path(count:count) = ' '
        count = count - 1
        go to 200
      endif
      end
c ======================================================================
      subroutine setpath(path,fname)
      implicit real * 8 (a - h, o - z)
      character*80 path,fname,newname
      integer*4 count,i
c
      count = 80
888   if ((count .gt. 0) .and. (path(count:count) .eq. ' ')) then
        count = count - 1
        go to 888
      endif
      count = count + 1
      i=1
      newname = path
998   if (fname(i:i) .eq. '"') then
        go to 999
      else if (i .eq. 80) then
	  go to 1000
	else
        i = i + 1
      endif
      go to 998
999   if (count .lt. 81  .and. i .lt. 81) then
        if (fname(i:i) .ne. '"') then
          newname(count:count) = fname(i:i)
          count = count + 1
        endif
        i = i + 1
        go to 999
      endif
	go to 1003
1000  i=2  ! JIG This is 2 to catch the | symbol?
1001  if (fname(i:i) .ne. ' ') then
	  go to 1002
	else
	  i = i + 1
	endif
      go to 1001
1002  if (count .lt. 81 .and. i .lt. 81) then
        if (fname(i:i) .ne. ' ') then
	    newname(count:count) = fname(i:i)
	    count = count + 1
	  endif
	  i = i + 1
	  go to 1002
	endif
1003  fname = newname
      end
c ======================================================================
      subroutine stopfile

      character*5 dummyvar
      character*10 stayopen
      common/emrl/stayopen

      if (stayopen.eq.'1') then
        write(*,*) ' '
        pause
      endif
      stop

      END
