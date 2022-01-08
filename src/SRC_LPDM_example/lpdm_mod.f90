module lpdm_mod

!--------------------------------------------!
! lagrangian particle dispersion model       !
! by Ji Nie, Li Chen, Z Kuang EPS, Harvard    !
!--------------------------------------------!
! > xygrids no more than 999x999
! > In each subdomain, number of particles no more than 9999999
! > The subdomain should not be larger than 214, otherwise the way to encode id needs to be changed. The 214 limit is because the id is a 4-byte integer, see more explanation in lpdm_sorting subroutine.
! > Pointers cannot be used for MPI
! > When run with a large number of particles, in lpdm_transfer, mpi buff limit of slave proc may cause model to freeze.
! > Now it is treated by separating large data chunk into small chunks. If model freezes in lpdm_transfer, try to decrease n_tran_max, and adjust n_tran

   use vars
   use tracers
   use microphysics

   implicit none
   include 'mpif.h'
   private

  !-------------------------------------------------------------------
  ! Particle Description
  ! note that SAM is in the moving frame with ug,vg. So don't add ug,vg to u and v
  type :: particle
    real(4) :: x,  ures             ! x: grid index
    real(4) :: y,  vres             !
    real(4) :: z,  wres             !
    real(4),dimension(3) :: uvw_b   ! uvw at time t-1
  end type


  type :: root_particle
    integer(4) :: xyz_out
  end type


  !-------------------- lpdm prm

  integer, parameter :: lpdm_top     = 120   !70   !64   !70   ! 85  ! top level for mirror reflection b.c.
  integer, parameter :: lpdm_num     = 1200 !1125 !120 !4480 !6400 !1120 !700  !1120 ! 4000  !900!960  ! particles in one column
  integer, parameter :: n_tran_max   = 80   !100  ! 300 is fine, 400 is not; maximum size of the array that can be transfered at one time in the LPDM transfer
  integer, parameter :: n_tran       = 350  !200  ! maximum number of loops in the transfer, n_tran*n_tran_max is the maximum total number of particles that can be transfered

  !-------------------
  ! define an array that holds all the particles in the subdomain. The factor of 2 is included to accomendate potential particle number fluctuations as particles move around.
  integer(4), parameter :: p_num_total= nx_gl*ny_gl*lpdm_num                          ! total number of par
  integer(4), dimension(int(p_num_total/nsubdomains_x/nsubdomains_y*2.0)) :: par_id   ! id of particles

  type (particle), dimension(int(p_num_total/nsubdomains_x/nsubdomains_y*2.0))      :: par_ay ! ay stands for array
  type (root_particle), dimension(int(p_num_total/nsubdomains_x/nsubdomains_y*2.0)) :: par_ay_dummy ! for MPI passing to par_ay_root only

  ! par_ay_root is only allocated in master processor
  integer(4), dimension(:),allocatable           :: par_id_root       ! root id array
  type (root_particle), dimension(:),allocatable :: par_ay_root       ! root array
  type (root_particle), dimension(:),allocatable :: par_ay_root_hold  ! root array

  integer(4) :: sub_pn                                                 ! actual number of par in current subdomain
  integer(4), dimension(nsubdomains_x*nsubdomains_y) :: sub_pnum       ! actual number of par in each subdomain

  !----------------- moving in/out
  integer(4) :: par_nin(8)    ! num of the pars moving in
  integer(4) :: par_nout(8)   ! num of the pars moving out

  ! define a array with a guess of length
  integer(4), dimension(8,int(p_num_total/nsubdomains_x/nsubdomains_y*1.0))      :: par_id_out ! moving out
  type (particle), dimension(8,int(p_num_total/nsubdomains_x/nsubdomains_y*1.0)) :: par_ay_out ! moving out

  !-----------------
  integer(4) :: inewtype,inewtype2   ! data type for mpi transfer

  public :: lpdm_ini,lpdm_update,lpdm_output,lpdm_top,lpdm_num

   contains



  !-----------------------------------------------------------
  subroutine lpdm_ini

    character *120     :: filename
    real, dimension(3) :: r             ! random number
    integer            :: ranki, rankj  ! rank position
    integer            :: i,j,k,m
    integer            :: kl
    type(particle)     :: par
    real               :: wrk

    if (masterproc) write (*,*)  'calling init in lpdm'

    if (modulo(p_num_total,lpdm_sub) .ne. 0) then  ! lpdm_sub is the number of output files
      print *, 'Error in lpdm_ini: lpdm_sub is not compatible with p_num_total'
      call task_abort()
    end if


    filename = './OUT_LPDM/'//trim(case)//'_'//trim(caseid)//'_sub1.xyz'
    open(14,file=filename,status='old',ACTION="read",err=87)
    write(*,*) 'Error in lpdm_ini: remove old lpdm output files first'
    call task_abort()
    87 continue


    call start_MPI    ! initial MPI

    rankj = rank/nsubdomains_x + 1
    ranki = rank+1 - (rankj-1)*nsubdomains_x

    if (masterproc) then
      allocate(par_id_root(p_num_total))
      allocate(par_ay_root(p_num_total))
      allocate(par_ay_root_hold(p_num_total))
    end if

    sub_pn = 0
    do i = 1,nx
      do j = 1,ny
        do m = 1,lpdm_num

          CALL RANDOM_NUMBER(r)

          par%x = nx*(ranki-1)+i+r(1)
          par%y = ny*(rankj-1)+j+r(2)

          wrk = r(3)*(presi(1)-presi(lpdm_top+1))+presi(lpdm_top+1)
          k = 1
          do while (wrk.lt.presi(k))
            k = k+1
          end do
          par%z = k-1+(wrk-presi(k-1))/(presi(k)-presi(k-1))

          call intp_uvw(par,1)

          par%uvw_b  = (/par%ures,par%vres,par%wres/)

          sub_pn = sub_pn+1

          par_id(sub_pn) = sub_pn+rank*10000000  ! use 10000000 (integer) instead of 1e7(real) to avoid precision loss
          par_ay(sub_pn) = par

        enddo
      enddo
    enddo

  end subroutine lpdm_ini


  !-------------------------------------------------------------
  subroutine lpdm_update()

    type(particle) :: par
    real      :: r
    integer*4 :: i,k

    real :: xx,yy,zz,zreal,zrealudt

    integer :: il,jl,kl,kb
    integer :: ranki,rankj

    rankj = rank/nsubdomains_x + 1
    ranki = rank+1 - (rankj-1)*nsubdomains_x

    ! Compute positions of particles at time t+1 - Adams-Bashforth, 2nd order
    do i = 1,sub_pn
      par = par_ay(i)

      k = floor(par%z)

      par%x   = par%x + ( 1.5*(par%ures) - 0.5*(par%uvw_b(1)) )*dt/dx
      par%y   = par%y + ( 1.5*(par%vres) - 0.5*(par%uvw_b(2)) )*dt/dy
      ! par%z   = par%z + ( 1.5*(par%wres) - 0.5*(par%uvw_b(3)) )*dt/(zi(k+1) - zi(k))

      ! Zeyuan 02.17.2021
      ! when vertical levels are not uniform, we don't directly
      ! update par%z, instead we first update actually height and then calculate
      ! which vertical level that height belongs to.
      zreal = zi(k) + (zi(k+1) - zi(k))*(par%z - k)
      zrealudt = zreal + ( 1.5*(par%wres) - 0.5*(par%uvw_b(3)) )*dt
      if (zrealudt > zreal) then
        do while (zrealudt > zi(k+1))
          k = k+1
        end do
      else
        do while (zrealudt < zi(k))
          k = k-1
        end do
      endif
      par%z = k + (zrealudt - zi(k))/(zi(k+1) - zi(k))

      call checkbnd(par)   ! check global boundary

      par_ay(i) = par
    enddo

    ! Move particles in/out of domain
    call lpdm_transfer()

    do i = 1,sub_pn

      par = par_ay(i)

      par%uvw_b = (/par%ures,par%vres,par%wres/)
      call intp_uvw(par,par_id(i))

      par_ay(i) = par
    enddo

  end subroutine lpdm_update


  !--------------------------------------------------------------
  subroutine intp_uvw(par,id)

    implicit none

    TYPE (particle),  intent(inout):: par
    integer*4, intent(in) :: id

    real :: xx,yy,zz
    real :: delx,dely,delz

    integer :: il,jl,kl
    integer :: ih,jh,kh
    integer :: ranki,rankj

    rankj = rank/nsubdomains_x + 1
    ranki = rank+1 - (rankj-1)*nsubdomains_x

    xx = par%x-nx*(ranki-1)
    yy = par%y-ny*(rankj-1)
    zz = par%z

    il = floor(xx)
    jl = floor(yy)
    kl = floor(zz)
    kl = max(1,kl)

    ih = il+1
    jh = jl+1
    kh = kl+1
    kh = min(nz_gl,kh)

    delx = xx-il
    dely = yy-jl
    delz = zz-kl

    !---------------- u
    par%ures = u(il,jl,kl)*(1.-delx)+u(ih,jl,kl)*delx

    !---------------- v
    par%vres = v(il,jl,kl)*(1.-dely)+v(il,jh,kl)*dely

    !---------------- w
    par%wres = ( rhow(kl)*w(il,jl,kl)*(1.-delz)+rhow(kh)*w(il,jl,kh)*delz )/rho(kl)
    !par%wres = w(il,jl,kl)*(1.-delz)+w(il,jl,kh)*delz

  end subroutine intp_uvw


  !-------------------------------------------------------------
  subroutine lpdm_output

    character *120 :: filename
    character *6   :: a
    integer*4      :: xstart,xend,nn,i
    integer*4      :: output_xyz(p_num_total)
    real*4         :: output_x(p_num_total)
    real*4         :: output_y(p_num_total)
    real*4         :: output_z(p_num_total)

    type(root_particle):: par

    if (modulo(nstep,lpdm_ndt) == 0) then

      call lpdm_reduce()

      if (masterproc) then
        print *, 'writing lpdm, nstep=', nstep,'   day=', day

        call lpdm_sorting()

        do i = 1,p_num_total

          par = par_ay_root(i)
          output_xyz(i) = par%xyz_out

        enddo

        do nn = 1,lpdm_sub

          write(a,'(I6)') nn
          xstart = 1+(nn-1)*(p_num_total/lpdm_sub)
          xend   = nn*(p_num_total/lpdm_sub)

          filename = './OUT_LPDM/'//trim(case)//'_'//trim(caseid)//'_sub'//trim(adjustl(a))//'.xyz'
          open (14,file=filename,status='unknown',ACTION="WRITE",form='unformatted',position='append')
          write(14) (output_xyz(i),i=xstart,xend)
          close(14)

        enddo

      endif

    endif

  end subroutine lpdm_output


  !-------------------------------------------------------------
  subroutine checkbnd(par)

    type (particle), intent(inout):: par

    if (par%x >= nx_gl+1) par%x = par%x-nx_gl
    if (par%y >= ny_gl+1) par%y = par%y-ny_gl
    if (par%x < 1)        par%x = par%x+nx_gl
    if (par%y < 1)        par%y = par%y+ny_gl

    if (par%z > nz_gl)      par%z = 2.*nz_gl-par%z
    if (par%z > lpdm_top+1) par%z = 2.*(lpdm_top+1)-par%z
    if (par%z < 1 )         par%z = 2 - par%z

  end subroutine checkbnd


  !-------------------------------------------------------------
  subroutine lpdm_sorting()

    ! Sort par_ay_root into ascending sequence using par_id
    ! The particles are labeled at the beginning of the run. This subroutine sorts them
    ! so they are in the same order as they had at the beginning.
    integer*4 :: i,id
    integer*4 :: pos1,pos2,pos

    par_ay_root_hold = par_ay_root

    do i = 1,p_num_total

      id=par_id_root(i)

      ! decode id
      ! pos1 is the index of subdomain the particle locates in, while pos2 is the index of particle within the current subdomain.

      pos1 = int(id/10000000) ! id is a 4-byte integer with the maximum possible value of 2147483647, so pos1 cannot be larger than 214
      pos2 = id-pos1*10000000
      pos  = pos1*(p_num_total/nsubdomains_x/nsubdomains_y)+pos2

      ! par_id_root(pos)=par_id_root(i)
      par_ay_root(pos)=par_ay_root_hold(i)

    enddo

  end subroutine lpdm_sorting


  !---------------------------------------------------!
  ! MPI                                               !
  !---------------------------------------------------!

  !----------------------------------
  subroutine start_MPI()
    ! create the date type for MPI data transform
    integer(4) :: iextent
    integer(4) :: iblock(8)
    integer(4) :: idisp(8)
    integer(4) :: itype(8)
    integer(4) :: ierr
    integer(4) :: iblock2(2)
    integer(4) :: idisp2(2)
    integer(4) :: itype2(2)

    ! Define first MPI derived type
    iblock = (/1,1,1,1,1,1,3,1/)

    ! Displacements
    idisp = (/0,4*1,4*2,4*3,4*4,4*5,4*8,4*9/)  ! in bytes

    ! Variable types are real
    itype = (/MPI_REAL,MPI_REAL,MPI_REAL,MPI_REAL,MPI_REAL,MPI_REAL,MPI_REAL,MPI_UB/)

    call MPI_TYPE_STRUCT(8,iblock,idisp,itype,inewtype,ierr)
    call MPI_TYPE_EXTENT(inewtype,iextent,ierr)
    call MPI_TYPE_COMMIT(inewtype,ierr)

    ! Define second MPI derived type
    iblock2 = (/1,1/)

    ! Displacements
    idisp2 = (/0,4*1/)  ! in bytes

    ! Variable types are integers
    itype2 = (/MPI_INTEGER4,MPI_UB/)
    call MPI_TYPE_STRUCT(2,iblock2,idisp2,itype2,inewtype2,ierr)
    call MPI_TYPE_EXTENT(inewtype2,iextent,ierr)
    call MPI_TYPE_COMMIT(inewtype2,ierr)

  end subroutine start_MPI


  !-------------------------------------------------------------
  subroutine lpdm_reduce
  ! gather all particle data, form a single array in master proc.

    integer(4) :: isend, irecv, ierr
    integer(4) :: istatus(MPI_STATUS_SIZE)
    integer(4) :: pos,i

    if (.not. masterproc ) then

      call MPI_ISEND(sub_pn,1,MPI_INTEGER4,0,4,MPI_COMM_WORLD,isend,ierr)
      call MPI_WAIT(isend,istatus,ierr)

    end if

    if (masterproc) then

      sub_pnum(1) = sub_pn
      do i = 2,nsubdomains_x*nsubdomains_y

        call MPI_IRECV(sub_pnum(i),1,MPI_INTEGER4,i-1,4,MPI_COMM_WORLD,irecv,ierr)
        call MPI_WAIT(irecv,istatus,ierr)

      end do

    end if

    if (masterproc .and. (sum(sub_pnum) .ne. p_num_total)  ) then

      print *,'Error in lpdm_reduce:total number of particles changed. bugs must be there.'
      call task_abort()

    endif

    if (masterproc) then

      do i = 1,sub_pn
        par_id_root(i) = par_id(i)

        !---------! |x|x|x|y|y|y|z|z|z  (grid index)
        par_ay_root(i)%xyz_out = floor(par_ay(i)%x)*1000000+floor(par_ay(i)%y)*1000+floor(par_ay(i)%z)
      enddo

    endif


    if (.not. masterproc) then

      do i = 1,sub_pn
        !---------! |x|x|x|y|y|y|z|z|z  (grid index)
        par_ay_dummy(i)%xyz_out = floor(par_ay(i)%x)*1000000+floor(par_ay(i)%y)*1000+floor(par_ay(i)%z)
      enddo


      call MPI_ISEND(par_ay_dummy,sub_pn,inewtype2,0,5,MPI_COMM_WORLD,isend,ierr)
      call MPI_WAIT(isend,istatus,ierr)
      call MPI_ISEND(par_id,sub_pn,MPI_INTEGER4,0,6,MPI_COMM_WORLD,isend,ierr)
      call MPI_WAIT(isend,istatus,ierr)

    endif

    if (masterproc) then

      pos = sub_pn + 1

      do i = 2,nsubdomains_x*nsubdomains_y

        call MPI_IRECV(par_ay_root(pos),sub_pnum(i),inewtype2,i-1,5,MPI_COMM_WORLD,irecv,ierr)
        call MPI_WAIT(irecv,istatus,ierr)
        call MPI_IRECV(par_id_root(pos),sub_pnum(i),MPI_INTEGER4,i-1,6,MPI_COMM_WORLD,irecv,ierr)
        call MPI_WAIT(irecv,istatus,ierr)
        pos = pos+sub_pnum(i)

      enddo

    endif

  end subroutine lpdm_reduce


  !-------------------------------------------------
  subroutine lpdm_transfer

    integer(4) :: isend, irecv, ierr
    integer(4) :: istatus(MPI_STATUS_SIZE)
    integer    :: ipos,jpos,irank
    integer    :: neighrank(8)
    integer*4  :: tmp,n
    integer*4  :: i1,i2,st,ed,j,i
    integer    :: iproc,ii
    integer*4,dimension(8,n_tran) :: par_nout_bf,par_nin_bf

    par_nout   = 0
    par_id_out = 0
    par_nin    = 0

    neighrank = (/ranknn,rankne,rankee,rankse,rankss,ranksw,rankww,ranknw/)

    ! prepare moving out particles
    n = 0
    do i = 1,sub_pn
      ipos = int(par_ay(i)%x-1)/nx+1
      jpos = int(par_ay(i)%y-1)/ny+1

      ! double periodic b.c.
      if (ipos > nsubdomains_x) ipos = ipos-nsubdomains_x
      if (jpos > nsubdomains_y) jpos = jpos-nsubdomains_y

      irank = (jpos-1)*nsubdomains_x+ipos-1
      if (irank .ne. rank) then   ! moving out

        tmp = 0

        ! the neighbors are clockwise numbered
        if (irank == ranknn) tmp = 1
        if (irank == rankne) tmp = 2
        if (irank == rankee) tmp = 3
        if (irank == rankse) tmp = 4
        if (irank == rankss) tmp = 5
        if (irank == ranksw) tmp = 6
        if (irank == rankww) tmp = 7
        if (irank == ranknw) tmp = 8

        if (irank < 0 .or. irank > (nsubdomains_x*nsubdomains_y-1).or.tmp == 0 ) then
          print *,'Error in lpdm_transfer:no neighbor found.', irank,par_ay(i)%x,par_ay(i)%y,ipos,jpos,rank,ranknn,rankne,rankee,rankse,rankss,ranksw,rankww,ranknw
          call task_abort()
        endif
        par_nout(tmp)=par_nout(tmp)+1

        if (par_nout(tmp) > int(p_num_total/nsubdomains_x/nsubdomains_y*1.0)) then
          print *,'Error in lpdm_transfer:par_ay_out array is outbounded. shoule change the size setting'
          call task_abort()
        endif

        par_ay_out(tmp,par_nout(tmp) ) = par_ay(i)
        par_id_out(tmp,par_nout(tmp) ) = par_id(i)

      else   ! not moving out, shrink the array

        n = n+1
        par_ay(n) = par_ay(i) !zk: this seems okay as n is alway <= i
        par_id(n) = par_id(i)
      endif
    enddo

    do iproc = 0,nsubdomains_x*nsubdomains_y-1
      if(rank.eq.iproc) then ! processor that's currently sending data
        do i = 1,8
          call MPI_ISEND(par_nout(i),1,MPI_INTEGER4,neighrank(i),1,MPI_COMM_WORLD,isend,ierr)
          call MPI_WAIT(isend,istatus,ierr)
        end do
      endif

      if(rank.ne.iproc) then ! processor that's not currently sending data
        do i = 1,8
          if(neighrank(i).eq.iproc) then !I'm a neighbor to the processor that's currently sending data
            call MPI_IRECV(par_nin(i),1,MPI_INTEGER4,neighrank(i),1,MPI_COMM_WORLD,irecv,ierr)
            call MPI_WAIT(irecv,istatus,ierr)
          endif
        end do
      endif

      call task_barrier()
    enddo

    sub_pn = n + sum(par_nin)
    n = n + 1

    if (sub_pn>int(p_num_total/nsubdomains_x/nsubdomains_y*2.0)) then
      print *,'Error in lpdm_transfer:par_ay array is outbounded. should change the size setting.'
      call task_abort()
    endif

    if (max(maxval(par_nout),maxval(par_nin))>n_tran_max*n_tran) then
      print *,'Error in lpdm_transfer:n_tran is too small'
      call task_abort()
    endif

    par_nout_bf = 0
    par_nin_bf  = 0
    do i = 1,8

      i1 = par_nout(i)/n_tran_max
      i2 = par_nout(i)-i1*n_tran_max

      par_nout_bf(i,1:i1) = n_tran_max
      par_nout_bf(i,i1+1) = i2

      i1 = par_nin(i)/n_tran_max
      i2 = par_nin(i)-i1*n_tran_max

      par_nin_bf(i,1:i1) = n_tran_max
      par_nin_bf(i,i1+1) = i2

    enddo

    do j = 1,n_tran    ! maximum loop
      do iproc = 0,nsubdomains_x*nsubdomains_y-1
        if(rank.eq.iproc) then ! work on current processor
          do i = 1,8
            if ( par_nout_bf(i,j) > 0 ) then

              st = (j-1)*n_tran_max+1
              ed = (j-1)*n_tran_max+par_nout_bf(i,j)
              call MPI_ISEND(par_ay_out(i,st:ed),par_nout_bf(i,j),inewtype,neighrank(i),2,MPI_COMM_WORLD,isend,ierr)
              call MPI_WAIT(isend,istatus,ierr)
              call MPI_ISEND(par_id_out(i,st:ed),par_nout_bf(i,j),MPI_INTEGER4,neighrank(i),3,MPI_COMM_WORLD,isend,ierr)
              call MPI_WAIT(isend,istatus,ierr)

            endif
          enddo
        endif

        if(rank.ne.iproc) then ! processor that's not currently sending data
          do i = 1,8
            if(neighrank(i).eq.iproc) then !I'm a neighbor to the processer that's currently sending data
              if (par_nin_bf(i,j) > 0) then

                call MPI_IRECV(par_ay(n),par_nin_bf(i,j),inewtype,neighrank(i),2,MPI_COMM_WORLD,irecv,ierr)
                call MPI_WAIT(irecv,istatus,ierr)
                call MPI_IRECV(par_id(n),par_nin_bf(i,j),MPI_INTEGER4,neighrank(i),3,MPI_COMM_WORLD,irecv,ierr)
                call MPI_WAIT(irecv,istatus,ierr)
                n = n + par_nin_bf(i,j)

              end if
            endif
          end do
        endif

        call task_barrier()

      enddo

    enddo

  end subroutine lpdm_transfer


end module lpdm_mod
