      subroutine fft991_crm(a,work,trigs,ifax,inc,jump,n,lot,isign)
      real a(*),work(*),trigs(*)
      integer ifax(*)
      nfax=ifax(1)
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
      igo=60
      go to 40
   30 continue
      call fft99a_crm(a,work,trigs,inc,jump,n,lot)
      igo=60
   40 continue
      ia=1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm_crm(a(ia),a(ia+inc),work(1),work(2),trigs, &
        ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm_crm(work(1),work(2),a(ia),a(ia+inc),trigs, &
         2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
      if (isign.eq.-1) go to 130
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=1
      do 100 l=1,lot
      i=ibase
      j=jbase
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
  110 continue
      ib=n*inc+1
      do 120 l=1,lot
      a(ib)=0.0
      a(ib+inc)=0.0
      ib=ib+jump
  120 continue
      go to 140
  130 continue
      call fft99b_crm(work,a,trigs,inc,jump,n,lot)
  140 continue
      return
      end





      subroutine fftfax_crm(n,ifax,trigs)
      integer ifax(13)
      real trigs(1)
      data mode /3/
      call fax_crm (ifax, n, mode)
      i = ifax(1)
      call fftrig_crm (trigs, n, mode)
      return
      end





      subroutine fax_crm(ifax,n,mode)
      dimension ifax(*)
      nn=n
      if (iabs(mode).eq.1) go to 10
      if (iabs(mode).eq.8) go to 10
      nn=n/2
      if ((nn+nn).eq.n) go to 10
      ifax(1)=-99
      return
   10 k=1
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
   50 l=5
      inc=2
   60 if (mod(nn,l).ne.0) go to 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      if (nn.eq.1) go to 80
      go to 60
   70 l=l+inc
      inc=6-inc
      go to 60
   80 ifax(1)=k-1
      nfax=ifax(1)
      if (nfax.eq.1) go to 110
      do 100 ii=2,nfax
      istop=nfax+2-ii
      do 90 i=2,istop
      if (ifax(i+1).ge.ifax(i)) go to 90
      item=ifax(i)
      ifax(i)=ifax(i+1)
      ifax(i+1)=item
   90 continue
  100 continue
  110 continue
      return
      end





      subroutine fftrig_crm(trigs,n,mode)
      real trigs(*), pi, del, angle
      pi=2.0*asin(1.0)
      imode=iabs(mode)
      nn=n
      if (imode.gt.1.and.imode.lt.6) nn=n/2
      del=(pi+pi)/float(nn)
      l=nn+nn
      do 10 i=1,l,2
      angle=0.5*float(i-1)*del
      trigs(i)=cos(angle)
      trigs(i+1)=sin(angle)
   10 continue
      if (imode.eq.1) return
      if (imode.eq.8) return
      del=0.5*del
      nh=(nn+1)/2
      l=nh+nh
      la=nn+nn
      do 20 i=1,l,2
      angle=0.5*float(i-1)*del
      trigs(la+i)=cos(angle)
      trigs(la+i+1)=sin(angle)
   20 continue
      if (imode.le.3) return
      del=0.5*del
      la=la+nn
      if (mode.eq.5) go to 40
      do 30 i=2,nn
      angle=float(i-1)*del
      trigs(la+i)=2.0*sin(angle)
   30 continue
      return
   40 continue
      del=0.5*del
      do 50 i=2,n
      angle=float(i-1)*del
      trigs(la+i)=sin(angle)
   50 continue
      return
      end










      subroutine fft99a_crm(a,work,trigs,inc,jump,n,lot)
      real a(*),work(*),trigs(*)
      real c,s
      nh=n/2
      nx=n+1
      ink=inc+inc
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
      do 10 l=1,lot
      work(ja)=a(ia)+a(ib)
      work(jb)=a(ia)-a(ib)
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   10 continue
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
      do 20 l=1,lot
      work(ja)=(a(ia)+a(ib))- &
         (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(jb)=(a(ia)+a(ib))+ &
         (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+ &
         (a(ia+inc)-a(ib+inc))
      work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))- &
         (a(ia+inc)-a(ib+inc))
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   20 continue
      iabase=iabase+ink
      ibbase=ibbase-ink
      jabase=jabase+2
      jbbase=jbbase-2
   30 continue
      if (iabase.ne.ibbase) go to 50
      ia=iabase
      ja=jabase
      do 40 l=1,lot
      work(ja)=2.0*a(ia)
      work(ja+1)=-2.0*a(ia+inc)
      ia=ia+jump
      ja=ja+nx
   40 continue
   50 continue
      return
      end





      subroutine fft99b_crm(work,a,trigs,inc,jump,n,lot)
      real work(*),a(*),trigs(*)
      real scale,c,s
      nh=n/2
      nx=n+1
      ink=inc+inc
      scale=1.0/float(n)
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
      do 10 l=1,lot
      a(ja)=scale*(work(ia)+work(ib))
      a(jb)=scale*(work(ia)-work(ib))
      a(ja+inc)=0.0
      a(jb+inc)=0.0
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   10 continue
      scale=0.5*scale
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
      do 20 l=1,lot
      a(ja)=scale*((work(ia)+work(ib)) &
        +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(jb)=scale*((work(ia)+work(ib)) &
        -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(ja+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
         +(work(ib+1)-work(ia+1)))
      a(jb+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1))) &
         -(work(ib+1)-work(ia+1)))
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   20 continue
      iabase=iabase+2
      ibbase=ibbase-2
      jabase=jabase+ink
      jbbase=jbbase-ink
   30 continue
      if (iabase.ne.ibbase) go to 50
      ia=iabase
      ja=jabase
      scale=2.0*scale
      do 40 l=1,lot
      a(ja)=scale*work(ia)
      a(ja+inc)=-scale*work(ia+1)
      ia=ia+nx
      ja=ja+jump
   40 continue
   50 continue
      return
      end



      subroutine vpassm_crm &
        (a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
      real a(*),b(*),c(*),d(*),trigs(*)
      real c1,c2,c3,c4,s1,s2,s3,s4
      real sin36/0.587785252292473/,cos36/0.809016994374947/, &
          sin72/0.951056516295154/,cos72/0.309016994374947/, &
          sin60/0.866025403784437/
      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo.gt.4) return
      go to (10,50,90,130),igo
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      do 20 l=1,la
      i=ibase
      j=jbase
      do 15 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
      i=i+inc3
      j=j+inc4
   15 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   20 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 40 k=la1,m,la
      kb=k+k-2
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      do 30 l=1,la
      i=ibase
      j=jbase
      do 25 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
      i=i+inc3
      j=j+inc4
   25 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   30 continue
      jbase=jbase+jump
   40 continue
      return
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      do 60 l=1,la
      i=ibase
      j=jbase
      do 55 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
      i=i+inc3
      j=j+inc4
   55 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   60 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 80 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      do 70 l=1,la
      i=ibase
      j=jbase
      do 65 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)= &
         c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
        -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)= &
         s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))) &
        +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)= &
         c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
        -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))) 
      d(jc+j)= &
         s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))) &
        +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
   65 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   70 continue
      jbase=jbase+jump
   80 continue
      return
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      do 100 l=1,la
      i=ibase
      j=jbase
      do 95 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
      i=i+inc3
      j=j+inc4
   95 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  100 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 120 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      do 110 l=1,la
      i=ibase
      j=jbase
      do 105 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)= &
         c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
        -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)= &
         s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
        +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)= &
         c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
        -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)= &
         s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
     	 +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)= &
         c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
        -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)= &
         s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
        +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  105 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  110 continue
      jbase=jbase+jump
  120 continue
      return
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
      do 140 l=1,la
      i=ibase
      j=jbase
      do 135 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
       -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
       +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
       +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
       -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
       -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
       +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
       +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
       -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  135 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  140 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 160 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      do 150 l=1,la
      i=ibase
      j=jbase
      do 145 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)= &
         c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
           -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) &
     	 -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
          +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(jb+j)= &
         s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
           -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) &
        +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
           +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(je+j)= &
         c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
           +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) &
        -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
           -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(je+j)= &
         s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i))) &
           +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))) &
        +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i))) &
           -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(jc+j)= &
         c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
           -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) &
        -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
     	    +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jc+j)= &
         s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
           -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) &
        +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
           +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      c(jd+j)= &
         c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
           +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) &
        -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
           -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jd+j)= &
         s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i))) &
           +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))) &
        +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i))) &
           -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      i=i+inc3
      j=j+inc4
  145 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  150 continue
      jbase=jbase+jump
  160 continue
      return
      end




        subroutine cosft_crm(a,work,trigs,ifax,inc,jump,n,lot,isign)
        implicit none
        real a(*),work(*),trigs(*)
        INTEGER ifax(*),inc,jump,n,lot,isign
!
! Cosine Fourier Transform on quarter-wavelength shifted data.
! The input parameters are the same as used in fft991.
! Coded by Marat Khairoutdinov based on the code "cosft2" from
! the "Numerical recipies" book.
!

        REAL(8) y(100000)
        REAL(8) sum,sum1,y1,y2,ytemp,coef
        REAL(8) theta,wi,wi1,wpi,wpr,wr,wr1,wtemp,PI
        REAL(8) wr0,wi0
        INTEGER i,ilot, jump1

        PI = acos(-1._8)

        theta=0.5_8*PI/n ! Initialize the recurrences.
        wr0=cos(theta)
        wi0=sin(theta)
        wpr=-2.0_8*wi0**2
        wpi=sin(2._8*theta)

!
! Forward Transform:
!


        if(isign.eq.-1)then

           do ilot = 1,lot

             jump1 = (ilot-1)*jump+1

             do i=1,n
               y(i) = a(jump1+(i-1)*inc)
             end do ! i

             wr1=wr0
             wi1=wi0
             do i=1,n/2
               y1=0.5*(y(i)+y(n-i+1))
               y2=wi1*(y(i)-y(n-i+1))
               y(i)=y1+y2
               y(n-i+1)=y1-y2
               wtemp=wr1 ! Carry out the recurrence.
               wr1=wr1*wpr-wi1*wpi+wr1
               wi1=wi1*wpr+wtemp*wpi+wi1
             enddo ! i

             do i=1,n
               a(jump1+(i-1)*inc) = y(i)
             end do ! i

           end do ! ilot


           call fft991_crm(a,work,trigs,ifax,inc,jump,n,lot,isign)

           do ilot = 1,lot

             jump1 = (ilot-1)*jump+1

             do i=3,n
               y(i) = a(jump1+(i-1)*inc)*dble(n)
             end do !
             y(1) = a(jump1)*float(n)
             y(2) = a(jump1+(n+1-1)*inc)*dble(n)


             wr=1.0d0
             wi=0.0d0
             do i=3,n,2 ! Even terms.
               wtemp=wr
               wr=wr*wpr-wi*wpi+wr
               wi=wi*wpr+wtemp*wpi+wi
               y1=y(i)*wr+y(i+1)*wi
               y2=-y(i+1)*wr+y(i)*wi
               y(i)=y1
               y(i+1)=y2
             enddo ! i
             sum=0.5*y(2) ! Initialize recurrence for odd terms
             do i=n,2,-2 ! Carry out recurrence for odd terms.
               sum1=sum
               sum=sum+y(i)
               y(i)=sum1
             enddo ! i

             coef = 2._8/dble(n)
             do i=1,n
               a(jump1+(i-1)*inc) = y(i)*coef
             end do ! i

           end do ! ilot

!
!   Inverse Transform:
!


        else if(isign.eq.1)then

           do ilot = 1,lot

             jump1 = (ilot-1)*jump+1

             do i=1,n
               y(i) = a(jump1+(i-1)*inc)*dble(n)/2.
             end do ! i

             ytemp=y(n)
             do i=n,4,-2 ! Form difference of odd terms.
               y(i)=y(i-2)-y(i)
             enddo  ! i
             y(2)=2.0*ytemp
             wr=1.0d0
             wi=0.0d0
             do i=3,n,2 ! Calculate Rk and Ik .
               wtemp=wr
               wr=wr*wpr-wi*wpi+wr
               wi=wi*wpr+wtemp*wpi+wi
               y1=y(i)*wr+y(i+1)*wi
               y2=y(i+1)*wr-y(i)*wi
               y(i)=y1
               y(i+1)=-y2
             enddo

             do i=3,n
               a(jump1+(i-1)*inc) = y(i)
             end do ! i
             a(jump1) = y(1)
             a(jump1+inc) = 0.
             a(jump1+(n+1-1)*inc) = y(2)
             a(jump1+(n+2-1)*inc) = 0.

           end do ! ilot


           call fft991_crm(a,work,trigs,ifax,inc,jump,n,lot,isign)

           do ilot = 1,lot

             jump1 = (ilot-1)*jump+1

             do i=1,n
               y(i) = a(jump1+(i-1)*inc)*0.5
             end do ! i

             wr1=wr0
             wi1=wi0
             do i=1,n/2 ! Invert auxiliary array.
               y1=y(i)+y(n-i+1)
               y2=(0.5/wi1)*(y(i)-y(n-i+1))
               y(i)=0.5*(y1+y2)
               y(n-i+1)=0.5*(y1-y2)
               wtemp=wr1
               wr1=wr1*wpr-wi1*wpi+wr1
               wi1=wi1*wpr+wtemp*wpi+wi1
             enddo

             coef = 2._8/dble(n)
             do i=1,n
               a(jump1+(i-1)*inc) = y(i) * coef
             end do ! i

           end do ! ilot

        endif! isign

        return
        END


