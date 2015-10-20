c-----------------------------------------------------------------------
c  Contains subrs for Cray inlining ( cft77 ... -i ...).
c  Some are Cray Fortran intrinsic functions, but this runs quicker.
c-----------------------------------------------------------------------
c
      subroutine gather (n, a, b, index)
 
c        Gathers indexed elements of array b into consecutive elements
c        of array a
 
      dimension a(*), b(*), index(*)
 
      do 10 i=1,n
        a(i) = b(index(i))
  10  continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine scatter (n, a, index, b)
 
c        Scatters consecutive elements of array b into indexed elements
c        of array a
 
      dimension a(*), index(*), b(*)
 
      do 10 i=1,n
        a(index(i)) = b(i)
  10  continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine reseti (ia, n, ival)
 
c        Sets integer array ia(n) to ival
 
      dimension ia(*)
 
      do 10 i=1,n
        ia(i) = ival
   10 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine resetr (a, n, val)
 
c        Sets real array a(n) to val
 
      dimension a(*)
 
      do 10 i=1,n
        a(i) = val
   10 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine scopy (n, a, inca, b, incb)
 
c        Copies array a to b, incrementing by inca and incb
c        respectively, up to a total length of n words of a.
 
      dimension a(*), b(*)
 
      ib = 1
      do 10 ia=1,n,inca
        b(ib) = a(ia)
        ib = ib + incb
   10 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      function sdot (n, a, inca, b, incb)
 
c        Forms dot product pf two vectors a and b, incrementing by
c        inca and incb respectively, up to a total length of n words
c        of a.
 
      dimension a(*), b(*)
 
      x = 0.
      ib = 1
      do 10 ia=1,n,inca
        x = x + a(ia)*b(ib)
        ib = ib + incb
   10 continue
      sdot = x
 
      return
      end
c
c-----------------------------------------------------------------------
c
      function ssum (n, a, inc)
 
c        Sums elements of array a, incrementing by inc,
c        up to a total length of n elements.
 
      dimension a(*)
 
      x = 0.
      do 10 i=1,n,inc
        x = x + a(i)
   10 continue
      ssum = x
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine zero (a, n)
 
c        Zeros n words starting at a(1)
 
      dimension a(*)
 
      do 10 i=1,n
        a(i) = 0.
   10 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine izero (ia, n)
 
c        Zeros n words starting at ia(1)
 
      dimension ia(*)
 
      do 10 i=1,n
        ia(i) = 0
   10 continue
 
      return
      end
 
c-----------------------------------------------------------------------
