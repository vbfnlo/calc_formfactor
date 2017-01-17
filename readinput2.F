c*************************************************************************  
c     This file contains subroutines for reading input files. 
c     The input files are in the format: key = value 
c*************************************************************************  
c   LIST OF ALL FUNCTIONS AND SUBROUTINES IN THIS FILE:
c
c     SUBROUTINE loadfile(filename)
c     SUBROUTINE closefile
c     SUBROUTINE read_Real(key, value, default_value)
c     SUBROUTINE read_Int(key, value, default_value)
c     SUBROUTINE read_IntList(key, N, values, default_values)
c     SUBROUTINE read_Logical(key, value, default_value)
c     SUBROUTINE read_String(key, value, default_value)
c     FUNCTION   GetInputPath()
c
c*************************************************************************  



c*************************************************************************  
      subroutine loadfile(filename,show_input)
c*************************************************************************  
c     Opens a textfile and reads the data in an internal buffer.
c     End-of-line comments, starting with ! or //, are removed.
c     Do only process one single file at a time.
c     DON"T FORGET TO CALL CLOSEFILE AFTERWARDS.
c*************************************************************************  
      implicit none

      character*(*) filename
     
#include "readinput.inc"
      
      integer i,le,k,j
      character*100 line
      logical show_input

      character*100 GetInputPath
      character*100 path
      external GetInputPath
      

      iunit=10
      usedefaults=.false.
      showinput=show_input
c add input file path to filename
      path=GetInputPath()
      path = path(1:len_trim(path))//"/"//filename
c Open file
      open (unit=iunit,
     &      file=path(1:len_trim(path)),
     &      status="old",err=77)

      fname = filename

c read text lines
      n_lines=0
 50   read(iunit,"(a100)",end=99,err=98) lines(n_lines+1)
      n_lines=n_lines+1
      if (n_lines.lt.maxlines) GOTO 50
 98   CONTINUE
      print *," Read error on file",filename
      print *," Using default values."
      usedefaults=.true.
 99   CONTINUE
      n_lines=n_lines+1

c remove comments
      do i=1,n_lines
         line=lines(i)
         if (INDEX(line,"!").ne.0) then
            line=line(1:INDEX(line,"!")-1)
         endif
         if (INDEX(line,"//").ne.0) then
            line=line(1:INDEX(line,"//")-1)
         endif
c upper case everything in front of a "="
         le = index(line,"=")
         if (le.gt.0) then
            do j=1,le-1
               k=ichar(line(j:j))
               if(k.ge.97.and.k.le.122) then
                  k=ichar(line(j:j))-32   
                  line(j:j)=char(k)
               endif
            enddo
         endif
         lines(i)=line
      enddo
cc     if (showinput) write(*,"(A,A,A)"), " reading input from file ",filename," : "
      
      return
 77   CONTINUE
      print *," Open error on file ",filename
      print *," Using default values."
      usedefaults=.true.

      return
      end

c*************************************************************************  
      subroutine closefile
c*************************************************************************  
c     Closes the current file.
c*************************************************************************  
      implicit none

#include "readinput.inc"
      
      close(iunit)
      return
      end

c*************************************************************************  
      SUBROUTINE read_Real(key, value, default_value)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"

      character*(*) key
      real*8 value, default_value
      character*100 line
      integer i,le,lk,lt
     
      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk))
     &       .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9) )) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value       
               if (showinput) write(*,"(T4,A,T25,A,G11.5)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      print *,"Warning: cannot read value for key = ",key,
     &      ". Using default value : ",default_value
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************  
      SUBROUTINE read_Int(key, value, default_value)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"

      character*(*) key
      
      integer value, default_value
      character*100 line
      integer i,le,lk,lt
     
      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk))
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value
               if (showinput) write(*,"(T4,A,T25,A,I5)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      print *,"Warning: cannot read value for key = ",key,
     &      ". Using default value : ",default_value
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************  
      SUBROUTINE read_IntList(key, N, values, default_values)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"

      character*(*) key
      
      integer N
      integer values(*), default_values(*)
      character*100 line
      integer i,j,le,lk,lt

      do i=1,N
         values(j)=default_values(j)
      enddo
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk))
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) (values(j),j=1,N)
               if (showinput) then 
c     maximum of 50 entries hardwired in format statement
                  write(*,"(T4,A,T25,A,50(I4))") key," = ",(values(j),j=1,N)
               endif
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      print *,"Warning: cannot read values for key = ",key,"."
      print *,"Using default values : ",(default_values(j),j=1,N)
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************  
      SUBROUTINE read_Logical(key, value, default_value)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"
      
      character*(*) key
      Logical value, default_value
      character*100 line
      integer i,le,lk,lt
     
      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk))
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then
               line=" "//line(le+1:len(line))
               read(line,*,ERR=10,END=20) value
               if (showinput) write(*,"(T4,A,T25,A,L5)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      print *,"Warning: cannot read value for key = ",key,
     &      ". Using default value : ",default_value
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END

c*************************************************************************  
      SUBROUTINE read_String(key, value, default_value)
c*************************************************************************  
      implicit none
     
#include "readinput.inc"

      character*(*) key
      character*(*) value, default_value
      character*100 line
      integer i,le,lk,lt
     
      value=default_value
      if (usedefaults) RETURN
      do i=1,n_lines
         line = lines(i)
         lk=INDEX(line,key)
         if (lk.ne.0) then
            le=INDEX(line,"=")
            lt = lk+LEN(key) 
            if (((le.ne.0).and.(le.gt.lk))
     &           .and.(line(lt:lt).eq." ".or.line(lt:lt).eq."=".or.line(lt:lt).eq.achar(9))) then           
               line=" "//line(le+1:len(line))
               read(line,*,ERR=20,END=20) value
               if (showinput) write(*,"(T4,A,T25,A,A)") key," = ",value
               RETURN
            endif
         endif
      enddo

 10   CONTINUE
      print *,"Warning: cannot read value for key = ",key,
     &      ". Using default value : ",default_value
      RETURN

 20   CONTINUE
      print *,"Input Error in line ",i," of ",fname
      print *,""
      GOTO 10

      RETURN
      END


c*************************************************************************  
      character*100 FUNCTION GetInputPath()
c*************************************************************************  
c     returns the path to input files
c     it is taken from the command line option input=my_input_path
c     or, if no command line option is given, from the 
c     environment variable VBF_INPUT_PATH
c*************************************************************************  
      implicit none
      character*100  s
      integer narg,i,j,k,lk,le,lt
      
c first, check for command line option
      narg = IARGC()
      if (narg.gt.0)then
         do i=1,narg
            call GETARG(i , s)
c upper case everything in front of a "="
            le = index(s,"=")
            if (le.gt.0) then
               do j=1,le-1
                  k=ichar(s(j:j))
                  if(k.ge.97.and.k.le.122) then
                     k=ichar(s(j:j))-32   
                     s(j:j)=char(k)
                  endif
               enddo
            endif
c look for input= statment
            lk=INDEX(s,"INPUT")
            if (lk.ne.0) then            
            lt = lk+LEN("INPUT") 
            if (((le.ne.0).and.(le.gt.lk))
     &           .and.(s(lt:lt).eq." ".or.s(lt:lt).eq."=".or.s(lt:lt).eq.achar(9))) then           
               s=s(le+1:len(s))
               GetInputPath = s
               RETURN
            endif
         endif
            
         enddo
      endif

c next, check for environment variable
c      call GETENV( "VBF_INPUT_PATH" , s )
c      if (len_trim(s).gt.0) then
c         GetInputPath = s
c         RETURN
c      endif

      GetInputPath = "."      
      RETURN
      END


