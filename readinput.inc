	integer Maxlines
        parameter (maxlines=200)
        character*100 lines(maxlines)
        integer n_lines
        integer iunit
	character*50 fname
	logical usedefaults,showinput
        COMMON/READINPUT2/iunit,n_lines,lines,usedefaults,showinput,fname
