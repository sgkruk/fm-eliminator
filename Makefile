ALL:	fm

CDD     =  /home/sgkruk/software/cddlib-094f/lib-src
GCCF1	= -Wall -Dgets=DONT_USE_GETS -Dlint $(DF) -g
GCCF2	= -Wshadow -Wpointer-arith -Wnested-externs -Winline
CFLAGS	=  $(GCCF1) $(GCCF2) -I$(CDD)
ALLL    = -lgsl -lgslcblas -lcdd -lm

TESTCUBE:		fm
		./fm  -i testcube.fm -o junk  -p x1
		diff -w -q junk testcube.out
		rm junk

lexer.o:	lexer.l parser.c fm.h
		flex -o lexer.c lexer.l 
		cc -c ${CFLAGS} -o lexer.o lexer.c

parser.o:	parser.y fm.h
		bison -d parser.y -o parser.c 
		cc -c ${CFLAGS} parser.c

fm.o:		fm.c
		cc -c ${CFLAGS} fm.c 

fm:		fm.o lexer.o parser.o 
		cc -o fm fm.o lexer.o parser.o  ${ALLL}
