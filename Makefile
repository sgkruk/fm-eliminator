ALL:		parser.o lexer.o fm testcdd projection

GCCF1	= -Wall -Dgets=DONT_USE_GETS -Dlint $(DF) -g
GCCF2	= -Wshadow -Wpointer-arith -Wnested-externs -Winline
CFLAGS	= -I/usr/local/include  -I/opt/local/include -g  $(GCCF1) $(GCCF2)
LOPT	= /opt/local/lib
LUSR	= /usr/local/lib
ALLL    = -lgsl -lgslcblas -lcdd -lm
TEST:		fm
		./fm -p y2 <test.fm

TEST0:		fm
		./fm -v -r x <test0.fm

TESTCUBE:		fm
		./fm  -r x -r z <testcube.fm

COVER:		cover.mod
			glpsol --model cover.mod --output cover.solution

testcdd:	t.c
			cc -o testcdd t.c ${ALLL}

projection:	projection.c
			cc -o projection projection.c ${ALLL}

lexer.o:	lexer.l parser.c fm.h
		flex -o lexer.c lexer.l 
		cc -c ${CFLAGS} -o lexer.o lexer.c

parser.o:	parser.y fm.h
		bison -d parser.y -o parser.c 
		cc -c ${CFLAGS} parser.c

fm.o:		fm.c
		cc -c ${CFLAGS} fm.c 

fm:		fm.o lexer.o parser.o 
		cc -o fm fm.o lexer.o parser.o  /usr/local/lib/libcdd.a ${ALLL}
