AR=ar rc
RANLIB=ranlib

OBJS = sqlite3.o

sqlite3.o: sqlite3.c
	${CC} -fPIC ${CFLAGS} -c sqlite3.c -o sqlite3.o

clean:
	-rm sqlite3.o libsqlite3.a

libsqlite3.a: $(OBJS)
	$(AR) $@ $(OBJS)
	-@ ($(RANLIB) $@ || true) >/dev/null 2>&1

all:libsqlite3.a
