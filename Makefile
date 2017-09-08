CPPFLAGS = -O2
CXX= g++
TARGET = lmprdf.out
OBJS = makeRDF.o main.o snapshot.o

all: $(TARGET)
$(TARGET): $(OBJS)
	g++  -o $(TARGET) $(OBJS)
lmprdf.out: makeRDF.o main.o snapshot.o
makeRDF.o: makeRDF.cpp makeRDF.h snapshot.h
main.o: main.cpp snapshot.h makeRDF.h
snapshot.o: snapshot.cpp snapshot.h


clean:
	rm -f *.bak
	rm -f *.map
	rm -f *.o
	rm $(TARGET)
	@echo "파일을 삭제했습니다."
