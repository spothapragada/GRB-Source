CC = g++
CFLAGS = -g -m32
LDFLAGS = -g -m32
binmodelfit_OBJECTS = binmodelfit.o
BinData_OBJECTS = BinData.o
FakeGRB_OBJECTS = FakeGRB.o
.SUFFIXES: .cc .cpp .c .o

All: FakeGRB BinData binmodelfit 

binmodelfit: $(binmodelfit_OBJECTS) 
	g++ -o binmodelfit $(binmodelfit_OBJECTS) $(CFLAGS) $(LDFLAGS) -lgsl -lcblas

BinData: $(BinData_OBJECTS)
	g++ -o BinData $(BinData_OBJECTS) $(CFLAGS) $(LDFLAGS) -lgsl

FakeGRB: $(FakeGRB_OBJECTS)
	g++ -o FakeGRB  $(FakeGRB_OBJECTS) $(CFLAGS) $(LDFLAGS) -lgsl  

clean:
	rm -f *.o

.c.o:
	$(CC) $(CFLAGS) -c $<

.cc.o:
	$(CC) $(CFLAGS) -c $<

.cpp.o:
	$(CC) $(CFLAGS) -c $<